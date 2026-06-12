# OPTIMJAM — résolution des amas denses (jams) dans FloeDyn

*Branche `optimjam2` — juin 2026. Référence : Q. Jouet (implémentation), avec assistance IA.*

---

## 1. L'idée en une phrase

> Dans un amas dense quasi-statique, presque tous les floes **veulent** bouger (le forçage les pousse)
> mais **ne peuvent pas** (la géométrie les bloque). OPTIMJAM identifie ces floes en observant leur
> **déplacement réel dans le temps**, les transforme temporairement en **ancres fixes** que le reste de
> l'amas voit comme des murs, résout le mouvement autour, et re-teste périodiquement leur mobilité.
> Les **forces**, elles, sont toujours calculées sur le système physique complet.

Activation : `--optim_jam 1 --jam_freeze 1` (les défauts correspondent au jeu de paramètres validé).
Avec `--optim_jam 0` (défaut), le comportement est **strictement identique à la baseline**.

---

## 2. Pourquoi : les deux pathologies des jams

Les simulations à amas denses (canal, silo, impasse ; des centaines de floes en contact) souffraient de
deux problèmes **distincts**, qu'il ne faut pas confondre :

**(A) Explosion du nombre de LCP** (coût par pas de temps). La stratégie historique des
« sous-graphes actifs » résout les contacts par petits paquets, en propageant l'énergie de proche en
proche jusqu'à convergence. Dans un amas dense, cela donne ~10⁴–10⁵ LCP par pas (mesuré : ~300 000
LCP/pas, ~22 s/pas, un run de 148 h à 1,16 milliard de LCP sans terminer). Cause racine : l'amas est
**hyperstatique** (plus de contacts que de degrés de liberté) → le réseau de forces est indéterminé →
les LCP sont dégénérés → le pivotage de Lemke cycle (d'où, historiquement, le bruit aléatoire `--bustle`
pour casser les égalités).

**(B) Effondrement du pas de temps adaptatif (« Zeno »)**. Le schéma résout les contacts, *puis*
applique la traînée pendant le déplacement : la traînée pousse des floes déjà en contact l'un vers
l'autre → interpénétration → rembobinage (`INTER`) → dt divisé par 2 → … → dt ~1e-7 →
`RECOVER STATES` (abandon du pas). Un meilleur solveur de contact **ne corrige pas** (B) : il faut
empêcher de bouger les floes qui, de toute façon, n'avanceront pas.

OPTIMJAM répond à (A) par un **solveur alternatif** (§3) et à (B) par le **gel temporel** (§4).

---

## 3. Le solveur Gauss-Seidel projeté (problème A)

`src/floe/lcp/solver/GS_solver.hpp`

C'est un **second** solveur de contact, à côté du Lemke historique (qui reste le solveur de référence
et le chemin par défaut). Mathématiquement, il résout **le même problème** que Lemke — les conditions
de complémentarité au niveau vitesse (Signorini : `0 ≤ pₙ ⊥ vitesse de séparation ≥ 0`, plus le
frottement de Coulomb), sur les **mêmes matrices** (l'assemblage `GraphLCP` est réutilisé : M, M⁻¹, J,
D, μ) — mais par une autre méthode :

- **Lemke** : méthode *directe* par pivotage. Exacte (sommet du polyèdre de complémentarité) mais
  fragile à la dégénérescence des systèmes hyperstatiques (cyclage).
- **Gauss-Seidel projeté** (PSOR / NSCD à la Jean-Moreau, cité dans l'article FloeDyn) : méthode
  *itérative* par relaxations successives contact par contact (projection Signorini sur le normal,
  puis projection sur l'intervalle de Coulomb `[−μpₙ, +μpₙ]`). Convergence linéaire vers *une*
  solution du même ensemble — **nativement tolérante à la dégénérescence** : là où Lemke doit choisir
  un sommet parmi un continuum, le PGS relaxe vers l'un d'eux sans état d'âme. En 2D, la facettisation
  du cône de Coulomb est exacte, donc les deux formulations ont le même ensemble de solutions.

Deux différences assumées avec le chemin Lemke :

1. **Restitution e = 0** (purement inélastique) : c'est la loi pertinente pour des contacts *tenus*.
   Les vraies collisions impulsives (restitution, Newton's cradle) restent routées sur Lemke.
2. **Résolution monolithique** : le GS résout la composante connexe entière comme *un seul* problème
   couplé — plus proche de la formulation théorique que la cascade séquentielle de petits LCP, dont on
   connaît le biais dissipatif (la cascade amortit plus qu'un solve global).

**Best-effort, dans la philosophie FloeDyn** : le solveur garde la meilleure solution rencontrée
(pénétration résiduelle minimale), avec tolérances *relatives à l'échelle* (les domaines FloeDyn
varient sur des ordres de grandeur), garde-fous `isfinite` et plafond de vitesse. Il ne « refuse » que
le déchet numérique avéré — auquel cas la composante part en **fallback Lemke** intégral (risque de
régression nul).

**Résultat (A)** : collisions ~2–9 ms/pas contre ~22 s/pas, sans bruit `bustle`.

---

## 4. Le gel temporel (problème B)

Boucle de gel dans `LCPManager::solve_contacts` (`src/floe/lcp/LCP_manager.hpp`) ; compteurs portés
par `KinematicFloe` (`src/floe/floes/kinematic_floe.hpp`, accesseurs `jam_*`) — volontairement **hors**
de `SpaceTimeState` : c'est du bookkeeping d'algorithme, pas de l'état physique.

### Pourquoi un critère temporel

Tous les critères *instantanés* ont été essayés et ont échoué : geler toute la composante (collages et
porte-à-faux non physiques), geler par vitesse post-solve (Zeno persiste), geler par déplacement prédit
`|v|·dt` (Zeno persiste). La raison de fond : un floe coincé dans un jam a une **vraie vitesse** — le
forçage le pousse, il *veut* avancer — et aucune mesure à l'instant t ne distingue « avance » de
« veut-avancer-mais-bloqué ». Le seul signal fiable est le **déplacement net effectif sur plusieurs
pas** : vitesse non nulle + progrès net nul = bloqué géométriquement.

### Le mécanisme

Chaque pas de temps, pour chaque floe d'une composante routée :

```
d = |position_réelle − position_de_référence|
si d > eps × diamètre   →  il bouge : référence ← position courante, compteur ← 0
sinon                   →  compteur += 1
                           si compteur ≥ N et pas un pas de sonde  →  GELÉ (set_jammed)
```

- **Gelé** = `dynamics_manager::move_floe` (`src/floe/dynamics/dynamics_manager.hpp:52`) saute
  l'intégration de **position** (la traînée, Coriolis etc. sont inchangés) ; et le solveur le traite en
  ancre (§5). Il est commité **au repos** (v = 0) : tenu immobile dans la structure ancrée — et la
  traînée ne peut pas accumuler une vitesse fantôme qui le transformerait en boulet de canon à la
  libération.
- **Jamais persistant** : `unjam_all_floes()` (`problem.hpp`, début de `step_solve`) remet tout le
  monde mobile à chaque pas ; le gel est re-décidé pas par pas. Une mauvaise décision ne vit qu'un pas.
- **Sondes (release)** : un floe gelé ne progresse pas, par construction — sans mécanisme dédié il
  resterait gelé à vie. Tous les `K` pas (décalage individuel par floe, pour éviter les relâchements
  collectifs synchronisés), un floe bloqué est **laissé mobile** pour prouver qu'il peut bouger.
  Encore coincé ? le solve cohérent lui donne v ≈ 0, la sonde ne coûte rien. Contrainte levée (arche
  effondrée, voisin parti) ? il bouge réellement → reset → libre. *C'est la sonde qui fait s'effondrer
  les piles physiquement.*

---

## 5. Le pipeline complet d'un pas de temps

```
step_solve (problem.hpp)
│
├─ 0. unjam_all_floes()                       — personne n'est gelé a priori
│
├─ 1. manage_collisions → LCPManager::solve_contacts (LCP_manager.hpp)
│     │   (rotation du cache warm-start, décompte du mode urgence)
│     │
│     └─ pour chaque composante connexe du graphe de contacts :
│        │
│        ├─ route_to_gs ?  =  ≥ min_contacts  ET  ancrée (contient un obstacle)
│        │                    ET  quasi-statique (approche max < rel_speed_max)
│        │   NON → chemin Lemke historique inchangé (chocs, petits clusters, floes libres)
│        │
│        ├─ a. DÉCISION DE GEL (temporelle, §4) — AVANT le solve
│        │
│        ├─ b. PASSE DYNAMIQUE (GS_solver, Mode::Dynamics)
│        │     gelés = ancres rigides (M⁻¹ = 0, v = 0) → les vitesses des mobiles sont
│        │     cohérentes avec ce qui bouge réellement (sinon : interpénétration au 1er
│        │     ordre contre les voisins épinglés → effondrement du dt).
│        │     Garde énergie : aucun floe ne sort avec plus d'énergie cinétique que la
│        │     composante entière n'en avait (analogue du test calcEc>1 de Lemke).
│        │     → commite les vitesses (gelés au repos). N'enregistre PAS d'impulsions.
│        │
│        ├─ c. PASSE FORCES (GS_solver, Mode::Forces)
│        │     masses VRAIES pour tous (seuls les vrais obstacles fixes), sur les vitesses
│        │     LIBRES d'avant résolution (la quantité de mouvement chargée par la traînée
│        │     est ce qui presse sur l'obstacle). Ne bouge personne.
│        │     → enregistre la CHAÎNE DE FORCES physique (modèle de fracture + visualisation).
│        │     (sauvegarde/restauration des vitesses autour, dans solve_contacts)
│        │
│        ├─ (b et c sont warm-startées : impulsions du pas précédent indexées par
│        │    (paire de floes, position du contact) → ~10-600 sweeps au lieu de ~5000-7000)
│        │
│        └─ échec de (b) → composante dégelée → fallback Lemke intégral
│
├─ 2. compute_time_step
│
├─ 3. safe_move_floe_group — move_floe saute la position des gelés ;
│     boucle INTER/rembobinage inchangée ; sur RECOVER : impression des paires en
│     interpénétration (diagnostic) + notify_recover (garde anti-cycle, §6)
│
└─ (fin de run) bilan #OPTIMJAM GS: commits, fallbacks, sweeps moyens/saturés par passe,
   clamps énergie, déclenchements d'urgence
```

Le principe directeur des passes b/c : **découpler la dynamique du diagnostic**. On s'autorise une
heuristique sur le *mouvement* (les ancres) parce qu'elle est bornée par construction (re-décidée
chaque pas, re-testée toutes les K pas) ; on ne s'autorise **aucune** approximation sur les *forces*
enregistrées (physique complète). L'histoire de la mise au point l'a confirmé : utiliser les
impulsions de la passe dynamique produisait des chaînes incohérentes et des forces divergentes (~1e21,
un floe pincé entre deux ancres infiniment rigides) — le problème masse-infinie déjà connu de Lemke
avec les obstacles.

---

## 6. Garde-fous (philosophie : « ne jamais bloquer »)

Tous calqués sur les garde-fous historiques (sortie de la boucle des sous-graphes actifs sans tout
résoudre, `RECOVER STATES`) : une solution imparfaite qui avance vaut mieux qu'une exacte
inatteignable.

| Garde | Déclencheur | Réponse | Où |
|---|---|---|---|
| Best-effort GS | solve non parfaitement convergé | commite la meilleure solution vue | `GS_solver.hpp` |
| Fallback Lemke | solution non-finie / divergente | composante dégelée, chemin Lemke | `solve_contacts` |
| Garde énergie | un floe sort avec plus d'EC que la composante n'en avait | rescaling chirurgical de ce floe | `GS_solver.hpp` (passe dyn.) |
| Diagnostic de paires | chaque `RECOVER` | imprime jusqu'à 5 paires en interpénétration (indice, `[frozen]`, position) | `detector.hpp` / `problem.hpp` |
| EMERGENCY anti-cycle | `RECOVER` répété **sans progrès du temps simulé** (cycle limite déterministe) | gèle 20 pas les floes à compteur ≥ 1 des composantes routées (les floes en chute libre ne sont **pas** épinglés) | `notify_recover`, `LCP_manager.hpp` |

Leçon empirique : même la version brutale de l'EMERGENCY (gel en bloc, abandonnée) ne *corrompait*
pas l'état — elle déformait transitoirement la trajectoire (bloc suspendu, fissure de cisaillement),
et la physique se ré-équilibrait dès le relâchement. Les gardes échangent de la fidélité transitoire
contre du progrès, jamais de la validité.

---

## 7. Paramètres

`--jam_params <min_contacts> <gs_max_iter> <rel_speed_max> <eps> <N> <K>`
**Défauts = jeu validé : `50 20000 0.5 3e-4 10 10`.**

| Paramètre | Rôle | Si ↑ | Si ↓ |
|---|---|---|---|
| `min_contacts` (50) | taille minimale d'une composante routée GS | on rate des jams moyens (restent Lemke) | overhead GS sur des cas que Lemke gère bien |
| `gs_max_iter` (20000) | plafond de sweeps (quasi gratuit grâce au warm-start) | meilleure convergence des chaînes de forces | chaînes tronquées → scintillement des impulsions |
| `rel_speed_max` (0.5 m/s) | porte quasi-statique du routage (vitesse d'*approche* aux contacts, le seul seuil dimensionné) | route des chocs que Lemke devrait traiter | les jams qui grincent restent sur Lemke |
| `eps` (3e-4) | seuil de progrès net, fraction du diamètre | gèle du mouvement lent légitime | le fluage passe sous le seuil → Zeno revient |
| `N` (10) | pas consécutifs sans progrès avant gel | plus fidèle, réagit plus lentement au blocage | gèle sur un creux transitoire |
| `K` (10) | période de sonde (release) | pack sur-figé, dégels tardifs | dégels réactifs (quasi gratuit depuis le solve cohérent) ; un peu plus de mobiles par pas |
| `--jam_freeze` (1) | 0 = GS seul, sans gel | — | pour comparer ; GS seul = Zeno garanti |
| `--jam_warmstart` (1) | 0 = cold start | — | pour mesurer le gain (~×10 sur le solve) |

**Budget d'erreur, pour les physiciens** : l'unique hypothèse de modélisation ajoutée est
*« un floe qui, sur N pas consécutifs, n'a pas progressé de plus de eps×diamètre est considéré
géométriquement bloqué ; on suspend l'intégration de sa position et on re-teste sa mobilité tous les
K pas. »* Tout le reste est soit exact (forces), soit du numérique standard (PGS vs pivotage sur le
même problème). Les trois nombres (eps, N, K) sont interprétables et fixent l'arbitrage
fidélité/vitesse. Pour resserrer vers la fidélité, ordre recommandé : K (10→5→3, réactivité du dégel),
puis N (10→15, geler plus tard), puis eps (3e-4→1e-4, prudence : trop bas = retour du Zeno) — un
curseur à la fois ; juges : vidéo naturelle, pas de retour des `RECOVER`, coût acceptable.

---

## 8. Conception des scènes (leçon importante)

Les derniers effondrements de dt observés n'étaient **pas** dans le solveur mais **fabriqués par la
géométrie des scènes** : murs en blocs carrés à 4 points séparés de petits interstices → le sommet du
polygone d'un cercle glissant plongeait dans la fente (~`(gap/2)·tan(π/n)`) puis **accrochait le coin
à 90° du bloc suivant** → interpénétration → tempêtes d'`INTER` (paires récidivistes mur/cercle dans
les logs). Règles, encodées dans `pack_creator.py` (`chamfered_block`, `blocks_along`,
`check_no_initial_overlap`) :

- coins de la face de contact **chanfreinés à 45°** (la jonction devient une rainure en V que le
  sommet remonte comme un contact incliné ordinaire) ;
- **interstice minuscule** (il ne sert qu'à ce que le détecteur ne voie pas de contact à t = 0 ;
  les obstacles ne bougent jamais) ;
- **faces subdivisées** (disques locaux serrés dans `OptimizedFloe` → distances de proximité moins
  pessimistes → dt moins bridé) ;
- offset `'normal'` pour les contours qui s'écartent de leur côté offset (coupes, entonnoirs),
  offset **fixe** pour les parois ondulées (flancs de blocs parallèles quelle que soit la courbure) ;
- les murs restent en **blocs séparés** : lecture des forces par bloc et décomposabilité MPI.

Scènes de validation : `in_calice` (coupe profonde, chaînes radiales — baseline A/B),
`in_entonnoir` (entonnoir borgne asymétrique 50°/70° — juge de la fidélité des arches),
`in_gosier` (parois sinusoïdales convergentes + rocher intrus — stress-test, arche sur le rocher).
Après correction des scènes : **zéro `RECOVER`** sur les runs.

---

## 9. Résultats (état au 2026-06-11)

- Problème A : collisions ~ms au lieu de ~22 s ; plus besoin de `--bustle`.
- Problème B : plus d'effondrement de dt sur les scènes corrigées ; simulations complètes
  (28 jours simulés en quelques heures) là où la baseline ne terminait pas.
- Chaînes de forces qualitativement correctes (arches en X ancrées base+parois, écrantage de type
  Janssen, cône déchargé sous la surface libre) ; forces croissant avec le forçage.
- Warm-start : ~×10 sur le solve à résidu égal. Garde énergie : plus d'« explosions ».

**Reste au backlog** : options CLI nommées + `max_iter` séparé dynamique/forces ; reset du tracking
jam et des caches warm-start sur RECOVER/fracture ; champs DIAG du log derrière un flag verbose ;
validation A/B quantitative vs Lemke (courbes d'énergie, débit de décharge) pour chiffrer le budget
d'erreur ; à terme, routage par régime physique (idée du directeur : sélection de solveur apprise —
désormais réaliste : deux solveurs et des métriques comparables existent).
