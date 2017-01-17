// IterLemke
// LCP solver based on original solver lemke.m from Paul L. Fackler & Mario J. Miranda
//
// Syntax:
// [z, w, info] = IterLemke(M, q)
// [z, w, info] = IterLemke(M, q, tol)
// [z, w, info] = IterLemke(M, q, tol, z0)
//
// Description:
// From a LCP system (M, q), it returns the solution z, the complementarity
// solution w and some diagnostic informations in the structure info :
// - info.sol_err : LCP error
// - info.num_err : numerical error (see Yamane & Nakamura, 2008)
// - info.min_num_error : minimal numerical error observed through the process
// - info.max_num_error : maximal numerical error observed through the process
// - info.error_code : error code that can take the following values :
//   - 0 <=> sucess;
//   - 1 <=> for the calculated numerical error, we cannot expect a better
//           solution;
//   - 2 <=> too much iterations without improvement of the solution.
//   - 3 <=> black hole : no way to have lemke.m return someting else than NaN
// - info.error_msg : error message
// - info.cnt : number of iteration (ie, the number of call to lemke.m)
// - info.nan_cnt : number of lemke that failed (that have return NaN)
// - info.lmk_iter : [min, mean, max] iterations in the lemke code
// - info.cpu : cpu time for the entire code
//
// Optional parameters:
// - tol : LCP error tolerance (1e-10 by default);
// - z0 : initial guess for the solution.

#ifndef FLOE_LCP_SOLVER_ITERLEMKE_HPP
#define FLOE_LCP_SOLVER_ITERLEMKE_HPP

#include "floe/lcp/lcp.hpp"

namespace floe { namespace lcp { namespace solver
{
#include "floe/lcp/lcp.hpp"

struct IL_info{
    double sol_err{Inf};
    double num_err{Inf};
    double min_num_err{Inf};
    double max_num_err{0};
    int error_code{0};
    std::string error_msg{'success'};
    int cnt{0};
    int nan_cnt{0};
    std::tuple<double, int, int> lmk_iter{{Inf, 0, 0}};
    double cpu{0};
}; 



int IterLemke(lcp_type& lcp, double tol=1e-10) //, decltype(lcp.z) z0)
{
    // MATLAB param : M, Q, tol, z0

    // Paramètres
    int best_timeout = 20; // Nombre d'itération max depuis la dernière meilleure solution
    double num_tol = 10; // On considère une solution bonne si son erreur est inférieure à  num_tol fois l'erreur numérique constatée
    // double lcp_tol = 1e-10; // Tolérance par défaut pour l'erreur LCP

    // Infos de diagnostique sur la solution trouvée
    IL_info info;

    // Solution initiale
    if (lcp.z == zero_vector<double>(lcp.dim))
    {
        lemke(lcp);
        //update info (useless ?)
    }
    
    // Erreur LCP initiale
    double sol_err = LCP_error(lcp);

    // Meilleure solution
    best = struct( ...
        'z', z, ...
        'w', M*z+q, ...
        'sol_err', sol_err, ...
        'num_err', NaN, ...
        'cnt', 0 ...
    );
    
    // Initialisation des pivots
    std::vector<bool> ind(lcp.dim, 0);
    
    // Tant que l'erreur LCP est supérieure à  la tolérance demandée 
    while (!info.error_code && sol_err > tol)
    {
        // à vrai dire, vu l'algo actuel, je ne sais même plus à  quoi sert
        // cette ligne ... mais ça fonctionne mieux que
        //   ind = z ~= 0;
        // Avant, les pivots s'enchaînaient et le xor permettait de suivre
        // l'évolution des pivots depuis le système originel, mais là  ...
        // ind = xor(ind, z ~= 0); // en C++
        for (std::size_t i = 0; i<ind.size(); ++i) ind[i] = (lcp.z(i) == 0);
        
        // Si la solution z était la bonne, M(ind,ind) serait inversible ...
        // Malheureusement, on en est loin, donc :
        auto inv_ind = getInvPrincipalSubMatrix( lcp.A );
        // On peut rajouter une tolérance sur le calcul du rang, genre
        // 1e-10 ... on le fait ?
        
        // Petite manip pour récupérer le champ de booléen correspondant
        tmp = find(ind);
        ind(:) = false;
        ind(tmp(inv_ind)) = true;
//         sum(ind)
        // On peut plutà´t faire ça en une ligne :
//         ind(ind) = inv_ind;
        
        // On pivote le système selon ces indices
        [M_tmp, q_tmp] = PivoteLCP( ind, ind, M, q );
        
//         spy([M_tmp,q_tmp(:,ones(1,10)),min(0,q_tmp(:,ones(1,10)))]);pause
        
        // Solution pour le système pivoté 
        [z_tmp, nan_cnt, iter] = lemke_nofail(M_tmp, q_tmp);
        w_tmp = M_tmp*z_tmp + q_tmp;
//         w_tmp = BestAb([M_tmp, q_tmp], [z_tmp ; 1]);

        // Mise à  jour des infos et verif du code d'erreur
        info = update(info, z, nan_cnt, iter);
        if info.error_code; break; end;

        // On reconstruit la solution équivalente pour le système originel
        z = z_tmp;
        z(ind) = w_tmp(ind);
        w = w_tmp;
        w(ind) = z_tmp(ind);
        
        // Erreur LCP et numérique
        sol_err = errLCP(z,M,q);
        num_err = norm(w - (M*z+q));
        info.min_num_err = min( info.min_num_err, num_err);
        info.max_num_err = max( info.max_num_err, num_err);
        
        // Mise à  jour de la meilleure solution
        if sol_err < best.sol_err
            best = struct( ...
                'z', z, 'w', w, 'cnt', info.cnt, ...
                'sol_err', sol_err, 'num_err', num_err ...
            );

            // On arrête là  si on considère que l'erreur numérique ne permettra
            // pas d'obtenir une meilleure solution
            if sol_err <= num_tol*num_err
                info.error_code = 1;
                info.error_msg = 'good solution compared to the numerical error';
                break;
            end
        end
        
        // On arrête si on n'a pas amélioré la solution depuis longtemps
        if info.cnt > best.cnt + best_timeout
            info.error_code = 2;
            info.error_msg = 'no improvement since long time';
            break;
        end
        
    }
    
    z = best.z;
    w = best.w;
    info.sol_err = best.sol_err;
    info.num_err = best.num_err;
    info.lmk_iter(2) = info.lmk_iter(2) / info.cnt;
    info.cpu = toc(info.cpu);
    
    warning(orig_state);
end
}

std::size_t getInvPrincipalSubMatrix( A , varargin ) // TODO types
{
    // A matrice carrée ...
    vector<bool> ind(A.size1(), 0);
    for (int i = 0; i<ind.size(); ++i){
        ind(i) = true;
        if rank( A(ind, ind) , varargin{2:end}) ~= sum(ind) ind[i] = false; // TODO C++ -> A(ind, ind) ? rank ?
    }
}
    

// Mise à  jour des infos après une résolution LCP
function info = update(info, z, nan_cnt, iter)
    info.cnt = info.cnt + 1;
    info.nan_cnt = info.nan_cnt + nan_cnt;
    info.lmk_iter(1) = min( info.lmk_iter(1), iter);
    info.lmk_iter(2) = info.lmk_iter(2) + iter;
    info.lmk_iter(3) = max( info.lmk_iter(3), iter);
    if any(isnan(z))
        info.error_code = 3;
        info.error_msg = 'NaN supermassive black hole !';
    end
end
        

// On relance le code lemke tant qu'il ne renvoie pas une solution correcte
function [z, nan_cnt, iter] = lemke_nofail(M, q)
    zer_tol = 1e-10;
    nan_cnt = 0;
    while zer_tol < 1
        [z,~,iter] = lemke(M, q, [], zer_tol);
        if ~any(isnan(z)); break; end;
        nan_cnt = nan_cnt + 1;
        zer_tol = 10*zer_tol;
    end;
end


}}} // namespace floe::lcp::solver

#endif // FLOE_LCP_SOLVER_ITERLEMKE_HPP