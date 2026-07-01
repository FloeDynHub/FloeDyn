# FloeDyn on CIMENT/Gricad clusters (kraken, dahu, ...)

Reproducible build + run via a pinned Guix profile. **Two scripts, two distinct jobs:**

| Script | When | What it does |
|--------|------|--------------|
| `guix_init.bash` | **ONCE** — first setup, or after changing dependencies | Builds the `floe` Guix profile through the pinned `channels.scm` + `manifest_floedyn.scm` (`guix time-machine … package`, `guix remove boost`), sets the build env, runs `waf configure`. **Mutates the profile → takes an exclusive lock.** Run it on a login node or on `dahu-workflow1` (long processes allowed). Then compile: `python3 ./waf build --target FLOE` (or FLOE_PBC / FLOE_MPI / FLOE_MPI_PBC). |
| `guix_env.bash` | **In every job** (OAR script) | Activates the already-built profile **read-only** (sets `LD_LIBRARY_PATH` for Boost 1.72, `PATH` for the profile's `mpirun`). **No profile mutation → no lock**, so many jobs can source it concurrently. |

## Why two scripts (do not merge them)
`guix package` / `guix remove` take an **exclusive lock** on the profile — even when nothing changes.
If every job runs `guix_init.bash`, launching a batch of jobs makes them all fight for that lock
(`error: profile … is locked by another process`), the profile setup fails, and the run crashes
(e.g. `H5::FileIException`, because HDF5 from the profile isn't available). So the per-job path must
never mutate the profile: build once with `guix_init.bash`, then run with `guix_env.bash`.

## Pinning
- `channels.scm` — exact Guix commit (determines package versions).
- `manifest_floedyn.scm` — package list, pinned to `gcc-toolchain@14` / `cgal@5` / `openmpi@4`
  (see the comments in that file for why).
- Boost 1.72 is hand-built outside Guix (in `$HOME`), used via `BOOST_ROOT` — the one non-reproducible
  piece for now.
