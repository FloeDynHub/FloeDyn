;; FloeDyn build dependencies. Exact versions come from the pinned channel (hpc/channels.scm),
;; resolved via `guix time-machine` in hpc/guix_init.bash. Versions are pinned to the known-good set:
;;   gcc-toolchain@14 (must match the libstdc++ ABI of the channel's libs, e.g. hdf5_cpp.so, which are
;;     built with the channel's recent gcc; a lower gcc fails to link on GLIBCXX_3.4.32 / CXXABI_1.3.15),
;;   cgal@5 (default 6.x breaks against Boost 1.72), openmpi@4 (openmpi 5 dropped the C++ MPI bindings).
;; boost is deliberately NOT listed: FloeDyn builds against a hand-built Boost 1.72 from $HOME (BOOST_ROOT
;; in guix_init.bash). Listing it and then `guix remove`-ing it made guix_init.bash churn 2 profile
;; generations per activation (manifest re-adds it, remove drops it); leaving it out keeps it idempotent.
(specifications->manifest
  '("gcc-toolchain@14"
    "openmpi@4"
    "python@3"
    "cgal@5"
    "gmp"
    "mpfr"
    "eigen"
    "hdf5"
    "cereal"
    "matio"
    ))
