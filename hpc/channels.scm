;; Pinned Guix state for reproducible FloeDyn builds (see hpc/guix_init.bash).
;; At this commit the manifest pins cgal@5 / gcc-toolchain@11 / openmpi@4, which match the known-good
;; toolchain: cgal 5.6.1 keeps the 5.x API (compatible with Boost 1.72; the default cgal 6.x breaks it),
;; and openmpi 4.x keeps the deprecated C++ MPI bindings FLOE_MPI relies on (openmpi 5 removed them).
;; Substitutes are available for this commit, so the build is fast (no source rebuild of an old Guix).
(list (channel
        (name 'guix)
        (url "https://git.guix.gnu.org/guix.git")
        (branch "master")
        (commit "9a62baae0e7de3ca774fd93c3038882aebc4eb84")
        (introduction
         (make-channel-introduction
          "9edb3f66fd807b096b48283debdcddccfea34bad"
          (openpgp-fingerprint
           "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
