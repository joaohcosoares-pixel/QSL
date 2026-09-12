# W-REVERSED numerical benchmark

One warm-up and at least five timed repeats per operation. Units: seconds.

| Sector | Operation | Ncell | Median | Min | Max |
|---|---|---:|---:|---:|---:|
| OLD | 1000 H(q) evaluations | — | 0.1179842 | 0.1145227 | 0.1239812 |
| OLD | 1000 eigvalsh 4x4 (H prebuilt) | — | 0.0100431 | 0.0095941 | 0.0103802 |
| OLD | physical bands 2001 (H build + eigvalsh) | — | 0.0073397 | 0.0073189 | 0.0086091 |
| OLD | physical bands 16001 (H build + eigvalsh) | — | 0.0603052 | 0.0600554 | 0.0613428 |
| OLD | integral 16001 (physical bands precomputed) | — | 0.0006715 | 0.0006532 | 0.0006850 |
| OLD | energy 16001 (H build + eigvalsh + integral) | — | 0.0614295 | 0.0607116 | 0.0626344 |
| OLD | finite A build + full eigvalsh(iA) | 5 | 0.0002802 | 0.0002743 | 0.0002994 |
| OLD | finite A build + full eigvalsh(iA) | 8 | 0.0004529 | 0.0004471 | 0.0005009 |
| OLD | finite A build + full eigvalsh(iA) | 17 | 0.0012045 | 0.0012000 | 0.0012316 |
| OLD | finite A build + full eigvalsh(iA) | 33 | 0.0034775 | 0.0034556 | 0.0036855 |
| OLD | finite A build + full eigvalsh(iA) | 65 | 0.0146036 | 0.0145863 | 0.0147135 |
| OLD | finite A build + full eigvalsh(iA) | 129 | 0.0868326 | 0.0819791 | 0.0879501 |
| W-REVERSED | 1000 H(q) evaluations | — | 0.1158681 | 0.1138051 | 0.1163873 |
| W-REVERSED | 1000 eigvalsh 4x4 (H prebuilt) | — | 0.0093319 | 0.0092718 | 0.0094667 |
| W-REVERSED | physical bands 2001 (H build + eigvalsh) | — | 0.0069356 | 0.0069133 | 0.0071629 |
| W-REVERSED | physical bands 16001 (H build + eigvalsh) | — | 0.0514835 | 0.0494139 | 0.0554720 |
| W-REVERSED | integral 16001 (physical bands precomputed) | — | 0.0004988 | 0.0004906 | 0.0005138 |
| W-REVERSED | energy 16001 (H build + eigvalsh + integral) | — | 0.0464317 | 0.0435691 | 0.0474290 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 5 | 0.0001826 | 0.0001787 | 0.0002104 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 8 | 0.0002930 | 0.0002895 | 0.0002994 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 17 | 0.0007887 | 0.0007862 | 0.0007962 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 33 | 0.0023817 | 0.0022508 | 0.0026958 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 65 | 0.0096551 | 0.0093807 | 0.0103786 |
| W-REVERSED | finite A build + full eigvalsh(iA) | 129 | 0.0521610 | 0.0490259 | 0.0550456 |

Finite-chain timings include building A and diagonalizing iA.

Thread environment: `{"MKL_NUM_THREADS": null, "OMP_NUM_THREADS": null, "OPENBLAS_NUM_THREADS": "1"}`.

Script SHA-256: `{"build_report_w_reversed.py": "e4ef45f1d3b91fd84725b2c125b8fc42bc80c284b9b41db158ad729924f0a2ac", "compare_old_vs_w_reversed.py": "2dbbffa1e9a1f04e332dc1fabd4aa37c5a331971a7c76c3c7db9626497956929", "draw_w_lattice.py": "d9715c4987f5d84ceeaec0610c80a0e5ff610f6104018863b99b2d3354b8d735", "quadrupolar_spin_liquid_w_reversed.py": "6bb17be928631f1a7837dc1c0124a8d37c213b012212b052688096f342993320", "symbolic_w_audit.py": "8c085c697ae8b026eb1d5ec4a0274e3a0b78d97ab4d8fd33580875534dd4ff71"}`.
