# T-REVERSED numerical benchmark

One warm-up and at least five timed repeats per operation. Units: seconds.

| Sector | Operation | Ncell | Median | Min | Max |
|---|---|---:|---:|---:|---:|
| OLD | 1000 H(q) evaluations | — | 0.0577331 | 0.0573916 | 0.0578806 |
| OLD | 1000 eigvalsh 4x4 (H prebuilt) | — | 0.0049187 | 0.0048699 | 0.0051318 |
| OLD | physical bands 2001 (H build + eigvalsh) | — | 0.0036639 | 0.0036311 | 0.0044772 |
| OLD | physical bands 16001 (H build + eigvalsh) | — | 0.0317081 | 0.0312469 | 0.0346434 |
| OLD | F(q) 16001 (H build + eigvalsh + negative sum) | — | 0.0325631 | 0.0316078 | 0.0327668 |
| OLD | integral 16001 (physical bands precomputed) | — | 0.0003367 | 0.0003277 | 0.0003561 |
| OLD | energy 16001 (H build + eigvalsh + integral) | — | 0.0327840 | 0.0323858 | 0.0333004 |
| OLD | finite A build + full eigvalsh(iA) | 5 | 0.0001513 | 0.0001379 | 0.0001639 |
| OLD | finite A build + full eigvalsh(iA) | 8 | 0.0002231 | 0.0002180 | 0.0002386 |
| OLD | finite A build + full eigvalsh(iA) | 17 | 0.0005969 | 0.0005845 | 0.0006117 |
| OLD | finite A build + full eigvalsh(iA) | 33 | 0.0017322 | 0.0017216 | 0.0017525 |
| OLD | finite A build + full eigvalsh(iA) | 65 | 0.0074139 | 0.0072350 | 0.0078351 |
| OLD | finite A build + full eigvalsh(iA) | 129 | 0.0457044 | 0.0448476 | 0.0462327 |
| T-REVERSED | 1000 H(q) evaluations | — | 0.0587762 | 0.0576297 | 0.0617663 |
| T-REVERSED | 1000 eigvalsh 4x4 (H prebuilt) | — | 0.0047704 | 0.0043694 | 0.0053609 |
| T-REVERSED | physical bands 2001 (H build + eigvalsh) | — | 0.0028568 | 0.0027802 | 0.0037831 |
| T-REVERSED | physical bands 16001 (H build + eigvalsh) | — | 0.0253836 | 0.0240054 | 0.0274032 |
| T-REVERSED | F(q) 16001 (H build + eigvalsh + negative sum) | — | 0.0245401 | 0.0241250 | 0.0259706 |
| T-REVERSED | integral 16001 (physical bands precomputed) | — | 0.0003365 | 0.0003287 | 0.0003688 |
| T-REVERSED | energy 16001 (H build + eigvalsh + integral) | — | 0.0250399 | 0.0245031 | 0.0254261 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 5 | 0.0001396 | 0.0001351 | 0.0001635 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 8 | 0.0002209 | 0.0002174 | 0.0002433 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 17 | 0.0005988 | 0.0005894 | 0.0006061 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 33 | 0.0017026 | 0.0016917 | 0.0017838 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 65 | 0.0073205 | 0.0072494 | 0.0074570 |
| T-REVERSED | finite A build + full eigvalsh(iA) | 129 | 0.0431558 | 0.0429746 | 0.0444735 |
| OLD + T-REVERSED | adversarial suite 1000 cases (both sectors, 500 mixed-sign cases) | — | 2.4161504 | 2.3493794 | 2.4804934 |

Finite-chain timings include building A and diagonalizing iA. The dedicated F(q) timing includes matrix construction and eigvalsh. The 1000-case randomized suite is repeated in full.

Thread environment: `{"MKL_NUM_THREADS": null, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}`.

Script SHA-256: `{"audit_explicit_majoranas.py": "49b0d8124dfba8d0fcc06e053503fc1f7f67b65b2c4cb95775ab3dd3014676bd", "build_report_t_reversed.py": "6bc37b4f070f52e9b54b6a23765c83dc1e4447c4a809ba8e0b5b6697810a96f0", "compare_old_vs_t_reversed.py": "495f507b4ed3cd123da99cc294d3ba452be7343adf2e008f2005ee1282bdbb5b", "draw_t_lattice.py": "53178a7b5cea5e4b054d2029f780a88b3f12b7b11784290aefdb55f4d364260b", "quadrupolar_spin_liquid_t_reversed.py": "a7b2ea8f54da3d352e47a9a8e996b598022481b5f1d752d11703db352e492351", "run_t_reversed.py": "92c993ff8db3566d09ab15146846c06fbe77c825b5409925636e0d8c27d3165f", "symbolic_t_audit.py": "d90f968dc8f4cce8dca4dd142da1a38b3463657aa2785a8179e8b552e3fd68e6"}`.
