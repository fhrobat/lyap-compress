# Results

Folder: `results/popeye_2026-11-10`. Machine and MATLAB version: see `machine.txt`.

matvecs of two-pass = matrix-vector products of both passes; residual = scaled residual norm computed a posteriori; ratio = time compress / time two-pass.

## 2D Laplacian

| N | k | matvecs compress | matvecs two-pass | matvecs restart | time compress | time two-pass | time restart | ratio | residual compress | residual two-pass |
|---|---|---|---|---|---|---|---|---|---|---|
| 179776 | 35 | 658 | 1316 | 7867 | 1.9 | 2.4 | 277.7 | 0.80 | 5.3e-07 | 5.3e-07 |
| 360000 | 38 | 936 | 1872 | - | 5.7 | 7.1 | - | 0.81 | 5.3e-07 | 5.3e-07 |
| 719104 | 41 | 1340 | 2680 | - | 19.6 | 23.7 | - | 0.83 | 4.1e-07 | 4.1e-07 |
| 1440000 | 44 | 1886 | 3772 | - | 62.9 | 72.4 | - | 0.87 | 5.9e-07 | 5.9e-07 |

<details><summary>All methods</summary>

| N | method | time | matvecs | residual | poles |
|---|---|---|---|---|---|
| 179776 | compress | 1.9 | 658 | 5.3e-07 | 35 |
| 179776 | two-pass | 2.4 | 1316 | 5.3e-07 | 35 |
| 179776 | restart | 277.7 | 7867 | 9.9e-07 | - |
| 179776 | cycle_direct | 1.9 | 0 | 1.6e-07 | 20 |
| 179776 | cycle_plain | 10.0 | 6487 | 1.6e-07 | 20 |
| 179776 | cycle_ichol | 7.1 | 504 | 1.6e-07 | 20 |
| 179776 | adaptive_direct | 3.5 | 0 | 7.7e-07 | 15 |
| 179776 | adaptive_plain | 7.7 | 4994 | 7.7e-07 | 15 |
| 179776 | adaptive_ichol | 5.6 | 377 | 7.7e-07 | 15 |
| 179776 | multishiftCG | 9.9 | 729 | 1.3e-08 | 35 |
| 179776 | multishiftCG_rem | 4.8 | 729 | 1.5e-08 | 35 |
| 179776 | extended-icholPCG | 29.2 | 2035 | 7.6e-07 | 48 |
| 179776 | extended-backslash | 3.1 | 23 | 5.8e-07 | 46 |
| 360000 | compress | 5.7 | 936 | 5.3e-07 | 38 |
| 360000 | two-pass | 7.1 | 1872 | 5.3e-07 | 38 |
| 360000 | cycle_direct | 3.9 | 0 | 5.2e-07 | 20 |
| 360000 | cycle_plain | 30.9 | 8710 | 5.2e-07 | 20 |
| 360000 | cycle_ichol | 18.3 | 649 | 5.2e-07 | 20 |
| 360000 | adaptive_direct | 8.0 | 0 | 7.7e-07 | 16 |
| 360000 | adaptive_plain | 26.9 | 7815 | 7.7e-07 | 16 |
| 360000 | adaptive_ichol | 16.4 | 564 | 7.7e-07 | 16 |
| 360000 | multishiftCG | 31.4 | 1036 | 1.2e-08 | 38 |
| 360000 | multishiftCG_rem | 14.4 | 1036 | 1.4e-08 | 38 |
| 360000 | extended-icholPCG | 80.5 | 3126 | 9.7e-07 | 54 |
| 360000 | extended-backslash | 8.0 | 27 | 6.8e-07 | 54 |
| 719104 | compress | 19.6 | 1340 | 4.1e-07 | 41 |
| 719104 | two-pass | 23.7 | 2680 | 4.1e-07 | 41 |
| 719104 | cycle_direct | 9.9 | 0 | 8.1e-08 | 24 |
| 719104 | cycle_plain | 114.9 | 14069 | 8.1e-08 | 24 |
| 719104 | cycle_ichol | 58.1 | 1007 | 8.1e-08 | 24 |
| 719104 | adaptive_direct | 20.6 | 0 | 9.9e-07 | 17 |
| 719104 | adaptive_plain | 86.9 | 10533 | 9.9e-07 | 17 |
| 719104 | adaptive_ichol | 44.4 | 738 | 9.9e-07 | 17 |
| 719104 | multishiftCG | 96.1 | 1469 | 1.3e-08 | 41 |
| 719104 | multishiftCG_rem | 43.0 | 1469 | 1.5e-08 | 41 |
| 719104 | extended-icholPCG | 264.7 | 5067 | 9.0e-07 | 64 |
| 719104 | extended-backslash | 20.8 | 31 | 7.2e-07 | 62 |
| 1440000 | compress | 62.9 | 1886 | 5.9e-07 | 44 |
| 1440000 | two-pass | 72.4 | 3772 | 5.9e-07 | 44 |
| 1440000 | cycle_direct | 21.5 | 0 | 2.5e-07 | 24 |
| 1440000 | cycle_plain | 340.5 | 19039 | 2.5e-07 | 24 |
| 1440000 | cycle_ichol | 155.7 | 1328 | 2.5e-07 | 24 |
| 1440000 | adaptive_direct | 48.8 | 0 | 8.5e-07 | 18 |
| 1440000 | adaptive_plain | 263.1 | 14897 | 8.5e-07 | 18 |
| 1440000 | adaptive_ichol | 121.3 | 1029 | 8.5e-07 | 18 |
| 1440000 | multishiftCG | 301.6 | 2087 | 1.3e-08 | 44 |
| 1440000 | multishiftCG_rem | 128.4 | 2087 | 1.5e-08 | 44 |
| 1440000 | extended-icholPCG | 878.2 | 8282 | 7.1e-07 | 76 |
| 1440000 | extended-backslash | 52.5 | 36 | 9.1e-07 | 72 |

</details>

## 3D Laplacian

| N | k | matvecs compress | matvecs two-pass | matvecs restart | time compress | time two-pass | time restart | ratio | residual compress | residual two-pass |
|---|---|---|---|---|---|---|---|---|---|---|
| 175616 | 33 | 172 | 344 | 218 | 0.5 | 0.7 | 14.1 | 0.66 | 1.6e-13 | 1.2e-13 |
| 357911 | 35 | 217 | 434 | - | 1.2 | 2.0 | - | 0.61 | 2.9e-13 | 2.9e-13 |
| 704969 | 37 | 299 | 598 | - | 3.7 | 5.7 | - | 0.65 | 8.2e-13 | 2.4e-13 |
| 1404928 | 39 | 365 | 730 | - | 10.3 | 14.8 | - | 0.70 | 1.2e-12 | 9.8e-13 |

<details><summary>All methods</summary>

| N | method | time | matvecs | residual | poles |
|---|---|---|---|---|---|
| 175616 | compress | 0.5 | 172 | 1.6e-13 | 33 |
| 175616 | two-pass | 0.7 | 344 | 1.2e-13 | 33 |
| 175616 | restart | 14.1 | 218 | 1.4e-12 | - |
| 175616 | cycle_plain | 5.2 | 2198 | 1.9e-13 | 24 |
| 175616 | cycle_ichol | 8.4 | 401 | 1.7e-13 | 24 |
| 175616 | adaptive_plain | 4.1 | 1773 | 3.7e-13 | 20 |
| 175616 | adaptive_ichol | 7.6 | 321 | 3.5e-13 | 20 |
| 175616 | extended-plainCG | 12.5 | 4206 | 1.0e-12 | 40 |
| 175616 | extended-icholPCG | 18.9 | 882 | 1.0e-12 | 42 |
| 357911 | compress | 1.2 | 217 | 2.9e-13 | 35 |
| 357911 | two-pass | 2.0 | 434 | 2.9e-13 | 35 |
| 357911 | cycle_plain | 13.9 | 2652 | 4.6e-13 | 24 |
| 357911 | cycle_ichol | 20.9 | 456 | 3.8e-13 | 24 |
| 357911 | adaptive_plain | 12.2 | 2224 | 7.7e-13 | 21 |
| 357911 | adaptive_ichol | 20.2 | 391 | 7.2e-13 | 21 |
| 357911 | extended-plainCG | 27.4 | 6201 | 1.1e-12 | 46 |
| 357911 | extended-icholPCG | 52.3 | 1240 | 9.4e-13 | 48 |
| 704969 | compress | 3.7 | 299 | 8.2e-13 | 37 |
| 704969 | two-pass | 5.7 | 598 | 2.4e-13 | 37 |
| 704969 | cycle_plain | 41.1 | 3687 | 5.1e-13 | 28 |
| 704969 | cycle_ichol | 56.2 | 618 | 1.7e-13 | 28 |
| 704969 | adaptive_plain | 31.5 | 2835 | 7.4e-13 | 23 |
| 704969 | adaptive_ichol | 50.1 | 483 | 3.4e-13 | 23 |
| 704969 | extended-plainCG | 92.2 | 8998 | 1.3e-12 | 54 |
| 704969 | extended-icholPCG | 142.6 | 1711 | 1.2e-12 | 54 |
| 1404928 | compress | 10.3 | 365 | 1.2e-12 | 39 |
| 1404928 | two-pass | 14.8 | 730 | 9.8e-13 | 39 |
| 1404928 | cycle_plain | 101.5 | 4446 | 1.3e-12 | 28 |
| 1404928 | cycle_ichol | 129.0 | 725 | 4.0e-13 | 28 |
| 1404928 | adaptive_plain | 84.0 | 3697 | 1.5e-12 | 24 |
| 1404928 | adaptive_ichol | 128.6 | 606 | 1.0e-12 | 24 |

</details>

## FEniCS Rail

| N | k | matvecs compress | matvecs two-pass | matvecs restart | time compress | time two-pass | time restart | ratio | residual compress | residual two-pass |
|---|---|---|---|---|---|---|---|---|---|---|
| 79841 | 30 | 532 | 1064 | - | 1.2 | 1.7 | - | 0.71 | 6.0e-05 | 6.0e-05 |
| 317377 | 33 | 596 | 1192 | - | 5.9 | 8.3 | - | 0.71 | 6.5e-05 | 6.5e-05 |
| 1265537 | 31 | 632 | 1264 | - | 27.3 | 45.6 | - | 0.60 | 6.4e-05 | 6.4e-05 |

<details><summary>All methods</summary>

| N | method | time | matvecs | residual | poles |
|---|---|---|---|---|---|
| 79841 | compress | 1.2 | 532 | 6.0e-05 | 30 |
| 79841 | two-pass | 1.7 | 1064 | 6.0e-05 | 30 |
| 79841 | cycle_plain_4_poles | 4.0 | 5117 | 3.0e-05 | 16 |
| 79841 | cycle_ichol_4_poles | 2.0 | 243 | 3.0e-05 | 16 |
| 79841 | cycle_plain_8_poles | 4.3 | 5478 | 8.6e-06 | 16 |
| 79841 | cycle_ichol_8_poles | 2.1 | 258 | 8.6e-06 | 16 |
| 79841 | cycle_plain_12_poles | 7.1 | 9280 | 3.8e-07 | 24 |
| 79841 | cycle_ichol_12_poles | 3.6 | 432 | 3.8e-07 | 24 |
| 317377 | compress | 5.9 | 596 | 6.5e-05 | 33 |
| 317377 | two-pass | 8.3 | 1192 | 6.5e-05 | 33 |
| 317377 | cycle_plain_4_poles | 26.3 | 6496 | 7.7e-05 | 16 |
| 317377 | cycle_ichol_4_poles | 9.9 | 305 | 7.7e-05 | 16 |
| 317377 | cycle_plain_8_poles | 28.0 | 6941 | 2.1e-05 | 16 |
| 317377 | cycle_ichol_8_poles | 10.8 | 323 | 2.1e-05 | 16 |
| 317377 | cycle_plain_12_poles | 48.2 | 12206 | 1.7e-06 | 24 |
| 317377 | cycle_ichol_12_poles | 18.5 | 559 | 1.7e-06 | 24 |
| 1265537 | compress | 27.3 | 632 | 6.4e-05 | 31 |
| 1265537 | two-pass | 45.6 | 1264 | 6.4e-05 | 31 |
| 1265537 | cycle_plain_4_poles | 235.2 | 10871 | 1.7e-05 | 20 |
| 1265537 | cycle_ichol_4_poles | 67.1 | 494 | 1.7e-05 | 20 |
| 1265537 | cycle_plain_8_poles | 193.8 | 9027 | 2.5e-05 | 16 |
| 1265537 | cycle_ichol_8_poles | 57.4 | 413 | 2.5e-05 | 16 |
| 1265537 | cycle_ichol_12_poles | 89.6 | 628 | 7.0e-06 | 24 |

</details>

## Heat equation

| N | k | matvecs compress | matvecs two-pass | matvecs restart | time compress | time two-pass | time restart | ratio | residual compress | residual two-pass |
|---|---|---|---|---|---|---|---|---|---|---|
| 146098 | 35 | 413 | 826 | - | 2.5 | 4.1 | - | 0.61 | 2.2e-07 | 2.2e-07 |
| 330043 | 37 | 569 | 1138 | - | 7.5 | 13.1 | - | 0.57 | 2.1e-07 | 2.1e-07 |
| 1059496 | 37 | 794 | 1588 | - | 38.7 | 68.1 | - | 0.57 | 2.9e-07 | 2.9e-07 |

<details><summary>All methods</summary>

| N | method | time | matvecs | residual | poles |
|---|---|---|---|---|---|
| 146098 | compress | 2.5 | 413 | 2.2e-07 | 35 |
| 146098 | two-pass | 4.1 | 826 | 2.2e-07 | 35 |
| 146098 | cycle_plain_4_poles | 17.3 | 5140 | 3.5e-07 | 24 |
| 146098 | cycle_ichol_4_poles | 11.3 | 503 | 3.5e-07 | 24 |
| 146098 | cycle_plain_8_poles | 17.5 | 5201 | 3.5e-08 | 24 |
| 146098 | cycle_ichol_8_poles | 12.1 | 514 | 3.5e-08 | 24 |
| 146098 | cycle_plain_12_poles | 17.4 | 5179 | 2.6e-08 | 24 |
| 146098 | cycle_ichol_12_poles | 12.7 | 513 | 2.6e-08 | 24 |
| 330043 | compress | 7.5 | 569 | 2.1e-07 | 37 |
| 330043 | two-pass | 13.1 | 1138 | 2.1e-07 | 37 |
| 330043 | cycle_plain_4_poles | 51.9 | 6039 | 7.4e-07 | 24 |
| 330043 | cycle_ichol_4_poles | 32.9 | 600 | 7.4e-07 | 24 |
| 330043 | cycle_plain_8_poles | 53.0 | 6088 | 7.1e-08 | 24 |
| 330043 | cycle_ichol_8_poles | 35.3 | 609 | 7.1e-08 | 24 |
| 330043 | cycle_plain_12_poles | 53.7 | 6114 | 4.6e-08 | 24 |
| 330043 | cycle_ichol_12_poles | 37.0 | 612 | 4.6e-08 | 24 |
| 1059496 | compress | 38.7 | 794 | 2.9e-07 | 37 |
| 1059496 | two-pass | 68.1 | 1588 | 2.9e-07 | 37 |
| 1059496 | cycle_plain_4_poles | 236.2 | 7171 | 8.6e-07 | 24 |
| 1059496 | cycle_ichol_4_poles | 139.6 | 769 | 8.6e-07 | 24 |
| 1059496 | cycle_plain_8_poles | 239.9 | 7199 | 9.6e-08 | 24 |
| 1059496 | cycle_ichol_8_poles | 149.1 | 778 | 9.6e-08 | 24 |
| 1059496 | cycle_plain_12_poles | 239.9 | 7258 | 6.6e-08 | 24 |
| 1059496 | cycle_ichol_12_poles | 155.5 | 789 | 6.6e-08 | 24 |

</details>

## Sensitivity to the estimate of the smallest eigenvalue

| i | variant | num_poles | matvecs | true_res |
|---|---|---|---|---|
| 1 | recomputed_k | 35 | 658 | 5.33e-07 |
| 1 | fixed_k | 35 | 658 | 5.33e-07 |
| 2 | recomputed_k | 32 | 669 | 3.12e-07 |
| 2 | fixed_k | 35 | 658 | 5.33e-07 |
| 3 | recomputed_k | 31 | 689 | 1.09e-07 |
| 3 | fixed_k | 35 | 658 | 5.33e-07 |
| 4 | recomputed_k | 30 | 709 | 4.16e-08 |
| 4 | fixed_k | 35 | 658 | 5.33e-07 |
| 5 | recomputed_k | 29 | 668 | 5.12e-07 |
| 5 | fixed_k | 35 | 658 | 5.33e-07 |
| 6 | recomputed_k | 28 | 686 | 2.72e-06 |
| 6 | fixed_k | 35 | 658 | 5.33e-07 |
| 7 | recomputed_k | 28 | 686 | 8.39e-06 |
| 7 | fixed_k | 35 | 658 | 5.34e-07 |
| 8 | recomputed_k | 27 | 704 | 6.3e-05 |
| 8 | fixed_k | 35 | 658 | 5.41e-07 |
| 9 | recomputed_k | 27 | 704 | 0.000145 |
| 9 | fixed_k | 35 | 658 | 5.9e-07 |
| 10 | recomputed_k | 26 | 655 | 0.000358 |
| 10 | fixed_k | 35 | 658 | 7.86e-07 |

## MATLAB runs

Exit status 0 = completed; 137 = killed (out of memory, limit in `machine.txt`); other = error, see `logs/`.

| script | size | exit status | duration (s) |
|---|---|---|---|
| test_2DLap | 1 | 0 | 379 |
| test_2DLap | 2 | 0 | 269 |
| test_2DLap | 3 | 0 | 829 |
| test_2DLap | 4 | 0 | 2499 |
| test_2DLap_sensitivity | all | 0 | 49 |
| test_3DLap | 1 | 0 | 84 |
| test_3DLap | 2 | 0 | 165 |
| test_3DLap | 3 | 0 | 446 |
| test_3DLap | 4 | 0 | 496 |
| test_fenics_rail | 1 | 0 | 38 |
| test_fenics_rail | 2 | 0 | 172 |
| test_fenics_rail | 3 | 0 | 767 |
| test_heat | 1 | 0 | 125 |
| test_heat | 2 | 0 | 340 |
| test_heat | 3 | 0 | 1451 |

