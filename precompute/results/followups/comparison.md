# Surrogate and release table comparison

All candidates use the same measured inputs and validation corpus. No candidate is promoted.

| Variant | Fit seconds | Table seconds | Epush dense | Dense max violation | Mixing Oklab p99 | Encoded channel p99 | Held-out max violation | Physical target reached |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| direct-lbfgsb | 218.65284419059753 | 573.8004369735718 | 3.51538205509085e-7 | 0.006013251202153513 | 0.03476567606022258 | 0.09121660558148359 | 0.006013248632715632 | false |
| gauge-lbfgsb | 107.96183681488037 | 554.6116240024567 | 2.9801466200771624e-9 | 0.002178178652680059 | 0.06170786570169339 | 0.1996545638816224 | 0.009105773063727503 | false |
| softplus-lbfgs | 838.4349110126495 | 711.8210561275482 | 1.2772871109866943e-9 | 0.0012253214083731364 | 0.06158925592175731 | 0.20105772291967738 | 0.012358246118034105 | false |
| softplus-lbfgsb | 132.36371302604675 | 696.0464770793915 | 6.971814233934321e-11 | 0.000346489020184405 | 0.06509229229645003 | 0.22953776384440544 | 0.004011866491919374 | false |

Matched runs stop at held-out gamut violation <= 1e-3, or at the common schedule/time budget. Gauge coordinates remove one common K/S scale per wavelength.


The softplus/L-BFGS-B row isolates the solver change; the direct/L-BFGS-B row then isolates coefficient coordinates. These are local, budget-limited fits. Layer errors are reported separately and are not additive.
