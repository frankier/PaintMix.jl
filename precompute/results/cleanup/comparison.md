# Surrogate and release table comparison

All candidates use the same measured inputs and validation corpus. No candidate is promoted.

| Variant | Fit seconds | Table seconds | Epush dense | Max gamut violation | Mixing Oklab p99 | Encoded channel p99 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| direct-lbfgsb | 239.84870505332947 | 736.4188771247864 | 3.047834748704887e-6 | 0.0163671968400976 | 0.0356442315992515 | 0.09819350180911951 |
| softplus-lbfgs | 891.7305769920349 | 712.8255181312561 | 7.132231810803058e-8 | 0.0019439165248735049 | 0.06036467227913475 | 0.19701668162747205 |
| softplus-lbfgsb | 137.76873302459717 | 714.387622833252 | 9.528900610884859e-8 | 0.0030667874613555138 | 0.06147403623196634 | 0.20983346397867889 |

The softplus/L-BFGS-B row isolates the solver change; the direct/L-BFGS-B row then isolates coefficient coordinates. These are local, budget-limited fits. Layer errors are reported separately and are not additive.
