# varycoef

This R package contains functions to work with Spatially Varying Coefficient Models. To install it, run 


`devtools::install_github("jakobdambon/varycoef")`


Be sure that the `devtools` package is already installed. Check out the vignette to find some examples:


`vignette("example", package = "varycoef")`


## Next steps 

| Version  | Functionality | Prognosed Roll-out |
|----------|---------------|------|
| 0.2.0    | Seperate fixed and random effects | 12th of April 2019 |
| 0.2.1    | enable tapering | 23rd of April 2019 |
| 0.2.2    | BUGFIX: Mulitple Observations at locations. New residuals and fitted methods for SVC_mle, SVC_mle with formula call, warning on call without control | 5th of June 2019 |
| 0.3.0    | add summary, printing and plotting functions for `SVC_mle` objects | TBD |
| 0.4.0    | parallelization using `optim.parallel` | TBD |
| 0.5.0    | SVC INLA | TBD |
| 0.6.0    | preparing for CRAN release | TBD |
