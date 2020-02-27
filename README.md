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
| 0.2.3    | add profile LL MLE, add functions to extract coefs and covariance parameters | 10th of July 2019 |
| 0.2.4    | add SVC_mle.control for SVC_mle object, update vignette | 10th of July 2019 |
| 0.2.5    | new handling of multiple observations | 10th of July 2019 |
| 0.2.6    | predictive variance | 5th of September 2019 |
| 0.2.7    | add summary, printing, plotting functions, and other methods for `SVC_mle` objects | 16th of September 2019 |
| 0.2.8    | Revisions for CRAN submission | 8th of October 2019 |
| 0.2.9    | Final CRAN version | 10th of October 2019 |
| 0.2.10   | parallelization using `optim.parallel`, citation info, orcID | 23rd of February 2020 |
| 0.2.11   | URL, Bugreports in DESCRIPTION | TBD |
