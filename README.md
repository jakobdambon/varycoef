# varycoef

## About

The R package `varycoef` is the software implementation of **Gaussian process-based spatially varying coefficient models** by [Dambon et al. (2021a)](https://www.sciencedirect.com/science/article/pii/S2211675320300646). It extends linear regression models such that the coefficients are depending on some coordinates in a $d$ dimensional space, i.e., the coefficient $b_j$ for a covariate $j$ is depending on a coordinate $s \in D \subset \mathbb{R}^d$ and therefore of the form $b_j(s)$. In most applications, these coordinates $s$ tend to be observation locations like longitude and latitude (see [Dambon et al. (2022)](https://sjes.springeropen.com/articles/10.1186/s41937-021-00080-2) as an example). However, the concept can be extend in the number of dimensions or, say, using observation time points to model time-varying coefficients.

The method relies on maximum likelihood estimation. It has been optimized to work with large data sets by applying covariance tapering [Furrer et al. (2020)](https://www.jstor.org/stable/27594195) if necessary and allows for a moderate number of spatially varying coefficients. The R package contains methods to estimate Gaussian process-based (spatially) varying coefficient models, (spatially) predict coefficients as well as the response, and variable selection methods. Latter are based on [Dambon et al. (2021b)](https://arxiv.org/abs/2101.01932).

## Getting Started

To install it, run

```
devtools::install_github("jakobdambon/varycoef")
```

for the latest version on this repository or download it from [CRAN](https://cran.r-project.org/web/packages/varycoef/index.html).


## Version History

| Version | Functionality                                                                                                                                        | Prognosed Roll-out     |
|----------------------|---------------------------------|-----------------|
| 0.3.3   | Removed some dependencies to `RandomFields`                                                                                                          | June 2022              |
| 0.3.2   | Parscale option in SVC_mle_control                                                                                                                   | 19th of July 2021      |
| 0.3.1   | New methods (summary output), pre JSS Submission                                                                                                     | 12th of May 2021       |
| 0.3.0   | Joint variable selection method available on CRAN                                                                                                    | 13th of January 2021   |
| 0.2.10  | parallelization using `optimParallel`, citation info, orcID                                                                                          | 23rd of February 2020  |
| 0.2.9   | Final CRAN version                                                                                                                                   | 10th of October 2019   |
| 0.2.8   | Revisions for CRAN submission                                                                                                                        | 8th of October 2019    |
| 0.2.7   | add summary, printing, plotting functions, and other methods for `SVC_mle` objects                                                                   | 16th of September 2019 |
| 0.2.6   | predictive variance                                                                                                                                  | 5th of September 2019  |
| 0.2.5   | new handling of multiple observations                                                                                                                | 10th of July 2019      |
| 0.2.4   | add SVC_mle.control for SVC_mle object, update vignette                                                                                              | 10th of July 2019      |
| 0.2.3   | add profile LL MLE, add functions to extract coefs and covariance parameters                                                                         | 10th of July 2019      |
| 0.2.2   | BUGFIX: Mulitple Observations at locations. New residuals and fitted methods for SVC_mle, SVC_mle with formula call, warning on call without control | 5th of June 2019       |
| 0.2.1   | enable tapering                                                                                                                                      | 23rd of April 2019     |
| 0.2.0   | Seperate fixed and random effects                                                                                                                    | 12th of April 2019     |
