# `varycoef`: An R package to Model Spatially Varying Coefficients using Gaussian Processes

## About

The R package `varycoef` is the software implementation of **Gaussian process-based spatially varying coefficient models** by [Dambon et al. (2021a)](https://www.sciencedirect.com/science/article/pii/S2211675320300646). It extends linear regression models such that the coefficients are depending on some coordinates in a `d` dimensional space, i.e., the coefficient `b_j` for a covariate `j` is depending on coordinates `s` and therefore of the form `b_j(s)`. These coefficients are modeled using Gaussian processes. In most applications, the coordinates `s` tend to be observation locations like longitude and latitude (see [Dambon et al. (2022)](https://sjes.springeropen.com/articles/10.1186/s41937-021-00080-2) as an example). However, the concept can be extended in the number of dimensions or, say, using observation time points to model time-varying coefficients.

The method relies on maximum likelihood estimation. It has been optimized to work with large data sets by applying covariance tapering by [Furrer et al. (2006)](https://www.jstor.org/stable/27594195) if necessary and allows for a moderate number of spatially varying coefficients. The R package contains methods to estimate Gaussian process-based (spatially) varying coefficient models, (spatially) predict coefficients as well as the response, and variable selection methods. Latter are based on [Dambon et al. (2021b)](https://arxiv.org/abs/2101.01932).

## Getting Started

To install it, run

```
devtools::install_github("jakobdambon/varycoef")
```

for the latest version on this repository or download it from [CRAN](https://cran.r-project.org/web/packages/varycoef/index.html).

## Model Assumptions

**Note: The exact definition of the model is given in [Dambon et al. (2021a)](https://www.sciencedirect.com/science/article/pii/S2211675320300646).**

### Linear Model

Let `y` be the response vector, let `X` be the covariate matrix, and let `xi` be the error term sampled from a zero-mean normal distribution with variance `s2`. Then the linear model is given by

```
y = Xb + xi
```

with coefficient vector `b`. The coefficients are also called *fixed effects*.

### Spatially Varying Coefficients

We now allow the coefficients to vary from their respective mean. Let `e(s)` contain the location-dependent differences. In a first step, with some slight abuse of notation, the linear model from above is extended by:

```
y = X(b + e(s)) + xi
```

Note that not all coefficients necessarily must be varying. Therefore, we introduce a second covariate matrix `W` to specify which coefficients are spatially varying. In the case where all coefficients are varying, we have `W = X`. Again with some abuse of notation, we have the SVC model:

```
y = Xb + We(s) + xi
```

### Gaussian Process-based Coefficients

We assume that the spatially varying coefficients are defined by Gaussian processes. That is, the deviations per covariate `e_k(s)` are defined as a zero-mean Gaussian process defined by some covariance function `c(d, par)` that models the spatial dependence between observations using the pairwise distances. Each coefficient is parameterized by a tuple of parameters `par` that consists of a range and variance. One main assumption of the model is that the individual coefficients per covariate, i.e., the individual Gaussian processes, are mutually independent and independent of the error.

### Connection to Mixed Effect Models

For a finite number of observations `n`, the model can be expressed as a so-called mixed effect model. That is, using the covariance functions of the Gaussian processes and the distance matrix of the observations, we can express each `e_k(s)` as a multivariate, zero-mean normal distribution. Another common name for these effects are random effects. Together with the fixed effects from an ordinary linear model, we receive a  so-called mixed effect model. The assumption of mutual independence between the Gaussian allows an easy construction of the joint covariance matrix `S_y` of the response `y`. 

## Examples

We will add links to vignettes and further examples here very soon.

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
