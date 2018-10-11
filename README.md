[![Github release](https://img.shields.io/github/release/haziqj/phd-thesis.svg)](https://github.com/haziqj/phd-thesis/releases)
[![pdf download](https://img.shields.io/badge/pdf-download-brightgreen.svg)](https://github.com/haziqj/phd-thesis/releases/download/v1.1/phd-thesis-final-f378773.pdf)
[![errata download](https://img.shields.io/badge/errata-download-brightgreen.svg)](https://github.com/haziqj/phd-thesis/releases/download/v1.1/errata.pdf)

# Regression modelling using priors depending on Fisher information covariance kernels (I-priors)

Haziq Jamil, 11 October 2018.

A thesis submitted to the Department of Statistics of the London School of Economics and Political Science for the degree of Doctor of Philosophy.

## Citation

Please cite using the following BibTeX entry:

```latex
@phdthesis{jamil2018phdthesis,
  title={Regression modelling using priors depending on Fisher information covariance kernels (I-priors)},
  author={Jamil, Haziq},
  school={London School of Economics and Political Science},
  year=2018,
  month=10
}
```

## Abstract

Regression analysis is undoubtedly an important tool to understand the relationship between one or more explanatory and independent variables of interest. 
In this thesis, we explore a novel methodology for fitting a wide range of parametric and nonparametric regression models, called the I-prior methodology [(Bergsma 2018)](https://arxiv.org/abs/1707.00274).

We assume that the regression function belongs to a reproducing kernel Hilbert or Kreĭn space of functions, and by doing so, allows us to utilise the convenient topologies of these vector spaces. 
This is important for the derivation of the Fisher information of the regression function, which might be infinite dimensional.
Based on the principle of maximum entropy, an I-prior is an objective Gaussian process prior for the regression function with covariance function proportional to its Fisher information. 

Our work focusses on the statistical methodology and computational aspects of fitting I-priors models. 
We examine a likelihood-based approach (direct optimisation and EM algorithm) for fitting I-prior models with normally distributed errors.
The culmination of this work is the `R` package [`iprior`](https://cran.r-project.org/package=iprior) which has been made publicly available on CRAN. 
The normal I-prior methodology is subsequently extended to fit categorical response models, achieved by "squashing" the regression functions through a probit sigmoid function.
Estimation of I-probit models, as we call it, proves challenging due to the intractable integral involved in computing the likelihood. 
We overcome this difficulty by way of variational approximations.
Finally, we turn to a fully Bayesian approach of variable selection using I-priors for linear models to tackle multicollinearity.

We illustrate the use of I-priors in various simulated and real-data examples. 
Our study advocates the I-prior methodology as being a simple, intuitive, and comparable alternative to similar leading state-of-the-art models. 

## Links

- My PhD website - [link](http://phd.haziqj.ml)
- Code for replciation - [link](http://myphdcode.haziqj.ml)
- Viva presentation - [link](https://haziqj.github.io/phd-thesis/)
- R/iprior - [GitHub](https://github.com/haziqj/iprior), [CRAN](https://cran.r-project.org/package=iprior)
- R/iprobit - [GitHub](https://github.com/haziqj/iprobit)
- R/ipriorBVS - [GitHub](https://github.com/haziqj/ipriorBVS)

------------------------------------------------------------------------

The copyright of this thesis rests with the author. 
Quotation from it is permitted, provided that full acknowledgement is made. 
This thesis may not be reproduced without my prior written consent.
I warrant that this authorisation does not, to the best of my belief, infringe the rights of any third party.