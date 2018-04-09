Regression analysis is undoubtedly an important tool to understand the relationship between one or more explanatory and independent variables of interest. 
In this thesis, we explore a novel methodology for fitting a wide-range of parametric and non-parametric regression models, called the I-prior methodology. 
We assume that the regression function belongs to a reproducing kernel Hilbert or Krein space of functions, and by doing so, it allows us to utilise the convenient topologies of these vector spaces. 
This is important for the derivation of the Fisher information of the regression function, which might be infinite-dimensional.
Based on the principle of maximum entropy, an I-prior is an objective Gaussian process prior for the regression function with covariance function proportional to its Fisher information. 
Our work focusses on the statistical methodology and computational aspects of fitting I-priors models. 
In the first part of the thesis, we examine a likelihood-based approach (direct optimisation and EM algorithm) for fitting I-prior models with normally distributed errors.
The culmination of this work is the R package 'iprior' which has been made publically available on CRAN.
In the second part, the normal I-prior methodology is extended to fit categorical response models, achieved by "squashing" the regression functions through a probit sigmoid function.
Estimation of I-probit models, as we call it, proves challenging due to the intractable integral involved in computing the likelihood.
We apply a fully Bayesian treatment with a variational approximation in order to obtain the required posterior distributions.
Finally, in the third part, we again turn to a fully Bayesian approach of variable selection for linear models using I-priors.
Our study advocates the I-prior methodology as being a simple, intuitive and comparable, and often times better, alternative to similar leading state-of-the-art models.
We illustrate the use of I-priors in various simulated and real-data examples.

