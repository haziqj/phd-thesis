# Regression modelling using priors with Fisher information covariance kernels (I-priors).

There are three main objectives of this project. The first is to outline an efficient computational method for estimating the parameters of an I-prior model. The second is to extend the I-prior methodology to categorical responses for classification and inference. The third is to explore the usefulness of I-priors towards variable selection for linear models.

# Introduction

We assume a general regression framework, where the responses are modelled as following a noisy version of a true regression function, which is to be estimated. We assume a normal distribution for the noise terms, having zero mean and precision Psi. This is the more general case, and for our purposes we make a further assumption of an isotropic precision matrix, ie independent error terms. The space of covariates can be unidimensional, multidimensional, or even functional covariates.

By choosing to work in a reproducing kernel Hilbert space of functions, we get very nice topologies and conveniences. In particular, this allows us to derive the Fisher information for our regression function. We derive an objective prior, one that is based on Jayne’s entropy maximisation principle, for our regression function (subject to certain constraints). This results in what we call the I-prior: a Gaussian distribution on the function with some mean chosen a priori, and covariance equal to the Fisher information. This has a nice intuition: much information regarding the function implies a larger variance, and hence a smaller influence of the prior mean and subsequently more influence of the data - and vice versa. 

Having defined such a prior, we are then interested in two things: the posterior distribution of our regression function; and the posterior predictive distribution for new data points. In each case, a point estimate can be taken to be the posterior mean, and its uncertainty quantified by the posterior variance.

# Advantages

The I-prior represents a unifying methodology for various regression models, including unidimensional or multidimensional smoothing models; random effects or multilevel models; longitudinal models; and even models with functional covariates. 

This is achieved by choosing the RKHS for which our regression function lives in. For a linear effect of the covariates, the canonical or linear RKHS is chosen. For a smooth effect of the covariates, the fractional Brownian motion (fBm) RKHS is chosen. The smoothness of the function in the fBm RKHS is controlled by the Hurst parameter. We find that the use of the fBm RKHS for smoothing gives superior results to other smooth RKHSs, such as squared exponential (SE) RKHS say. FBm I-prior paths are smooth enough to work, but not as smooth as SE RKHS say, which may tend to oversmooth the effects. For a categorical effect of covariates, the Pearson RKHS is chosen. Such covariates could, besides being truly categorical in nature, can provide information such as groupings/levels/time in random effects/multilevel/longitudinal models.

There is a one-to-one mapping between the set of positive definite functions (kernels) and the set of RKHSs. Conveniently, this means that the choice of kernel gives rise to a specific RKHS, and furthermore, kernel functions can be combined to achieve a combination of effects, and we are guaranteed a RKHS for our regression function. Note that the scaling of the RKHS towards the function may be arbitrary, and as such we typically multiply the kernel function with a scale parameter, which is to be estimated.

As we will see, I-prior models are straightforward to estimate and provides an easy method for model comparison and inference through the RKHS scale parameters. We also find that the I-prior often gives better prediction for new data compared to other methods such as Tikhonov regularisation.

# Estimation

In this continuous-response case, the assumed normality of the relevant distributions mean that the posterior densities can be derived in closed form. There are several model (hyper-)parameters, such as the error precision, RKHS scale parameters, and any other parameters relating to the kernel, which need to be estimated. These may be estimated in several different ways.

By integrating out the prior from the joint distribution of the responses and I-prior, we arrive at the marginal likelihood. This can then be maximised with respect to the model parameters. This method of maximising the (marginal) likelihood is also known as an empirical-Bayes approach. A Newton-based approach can be employed for this maximisation, and for very simple problems (ie small number of scale parameters) this might well be the fastest approach. However, as is well know, Newton methods are quite sensitive to starting values, and we have often encountered sub-optimal results due to the presence of multiple local optima of I-prior marginal likelihoods.

The second approach would be employing the iterative Expectation-Maximisation (EM) algorithm. The I-prior model is unique, in that it allows for a convenient reparameterisation of the prior function into isotropic, normal random effects. As such, the resulting EM algorithm is very straightforward to employ - and in fact, for most I-prior models, the M-step for the model parameters can be found in closed form, which makes for an efficient updating cycle. The EM algorithm is implemented in the R package iprior, which has been published to CRAN and is also available on GitHub.

Finally, one could also attempt a fully Bayesian treatment of the model, in which prior distributions are assigned to the model parameters. Markov Chain Monte Carlo (MCMC) methods such as Gibbs sampling can be employed to obtain posterior estimates for the parameters. In our experience, the Gibbs sampling can suffer from severe autocorrelation in the posterior samples, especially when there are a lot of scale parameters to estimate. A better approach to obtaining posterior samples is Hamiltonian Monte Carlo, which makes use of Hamiltonian dynamics to efficiently map out the posterior state space and can avoid autocorrelated samples.

# Computational Hurdles

Regardless of the estimation procedure, computational complexity is dominated by the inversion of a n-by-n covariance matrix. In the case of Newton-based approaches, this needs to be evaluated at each Newtons step. In the case of the EM algorithm, every update cycle also involves such an inversion. For stochastic MCMC sampling methods, each sample would also involve this inversion operation.

I-priors, while being philosophically different from Gaussian process priors, do share the same computational hurdle. As such, several methods exist in the machine learning literature to overcome this issue. Amongst others, is a method to approximate the covariance kernel by a low-rank matrix, so that the most expensive operation of inverting a n-by-n matrix is greatly reduced. Our approach is to apply the Nystrom method of low-rank matrix approximation, and we find that this works reasonably well for certain types of RKHS. In general, one could also predetermine whether or not the Nystrom method works by inspecting the rate of decay of the eigenvalues of the covariance kernel.

Another computational hurdle is to ensure numerical stability. We find that due to the structure of the marginal covariance of the responses, numerical instabilities can and are likely to occur - which may give rise to embarassments such as negative covariances. We employ a stable eigendecomposition regime which allows us to efficiently calculate matrix squares and inverses by making use of the spectral theorem.

# Categorical Responses

Suppose now we are interested in a regression model where the responses are categorical. Assume a categorical distribution on the responses with certain probabilities for each class and for each observation. This is of course a generalisation of the Bernoulli distribution to more than two possible outcomes. The question is how can we relate the effect of the covariates through the function, which has unrestricted range, to the responses, which may only take one of m several outcomes? In the spirit of generalised linear models, we answer this by making use of an appropriate link function, and our case, the probit link function. In the binary case, this amounts to squashing our regression function through the (inverse) probit link function in order to model probabilities which are between zero and one. This idea is then extended to the multinoulli case, giving rise to a multinomial probit I-prior model, which we call I-probit.

The main issue with estimation now is that because our responses no longer follow a Gaussian distribution, the relevant marginal distribution, on which the posterior depends, can no longer be found in closed form. The integral required to perform the calculation is intractable, and the focus now is on methods to adequately approximate the integral.

In the Bayesian literature, the Laplace approximation amounts to approximating the posterior distribution with a normal distribution centred around the mode of the integrand. Additionally, the covariance matrix is equal to the inverse (negative) Hessian. Having approximated the posterior by a Gaussian distribution, one could then proceed to find the marginal easily, which is then maximised. Due to the Newton step in the Laplace step, the whole procedure scales cubicly with both the sample size and the number of outcomes, which makes it undesirable to implement.

MCMC methods such as Gibbs sampling or Hamiltonian Monte Carlo can also provide a stochastic approximation to the integral. Unfortunately, the difficulties faced in the continuous case for MCMC methods also present themselves in the categorical case.

Deterministic approaches such as quadrature methods prove unfeasible. Quadrature methods scale exponentially with the variables of integration, in our case, is the sample size.

# Variational Approximation

We consider a type of approximation based on minimising the Kullback-Leibler (KL) divergence from the approximating density to the true posterior density. This is done without making any distributional assumptions, only that the our approximating density factorises over its components (ie an independence assumption) - this is known as the mean-field approximation, which has its roots in the physics literature. As an aside, the term 'variational' stems from the fact that a minimisation of a functional, rather than a function, is involved, and this requires the calculus of variations.

By working in a fully Bayesian setting, we append the model parameters to the list of unknowns in which to estimate, and employ the variational approximation to find a suitable approximation to the required posterior density. The result is an iterative algorithm, similar to the EM. 

As this variational EM works harmoniously with exponential family distributions, the probit link is much preferred over the logit. Most of the EM updating cycles which could be found in closed form are also applicable in the variational EM. Unlike Gaussian process priors, the variational method does not typically result in closed-form updates for the RKHS scale parameters. In such cases, an additional step such as importance sampling is required, which arguably reduces efficiency of the whole variational scheme.

The variational EM is implemented in the R package iprobit. This has been shown to work well for several toy examples as well as real world applications. In the binary case, the I-prior outperforms other popular classication methods including k-nearest neighbours, support vector machines, and Gaussian process classification.

# Variable Selection for Linear Models

As mentioned earlier, model selection can easily be done by comparing likelihoods (empirical Bayes factors). However, with a large number of variables, these pairwise comparisons quickly become unfeasible to perform. We suggest a fully Bayesian approach to estimating posterior model probabilities, and selecting models based on these quantities. For example, one may choose the model which gives the largest posterior model probability (maximum probability model). These are done using MCMC methods (Gibbs sampling).

We restrict ourselves to linear models only. We can easily derive an equivalent I-prior representation by working in the feature space of the betas (linear effects). As a side note, if the dimensions of the linear effects is much, much less than the sample size, then it is worth working in this representation.

We believe the I-prior performs superiorly in cases where there is multicollinearity. This is evidenced by the simulation results that we conducted on a 100-variate experiment, and also in the real-data examples comparing our method with others such as greedy selection, g-priors, and regularisation (ridge and Lasso).

# Conclusions

In this project, we aimed to provide an efficient computational method for estimating I-prior models. The bulk of the computational complexity is proposed to be lessened by way of low-rank matrix approximations, or conversion to a linear feature space if appropriate. The advantages of the original I-prior methodology transfer nicely over to the categorical responses case as well, where estimation was done using variational inference. Finally, the superiority of I-priors in variable selection for linear models is seen, especially in cases where multicollinearity exists between the predictors.
