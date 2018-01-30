library(quanteda)

# Load my text
my.text <- texts(readtext::readtext("wordcloud.md"))
names(my.text) <- "README"
summary(my.text)

# Reprocess
my.text.lower <- char_tolower(my.text)
my.text.word <- as.character(tokens(my.text.lower, remove_punct = TRUE))
my.text.wordp <- as.character(tokens(my.text.lower, remove_punct = FALSE))

my.dict <- dictionary(list(
  "I-prior"       = c("I-prior", "I-priors", "I-probit"),
  "Gaussian"      = c("normal*", "Gaussian", "process"),
  "function"      = "function*",
  "model"         = "model*",
  "posterior"     = c("posterior*", "mean*"),
  "MCMC"          = c("markov", "chain", "monte", "carlo", "mcmc", "stochastic"),
  "computation"   = c("computation*", "algorithm", "procedure", "numerical", "stab*"),
  "approximation" = "approxim*",
  "method"        = "method*",
  "RKHS"          = c("rkhs", "space", "rkhss", "reproducing", "hilbert", "space"),
  "Bayes"         = c("bayes*", "prior*"),
  "estimation"    = "estimat*",
  "method"        = c("method*", "approach*", "framework"),
  "Fisher-information" = c("fisher", "information"),
  "inverse"       = "invers*",
  "variable-selection" = c("selection", "chosen", "choos*", "variable*", "100-variate"),
  "Newton"        = c("newton*", "maxim*", "newton-based", "*optima*", "local", "hessian"),
  "exponential-family"   = c("exponential*", "se", "family"),
  "distributions" = c("distribution*", "densit*", "joint"),
  "probability"   = c("probabilit*"),
  "smooth"        = "*smooth*",
  "update"        = c("updat*", "cycle*"),
  "closed-form"   = c("closed", "closed-form", "form"),
  "likelihood"    = "likelihood*",
  "advantage"     = c("advantage*", "convenien*"),
  "kernel"        = "kernel*",
  "efficient"     = "efficient*",
  "fBm"           = c("fbm", "fractional", "brownian", "motion", "hurst"),
  "intractable"   = c("intractable", "unfeasible"),
  "classification" = c("class*", "categorical", "multino*", "outcome*", "m"),
  "inference"     = c("inference", "compar*"),
  "linear"        = c("linear", "canonical"),
  "matrix"        = c("matrix", "inver*"),
  "independent"   = c("independent", "isotropic"),
  "entropy"       = c("entropy", "jayne*", "objective", "principle"),
  "variance"      = c("variance*", "covariance*"),
  "regularisation" = c("regularisation", "tikhonov", "lasso", "ridge"),
  "variational"   = c("variation*", "calculus", "mean-field", "minimis*"),
  "probit"        = c("probit", "link"),
  "R"             = c("iprior", "iprobit", "package", "r", "github", "cran"),
  "Nyström"       = c("Nystrom", "low-rank"),
  "parameters"    = c("*parameter*", "hyper"),
  "machine-learning" = c("machine*", "learning", "literature"),
  "multicollinearity" = c("multicollinearity", "autocorr*"),
  "efficient"     = c("efficien*"),
  "Kullback-Leibler" = c("kullback", "leibler", "kullback-leibler", "kl", "divergence"),
  "eigen"         = c("eigen*", "spectral", "theorem"),
  "EM"            = c("em", "expectation-maximisation", "m-step", "e-step"),
  "Hamiltonian-MC" = c("hamilton*", "dynamic*"),
  "Gibbs"         = c("gibbs", "sampl*"),
  "Laplace"       = "laplace",
  "GLM"           = "general*",
  "random-effects" = c("random", "effects"),
  "feature-space" = "feature",
  "Bernoulli"     = c("bernoulli", "binary"),
  "empirical-Bayes" = c("empirical", "*factor*", "empirical-bayes"),
  "scale"         = "scal*",
  "simple"        = c("simple", "eas*"),
  "better"        = c("better", "superior*"),
  "integration"   = "integr*",
  "multidimensional" = "*dimensional",
  "fast"          = c("fast*", "quick*"),
  "longitudinal"  = c("longitudinal", "time"),
  "simulations"   = c("simulation*", "experiment*", "toy", "example*"),
  "iterative"     = c("iterative", "scheme"),
  "precision"     = c("precision", "psi"),
  "real-data"     = c("real*", "data", "real-data", "*world", "real-world")
), tolower = FALSE)

remove.this <- c(
  "can", "case", "also", "well", "may", "find", "one", "number", "due", "ie", "step",
  "allow*", "now", "make*", "done", "several", "result*", "zero", "second",
  "true", "gives", "using", "found", "use", "towards", "two", "us",
  "nice", "project", "work*", "say", "terms", "assum*", "fully", "incl",
  "requir*", "incl*", "interest*", "different", "relevant", "often", "call",
  "employed", "typical*", "much", "however", "even", "note", "final*", "implement*",
  "derive", "based", "rise", "provide*", "certain", "other*", "hurdle*",
  "operation", "cycle", "known", "longer", "amount*", "especially", "perform",
  "whole", "issue", "set", "main", "cases", "three", "first", "n-by-n", "appropriate",
  "negative", "representation", "making", "size", "employ", "influence",
  "fact", "equal", "third"
)

# Wordcloud
my.text.dfm <- dfm(my.text, remove_punct = TRUE, verbose = FALSE,
                   remove = c(stopwords(), remove.this))
my.text.dfm <- dfm_lookup(my.text.dfm, my.dict, exclusive = FALSE, capkeys = FALSE)
topfeatures(my.text.dfm, 100)
textplot_wordcloud(my.text.dfm, colors = iprior::ggColPal(6), max.words = Inf,
                   min.freq = 2, scale = c(5.7, 1.7), rot.per = 0.2, random.color = FALSE)

textplot_wordcloud2 <- function(x, ...) wordcloud2::wordcloud2(data.frame(word = featnames(x), freq = colSums(x)), ...)
textplot_wordcloud2(my.text.dfm, fontFamily = "Helvetica", color = iprior::ggColPal(6),
                    minRotation = 0, maxRotation = pi/2, rotateRatio = 1)
