library(cmdstanr)
library(matrixStats)

#### Facilitator functions ####
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
COMP_lpmf <- function(k, theta) k * log(theta[1]) - theta[2] * lfactorial(k)
TrueValue <- function(mu, nu){
  if (nu == 1) mu else if (nu == 2) log(besselI(2*sqrt(mu), nu = 0)) else
    logSumExp(sort(COMP_lpmf(k = 0:1e5, c(mu, nu))))}

# system("rm COMP_normalizing")
model <- cmdstan_model("stan/COMP_normalizing.stan",
                       include_paths = "stan")

mu <- 5; nu <- 1/2
epsilon <- 1e-16; M <- 1e5

cdata <- list(Mu = mu, Nu = nu, Epsilon = epsilon, maxIter = M,
              TrueValue = TrueValue(mu, nu) )

comparison <- stanfit(
  model$sample(
    data = cdata,
    iter_warmup = 0,
    iter_sampling = 1,
    fixed_param = TRUE
  )
)

( ext <- rstan::extract(comparison) )

Rmpfr::mpfr(cdata$TrueValue, 1000)
Rmpfr::mpfr(ext$truth_1, 1000)
Rmpfr::mpfr(ext$truth_2, 1000)
Rmpfr::mpfr(ext$difference, 1000)
Rmpfr::mpfr(ext$difference2, 1000)
