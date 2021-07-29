library(cmdstanr)
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())

model <- cmdstanr::cmdstan_model("COMP_normalizing.stan",
                                 include_paths = ".")

mu <- 2; nu <- 0.5
epsilon <- 1e-16; M <- 1e5

COMP_lpmf <- function(k, theta){
  mu <- theta[1]
  nu <- theta[2]
  return(
    k * log(mu) - nu*lfactorial(k)  
  )
}

comparison <- stanfit(
  model$sample(
    data = list(Mu = mu, Nu = nu, Epsilon = epsilon, maxIter = M,
                TrueValue = matrixStats::logSumExp(
                  sort(
                    COMP_lpmf(k = 0:M, theta = c(mu, nu))
                  )
                )),
    iter_warmup = 0,
    iter_sampling = 1,
    fixed_param = TRUE
  )
)

rstan::extract(comparison)
