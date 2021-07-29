model <- cmdstanr::cmdstan_model("COMP_normalizing.stan")

mu <- 2; nu <- 0.5
epsilon <- 1e-16; M <- 1e5

opt_adaptive <- adaptive_impl$sample(
  data = list(Mu = mu, Nu = nu, Epsilon = epsilon, maxIter = M,
              TrueValue = matrixStats::logSumExp(
                sort(
                  0:M |> {\(k) k * log(mu) - nu * lfactorial(k)}()
                )
              ))
)