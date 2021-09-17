// Sum-To-Threshold inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
real[] infiniteSumToThreshold(real[] p, int[] ns, real epsilon, int maxIter, int n0) {
  vector[maxIter + 1] storeVal;
  real leps = log(epsilon);
  int n = 2;
  int n0_ = n0;

  // Setting up first iterations
  storeVal[1] = logFunction(n0_, p, ns);
  n0_ += 1;
  storeVal[2] = logFunction(n0_, p, ns);
  
  // Find the maximum
  while (storeVal[n] > storeVal[n - 1]) {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction(n0_, p, ns);
    if (n >= maxIter) return({log_sum_exp(storeVal[:n]), 1. * n}); // Return if maxIter is reached
  }
  // Start testing convergence after the maximum
  while ( storeVal[n] > leps )  {
    n0_ += 1;
    n += 1;
    storeVal[n] = logFunction(n0_, p, ns);
    if (n >= maxIter) break;
  }
  return {log_sum_exp(sort_asc(storeVal[:n])), 1. * n};
}
