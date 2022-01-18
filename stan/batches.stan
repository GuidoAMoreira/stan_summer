// c-folding inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
real[] infiniteBatches(real[] p, int c, int N_start, real epsilon, int maxIter, int n0) {
  int N_inc = c * N_start;
  real storeVal[N_inc + 1];
  real leps = log(epsilon);
  int n0_ = n0;
  int nextCheckPoint = N_start - 1;
  int lastCheckPoint = nextCheckPoint + 1;
  real summedIncrement;
  int n = N_inc;
  int pos = 1;
  real check;
  // Setting up first batches
  for (i in 1:N_inc) {
    storeVal[i + 1] = logFunction(n0_, p);
    n0_ += 1;
  }
  storeVal[1] = log_sum_exp(sort_asc(storeVal[:(N_start + 1)]));
  summedIncrement = log_sum_exp(sort_asc(storeVal[(N_start + 2):]));

  // Find the maximum
  while (n < maxIter && (
    summedIncrement > leps ||
    storeVal[N_inc + 1] - storeVal[N_inc] >
      -log1p_exp(storeVal[N_inc + 1] - summedIncrement) ||
    storeVal[N_inc + 1] > storeVal[N_inc] ||
    is_inf(storeVal[N_inc + 1])
  )) {
    storeVal[1] = log_sum_exp(storeVal[1], summedIncrement);
    for (i in 2:(N_inc + 1)) {
      storeVal[i] = logFunction(n0_, p);
      n0_ += 1;
    }
    summedIncrement = log_sum_exp(sort_asc(storeVal[2:]));
    n += N_inc;
    if (n >= maxIter) break; // Return if maxIter is reached
  }
  return {log_sum_exp(storeVal[1], summedIncrement), 1. * n};
}
