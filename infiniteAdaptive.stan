// Internal use
int adaptiveConvergenceCheck(real oldT, real newT, real lepsilon, real log1mL) {
  if (!log1mL) return (oldT - newT < log1p_exp(newT - lepsilon));
  
  real logZ = newT + oldT - log_diff_exp(oldT, newT);
  real ls = newT - log1mL;
  
  if (logZ > ls){
    return (log_diff_exp(logZ, ls) >= lepsilon);
  }else{
    return (log_diff_exp(ls, logZ) >= lepsilon);
  }
}

// Adaptive inifinite sum algorithm
// Requires definition of logFunction with two arguments:
// int k and real[] parameters
real[] infiniteAdaptive(real[] p, real epsilon, int maxIter, real logL, int n0) {
  vector[maxIter+2] storeVal;
  real leps = log(epsilon) + log2();
  int n = 1;
  int n0_ = n0;
  real old_term;
  real new_term;
  real log1mL = log_diff_exp(0, logL);
  
  // Setting up first iterations
  old_term = logFunction(n0_, p);
  n0_ += 1;
  new_term = logFunction(n0_, p);
  n0_ += 1;
  storeVal[1] = old_term;
  
  // Find the maximum
  while (new_term >= old_term) {
    old_term = new_term;
    new_term = logFunction(n0_, p);
    n0_ += 1;
    n += 1;
    storeVal[n] = old_term;
    if (n == maxIter) return({log_sum_exp(storeVal[:n]), n}); // Return if maxIter is reached
  }
  // Start testing convergence after the maximum
  while ( adaptiveConvergenceCheck(old_term, new_term, leps, log1mL) == 0 )  {
    old_term = new_term;
    new_term = logFunction(n0_, p);
    n0_ += 1;
    n += 1;
    storeVal[n] = old_term;
    if (n == maxIter) break;
  }
  //storeVal[n+1] = new_term;
  storeVal[n+1] = new_term - log_diff_exp(0, new_term - old_term) - log2();
  storeVal[n+2] = new_term - log2();
  return {log_sum_exp(storeVal[:(n+2)]), n};
}
