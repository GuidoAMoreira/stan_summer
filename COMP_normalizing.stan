functions{
  real logFunction(int k, real[] p){
    return k * p[1] - p[2] * lgamma(k + 1);
  }
  
  // logFunction must be defined beforehand
  #include infiniteAdaptive.stan
}

data{
  real<lower=0> Mu;
  real<lower=0> Nu;
  real<lower=0> Epsilon;
  real TrueValue;
  int<lower=0> maxIter;
}

transformed data {
  real lMu = log(Mu);
  real params[2];
  
  params[1] = lMu;
  params[2] = Nu;
}

generated quantities {
  real estimated[2] = infiniteAdaptive(params, Epsilon, maxIter, negative_infinity(), 0);
  real difference;
  
  if (TrueValue > estimated[1])
    difference = exp(log_diff_exp(TrueValue, estimated[1]));
  else
    difference = exp(log_diff_exp(estimated[1], TrueValue));
}
