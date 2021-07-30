functions{
  real logFunction(int k, real[] p){
    return k * p[1] - p[2] * lgamma(k + 1);
  }
  real true_value(real[] p, int M){
    if(p[2] == 1){
      return(exp(p[1]));
    }else{
      if(p[2] == 2){
        return(
          log(modified_bessel_first_kind(0, 2*sqrt(exp(p[1]))))
          );
      }else{
        real lps[M];
        for(i in 1:M) lps[i] = logFunction(i-1, p);
        return( log_sum_exp(sort_asc(lps)) );
      }
    } 
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
  real params[2] = {log(Mu), Nu};
  real TV = true_value(params, maxIter);
}

generated quantities {
  real estimated[2] = infiniteAdaptive(params, Epsilon, maxIter, negative_infinity(), 0);
  real difference;
  real difference2;
  real truth_1 = TrueValue;
  real truth_2 = TV;
  
  if (TrueValue > estimated[1])
  difference = exp(log_diff_exp(TrueValue, estimated[1]));
  else
  difference = exp(log_diff_exp(estimated[1], TrueValue));
  print(TV);
  if (TV > estimated[1])
  difference2 = exp(log_diff_exp(TV, estimated[1]));
  else
  difference2 = exp(log_diff_exp(estimated[1], TV));
}
