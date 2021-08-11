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
  real difference_check(real v1, real v2){
    if (v1 > v2)
      return(exp(log_diff_exp(v1, v2)));
    else
      return(exp(log_diff_exp(v2, v1)));
  }
  // logFunction must be defined beforehand
  #include infiniteAdaptive.stan
  #include infiniteSumToThreshold.stan
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
  real estimatedAdaptive[2] = infiniteAdaptive(params, Epsilon, maxIter, negative_infinity(), 0);
  real estimatedSumToThreshold[2] = infiniteSumToThreshold(params, Epsilon, maxIter, 0);
  real differenceAdaptive;
  real differenceSumToThreshold;
  real differenceAdaptive2;
  real differenceSumToThreshold2;
  real truth_1 = TrueValue;
  real truth_2 = TV;
  
  differenceAdaptive = difference_check(TrueValue, estimatedAdaptive[1]);
  differenceSumToThreshold = difference_check(TrueValue, estimatedSumToThreshold[1]);
  differenceAdaptive2 = difference_check(TV, estimatedAdaptive[1]);
  differenceSumToThreshold2 = difference_check(TV, estimatedSumToThreshold[1]);
}
