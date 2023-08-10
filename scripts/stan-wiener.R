stan_wiener <- function() {
  
  stan_wiener_base <- write_stan_file("
data {
  int n_data_1;
  int n_data_0;
  array[n_data_1] int rt_1;
  array[n_data_0] int rt_0;
}

parameters {
  real<lower=0> alpha;
  real<lower=0> tau;
  real<lower=0,upper=1> beta;
  real delta;
}

model {
    rt_0 ~ wiener(alpha, tau, beta, delta);

    rt_1 ~ wiener(alpha, tau, 1 - beta, -delta);
    
    alpha ~ normal(1, .01);
    tau ~ normal(.2, .005);
    beta ~ normal(.5, .005);
    delta ~ normal(.2, .005);

}

")
  return(stan_wiener_base)
}
