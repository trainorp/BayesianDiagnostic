functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=2> ncat;  // number of categories 
  int<lower=1> K_muType1MI;  // number of population-level effects 
  matrix[N, K_muType1MI] X_muType1MI;  // population-level design matrix 
  int<lower=1> K_muType2MI;  // number of population-level effects 
  matrix[N, K_muType2MI] X_muType2MI;  // population-level design matrix 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc_muType1MI = K_muType1MI - 1; 
  matrix[N, K_muType1MI - 1] Xc_muType1MI;  // centered version of X_muType1MI 
  vector[K_muType1MI - 1] means_X_muType1MI;  // column means of X_muType1MI before centering 
  int Kc_muType2MI = K_muType2MI - 1; 
  matrix[N, K_muType2MI - 1] Xc_muType2MI;  // centered version of X_muType2MI 
  vector[K_muType2MI - 1] means_X_muType2MI;  // column means of X_muType2MI before centering 
  for (i in 2:K_muType1MI) { 
    means_X_muType1MI[i - 1] = mean(X_muType1MI[, i]); 
    Xc_muType1MI[, i - 1] = X_muType1MI[, i] - means_X_muType1MI[i - 1]; 
  } 
  for (i in 2:K_muType2MI) { 
    means_X_muType2MI[i - 1] = mean(X_muType2MI[, i]); 
    Xc_muType2MI[, i - 1] = X_muType2MI[, i] - means_X_muType2MI[i - 1]; 
  } 
} 
parameters { 
  vector[Kc_muType1MI] b_muType1MI;  // population-level effects 
  real temp_muType1MI_Intercept;  // temporary intercept 
  vector[Kc_muType2MI] b_muType2MI;  // population-level effects 
  real temp_muType2MI_Intercept;  // temporary intercept
  vector[Kc_muType1MI] lambda;
  real<lower=0> tau; 
} 
transformed parameters { 
} 
model { 
  vector[N] muType1MI = Xc_muType1MI * b_muType1MI + temp_muType1MI_Intercept; 
  vector[N] muType2MI = Xc_muType2MI * b_muType2MI + temp_muType2MI_Intercept; 
  // linear predictor matrix 
  vector[ncat] mu[N]; 
  for (n in 1:N) { 
    mu[n] = [0, muType1MI[n], muType2MI[n]]';
  } 
  // priors including all constants
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, 1); 
  b_muType1MI ~ normal(0, 1); 
  temp_muType1MI_Intercept ~ normal(0, 4); 
  b_muType2MI ~ normal(0, 1); 
  temp_muType2MI_Intercept ~ normal(0, 4); 
  // likelihood including all constants 
  if (!prior_only) { 
    for (n in 1:N) { 
      target += categorical_logit_lpmf(Y[n] | mu[n]); 
    } 
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_muType1MI_Intercept = temp_muType1MI_Intercept - dot_product(means_X_muType1MI, b_muType1MI); 
  // actual population-level intercept 
  real b_muType2MI_Intercept = temp_muType2MI_Intercept - dot_product(means_X_muType2MI, b_muType2MI); 
} 

