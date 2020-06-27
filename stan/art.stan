functions {
}
data {
  int<lower=1> I;
  int<lower=1> N;
  int<lower=1> L;
  int<lower=1> G;
  int<lower=1> stuid[L];
  int<lower=1> item[L];
  int<lower=1> seg_g[L];
  real<lower=machine_precision()> len[L];
  int<lower=0,upper=1> status_T[L];
  int<lower=0,upper=1> status_F[L];
}
transformed data {
  real<lower=machine_precision()> sa = 1;
  real<lower=machine_precision()> sb = 1;
  real<lower=machine_precision()> la = 0.001;
  real<lower=machine_precision()> lb = 0.001;
}
// transformed data {
//   // STAN takes row vectors for matrix initialization.
//   // hyperparameter and transformed data, etc.
//   int<lower=0,upper=1> u[n,p,p];
//   real<lower=machine_precision()> a = 0.001;
//   real<lower=machine_precision()> b = 0.001;
//   for (k in 1:n) {
//     for (i in 1:p) {
//       for (j in 1:p) {
//         if (i != j)
//           u[k,i,j] = Y[k,i] * Y[k,j];
//         else
//           u[k,i,j] = 0; }
//     }
//   }
// }
parameters {
  real theta[N,2];
  real beta[I,2];
  real<lower=machine_precision()> gamma[2];
  real<lower=machine_precision()> sigma;
  matrix[N,2] z1;
  matrix[N,2] z2;
  matrix[I,2] w;
  real<lower=machine_precision()> lambda[I,G,2];
}
model {
  for (l in 1:L) {
    status_T[l] ~ poisson(lambda[item[l],seg_g[l],1] * len[l] * exp(beta[item[l],1] + theta[stuid[l],1] - gamma[1] * distance(row(z1,stuid[l]), row(w, item[l]))));
    status_F[l] ~ poisson(lambda[item[l],seg_g[l],2] * len[l] * exp(beta[item[l],2] + theta[stuid[l],2] - gamma[2] * distance(row(z2,stuid[l]), row(w, item[l]))));
  }
  for (c in 1:2) {
    gamma[c] ~ lognormal(0, 1);
      for (k in 1:N) {
        z1[k,c] ~ normal(0,1); // normal(mu, sigma)
        z2[k,c] ~ normal(0,1);
        theta[k,c] ~ normal(0,sigma);
          }
    for (i in 1:I) {
      for(g in 1:G) {
        lambda[i,g,c] ~ gamma(la,lb);
          }
    }
  }
  for (i in 1:I) {
    w[i] ~ normal(0,1);
  }
  target += inv_gamma_lpdf(sigma^2 | sa,sb);
}
