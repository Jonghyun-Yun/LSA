functions {
  // real lambda_g(int g, int y, row_vector lambda1, row_vector lambda2) {
  //   real out;
  //   if (y == 1) out = lambda1[g];
  //   else if (y == 0) out = lambda2[g];
  //   else out = 1;
  //   return log(out);
  // }
  // real vcum_hazard(row_vector lambda, row_vector len, int[] g, row_vector H) {
  //   int N = size(g);
  //   real out = 0;
  //   for (n in 1:N){
  //     if (g[n] > 1) out += sum(lambda[1:(g[n]-1)] .* len[1:(g[n]-1)]) + lambda[g[n]] * H[n];
  //     else out += lambda[g[n]] * H[n];
  //   }
  //   return out;
  // }
  // real log_hazard(int y, row_vector beta, row_vector theta, vector gamma, row_vector z1, row_vector z2, row_vector w) {
  //   real out;
  //   if (y == 1) out = beta[1] + theta[1] - gamma[1] * distance(z1, w);
  //   else if (y == 0) out = beta[2] + theta[2] - gamma[2] * distance(z2, w);
  //   else out = 0;
  //   return out;
  // }
  real vlog_hazard(int[] y, row_vector beta, matrix theta,
                   vector gamma, matrix z1, matrix z2, row_vector w) {
    int N = size(y);
    real out = 0;
    for (n in 1:N){
      if (y[n] == 1) out += beta[1] + theta[n,1] - gamma[1] * distance(z1[n], w);
      else if (y[n] == 0) out += beta[2] + theta[n,2] - gamma[2] * distance(z2[n], w);
    }
    return out;
  }
  real cum_hazard(row_vector lambda, row_vector len, int g, real H) {
    real out;
    if (g > 1) out = sum(lambda[1:(g-1)] .* len[1:(g-1)]) + lambda[g] * H;
    else out = lambda[g] * H;
    return out;
  }
  real cum_vhazard(row_vector beta, matrix theta,
                   vector gamma, matrix z1, matrix z2, row_vector w,
                   row_vector lambda1, row_vector lambda2,
                   row_vector len, int[] g, row_vector H) {
    int N = size(g);
    real out = 0;
    for (n in 1:N){
      out += cum_hazard(lambda1, len, g[n], H[n]) * exp(beta[1] + theta[n,1] - gamma[1] * distance(z1[n], w));
      out += cum_hazard(lambda2, len, g[n], H[n]) * exp(beta[2] + theta[n,2] - gamma[2] * distance(z2[n], w));
    }
    return out;
  }
  real vlambda_g(int[] g, int[] y, row_vector lambda1, row_vector lambda2) {
    int N = size(y);
    real out = 0;
    for (n in 1:N){
      if (y[n] == 1) out += log(lambda1[g[n]]);
      else if (y[n] == 0) out += log(lambda2[g[n]]);
    }
    return out;
  }
  real part_sum(int[] mi,
                int start, int end,
                int[] mseg,
                row_vector lambda1, row_vector lambda2,
                row_vector beta, matrix theta, vector gamma,
                matrix z1, matrix z2, row_vector w,
                row_vector mlen, row_vector mh
                ) {
    return vlambda_g(mseg[start:end], mi, lambda1, lambda2) +
      vlog_hazard(mi, beta, theta[start:end], gamma, z1[start:end], z2[start:end], w) +
      cum_vhazard(beta, theta[start:end], gamma, z1[start:end], z2[start:end],
      w, lambda1, lambda2, mlen, mseg[start:end], mh[start:end]);
  }
}
data {
  int I;
  int N;
  int G;
  int<lower=1,upper=G> mseg[I,N];
  row_vector<lower=0>[G-1] mlen;
  int<lower=0,upper=2> mi[I,N];
  matrix<lower=0>[I,N] mt;
  matrix<lower=0>[I,N] mh;
}
transformed data {
  // STAN takes row vectors for matrix initialization.
  // hyperparameter and transformed data, etc.
  real<lower=machine_precision()> sa = 1;
  real<lower=machine_precision()> sb = 1;
  real<lower=machine_precision()> la = 0.001;
  real<lower=machine_precision()> lb = 0.001;
}
parameters {
  matrix[N,2] theta;
  matrix[I,2] beta;
  real<lower=machine_precision()> sigma;
  matrix[N,2] z1;
  matrix[N,2] z2;
  matrix[I,2] w;
  vector<lower=machine_precision()>[2] gamma;
  matrix<lower=machine_precision()>[I,G] lambda1;
  matrix<lower=machine_precision()>[I,G] lambda2;
}
model {
  for (i in 1:I) {
    // for (k in 1:N) {
    //   target += lambda_g(mseg[i,k], mi[i,k], lambda1[i], lambda2[i]) +
    //     log_hazard(mi[i,k], beta[i], theta[k], gamma, z1[k], z2[k], w[i]) -
    //     cum_hazard(lambda1[i], mlen, mseg[i,k], mh[i,k]) *
    //     exp(log_hazard(1, beta[i], theta[k], gamma, z1[k], z2[k], w[i])) -
    //     cum_hazard(lambda2[i], mlen, mseg[i,k], mh[i,k]) *
    //     exp(log_hazard(2, beta[i], theta[k], gamma, z1[k], z2[k], w[i]));
    // }

    int grainsize = 1;
    target += reduce_sum(part_sum, mi[i], grainsize,
                         mseg[i],
                         lambda1[i], lambda2[i],
                         beta[i], theta, gamma,
                         z1, z2, w[i], mlen, mh[i]);

    // target += vlambda_g(mseg[i], mi[i], lambda1[i], lambda2[i]) +
    //   vlog_hazard(mi[i], beta[i], theta, gamma, z1, z2, w[i]) +
    //   cum_vhazard(beta[i], theta, gamma, z1, z2,
    //   w[i], lambda1[i], lambda2[i], mlen, mseg[i], mh[i]);

  }
  gamma ~ lognormal(0, 1);
  to_vector(z1) ~ std_normal(); // normal(mu, sigma)
  to_vector(z2) ~ std_normal(); // normal(mu, sigma)
  to_vector(theta) ~ normal(0,sigma);
  to_vector(beta) ~ normal(0,1);
  to_vector(lambda1) ~ gamma(la,lb);
  to_vector(lambda2) ~ gamma(la,lb);
  to_vector(w) ~ normal(0,1);
  target += inv_gamma_lpdf(sigma^2 | sa,sb);
}
