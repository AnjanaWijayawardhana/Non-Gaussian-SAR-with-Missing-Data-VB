data {
  int<lower=0> n;                // Number of data points
  int<lower=0> p;                // Number of predictors+1, regression 
  int<lower=0> q;                // Number logistic paras
  int<lower=0> no;
  int<lower=0> nu;
  
  matrix[n, (q-1)] xm;
  matrix[n, p] x;                // Predictor matrix
  matrix[nu, p] xu;
  matrix[no, p] xo;
  vector[no] yo;                 // Response vector
  int<lower=0, upper=1> m[n];
  
  matrix[n, n] w;
  matrix[n, n] I;
  matrix[n, n] wpluwt;
  matrix[n, n] wtw;
  
  vector[p] mu_beta;
  matrix[p, p] sigma2_beta;
  vector[q] mu_psi;
  matrix[q, q] sigma2_psi;
}

parameters {
  vector[p] beta;
  real omegadash; //s2
  real rhodash; // rho
  real gamadash; // gama
  vector[q] psi;
  vector[nu] yu; //latent (yu) treats as parameters, 
  //so should be defined as a parameter
}



transformed parameters {
  real sigma2;  // Standard deviation
  real rho;     // Correlation
  real gama;    // Transformed parameter
  
  matrix[n, n] M;
  vector[n] z;
  
  // Create a combined vector 'y' (observed + latent values)
  vector[n] y;
  for (i in 1:no) {
    y[i] = yo[i];  // Copy observed data
  }
  for (i in 1:nu) {
    y[no + i] = yu[i];  // Append latent data
  }
  
  // Compute transformations
  sigma2 = exp(omegadash); // Exponential transformation
  rho = (exp(rhodash) - 1) / (exp(rhodash) + 1); // Transform rhodash to rho
  gama = 2 * exp(gamadash) / (exp(gamadash) + 1); // Transform gamadash to gama
  
  // Compute matrix M
  M = I - rho * wpluwt + rho * rho * wtw;
  
  // Compute z using conditional operations without if_else
  for (i in 1:n) {
    // Use element-wise condition on y[i]
    if (y[i] > 0) {
      z[i] = (pow(y[i]+1,gama)-1)/gama; // Positive values
    } else {
      z[i] = -(pow(-y[i]+1,2-gama)-1)/(2-gama); // Non-positive values
    }
  }
}


model {
  // Jacobian adjustment for the transformation
  real log_product_jacob=0;
  
  // Priors
  target += multi_normal_lpdf(beta | mu_beta, sigma2_beta);  // Prior for regression coefficients
  target += normal_lpdf(omegadash | 0, 10000);              // Prior for omegadash
  target += normal_lpdf(rhodash | 0, 10000);                // Prior for rhodash
  target += normal_lpdf(gamadash | 0, 10000);               // Prior for gamadash
  target += multi_normal_lpdf(psi | mu_psi, sigma2_psi);    // Prior for logistic parameters
  
  
  
  for (i in 1:n) {
    // Use element-wise condition on y[i]
    if (y[i] > 0) {
      log_product_jacob=log_product_jacob+(gama-1)*log(y[i]+1);
    } else {
      log_product_jacob=log_product_jacob+(1-gama)*log(-y[i]+1); // Non-positive values
    }
  }
  
  target += log_product_jacob;
  
  // Likelihood
  target += -0.5 * n * log(2 * pi())
  - 0.5 * n * log(sigma2)
  + 0.5 * 2 * sum(log(diagonal(cholesky_decompose(M))))
  - 0.5 * (transpose(z - x * beta) * M * (z - x * beta)) / sigma2;
  
  // Missing data model
  target += bernoulli_logit_lpmf(m | xm * psi[1:(q - 1)] + y * psi[q]);
  
}
