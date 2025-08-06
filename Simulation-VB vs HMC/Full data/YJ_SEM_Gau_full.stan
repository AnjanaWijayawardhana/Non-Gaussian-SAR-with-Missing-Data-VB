data {
  int<lower=0> n;                // Number of data points
  int<lower=0> p;                // Number of predictors+1, regression 


  

  matrix[n, p] x;                // Predictor matrix
  vector[n]   y;     // Response vector
  
  matrix[n, n] w;
  matrix[n, n] I;
  matrix[n, n] wpluwt;
  matrix[n, n] wtw;
  
  vector[p] mu_beta;
  matrix[p, p] sigma2_beta;
}

parameters {
  vector[p] beta;
  real omegadash; //s2
  real rhodash; // rho
  real gamadash; // gama
}



transformed parameters {
  real sigma2;  // Standard deviation
  real rho;     // Correlation
  real gama;    // Transformed parameter
  
  matrix[n, n] M;
  vector[n] z;
  
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
  
  
}
