data{
  int NG; //number of genes
  int NC; //number of cells
  int NR; //number of replicates
  matrix[NG,NR] Y; //response signal
  matrix[NG,NC] X; //signatures
  real tol; //tolerance in data likelihood
}
transformed data{
  //matrix[NR,NR] YtY = Y' * Y;
  matrix[NR,NC] YtX = Y' * X;
  matrix[NC,NC] XtX = X' * X;
}
parameters{
  matrix[NC,NR] bet; //beta coefficients
}
transformed parameters{
  matrix[NC,NR] W = exp(bet); //unnormalized weights
}
model{
  //priors
  to_vector(bet) ~ normal(0,5);
  
  // Y[i,j] ~ normal((XW)[i,j],tol)
  // LL = -sum( (Y[i,j]-(XW)[i,j])^2/(2*tol^2) + log(tol) )
  // = -sum(Y[i,j]^2 + (XW)[i,j]^2 - 2*Y[i,j]*(XW)[i,j] )/(2*tol^2) - sum(log(tol))
  // = sum(Y[i,j]*(XW)[i,j])/tol^2 - sum((XW)[i,j]^2)/(2*tol^2)
  //
  // data likelihood (may need to...)
  //to_vector(Y) ~ normal( to_vector(X * W), tol);
  //NOTE ignoring log(tol) term as not doing inference on this
  target += trace(YtX * W) / (tol*tol); //
  target += -trace_quad_form( XtX, W ) / (2*tol*tol);
}
generated quantities{
  // TODO unclear if outputting proportions at replicate level
  // or introducing hierarchical model for betas and outputting these
  // NOTE in transofming need to use max trick - may be automaticall available in stan
  // https://mc-stan.org/docs/stan-users-guide/floating-point.html#log-sum-of-exponentials
  matrix[NC,NR] P;
  for(j in 1:NR) P[:,j] = softmax(bet[:,j]);
}
