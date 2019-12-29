
#' Subgroup detection in setting of regression with ADMM algorithm
#'
#' The model is \deqn{  Y_i = X_i \beta + Z_i \theta_i.} This ADMM algorithm
#' can detect the subgroup effects of \eqn{\theta_i} and estimate the regression
#' coefficient \eqn{\beta} simultaneously.
#'
#'@param X a n*p matrix, \code{p}-dimensional covariate.
#'@param Y a n*1 vector, response.
#'@param Z0 a n*q matrix, treatment indicator.
#'@param peanlty \code{"SCAD"} or \code{"MCP"}.
#'@param ga a parameter for penalty function, default 3.7.
#'@param nu a tuning parameter for augmented Largranian objection function,
#' the default value is 1.
#'@param lam_vec a series of candidate tuning parameter \code{lambda} for penalty function,
#' which is a vector.
#'@param maxiter The maximum number of iterations. Defaults to 1000.
#'@param ep the absolute convergence tolerance.
#'@return \itemize{
#'    \item \code{BIC} the values of BIC for each \code{lam_vec}.
#'    \item \code{lam_opt} the optimal \code{lambda} that have.
#'    the mimimal \code{BIC}.
#'    \item \code{beta_Est_lam} the regression coefficients matrix
#'    of covariates for different value of \code{lam_vec}, each column denotes a
#'    value of \code{lam_vec}.
#'    \item \code{theta_Est_lam}  the regression coefficients matrix of treatment indicator
#'    for each individuals, each column denotes a different value of \code{lam_vec}.
#'    \item \code{theta_est_opt} the regression coefficients of treatment indicator for
#'    optimal \code{lambda}.
#'    \item \code{beta_est_opt} the regression coefficients of covariates for optimal \code{lambda}.
#'}
#'@examples
#'set.seed(1)
#'n = 200; p = 2
#'x = matrix(rnorm(n*p, mean = 2), nrow = n, ncol = p) # covariates
#'z = sample(c(0,1), n, T)    # treatment indicator
#'beta = c(-2, -3)            # regression coefficients of covariates
#'mu = c(-2, 2)               # regression coefficients of treatment indicator
#'
#'y = x %*% beta + z * mu[1] + rnorm(n, sd = 0.5)
#'index = sample(c(0,1), n, T)
#'y[index == 1] = y[index == 1] + z[index == 1] * (mu[2] - mu[1])
#'lam_vec = seq(0,2,0.2)
#'result = ADMM.q1(y, x, z, penalty = 'SCAD', lam_vec = lam_vec)
#'result$beta_est_opt
#'unique( round(result$theta_est_opt, 1) )
#'
#'# solution paths
#'py_mu = result$theta_Est_lam
#'plot(x=lam_vec, y=py_mu[1,], xlab=expression(lambda),
#'    ylab=expression(paste("Estimated value of  ", vartheta, sep=" ")),
#'    xlim=c(0,2), ylim=c(min(py_mu),max(py_mu)),
#'    main="Solution path for domain I", type='l')
#'for(j in 1:nrow(py_mu) ){
#'    lines(lam_vec, py_mu[j,], lty=j, col=j+1,lwd=2 )
#'}
#'abline(v=result$lam_opt, lty=3, col="lightblue",lwd =3)
#'


ADMM <- function(Y, X, Z0,  penalty, ga = 3.7, nu = 1,
                 lam_vec, maxiter=1000, ep=5e-4 )
{
  # Y: n x 1;
  # X: n x p;
  # Z0: n x q;
  # Z: n x nq;
  # penalty = "SCAD"; ga = 3.7; nu = 1;

  n = nrow(X);
  p = ncol(X);            # Dimension of Covariates.
  q = 1;                  # Dimension of Z

  n_lam = length(lam_vec);

  # Z = matrix(0, nrow=n, ncol=n*q); # Z <- Matrix(Z, sparse = T)

  # -----------------------------------------
  Z = matrix(0, nrow=n, ncol=n*q)
  # for(i in 1:n){
  #   Z[i, seq(1,q) + q*(i-1)] = Z0[i];
  # }
  for(j in 1:q){
    Z[seq_len(n), (seq_len(n)-1) * q + j]  <- diag(Z0)
  }
  Z <- Matrix::Matrix(Z, sparse = T)
  # ----------------------------------------

  if( penalty == "SCAD"){
    # SCAD penalty;

    penalty.fun <- function(eta, lambda, nu, ga){
      # eta_norm = apply( eta, 2, function(x){ sqrt(sum(x^2)) } );
      # eta_norm = col_norms(eta)
      # eta_norm = sqrt( colSums(eta * eta) )
      eta_norm = sqrt( .colSums(eta * eta, m = nrow(eta), n = ncol(eta)) )


      Indicator_1 = eta_norm <= (lambda + lambda/nu);
      # Indicator_2 = eta_norm <= (ga * lambda)  &
      #   eta_norm > (lambda + lambda/nu);
      Indicator_2 = dplyr::between(eta_norm, lambda + lambda/nu, ga * lambda)

      if( sum(Indicator_1) > 0 ){
        val = eta[, Indicator_1]
        eta[, Indicator_1] = sign(val) * matrix( pmax(0, abs(val) -
                                                        lambda/nu), nrow= q);
      }

      if( sum(Indicator_2) > 0 ){
        val = eta[, Indicator_2];
        eta[, Indicator_2] = sign(val) * matrix( pmax(0, abs(val) - ga * lambda /(nu*(ga-1)) ), nrow= q)/(1- 1/((ga-1)*nu));
      }

      return( as.vector(eta) )
    }
  } else{
    # MCP Penalty;

    penalty.fun <- function(eta, lambda, nu, ga)
    {
      # eta_norm = apply( eta, 2, function(x){ sqrt(sum(x^2)) } );
      # eta_norm = col_norms(eta)
      # eta_norm = sqrt( colSums(eta * eta) )
      eta_norm = sqrt( .colSums(eta * eta, m = nrow(eta), n = ncol(eta)) )

      if( sum( eta_norm <= (ga * lambda) ) > 0 )
      {
        val = eta[, eta_norm <= (ga * lambda)];
        eta[, eta_norm <= (ga * lambda) ] = sign(val) * matrix(pmax(0, abs(val) - lambda/nu), nrow= q)/( 1 - 1/(ga * nu) );
      }
      return( as.vector(eta) )
    }

  }

  # ----------------------- #
  #     ADMM Algorithm      #
  # ----------------------- #

  # -----------------------------------------------------------------------------
  # del_mat = big.matrix(nrow=n*(n-1)/2, ncol=n, init = 0)
  I_n = Matrix::Matrix( diag(n), sparse = T)
  del_mat = Matrix::Matrix(0, nrow=n*(n-1)/2, ncol=n, sparse = T)
  # tmp <- expand.grid(1:n, 1:n) %>% dplyr::filter(Var1 < Var2) %>% dplyr::select(Var2, Var1)
  tmp <- expand.grid(1:n, 1:n)
  tmp <- dplyr::filter(tmp, Var1 < Var2)
  tmp <- dplyr::select(tmp, Var2, Var1)

  names(tmp) <- c('j','i')
  del_mat <- Matrix::t(I_n[, tmp$i] - I_n[, tmp$j])

  # -----------------------------------------------------------------------------


  nu_del_del = nu*( n*diag(n*q) - kronecker(rep(1,n), diag(q)) %*% t( kronecker(rep(1,n), diag(q)) ) );

  # Qx = I_n - X %*% solve(t(X) %*% X) %*% t(X);
  Qx = I_n - X %*% solve( crossprod(X)) %*% t(X)

  # D1 = solve( t(Z) %*% Qx %*% Z + nu_del_del );   ## np x nq;
  # D2 = D1 %*% t(Z) %*% Qx %*% Y;
  ## np x 1;

  tmp <- Matrix::crossprod(Z, Qx)
  D1 = Matrix::solve( tmp %*% Z + nu_del_del );
  D2 = D1 %*% tmp %*% Y

  # V0 = solve(t(X) %*% X) %*% t(X);                ## q x n;
  V0 = solve( crossprod(X)) %*% t(X)
  V1 = V0 %*% Y;                                  ## q x 1;
  V2 = V0 %*% Z;                                  ## q x nq;


  # cat('-------- loop begins ------------\n ')

  # --------------  Initial Value Begin  ------------- #

  coeff_lm = as.numeric( lm(Y~X+Z0-1)$coefficients );
  Y_tild = Y - X %*% (coeff_lm[1:p]);
  mu_00 = Matrix::solve( tmp %*% Z + 0.001 * nu_del_del/nu ) %*% tmp %*% Y;


  K00 = floor(sqrt(n));
  mu_med = apply(matrix(mu_00, nrow=q, ncol=n), 2, median);
  q_x = quantile( x = sort(mu_med), prob = seq(1, K00-1)/K00 );

  mu_00 = matrix(0, nrow = n, ncol = q);

  for( i in 1:K00){
    if(i==1){
      ind_coef = mu_med <= q_x[1];
    } else if( i > 1 & i <= (K00-1) ){
      ind_coef = (mu_med <= q_x[i]) & (mu_med > q_x[i-1]);
    } else {
      ind_coef = mu_med > q_x[K00-1];
    }

    if( sum(abs(Z0[ind_coef])) == 0)
    {
      mu_00[ind_coef,] = 0;
    } else
    {
      mu_00[ind_coef,] = lm( Y_tild[ind_coef] ~ Z0[ind_coef] - 1 )$coefficients;
    }
  }

  mu_Mat = - matrix( rep(mu_00, n), nrow=n, ncol=n, byrow = TRUE ) + c(mu_00);
  theta_mu_00 = mu_Mat[upper.tri(mu_Mat)];

  # --------------  Initial Value End  ---------------- #
  # --------------------------------------------------- #
  # ------------  BIC & Estimation Begin  ------------- #

  BIC_lam = K_Est_lam = rep(0, n_lam);
  mu_Est_lam = matrix(0, nrow=n*q, ncol=n_lam);
  beta_Est_lam = matrix(0, nrow=p, ncol=n_lam);
  # theta_Est_lam = matrix(0, nrow=n*(n-1)*q/2, ncol=n_lam);
  # Iter_Conv_lam = matrix(0, nrow=2, ncol=n_lam);

  for(vi in 1:n_lam){
    lambda = lam_vec[vi];

    theta_new = theta_0 = theta_mu_00;
    ups_new = ups_0 = rep(0, n*(n-1)*q/2);
    Converge_ind=0;  iter=1;  error=1;  # maxiter=700;  ep=0.005;

    while( iter <= maxiter & error > ep ){
      # ------------------ #
      # update  mu,  beta  #
      # ------------------ #

      D3 <- nu * Matrix::Matrix(theta_0 - ups_0/nu, nrow=q, ncol=n*(n-1)/2, sparse = T) %*% del_mat

      # Update 'mu'
      mu_new =  D2 + D1 %*% Matrix::t(D3)
      mu_new = mu_new[,1]
      # ---------------- #
      #  update 'theta'  #
      # ---------------- #
      mu_Mat = mu_new - matrix( rep(mu_new, n), nrow=n, ncol=n, byrow = TRUE )
      theta_0 = mu_Mat[upper.tri(mu_Mat)];

      eta = matrix( theta_0 + ups_0/nu, nrow= q, ncol= n*(n-1)/2, byrow=FALSE );
      theta_new = penalty.fun(eta = eta, lambda=lambda, nu=nu, ga=ga);

      # ---------------- #
      # update 'upsilon' #
      # ---------------- #
      ups_new = ups_0 + nu * (theta_0 - theta_new);

      error = sqrt(sum(theta_0 - theta_new)^2);

      # cat('\r lambda = ', lam_vec[vi],' ------ iter = ', iter, '  error = ',error,' -------')

      if( error < ep ){
        Converge_ind = 1;
        break;
      }else{
        theta_0 = theta_new;
        ups_0 = ups_new;
        iter = iter + 1;
      }

    }

    # cat('\n')
    # -- Update 'beta' -- #
    beta_new = (V1 - V2 %*% mu_new )[,1]

    # -- Sort Results -- #

    MU1 = matrix(round(mu_new, digit=0), nrow=q, ncol=n);
    Table = table(MU1);

    beta_Est_lam[,vi] = beta_new;
    mu_Est_lam[,vi] = mu_new;
    # theta_Est_lam[,vi] = theta_new;
    # Iter_Conv_lam[,vi] = c(iter, Converge_ind);

    K_Est_lam[vi] = sum(Table!=0);
    # length( unique( round(mu_new, digits=1) ) );
    # median( apply(MU1, 1, function(x){length(unique(x))}) );

    BIC_lam[vi] = log(  Matrix::mean( ( Y - X %*% beta_new - Z %*% mu_new )^2 ) ) +
      (n^{1/2.7})*log(n*q+p) * log(n)/(n) * (K_Est_lam[vi] * q + p);

    # cat('\n')
  }

  # ------------  BIC & Estimation End  --------------- #
  # --------------------------------------------------- #

  ind_opt = which.min(BIC_lam);
  lam_opt = lam_vec[ind_opt];

  beta_new = beta_Est_lam[, ind_opt];
  mu_new = mu_Est_lam[, ind_opt];

  return(list( BIC = BIC_lam, lam_opt = lam_opt,
               beta_Est_lam = beta_Est_lam, mu_Est_lam = mu_Est_lam,
               mu_est_opt = mu_new, beta_est_opt = beta_new))

}
