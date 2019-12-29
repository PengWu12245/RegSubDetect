
#' Subgroup detection in setting of regression with ADMM algorithm
#'
#' The model is \deqn{  Y_i = X_i \beta + Z_i \theta_i.} This ADMM algorithm
#' can detect the subgroup effects of \eqn{\theta_i} and estimate the regression
#' coefficient \eqn{\beta} simultaneously.
#'
#'@param X a n*p matrix, \code{p}-dimensional covariate.
#'@param Y a n*1 vector, response.
#'@param Z0 a n*1 matrix, treatment indicator.
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
#'@seealso \code{\link{ADMM}}
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

ADMM.q1 <- function(Y, X, Z0,  penalty, lam_vec,
                    ga = 3.7, nu = 1, maxiter=1000, ep=5e-4 ){
  # Y: n x 1;
  # X: n x p;
  # Z0: n x 1;
  # Z: n x n;  = diag(Z)
  # model: Y = X beta + Z theta    note: here is Z, not Z0
  # nu: tuning parameter for augmented Largrange function

  n = nrow(X);
  p = ncol(X);            # Dimension of Covariates.
  n_lam = length(lam_vec);
  Z = diag(Z0)

  # ---------- 1. define updating formula of penalty functions -------------
  if( penalty == "SCAD"){
    # SCAD penalty;

    penalty.fun <- function(eta, lambda, nu, ga){

      eta_norm = abs(eta)
      Indicator_1 = eta_norm <= (lambda + lambda/nu)
      Indicator_2 = dplyr::between(eta_norm, lambda + lambda/nu, ga * lambda)

      if( sum(Indicator_1) > 0 ){
        val = eta[Indicator_1]
        eta[Indicator_1] = sign(val) * pmax(0, abs(val) -lambda/nu)
      }

      if( sum(Indicator_2) > 0 ){
        val = eta[Indicator_2];
        eta[Indicator_2] = sign(val) * pmax(0, abs(val) - ga * lambda /(nu*(ga-1)) )/(1- 1/((ga-1)*nu));
      }

      return( eta )
    }
  } else{
    # MCP Penalty;

    penalty.fun <- function(eta, lambda, nu, ga)
    {
      eta_norm = abs(eta)
      ind = eta_norm <= (ga * lambda)
      if( sum( ind ) > 0 )
      {
        val = eta[ind];
        eta[ind] = sign(val) * pmax(0, abs(val) - lambda/nu)/( 1 - 1/(ga * nu) );
      }
      return( eta )
    }

  }


  # -------------- 2. preparations for ADMM algorithm  -----------------

  I_n = Matrix::Matrix( diag(n), sparse = T)
  A = Matrix::Matrix(0, nrow=n*(n-1)/2, ncol=n, sparse = T)
  tmp <- expand.grid(1:n, 1:n)
  tmp <- dplyr::filter(tmp, Var1 < Var2)
  tmp <- dplyr::select(tmp, Var2, Var1)
  names(tmp) <- c('j','i')
  A <- Matrix::t(I_n[, tmp$i] - I_n[, tmp$j])
  # A also will be used in updating theta

  nu_A_A = nu*( n*diag(n) - kronecker(rep(1,n), diag(1)) %*%
                  t( kronecker(rep(1,n), diag(1)) ) );

  D = I_n - X %*% solve( crossprod(X)) %*% t(X)
  tmp <- Matrix::crossprod(Z, D)
  D1 = Matrix::solve( tmp %*% Z + nu_A_A );
  D2 = D1 %*% tmp %*% Y
  # D1, D2 also will be used in updating theta

  V0 = solve( crossprod(X)) %*% t(X)
  V1 = V0 %*% Y                                 ## p x 1;
  V2 = V0 %*% Z                                 ## p x n;
  # V1, V2 will be used in updating beta


  # ----------------- 3. ADMM algorithm  ----------------------

  # ---------- 3.1 Initial Value of eta and ups   -------------

  coeff_lm = as.numeric( lm(Y~X+Z0-1)$coefficients );
  Y_tild = Y - X %*% (coeff_lm[1:p]);
  theta_0 = Matrix::solve( tmp %*% Z + 0.001 * nu_A_A/nu ) %*% tmp %*% Y;

  K00 = floor(sqrt(n));
  theta_med = as.numeric(theta_0) # apply(matrix(theta_0, nrow=q, ncol=n), 2, median);
  q_x = quantile( x = theta_med, prob = seq(1, K00-1)/K00 );

  theta_0 = matrix(0, nrow = n, ncol = 1);

  for( i in 1:K00){
    if(i==1){
      ind_coef = theta_med <= q_x[1];
    } else if( i > 1 & i <= (K00-1) ){
      ind_coef = (theta_med <= q_x[i]) & (theta_med > q_x[i-1]);
    } else {
      ind_coef = theta_med > q_x[K00-1];
    }

    if( sum(abs(Z0[ind_coef])) == 0)
    {
      theta_0[ind_coef,] = 0;
    } else
    {
      theta_0[ind_coef,] = lm( Y_tild[ind_coef] ~ Z0[ind_coef] - 1 )$coefficients;
    }
  }

  eta_0 = as.numeric(A %*% theta_0)
  ups_0 = rep(0, n*(n-1)/2)


  # ------------   Estimation Begin  ------------- #

  BIC_lam = K_Est_lam = rep(0, n_lam);
  theta_Est_lam = matrix(0, nrow=n, ncol=n_lam);
  beta_Est_lam = matrix(0, nrow=p, ncol=n_lam);

  for(vi in 1:n_lam){
    lambda = lam_vec[vi];

    eta_last = eta_0#= as.numeric(A %*% theta_0)
    ups_last = ups_0#= rep(0, n*(n-1)/2)
    Converge_ind=0;  iter=1;  error=1;  # maxiter=700;  ep=0.005;

    while( iter <= maxiter & error > ep ){

      # ---------------  3.2 Update 'theta'  --------------------

      D3 <- nu * Matrix::Matrix(eta_last - ups_last/nu, nrow=1, ncol=n*(n-1)/2,
                                sparse = T) %*% A
      theta_new =  D2 + D1 %*% Matrix::t(D3)
      theta_new = theta_new[,1]

      # ---------------- 3.3 Updata 'eta'   ---------------------
      eta_tmp = as.numeric(A %*% theta_new)
      eta = eta_tmp + ups_last/nu
      eta_new = penalty.fun(eta = eta, lambda=lambda, nu=nu, ga=ga);

      # --------- 3.4 update lagrange multipliers 'ups' ---------
      tmp = eta_tmp - eta_new
      ups_new = ups_last + nu * tmp;
      error = sqrt(sum(tmp)^2);
      # error = sqrt(sum(tmp*tmp))
      # cat('\r lambda = ', lam_vec[vi],' ------ iter = ', iter, '  error = ',error,' -------')

      if( error < ep ){
        Converge_ind = 1;
        break;
      }else{
        eta_last = eta_new;
        ups_last = ups_new;
        iter = iter + 1;
      }

    }

    # ------------- 3.5 Update 'beta' --------------------------
    beta_new = (V1 - V2 %*% theta_new )[,1]

    # ------------- 3.6 Sort Results --------------------------

    MU1 = matrix(round(theta_new, digit=0), nrow=1, ncol=n);
    Table = table(MU1);

    beta_Est_lam[,vi] = beta_new;
    theta_Est_lam[,vi] = theta_new;

    K_Est_lam[vi] = sum(Table!=0);

    BIC_lam[vi] = log(  Matrix::mean( ( Y - X %*% beta_new - Z %*% theta_new )^2 ) ) +
      (n^{1/2.7})*log(n+p) * log(n)/(n) * (K_Est_lam[vi]  + p);

    # cat('\n')
  }

  # ------------  4. BIC & Estimation End  ---------------

  ind_opt = which.min(BIC_lam);
  lam_opt = lam_vec[ind_opt];

  beta_new = beta_Est_lam[, ind_opt];
  theta_new = theta_Est_lam[, ind_opt];

  return(list( BIC = BIC_lam, lam_opt = lam_opt,
               beta_Est_lam = beta_Est_lam, theta_Est_lam = theta_Est_lam,
               theta_est_opt = theta_new, beta_est_opt = beta_new))

}



