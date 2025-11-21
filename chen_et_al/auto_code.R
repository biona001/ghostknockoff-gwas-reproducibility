library("knockoff")
library("expm")
library("CVXR")
library("glmnet")
library("corpcor")
library("MASS")
library("susieR")
library("Rfast")

n = 600
p = 200
n0 = 30
amplitude = 4
rho = 0.5

# standardize X
scale_X <- function(X){
  n = nrow(X)
  X = X - outer(rep(1,n),colMeans(X))
  D = diag(sqrt(n/apply(X,2,function(t) sum(t^2))))
  return (X%*%D)
}

# compute lambda for GK-sqrtlasso
get_lambda_sqrtlasso_1 <- function(X){
  n = nrow(X)
  eps = rnorm(n)
  return (max(abs(t(X)%*%eps))/sqrt(sum(eps^2)))
}

get_lambda_sqrtlasso <- function(X, num=200){
  mean(replicate(num, get_lambda_sqrtlasso_1(X)))
}

# noise std estimation
get_sigma <- function(Z, y_norm, Sigma_inv, n){
  p = nrow(Z)
  return (sqrt(max((p+n+1)*(y_norm^2)/(n^2+n) - t(Z)%*%Sigma_inv%*%Z/(n^2+n),0)))
}

# compute lambda for GK-pseudolasso (lasso-min)
get_lambda_ghost_1 <- function(sigma, chol_Sigma_0, n){
  p = nrow(chol_Sigma_0)
  X = matrix(rnorm(n*p),n,p)%*%chol_Sigma_0
  eps = rnorm(n)
  return (max(abs(t(X)%*%eps))*sigma/n)
}

get_lambda_ghost <- function(sigma, chol_Sigma_0, n, num=200){
  mean(replicate(num, get_lambda_ghost_1(sigma, chol_Sigma_0, n)))
}

# create Gaussian knockoffs 
create_gaussian_s <- function(X,Sigma,s,chol_V){
  n = nrow(X)
  p = ncol(X)
  D = diag(s)
  E = matrix(rnorm(p*n), n, p)
  P = diag(p) - solve(Sigma, D)
  V = 2*D - D%*%solve(Sigma, D)
  return(X%*%P + E%*%chol_V)
}

# GK-sqrtlasso
knockoff_sqrt_lasso <- function(X,X_tilde,y,fdr=0.2,offset=1,kappa=0.3){
  n = nrow(X)
  p = ncol(X)
  beta = Variable(2*p)
  lambda0 = kappa*get_lambda_sqrtlasso(cbind(X,X_tilde))
  obj = cvxr_norm(y - cbind(X,X_tilde)%*%beta, 2) + lambda0*cvxr_norm(beta, 1)
  prob = Problem(Minimize(obj))
  result = solve(prob)
  beta_sqrt_lasso = result$getValue(beta)
  W = abs(beta_sqrt_lasso[1:p]) - abs(beta_sqrt_lasso[(p+1):(2*p)])
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  return(list(statistic=W,threshold=t,selected=selected))
}

# GK-pseudolasso with tuning parameter chosen by lasso-min 
knockoff_ghost_lasso <- function(Z,y_norm,Z_tilde,n,Sigma_0,Sigma_inv,chol_Sigma_0,fdr=0.2,offset=1,kappa=0.6){
  p = nrow(Z)
  sigma = get_sigma(Z, y_norm, Sigma_inv, n)
  beta = Variable(2*p)
  lambda0 = kappa*get_lambda_ghost(sigma=sigma, chol_Sigma_0=chol_Sigma_0, n=n)
  obj = 1/2*quad_form(beta, Sigma_0 + 1e-2*diag(2*p)) - matrix(c(as.vector(Z),as.vector(Z_tilde)),1)%*%beta/n + lambda0*cvxr_norm(beta, 1)
  prob = Problem(Minimize(obj))
  result = solve(prob)
  beta_ghost = result$getValue(beta)
  W = abs(beta_ghost[1:p]) - abs(beta_ghost[(p+1):(2*p)])
  W[abs(W) <= 1e-10] = 0
  return(W)
}

# generate training and validation summary statistics (pseudo-sum)
gen_pseudo_summary_stats <- function(r,r_tilde,chol_Sigma_0,n,train_valid_ratio = 4){
  n_valid = n%/%(1+train_valid_ratio)
  n_train = n - n_valid
  p = length(r)
  R = chol_Sigma_0%*%rnorm(2*p)
  r_comb = c(r, r_tilde)
  r_comb_train = r_comb + sqrt(n_valid/n/n_train)*R
  r_comb_valid = 1/n_valid*(n*r_comb - n_train*r_comb_train)
  return(list(r_comb_train = r_comb_train,  r_comb_valid = r_comb_valid))
}

# get lambda (pseudo-sum)
pseudo_summary_stats_get_lambda <- function(lambda,r,r_tilde,Sigma_0,chol_Sigma_0,n,train_valid_ratio = 4){
  summary_stats_p = gen_pseudo_summary_stats(r,r_tilde,chol_Sigma_0,n,train_valid_ratio = 4)
  r_comb_train = summary_stats_p$r_comb_train
  r_comb_valid = summary_stats_p$r_comb_valid
  pseudo_cor = rep(0, length(lambda))
  for (i in 1:length(lambda)){
    beta = Variable(2*p)
    obj = 1/2*quad_form(beta, Sigma_0 + 1e-2*diag(2*p)) - sum(beta*r_comb_train) + lambda[i]*cvxr_norm(beta, 1)
    prob = Problem(Minimize(obj))
    result = solve(prob)
    beta_pseudo = result$getValue(beta)
    if (max(abs(beta_pseudo)) < 1e-10){
      pseudo_cor[i] = 0
    }
    else{
      pseudo_cor[i] = sum(beta_pseudo*r_comb_valid)/sqrt((t(beta_pseudo)%*%Sigma_0%*%beta_pseudo)[1])
    }
  }
  return(lambda[which.max(pseudo_cor)])
}

# GK-pseudolasso with tuning parameter chosen by pseudo-sum 
knockoff_pseudo_summary_stats <- function(r,r_tilde,Sigma_0,chol_Sigma_0,n){
  p = length(r)
  beta = Variable(2*p)
  r_comb = c(r, r_tilde)
  lambda_max = max(abs(r_comb))
  lambda_seq = exp(seq(log(lambda_max/1000),log(lambda_max),length.out=100))
  lambda_opt = pseudo_summary_stats_get_lambda(lambda=lambda_seq,r,r_tilde,Sigma_0,chol_Sigma_0,n,train_valid_ratio = 4)
  obj = 1/2*quad_form(beta, Sigma_0 + 1e-2*diag(2*p)) - sum(r_comb*beta) + lambda_opt*cvxr_norm(beta, 1)
  prob = Problem(Minimize(obj))
  result = solve(prob)
  beta_ghost = result$getValue(beta)
  W = abs(beta_ghost[1:p]) - abs(beta_ghost[(p+1):(2*p)])
  W[abs(W) <= 1e-10] = 0
  return(W)
}


set.seed(12345)
rand = sample(1:p,n0) 
rand_sign = sample(c(-1,1), size=p, replace=TRUE)

# AR(1) covariates
simulation <- function(rho,k = 200){
  Covariance = toeplitz(rho^(0:(p-1)))
  Sigma = Covariance
  Sigma_inv = solve(Sigma)
  s = create.solve_sdp(Sigma)
  Sigma_0 = matrix(0,2*p,2*p)
  Sigma_0[1:p,1:p] = Sigma
  Sigma_0[1:p,(p+1):(2*p)] = Sigma - diag(s)
  Sigma_0[(p+1):(2*p),1:p] = Sigma - diag(s)
  Sigma_0[(p+1):(2*p),(p+1):(2*p)] = Sigma
  Sigma_0_inv = solve(Sigma_0)
  chol_Sigma_0 = chol(Sigma_0)
  D = diag(s)
  P = diag(p) - solve(Sigma,D)
  V = 2*D - D%*%solve(Sigma,D)
  chol_V = chol(V)
  
  help_fun <- function(rho){
    X = matrix(rnorm(n*p),n,p)%*%chol(Covariance)
    beta = rep(0,p)
    beta[rand] = amplitude
    beta = beta*rand_sign
    y = X%*%beta + sqrt(n)*rnorm(n) # generate y via linear model 
    y = (y - mean(y))/sd(y)
    
    X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
    X_comb = cbind(X,X_tilde)
    
    Z = t(X)%*%y
    y_norm = sqrt(sum(y^2))
    Z_tilde = t(P)%*%Z + y_norm*mvrnorm(1,rep(0,p),V)
    
    # GK-sqrtlasso
    result_sqrtlasso = knockoff_sqrt_lasso(X,X_tilde,y,fdr=0.2)
    power_sqrtlasso = length(intersect(result_sqrtlasso$selected,rand))/n0
    fdp_sqrtlasso = length(setdiff(result_sqrtlasso$selected,rand))/max(length(result_sqrtlasso$selected),1)
    
    # susie_suff_fit = susie_suff_stat(XtX = t(X_comb)%*%X_comb, Xty = t(X_comb)%*%y,yty = sum(y^2),n = n)
    # coef_susie_suff = susie_suff_fit$pip
    # W_susie_suff = abs(coef_susie_suff[1:p]) - abs(coef_susie_suff[(p+1):(2*p)])
    # t_susie_suff = knockoff.threshold(W_susie_suff, fdr=0.2)
    # result_susie_suff = sort(which(W_susie_suff >= t_susie_suff))
    # power_susie_suff = length(intersect(result_susie_suff,rand))/n0
    # fdp_susie_suff = length(setdiff(result_susie_suff,rand))/max(length(result_susie_suff),1)
    
    # GK-lassomax
    W_lassomax = stat.glmnet_lambdadiff(X,X_tilde,y)
    t_lassomax = knockoff.threshold(W_lassomax, fdr=0.2)
    result_lassomax = sort(which(W_lassomax >= t_lassomax))
    power_lassomax = length(intersect(result_lassomax,rand))/n0
    fdp_lassomax = length(setdiff(result_lassomax,rand))/max(length(result_lassomax),1)
    
    # KF-lassocv
    W_cvlasso = stat.lasso_coefdiff(X,X_tilde,y)
    t_cvlasso = knockoff.threshold(W_cvlasso, fdr=0.2)
    result_cvlasso = sort(which(W_cvlasso >= t_cvlasso))
    power_cvlasso = length(intersect(result_cvlasso,rand))/n0
    fdp_cvlasso = length(setdiff(result_cvlasso,rand))/max(length(result_cvlasso),1)
    
    # GK-marginal
    W_mar = abs(Z) - abs(Z_tilde)
    t_mar = knockoff.threshold(W_mar, fdr=0.2)
    result_mar = sort(which(W_mar >= t_mar))
    power_mar = length(intersect(result_mar,rand))/n0
    fdp_mar = length(setdiff(result_mar,rand))/max(length(result_mar),1)
    
    # GK-pseudolasso (lasso-min)
    W_ghost = knockoff_ghost_lasso(Z,y_norm,Z_tilde,n,Sigma_0,Sigma_inv,chol_Sigma_0,fdr=0.2,offset=1)
    t_ghost = knockoff.threshold(W_ghost, fdr=0.2)
    result_ghost = sort(which(W_ghost >= t_ghost))
    power_ghost = length(intersect(result_ghost,rand))/n0
    fdp_ghost = length(setdiff(result_ghost,rand))/max(length(result_ghost),1)
    
    # GK-susie-rss
    susie_fit = susie_rss(z = rbind(Z,Z_tilde)/sqrt(n)/sqrt(mean(y^2)-(rbind(Z,Z_tilde)/n)^2), R=Sigma_0, n=n,var_y=var(y))
    susie_coef = susie_fit$pip
    W_susie = abs(susie_coef[1:p]) - abs(susie_coef[(p+1):(2*p)])
    t_susie = knockoff.threshold(W_susie, fdr=0.2)
    result_susie = sort(which(W_susie >= t_susie))
    power_susie = length(intersect(result_susie,rand))/n0
    fdp_susie = length(setdiff(result_susie,rand))/max(length(result_susie),1)
    
    r = as.vector(Z)/n
    r_tilde = as.vector(Z_tilde)/n
    
    # GK-pseudolasso (pseudo-sum)
    W_pseudo = knockoff_pseudo_summary_stats(r,r_tilde,Sigma_0,chol_Sigma_0,n)
    t_pseudo = knockoff.threshold(W_pseudo, fdr=0.2)
    result_pseudo = sort(which(W_pseudo >= t_pseudo))
    power_pseudo = length(intersect(result_pseudo,rand))/n0
    fdp_pseudo = length(setdiff(result_pseudo,rand))/max(length(result_pseudo),1)
    
    return(c(power_sqrtlasso,fdp_sqrtlasso,power_lassomax,fdp_lassomax,power_cvlasso,fdp_cvlasso,power_mar,fdp_mar,power_ghost,fdp_ghost,power_susie,fdp_susie,power_pseudo,fdp_pseudo))
    
  }
  result = sapply(rep(rho,k),help_fun)
  return(c(apply(result,1,mean), apply(result,1,sd)/sqrt(k)))
}

result = simulation(rho)
print(result)
