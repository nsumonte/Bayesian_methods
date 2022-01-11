library(invgamma)
library(MASS)
library(mvtnorm)
### funcion que calcula la merginal posteriori
marginal_posteriori <- function(X,y,b_0,B_0,sigma_0,v_0,vect_gamma){
  set.seed(1)
  sigma_gen_priori <- rinvchisq(1,v_0,sigma_0)
  B_0_gamma <- t(B_0)[vect_gamma == 1, DROP = TRUE]
  B_0_gamma1 <- t(B_0_gamma)[vect_gamma == 1, DROP = TRUE]
  b_0_gamma <- matrix(b_0)[vect_gamma == 1, DROP = TRUE]
  X_gamma <- t(X)[vect_gamma == 1, DROP = TRUE]
  X_gamma1 <- t(X_gamma)
  sigma <- sigma_gen_priori*B_0_gamma1
  betas_priori <- mvrnorm(1,b_0_gamma,sigma)
  y_posteriori <- mvrnorm(1,t(X_gamma)%*%(betas_priori), sigma_gen_priori*diag(nrow(t(X_gamma))))
  n <- dim(X_gamma1)[1]
  B_n <- solve(t(X_gamma1)%*%X_gamma1 + solve(B_0_gamma1))
  b_n <- solve(t(X_gamma1)%*%X_gamma1 + solve(B_0_gamma1))%*%(t(X_gamma1)%*%y + solve(B_0_gamma1)%*%b_0_gamma)
  v_n <- v_0 + n
  sigma_n <- (v_0*sigma_0 + t(b_0_gamma)%*%solve(B_0_gamma1)%*%b_0_gamma + t(y)%*%y - t(b_n)%*%solve(B_n)%*%b_n)/v_n
  sigma_gen_posteriori <- rinvchisq(1,v_n,sigma_n)
  betas_posteriori <- mvrnorm(1, b_n, c(sigma_n)*B_n)
  densidades_beta_prioris <- dmvnorm(betas_priori, b_0_gamma, sigma)
  densidades_beta_posterioris <- dmvnorm(betas_posteriori, b_n, c(sigma_n)*B_n)
  densidades_y_posterioris <- dmvnorm(y_posteriori, t(X_gamma)%*%(betas_priori) ,sigma_gen_priori*diag(nrow(t(X_gamma))))
  densidad_sigma_priori <- dinvchisq(sigma_gen_priori,v_0,sigma_0)
  densidad_sigma_posteriori <- dinvchisq(sigma_gen_posteriori,v_n,sigma_n)
  y_marginal <- log(densidades_y_posterioris) + log(densidades_beta_prioris*densidad_sigma_priori) - log(densidad_sigma_posteriori*densidades_beta_posterioris)
  return(y_marginal)
}

### datos de swiss
y = swiss[,1]
agricultura = swiss[,2]
examen = swiss[,3]
educacion = swiss[,4]
catolico = swiss[,5]
mortalidad_infantil = swiss[,6]
X = cbind(agricultura, examen, educacion , catolico, mortalidad_infantil)
sigma_0 = try(summary(lm(y~1+X))$sigma^2,silent=TRUE)
B_0 = solve(t(X)%*%X)*sigma_0
b_0 = solve(t(X)%*%X)%*%t(X)%*%y
v_0 = 1
### simulamos las marginales para datos pedidos.
vect_gamma = c(1,1,1,1,1)
arg <- marginal_posteriori(X,y,b_0,B_0,sigma_0,v_0,vect_gamma)
arg
vect_gamma1 = c(1,0,1,1,1)
arg1 <- marginal_posteriori(X,y,b_0,B_0,sigma_0,v_0,vect_gamma1)
arg1
vect_gamma2 = c(0,0,1,1,1)
arg2 <- marginal_posteriori(X,y,b_0,B_0,sigma_0,v_0,vect_gamma2)
arg2
vect_gamma3 = c(0,1,1,1,1)
arg3 <- marginal_posteriori(X,y,b_0,B_0,sigma_0,v_0,vect_gamma3)
arg3

