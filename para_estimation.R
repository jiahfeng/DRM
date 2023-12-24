#######################################################
# para_estimation.R
# Jiahui Feng
# update Dec. 23, 2023
# 
# Simulation study of parameter estimation under the density ratio
# change point model
#
########################################################


library(MASS)       

# Generates random numbers from a truncated normal distribution
rtruncnorm <- function(n, mu, sigma, low, high) {
  # n: sample size
  # mu: mean
  # sigma: standard deviation
  # low: lower bound
  # high: higher bound
  p_low <- pnorm(low, mu, sigma)
  p_high <- pnorm(high, mu, sigma)
  qnorm(runif(n, p_low, p_high), mu, sigma)
}


# likelihood function of density-ration model
lik <- function(par, eta, Y, Z){
  # par: regression parameter 
  # eta: change point
  # Y: class indicator
  # Z: observed covariate
  gamma <- par[1]
  alpha <- par[2]
  -sum(log(n0 + n1 * exp(gamma + alpha * ifelse(Z-eta>=0, Z-eta, 0)))) + sum(Y * (gamma + alpha * ifelse(Z-eta>=0, Z-eta, 0)))
}



# Generates random numbers from a truncated normal distribution
rtruncnorm <- function(n, mu, sigma, low, high) {
  # n: sample size
  # mu: mean
  # sigma: standard deviation
  # low: lower bound
  # high: higher bound
  p_low <- pnorm(low, mu, sigma)
  p_high <- pnorm(high, mu, sigma)
  qnorm(runif(n, p_low, p_high), mu, sigma)
}


# likelihood function of density-ration model
lik <- function(par, eta, Y, Z){
  # par: regression parameter 
  # eta: change point
  # Y: class indicator
  # Z: observed covariate
  gamma <- par[1]
  alpha <- par[2]
  -sum(log(n0 + n1 * exp(gamma + alpha * ifelse(Z-eta>=0, Z-eta, 0)))) + sum(Y * (gamma + alpha * ifelse(Z-eta>=0, Z-eta, 0)))
}



# N0: sample size for group 0
# N1: sample size for group 1
# teta: value of change point
for(N0 in c(500,1000)){
  for(N1 in c(500,1000)){
    for(teta in c(-0.5,0,0.5)){
      set.seed(11111)
      # simulation setting
      n0 <- N0
      n1 <- N1
      alpha <- 1
      true.eta <- teta
      gamma <- -log(exp(0.5*alpha^2-alpha*true.eta)*(1-pnorm(true.eta-alpha))+pnorm(true.eta))
      # generate data to set an empirical interval for change point candidates
      pU <-pnorm(q = true.eta)/((exp(0.5*alpha^2-alpha*true.eta))*(1-pnorm(true.eta-alpha))+pnorm(q = true.eta))
      temp_Z0 <- rnorm(n0*1000)
      temp_U <- rbinom(n1*1000, size = 1, prob = pU)
      temp_Z1 <- ifelse(temp_U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                        rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
      temp_Z <- c(temp_Z0, temp_Z1)
      space_all <- seq(quantile(temp_Z,0.1), quantile(temp_Z, 0.9), length.out=50)
      
      # main part of simulation
      sim <- c()
      for(i in 1:10000){
        set.seed(i)
        # data generation 
        Z0 <- rnorm(n0, 0, 1)
        U <-rbinom(n1, size = 1,prob = pU)
        Z1 <- ifelse(U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                     rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
        
        dat <- data.frame(Y=c(rep(0,n0),rep(1,n1)),Z=c(Z0,Z1))
        Y <- dat$Y
        Z <- dat$Z
        n <- n0+n1
        
        # parameter estimation given the value of eta
        space <- space_all[which(space_all < max(Z))]
        store <- c()
        for(eta in space){
          gamma.d <- 0
          alpha.d <- 0
          S <- 0
          repeat{
            S <- S+1
            gamma.store <- gamma.d; alpha.store <- alpha.d; 
            X_eta <- ifelse(Z-eta>=0,Z-eta,0)
            r_zeta <- exp(gamma.d+alpha.d*X_eta)
            n1e <- n1*r_zeta/(n0+n1*r_zeta)
            n2e <- n0*n1*r_zeta/(n0+n1*r_zeta)^2
            dl_g <- -sum(n1e-Y)
            dl_a <- -sum(n1e*X_eta-Y*X_eta)
            dl2_g <- -sum(n2e)
            dl2_ga <- -sum(n2e*X_eta)
            dl2_a <- -sum(n2e*X_eta^2)
            
            update <- c(gamma.d, alpha.d)-ginv(matrix(c(dl2_g,dl2_ga,dl2_ga,dl2_a),nrow=2,ncol=2))%*%c(dl_g,dl_a)
            gamma.d <- update[1]
            alpha.d <- max(-5,min(5,update[2]))
            dist <- max(abs(c(gamma.d-gamma.store,alpha.d-alpha.store)))
            if(dist<10^{-4}|S==50){break}
          }
          store <- rbind(store,c(eta,gamma.d,alpha.d,lik(par = c(gamma.d,alpha.d),eta = eta,Y = Y,Z = Z)))
        }
        
        op <- which.max(store[,4])
        eta.star <- store[op,1]
        gamma.star <- store[op,2]
        alpha.star <- store[op,3]
        
        X_eta.star <- ifelse(Z-eta.star>=0,Z-eta.star,0)
        r_zeta.star <- exp(gamma.star+alpha.star*ifelse(Z-eta.star>=0,Z-eta.star,0))
        n2e.star <-  n0*n1*r_zeta.star/(n0+n1*r_zeta.star)^2
        dl2_g.star <- -sum(n2e.star)
        dl2_ga.star <- -sum(n2e.star*X_eta.star)
        dl2_a.star <- -sum(n2e.star*X_eta.star^2)
        
        V2_g.star <- -dl2_g.star/n - dl2_g.star^2/(n0*n1)
        V2_ga.star <- -dl2_ga.star/n - dl2_g.star*dl2_ga.star/(n0*n1)
        V2_a.star <- -dl2_a.star/n - dl2_ga.star^2/(n0*n1)
        V <- matrix(c(V2_g.star, V2_ga.star, V2_ga.star, V2_a.star), nrow = 2, ncol = 2)
        
        D <- matrix(c(dl2_g.star,dl2_ga.star,dl2_ga.star,dl2_a.star),nrow=2,ncol=2) / n
        
        var_est <- ginv(D)%*%V%*%ginv(D)/n
        
        sim <- rbind(sim, c(eta.star-true.eta, gamma.star-gamma, alpha.star-alpha, var_est[1,1], var_est[2,2]))
      }
      sim.df <- data.frame(sim) 
      
      names(sim.df) <- c("eta", "gamma", "alpha", "var_gamma", "var_alpha")
      write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_eta=",teta,"_estimation.csv"))
    }
  }
}
