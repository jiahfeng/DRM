#######################################################
# Wn_montecarlo.R
# Jiahui Feng
# update Dec. 23, 2023
# 
# Score test for the existence of change point using Monte Carlo method
#
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



## Under the null hypothesis without change point ##
for(N0 in c(50,100)){
  for(N1 in c(50,100)){
    set.seed(11111)
    # inital interval for selecting change point
    space <- seq(qnorm(0.1), qnorm(0.9), length.out = 50)
    len_eta <- length(space)
    Z_norm <- mvrnorm(500000, mu = c(rep(0, len_eta)), Sigma = diag(len_eta))
    
    sim <- c()
    for(i in 1:10000){
      set.seed(i)
      n0 <- N0
      n1 <- N1
      # data generation
      Z0 <- rnorm(n0, 0, 1)
      Z1 <- rnorm(n1, 0, 1)
      dat <- data.frame(Y=c(rep(0,n0),rep(1,n1)),Z=c(Z0,Z1))
      Y <- dat$Y
      Z <- dat$Z
      n <- n0+n1
      
      space_temp <- space[which(space<max(Z))]
      len_eta_temp <- ifelse(length(space)==length(space_temp), len_eta, length(space_temp))
      Z_norm_temp <- if(len_eta==len_eta_temp) Z_norm else mvrnorm(500000, mu = c(rep(0, len_eta_temp)), Sigma = diag(len_eta_temp))
      
      # test statistic given the value of eta
      store <- c()
      for (eta in space_temp) {
        gamma.d <- 0
        alpha.d <- 0
        X_eta <- ifelse(Z-eta>=0,Z-eta,0)
        r_zeta <- exp(gamma.d+alpha.d*X_eta)
        n1e <- n1*r_zeta/(n0+n1*r_zeta)
        n2e <- n0*n1*r_zeta/(n0+n1*r_zeta)^2
        dl_g <- -sum(n1e-Y)
        dl_a <- -sum(n1e*X_eta-Y*X_eta)
        dl2_g <- -sum(n2e)
        dl2_ga <- -sum(n2e*X_eta)
        dl2_a <- -sum(n2e*X_eta^2)
        TMI <- -matrix(c(dl2_g,dl2_ga,dl2_ga,dl2_a),nrow=2,ncol=2) / n
        variance <- 1 / ginv(TMI)[2,2]
        stat_1 <- dl_a / sqrt(n)
        stat_2 <- stat_1^2/variance
        store <- rbind(store, c(stat_1, stat_2, variance))
      }
      
      Wn <- max(abs(store[,1])) #score statistic
      Wn_star_1 <- max(store[,2]) #normalized score statistic
      
      # variance matrix
      Sigma <- matrix(NA, len_eta_temp, len_eta_temp)
      for (j in 1:len_eta_temp) {
        for (k in 1:len_eta_temp) {
          max_data_1 <- ifelse(Z-space_temp[j]>=0,Z-space_temp[j],0)
          max_data_2 <- ifelse(Z-space_temp[k]>=0,Z-space_temp[k],0)
          D_alphaalpha <- mean(max_data_1 * max_data_2)
          D_gammaalpha_1 <- mean(max_data_1)
          D_gammaalpha_2 <- mean(max_data_2)
          Sigma[j,k] <- (n0/n) * (n1/n) * (D_alphaalpha - D_gammaalpha_1 * D_gammaalpha_2)
        }
      }
      
      
      # maximal score test
      temp_Wn <- eigen(Sigma)
      eigen_value_Wn <- temp_Wn$values
      eigen_value_Wn[which(eigen_value_Wn<0)] <- 0
      eigen_matrix_Wn <- diag(eigen_value_Wn)
      eigen_vector_Wn <- temp_Wn$vectors
      sd_Mat_Wn <- eigen_vector_Wn %*% sqrt(eigen_matrix_Wn) %*% t(eigen_vector_Wn)
      z_Wn <- Z_norm_temp %*% sd_Mat_Wn
      z_max_Wn <- apply(abs(z_Wn), 1, max)
      pv_max_Wn <- mean(z_max_Wn > Wn) 
      
      
      # maximal normalized score test
      Balance <- diag(1/sqrt(diag(Sigma)))
      Cov_Mat_Wn_star <- Balance %*% Sigma %*% Balance
      temp_Wn_star <- eigen(Cov_Mat_Wn_star)
      eigen_value_Wn_star <- temp_Wn_star$values
      eigen_value_Wn_star[which(eigen_value_Wn_star<0)] <- 0
      eigen_matrix_Wn_star <- diag(eigen_value_Wn_star)
      eigen_vector_Wn_star <- temp_Wn_star$vectors
      sd_Mat_Wn_star <- eigen_vector_Wn_star %*% sqrt(eigen_matrix_Wn_star) %*% t(eigen_vector_Wn_star)
      z_Wn_star <- Z_norm_temp %*% sd_Mat_Wn_star
      z_max_Wn_star <- apply(abs(z_Wn_star), 1, max)
      pv_max_Wn_star_1 <- mean(z_max_Wn_star > sqrt(Wn_star_1))
  
      sim <- rbind(sim, c(Wn, Wn_star_1, pv_max_Wn, pv_max_Wn_star_1))
    }
    sim.df<-data.frame(sim)
    
    names(sim.df)<-c("score_1","score_2","pv_1","pv_2")
    write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_H0_montecarlo.csv"))
  }
}



## Under the alternative hypothesis with change point ##
for(N0 in c(50,100)){
  for(N1 in c(50,100)){
    for(teta in c(-0.5,0,0.5)){
      set.seed(11111)
      n0 <- N0
      n1 <- N1
      alpha <- 0.5
      true.eta <- teta
      gamma <- -log(exp(0.5*alpha^2-alpha*true.eta)*(1-pnorm(true.eta-alpha))+pnorm(true.eta))
      # generate data to set an empirical interval for change point candidates
      pU <- pnorm(q = true.eta)/((exp(0.5*alpha^2-alpha*true.eta))*(1-pnorm(true.eta-alpha))+pnorm(q = true.eta))
      temp_Z0 <- rnorm(n0*1000)
      temp_U <- rbinom(n1*1000, size = 1, prob = pU)
      temp_Z1 <- ifelse(temp_U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                        rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
      temp_Z <- c(temp_Z0, temp_Z1)
      space <- seq(quantile(temp_Z,0.1), quantile(temp_Z, 0.9), length.out=50)
      len_eta <- length(space)
      Z_norm <- mvrnorm(500000, mu = c(rep(0, len_eta)), Sigma = diag(len_eta))
      
      sim <- c()
      for(i in 1:10000){
        set.seed(i)
        # data generation
        Z0 <- rnorm(n0, 0, 1)
        U <- rbinom(n1, size = 1,prob = pU)
        Z1 <- ifelse(U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                      rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
        
        dat <- data.frame(Y=c(rep(0,n0),rep(1,n1)),Z=c(Z0,Z1))
        Y <- dat$Y
        Z <- dat$Z
        n <- n0+n1
        
        space_temp <- space[which(space<max(Z))]
        len_eta_temp <- ifelse(length(space)==length(space_temp), len_eta, length(space_temp))
        Z_norm_temp <- if(len_eta==len_eta_temp) Z_norm else mvrnorm(500000, mu = c(rep(0, len_eta_temp)), Sigma = diag(len_eta_temp))
        
        # test statistic
        store <- c()
        for (eta in space_temp) {
          gamma.d <- 0
          alpha.d <- 0
          X_eta <- ifelse(Z-eta>=0,Z-eta,0)
          r_zeta <- exp(gamma.d+alpha.d*X_eta)
          n1e <- n1*r_zeta/(n0+n1*r_zeta)
          n2e <- n0*n1*r_zeta/(n0+n1*r_zeta)^2
          dl_g <- -sum(n1e-Y)
          dl_a <- -sum(n1e*X_eta-Y*X_eta)
          dl2_g <- -sum(n2e)
          dl2_ga <- -sum(n2e*X_eta)
          dl2_a <- -sum(n2e*X_eta^2)
          TMI <- -matrix(c(dl2_g,dl2_ga,dl2_ga,dl2_a),nrow=2,ncol=2) / n
          variance <- 1 / ginv(TMI)[2,2]
          stat_1 <- dl_a / sqrt(n)
          stat_2 <- stat_1^2/variance
          store <- rbind(store, c(stat_1, stat_2, variance))
        }
        
        Wn <- max(abs(store[,1])) #score statistic
        Wn_star_1 <- max(store[,2]) #normalized score statistic
        
        Sigma <- matrix(NA, len_eta_temp, len_eta_temp)
        for (j in 1:len_eta_temp) {
          for (k in 1:len_eta_temp) {
            max_data_1 <- ifelse(Z-space_temp[j]>=0,Z-space_temp[j],0)
            max_data_2 <- ifelse(Z-space_temp[k]>=0,Z-space_temp[k],0)
            D_alphaalpha <- mean(max_data_1 * max_data_2)
            D_gammaalpha_1 <- mean(max_data_1)
            D_gammaalpha_2 <- mean(max_data_2)
            Sigma[j,k] <- (n0/n) * (n1/n) * (D_alphaalpha - D_gammaalpha_1 * D_gammaalpha_2)
          }
        }
        
        # maximal score test 
        temp_Wn <- eigen(Sigma)
        eigen_value_Wn <- temp_Wn$values
        eigen_value_Wn[which(eigen_value_Wn<0)] <- 0
        eigen_matrix_Wn <- diag(eigen_value_Wn)
        eigen_vector_Wn <- temp_Wn$vectors
        sd_Mat_Wn <- eigen_vector_Wn %*% sqrt(eigen_matrix_Wn) %*% t(eigen_vector_Wn)
        z_Wn <- Z_norm_temp %*% sd_Mat_Wn
        z_max_Wn <- apply(abs(z_Wn), 1, max)
        pv_max_Wn <- mean(z_max_Wn > Wn)
        
        # maximal normalized score test 
        Balance <- diag(1/sqrt(diag(Sigma)))
        Cov_Mat_Wn_star <- Balance %*% Sigma %*% Balance
        temp_Wn_star <- eigen(Cov_Mat_Wn_star)
        eigen_value_Wn_star <- temp_Wn_star$values
        eigen_value_Wn_star[which(eigen_value_Wn_star<0)] <- 0
        eigen_matrix_Wn_star <- diag(eigen_value_Wn_star)
        eigen_vector_Wn_star <- temp_Wn_star$vectors
        sd_Mat_Wn_star <- eigen_vector_Wn_star %*% sqrt(eigen_matrix_Wn_star) %*% t(eigen_vector_Wn_star)
        z_Wn_star <- Z_norm_temp %*% sd_Mat_Wn_star
        z_max_Wn_star <- apply(abs(z_Wn_star), 1, max)
        pv_max_Wn_star_1 <- mean(z_max_Wn_star > sqrt(Wn_star_1))
        
        
        sim <- rbind(sim, c(Wn, Wn_star_1, pv_max_Wn, pv_max_Wn_star_1))
      }
      sim.df <- data.frame(sim)
      
      names(sim.df) <- c("score_1","score_2","pv_1","pv_2")
      write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_eta=",teta,"_H1_montecarlo.csv"))
    }
  }
}
