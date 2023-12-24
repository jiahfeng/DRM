#######################################################
# Wn_bootstrap.R
# Jiahui Feng
# update Dec. 23, 2023
# 
# Score test for the existence of change point using bootstrap method
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
    sim <- c()
    for (i in 1:10000) {
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
      
      # inital interval for selecting change point
      space <- seq(qnorm(0.1), qnorm(0.9), length.out = 50)
      
      # test statistic
      store <- c()
      for (eta in space) {
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
      
      
      # bootstarp
      Tn.bts.1 <- c()
      Tn.bts.2 <- c()
      B <- 1000
      for(j in c(1:B)){
        set.seed(j)
        boot.Z <- sample(x = Z,size = n,replace = TRUE) # bootstrap
        
        store.bts <- c()
        for (eta in space) {
          gamma.d <- 0
          alpha.d <- 0
          X_eta <- ifelse(boot.Z-eta>=0,boot.Z-eta,0)
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
          store.bts <- rbind(store.bts, c(stat_1, stat_2, variance))
        }
        
        bts.Wn <- max(abs(store.bts[,1])) #score statistic
        bts.Wn_star_1 <- max(store.bts[,2]) #normalized score statistic
        
        Tn.bts.1 <- c(Tn.bts.1, bts.Wn)
        Tn.bts.2 <- c(Tn.bts.2, bts.Wn_star_1)
      }
      
      sim <- rbind(sim, c(Wn, Wn_star_1, quantile(Tn.bts.1,c(0.9, 0.95, 0.99)), ifelse(Wn>quantile(Tn.bts.1,c(0.9, 0.95, 0.99)), 1, 0),
                          quantile(Tn.bts.2,c(0.9, 0.95, 0.99)), ifelse(Wn_star_1>quantile(Tn.bts.2,c(0.9, 0.95, 0.99)), 1, 0)))
    }
    sim.df <- data.frame(sim)
    
    names(sim.df) <- c("score_1","score_2","Tn.bts90","Tn.bts95","Tn.bts99","Tn.rej90","Tn.rej95","Tn.rej99",
                     "Tn.bts.star90","Tn.bts.star95","Tn.bts.star99","Tn.rej.star90","Tn.rej.star95","Tn.rej.star99")
    write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_H0_bootstrap.csv"))
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
      pU <-pnorm(q = true.eta)/((exp(0.5*alpha^2-alpha*true.eta))*(1-pnorm(true.eta-alpha))+pnorm(q = true.eta))
      temp_Z0 <- rnorm(n0*1000)
      temp_U <- rbinom(n1*1000, size = 1, prob = pU)
      temp_Z1 <- ifelse(temp_U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                        rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
      temp_Z <- c(temp_Z0, temp_Z1)
      space <- seq(quantile(temp_Z,0.1), quantile(temp_Z, 0.9), length.out=50)
      
      sim <- c()
      for (i in 1:10000) {
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
        
        # test statistic
        store <- c()
        for (eta in space) {
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
        
        # bootstrap
        Tn.bts.1 <- c()
        Tn.bts.2 <- c()
        B <- 1000
        for(j in c(1:B)){
          set.seed(j)
          boot.Z <- sample(x = Z,size = n,replace = TRUE)
          
          store.bts <- c()
          for (eta in space) {
            gamma.d <- 0
            alpha.d <- 0
            X_eta <- ifelse(boot.Z-eta>=0,boot.Z-eta,0)
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
            store.bts <- rbind(store.bts, c(stat_1, stat_2, variance))
          }
          
          bts.Wn <- max(abs(store.bts[,1])) #score statistic
          bts.Wn_star_1 <- max(store.bts[,2]) #normalized score statistic
          
          Tn.bts.1 <- c(Tn.bts.1, bts.Wn)
          Tn.bts.2 <- c(Tn.bts.2, bts.Wn_star_1)
        }
        
        sim <- rbind(sim, c(Wn, Wn_star_1, quantile(Tn.bts.1,c(0.9, 0.95, 0.99)), ifelse(Wn>quantile(Tn.bts.1,c(0.9, 0.95, 0.99)), 1, 0),
                            quantile(Tn.bts.2,c(0.9, 0.95, 0.99)), ifelse(Wn_star_1>quantile(Tn.bts.2,c(0.9, 0.95, 0.99)), 1, 0)))
      }
      sim.df <- data.frame(sim)
      
      names(sim.df) <- c("score_1","score_2","Tn.bts90","Tn.bts95","Tn.bts99","Tn.rej90","Tn.rej95","Tn.rej99",
                       "Tn.bts.star90","Tn.bts.star95","Tn.bts.star99","Tn.rej.star90","Tn.rej.star95","Tn.rej.star99")
      write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_eta=",teta,"_H1_bootstrap.csv"))
    }
  }
}


