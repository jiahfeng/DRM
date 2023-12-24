#######################################################
# CI_changepoint.R
# Jiahui Feng
# update Dec. 23, 2023
# 
# simulation study to illustrate the confidence interval 
# of the change point using m-out-of-n bootstrap
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



# N0: sample size for group 0
# N1: sample size for group 1
# teta: value of change point
for(N0 in c(500,1000)){
  for(N1 in c(500,1000)){
    for(teta in c(-0.5,0,0.5)){
      set.seed(99999)
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
      space <- seq(quantile(temp_Z, 0.1), quantile(temp_Z, 0.9), length.out=50)
      
      # main part of simulation
      sim <- c()
      for (i in 1:10000) {
        set.seed(i)
        # data generation
        Z0 <- rnorm(n0, 0, 1)
        U <- rbinom(n1, size = 1,prob = pU)
        Z1 <- ifelse(U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                      rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
        dat <- data.frame(Y = c(rep(0, n0), rep(1, n1)), Z = c(Z0, Z1))
        Y <- dat$Y
        Z <- dat$Z
        n <- n0 + n1
        
        # parameter estimation given the value of eta
        store <- c()
        for(eta in space){
          gamma.d <- 0
          alpha.d <- 0
          S <- 0
          repeat{
            S <- S + 1
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
          store <- rbind(store, c(eta,gamma.d,alpha.d,lik(par = c(gamma.d,alpha.d),eta = eta,Y = Y,Z = Z)))
        }
        
        # select the set of parameters with maximum likelihood function value
        op <- which.max(store[,4])
        eta.star <- store[op,1]
        gamma.star <- store[op,2]
        alpha.star <- store[op,3]
        r_zeta.star <- exp(gamma.star+alpha.star*ifelse(Z-eta.star>=0,Z-eta.star,0))
        
        
        # m-out-of-n bootstrap
        B <- 200
        r <- 0.75
        Ftilde.bts <- c()
        result.bts <- list()
        for (j in 1:5) {
          m <- as.integer(n * r^j)
          eta.bts <- c()
          for (k in 1:B) {
            set.seed(i+k)
            dat.bts <- dat[sample(n, size = m, replace = TRUE),] # bootstap sample
            Z0.bts <- dat.bts[which(dat.bts$Y==0), 2]
            Z1.bts <- dat.bts[which(dat.bts$Y==1), 2]
            
            n0.bts <- length(Z0.bts)
            n1.bts <- length(Z1.bts)
            dat.bts <- data.frame(Y=c(rep(0, n0.bts),rep(1, n1.bts)),Z=c(Z0.bts,Z1.bts))
            Y.bts <- dat.bts$Y
            Z.bts <- dat.bts$Z
            
            store.bts <- c()
            for(eta in space){
              gamma.d <- 0
              alpha.d <- 0
              S <- 0
              repeat{
                S <- S+1
                gamma.store <- gamma.d; alpha.store <- alpha.d; 
                X_eta <- ifelse(Z.bts-eta>=0,Z.bts-eta,0)
                r_zeta <- exp(gamma.d+alpha.d*X_eta)
                n1e <- n1.bts*r_zeta/(n0.bts+n1.bts*r_zeta)
                n2e <- n0.bts*n1.bts*r_zeta/(n0.bts+n1.bts*r_zeta)^2
                dl_g <- -sum(n1e-Y.bts)
                dl_a <- -sum(n1e*X_eta-Y.bts*X_eta)
                dl2_g <- -sum(n2e)
                dl2_ga <- -sum(n2e*X_eta)
                dl2_a <- -sum(n2e*X_eta^2)
                
                update <- c(gamma.d, alpha.d)-ginv(matrix(c(dl2_g,dl2_ga,dl2_ga,dl2_a),nrow=2,ncol=2))%*%c(dl_g,dl_a)
                gamma.d <- update[1]
                alpha.d <- max(-5,min(5,update[2]))
                dist <- max(abs(c(gamma.d-gamma.store,alpha.d-alpha.store)))
                if(dist<10^{-4}|S==50){break}
              }
              store.bts <- rbind(store.bts,c(eta,gamma.d,alpha.d,lik(par = c(gamma.d,alpha.d),eta = eta,Y = Y.bts,Z = Z.bts)))
            }
            op <- which.max(store.bts[,4])
            eta.bts <- c(eta.bts, store.bts[op,1])
          }
          # Empirical bootstrap distribution for the change point estimator
          temp <- sapply(Z, function(x){
            mean(sqrt(m)*(eta.bts - eta.star) <= x)
          })
          
          Ftilde.bts <- cbind(Ftilde.bts, temp)
          
          result.bts[[j]] <- eta.bts
        }
        
        
        ks_distance <- c()
        for (l in c(1:4)) {
          ks_distance <- c(ks_distance, max(abs(Ftilde.bts[,l] - Ftilde.bts[, (l+1)])))
        }
        
        # select m with minimum KS distance
        index <- which.min(ks_distance)
        eta.bts.m <- result.bts[[index]]
        m.star <- as.integer(n * r^index)
        
        
        rate_1 <- sqrt(m.star/n)
        rate_2 <- m.star/n
        
        Q_left <- quantile(eta.bts.m, 0.025) - eta.star
        Q_right <- quantile(eta.bts.m, 0.975) - eta.star
        
        sim <- rbind(sim, c(m.star, eta.star, Q_left,Q_right, 
                 ifelse(true.eta>eta.star+Q_left&true.eta<eta.star+Q_right, 1, 0),
                 ifelse(true.eta>eta.star+rate_1*Q_left&true.eta<eta.star+rate_1*Q_right, 1, 0),
                 ifelse(true.eta>eta.star+rate_2*Q_left&true.eta<eta.star+rate_2*Q_right, 1, 0)))
        
      }
      sim.df <- data.frame(sim) 
      names(sim.df) <- c("m", "eta", "CI_left", "CI_right", "standard", "sqrt_movern", "movern")
      write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_eta=",teta,"_CI.csv"))
    }
  }
}

