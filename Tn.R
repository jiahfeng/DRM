#######################################################
# Tn.R
# Jiahui Feng
# update Dec. 23, 2023
# 
# Goodness-of-fit test to validate the model assumption
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


## under the null hypothesis
for(N0 in c(200,400)){
  for(N1 in c(200,400)){
    for(teta in c(-0.5,0,0.5)){
      for(alpha in c(1,-1)){ # for alpha = 0, just set teta = 0
        sim <- c()
        for (i in 1:1000) {
          set.seed(i)
          n0 <- N0
          n1 <- N1
          # data generation
          true.eta <- teta
          gamma <- -log(exp(0.5*alpha^2-alpha*true.eta)*(1-pnorm(true.eta-alpha))+pnorm(true.eta))
          Z0 <- rnorm(n0, 0, 1)
          pU <- pnorm(q = true.eta)/((exp(0.5*alpha^2-alpha*true.eta))*(1-pnorm(true.eta-alpha))+pnorm(q = true.eta))
          U <- rbinom(n1, size = 1,prob = pU)
          Z1 <- ifelse(U==0, rtruncnorm(n = n1,mu = alpha,sigma = 1,low = true.eta,high = Inf), 
                        rtruncnorm(n = n1,mu = 0,sigma = 1,low = -Inf,high = true.eta))
          dat <- data.frame(Y=c(rep(0,n0),rep(1,n1)),Z=c(Z0,Z1))
          Y <- dat$Y
          Z <- dat$Z
          n <- n0+n1
          
          # candidates for the change point eta
          space <- seq(-1.5,1.5,0.1)
          
          # estimation given the value of eta
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
          # select the set of parameters with maximum likelihood function value
          op <- which.max(store[,4])
          eta.star <- store[op,1]
          gamma.star <- store[op,2]
          alpha.star <- store[op,3]
          r_zeta.star <- exp(gamma.star+alpha.star*ifelse(Z-eta.star>=0,Z-eta.star,0))
          
          Fx <- sapply(sort(Z), function(o){sum(ifelse(o>Z,1,0)/(n0+n1*r_zeta.star))}) # estimated distribution function of F
          Gx <- sapply(sort(Z), function(o){sum(r_zeta.star*ifelse(o>Z,1,0)/(n0+n1*r_zeta.star))}) # estimated distribution function of G
          EFx <- sapply(sort(Z), function(o){sum(ifelse(o>Z0,1,0)/n0)}) # empirical distribution function of F
          EGx <- sapply(sort(Z), function(o){sum(ifelse(o>Z1,1,0)/n1)}) # empirical distribution function of G
          
          Tn <- sqrt(n)*max(abs(Fx-EFx)) # test statistic
          
          
          # bootstrap
          Tn.bts <- c()
          B <- 2000
          for(j in c(1:B)){
            set.seed(i+j)
            # generate bootstrap sample based on the estimated distribution functions
            U0 <- runif(n0)
            U1 <- runif(n1)
            Z0.bts <- sapply(c(1:n0),function(o){sort(Z)[which.min(abs(U0[o]-Fx))]})
            Z1.bts <- sapply(c(1:n1),function(o){sort(Z)[which.min(abs(U1[o]-Gx))]})
            dat.bts <- data.frame(Y=c(rep(0,n0),rep(1,n1)), Z=c(Z0.bts,Z1.bts))
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
                n1e <- n1*r_zeta/(n0+n1*r_zeta)
                n2e <- n0*n1*r_zeta/(n0+n1*r_zeta)^2
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
            r_zeta.star <- exp(store.bts[op,2]+store.bts[op,3]*ifelse(Z.bts-store.bts[op,1]>=0,Z.bts-store.bts[op,1],0))
            Fx.bts <- sapply(sort(Z.bts), function(o){sum(ifelse(o>Z.bts,1,0)/(n0+n1*r_zeta.star))})
            EFx.bts <- sapply(sort(Z.bts), function(o){sum(ifelse(o>Z0.bts,1,0)/n0)})
            
            Tn.bts <- c(Tn.bts,sqrt(n)*max(abs(Fx.bts-EFx.bts))) # test statistic of the bootstrap sample
          }
          sim <- rbind(sim, c(eta.star,gamma.star,alpha.star,Tn,
                              quantile(Tn.bts,c(0.9,0.95,0.99)),ifelse(Tn>=quantile(Tn.bts,c(0.9,0.95,0.99)),1,0)))
        }
        sim.df <- data.frame(sim)
        
        names(sim.df) <- c("eta","gamma","alpha","Tn", "Tn.bts90","Tn.bts95","Tn.bts99","Tn.rej90","Tn.rej95","Tn.rej99")
        write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_eta=",teta,"_alpha=",alpha,"_H0.csv"))
      }
    }
  }
}



## under the alternative hypothesis of not being a two-sample density ratio model ##
for(N0 in c(200,400)){
  for(N1 in c(200,400)){
    for (model in c(1,2,3)) {
      sim <- c()
      for (i in 1:1000) {
        set.seed(i)
        n0 <- N0
        n1 <- N1
        # data generation
        Z0 <- rnorm(n0, 0, 1)
        
        if (model == 1){
          beta <- -0.5
          Z1 <- rnorm(n1, 0, sqrt(1/(1-2*beta))) ## Model 1: g(x) = exp(gamma - 0.5*x^2) * f(x)
        }
        
        # rejection sampling
        if (model == 2){
          x_lower <- -3
          x_upper <- 3
          g <- function(x, gamma) exp(gamma + 2*sin(x)) * dnorm(x) ## Model 2: g(x) = exp(2 + sin(x)) * f(x) 
          # cdf of g
          intg <- function(gamma) {
            integrate(g, lower = x_lower, upper = x_upper, gamma=gamma)
          }
          objective <- function(gamma, value) {
            abs(intg(gamma)$value - value)
          }
          
          gammahat <- optimize(objective, c(-10, 10), value=1)$minimum # find the value of gamma such that the integration of g(x) is 1 
          g_max <- optimize(g, interval=c(x_lower, x_upper), maximum=TRUE, gamma = gammahat)$objective
          
          accept = NULL
          accept.count = 0
          
          while (accept.count < n1) {
            X = runif(1, min=x_lower, max=x_upper) #uniform distribution
            Y = runif(1, min=0, max=g_max)
            if (Y <= g(X, gammahat)){
              accept.count = accept.count+1
              accept[accept.count] = X
            } 
          }
          Z1   <-accept
        }
        
        
        if (model == 3){
          beta <- -0.5
          true.eta <- 0
          alpha <- 1
          # gamma<- -log(exp(0.5*alpha^2-alpha*true.eta)*(1-pnorm(true.eta-alpha))+pnorm(true.eta))
          pU <- pnorm(sqrt(1-2*beta)*true.eta) / (exp(alpha^2/(2*(1-2*beta)) - alpha*true.eta) * 
                                                    (1 - pnorm(sqrt(1-2*beta)*(true.eta-alpha))) + pnorm(sqrt(1-2*beta)*true.eta)) 
          U <-rbinom(n1, size = 1,prob = pU)
          Z1 <-ifelse(U==0, rtruncnorm(n = n1,mu = alpha/(1-2*beta),sigma = sqrt(1/(1-2*beta)),low = true.eta,high = Inf), 
                      rtruncnorm(n = n1,mu = 0,sigma = sqrt(1/(1-2*beta)),low = -Inf,high = true.eta))
          ## Model 3: g(x) = exp(gamma + max(x,0) - 0.5*x^2) * f(x)
        }
        
        
        dat <- data.frame(Y=c(rep(0,n0),rep(1,n1)),Z=c(Z0,Z1))
        Y <- dat$Y
        Z <- dat$Z
        n <- n0+n1
        
        # candidates for the change point eta
        space <- seq(-1.5,1.5,0.1)
        
        
        # estimation given the value of eta
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
        
        # select the set of parameters with maximum likelihood function value
        op <- which.max(store[,4])
        eta.star <- store[op,1]
        gamma.star <- store[op,2]
        alpha.star <- store[op,3]
        r_zeta.star <- exp(gamma.star+alpha.star*ifelse(Z-eta.star>=0,Z-eta.star,0))
        
        Fx <- sapply(sort(Z), function(o){sum(ifelse(o>Z,1,0)/(n0+n1*r_zeta.star))}) # estimated distribution function of F
        Gx <- sapply(sort(Z), function(o){sum(r_zeta.star*ifelse(o>Z,1,0)/(n0+n1*r_zeta.star))}) # estimated distribution function of G
        EFx <- sapply(sort(Z), function(o){sum(ifelse(o>Z0,1,0)/n0)}) # empirical distribution function of F
        EGx <- sapply(sort(Z), function(o){sum(ifelse(o>Z1,1,0)/n1)}) # empirical distribution function of G
        
        Tn <- sqrt(n)*max(abs(Fx-EFx)) # test statistic
        
        # Bootstrap
        Tn.bts <- c()
        B <- 2000
        for(j in c(1:B)){
          set.seed(i+j)
          # generate bootstrap sample based on the estimated distribution functions
          U0 <- runif(n0)
          U1 <- runif(n1)
          Z0.bts <- sapply(c(1:n0),function(o){sort(Z)[which.min(abs(U0[o]-Fx))]})
          Z1.bts <- sapply(c(1:n1),function(o){sort(Z)[which.min(abs(U1[o]-Gx))]})
          dat.bts <- data.frame(Y=c(rep(0,n0),rep(1,n1)), Z=c(Z0.bts,Z1.bts))
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
              n1e <- n1*r_zeta/(n0+n1*r_zeta)
              n2e <- n0*n1*r_zeta/(n0+n1*r_zeta)^2
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
          r_zeta.star <- exp(store.bts[op,2]+store.bts[op,3]*ifelse(Z.bts-store.bts[op,1]>=0,Z.bts-store.bts[op,1],0))
          Fx.bts <- sapply(sort(Z.bts), function(o){sum(ifelse(o>Z.bts,1,0)/(n0+n1*r_zeta.star))})
          EFx.bts <- sapply(sort(Z.bts), function(o){sum(ifelse(o>Z0.bts,1,0)/n0)})
          
          Tn.bts <- c(Tn.bts,sqrt(n)*max(abs(Fx.bts-EFx.bts))) # test statistic of the bootstrap sample
        }
        sim <- rbind(sim, c(eta.star,gamma.star,alpha.star,Tn,quantile(Tn.bts,c(0.9,0.95,0.99)),
                            ifelse(Tn>=quantile(Tn.bts,c(0.9,0.95,0.99)),1,0)))
      }
      sim.df <- data.frame(sim)
      
      names(sim.df) <- c("eta","gamma","alpha","Tn", "Tn.bts90","Tn.bts95","Tn.bts99","Tn.rej90","Tn.rej95","Tn.rej99")
      write.csv(sim.df, paste0("n0=",N0,"_n1=",N1,"_beta=",beta, "_model=", model, "_H1.csv"))
    }
  }
}

