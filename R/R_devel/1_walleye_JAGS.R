# rm(list=ls())
library(R2jags)
library(lme4)
library(MCMCpack)
library(arm)
library(lattice)
library(PerformanceAnalytics)



# Read in data
dat <- read.csv('walleye_recruitment_GDD.csv')
dim(dat)
head(dat)

length(unique(dat$WBIC))





# Rename and scale covarites
dat$y <- dat$recruitment.binary
# Log-transform and then standardize GDD
dat$x1 <- log(dat$GDD_wtr_5c) 
dat$x <- (dat$x1 - mean(dat$x1))/sd(dat$x1)
# summary(dat$x)
# dat$area <- as.numeric(scale(log(dat$area.hectares)))
# # hist(dat$area)
# dat$cond <- as.numeric(scale(log(dat$Conductance)))
# dat$bass <- dat$mean.lmb.standardized # Should probably log-transform and then standardize
dat$lake <- dat$Lake
# dat$lat <- as.numeric(scale(dat$latitude))
dat$lakenum <- as.numeric(as.factor(as.numeric(dat$WBIC)))
summary(dat)


# Remove lakes with missing data
dat <- dat[complete.cases(dat),]
length(unique(dat$lakenum))

# renumber lakenum
dat$lakenum <- as.numeric(as.factor(as.numeric(dat$lakenum)))



#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1+exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta[group[i]] * x[i]	# linear predictor    
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
    BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta ~ dnorm(0, 0.0001)
    gamma0b ~ dnorm(0, 0.0001)
    gamma0b2 ~ dnorm(0, 0.0001)
    gamma0b3 ~ dnorm(0, 0.0001)

    gamma1b ~ dnorm(0, 0.0001)
    gamma1b2 ~ dnorm(0, 0.0001)
    gamma1b3 ~ dnorm(0, 0.0001)
 
    
    
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()



# Number of parameters
K <- 2

# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)

W <- diag(K)

# Level-2 covariates
# lake area
area <- as.numeric(by(dat$area.hectares, dat$lakenum, mean)) 
area <- as.numeric(scale(log(area)))
hist(area)

# Conductivity
cond <- as.numeric(by(dat$Conductance, dat$lakenum, mean)) 
cond <- as.numeric(scale(log(cond)))
hist(cond)

# Latitude
lat <- as.numeric(by(dat$lat, dat$lakenum, mean)) 
hist(lat)

# GDD
gdd <- as.numeric(by(dat$GDD_wtr_5c, dat$lakenum, mean)) 
gdd2 <- as.numeric(scale(log(gdd)))
mat1 <- cbind(area, cond, lat, gdd2) # gdd and lat correlated: r = -0.75
cor(mat1)

# Number of lakes
J <- length(unique(dat$lakenum))


# Load data
data <- list(y = dat$y, group = dat$lakenum, n = dim(dat)[1], J = J,
             x=dat$x, K=K, W=W, z1 = area, z2=cond, z3=gdd2 )


# Get data in form for lmer
dat$cond <- cond[dat$lakenum]
dat$area <- area[dat$lakenum]
dat$gdd2 <- gdd2[dat$lakenum]


# Fit model using glmer
m1 <- glmer(y ~ 1 + x + area + x:area + cond + x:cond + (1+x|lakenum), 
            control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)),
            family = binomial, data=dat)
summary(m1)


# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), mu.beta=rnorm(1), 
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma0b = rnorm(1), gamma1b=rnorm(1),gamma0b2=rnorm(1),gamma0b3=rnorm(1),
        gamma1b2=rnorm(1),gamma1b3=rnorm(1) )
}


# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B","gamma0b","gamma1b","gamma0b2","gamma0b3",
                "gamma1b2","gamma1b3")


# MCMC settings
ni <- 20000
nt <- 2
nb <- 10000
nc <- 3

# use jagsUI and parallel processing for this one.
start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out1, dig = 3)

# Find which parameters, if any, have Rhat > 1.1
which(out1$BUGSoutput$summary[, c("Rhat")] > 1.1)

