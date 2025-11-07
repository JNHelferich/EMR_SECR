##Simulating some EMR data and using the discretized emr bayes code

#Load packages
library(nimble)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(sp)
library(terra)
library(raster)
library(NLMR)
library(landscapetools)
library(basicMCMCplots)
library(oSCR)
library(abind)




study_bound <- vect("source/studyarea_proj.shp")
water <- vect("source/water_project.shp")

tr_ras <- rast(study_bound, res = 100)
values(tr_ras) <- 1
tr_ras <- mask(tr_ras, study_bound)
traps <- as.data.frame(tr_ras, xy = TRUE)
traps.vect <- vect(traps, geom = c("x", "y"), crs = crs(study_bound))

plot(study_bound)
plot(traps.vect, add = TRUE)

traps <- traps[,-3]

n.traps <- nrow(traps)
#Working code, but trying to do multi-session

########################## Set Up Trap Array and Hab Cov ##########################


#makebuffer
traps_buffer <- buffer(traps.vect, width = 200)
traps_buffer <- terra::aggregate(traps_buffer)

#Let's get this down to just east and a reasonable size
template_raster <- rast(extent = ext(traps_buffer),  # Use the same extent
                        resolution = 25,          # Set the resolution to 500m
                        crs = crs(study_bound))   # Use the same coordinate reference system


water.dist <- terra::distance(template_raster, water)
water.dist <- mask(water.dist, traps_buffer)


#create dataset with grid and covariate (density between 2 and 3 meters)
habcov<- as.data.frame(water.dist$lyr.1, xy = TRUE)
colnames(habcov) <- c("x", "y", "wtr")
SS_size <- nrow(habcov)

#Scale waterdist
habcov$wtr <- scale(habcov$wtr)
habcov$wtr <- as.numeric(habcov$wtr)

#scale x, y coords
traps$xscale <- (traps$x -mean(c(max(habcov$x),min(habcov$x))))/1000
traps$yscale <- (traps$y - mean(c(max(habcov$y),min(habcov$y))))/1000

habcov$xscale <- (habcov$x - mean(c(max(habcov$x),min(habcov$x))))/1000
habcov$yscale <- (habcov$y - mean(c(max(habcov$y),min(habcov$y))))/1000

pixelarea = 0.2^2

########################## Simulate data ##########################

###Set seed and set up simulation variables
set.seed(121)

#Making dataset for the spring
#how many we want to simulate
sexratio <- .5
sigma.sp <- c(.2, .3)#sigma for females and males
lambda0.sp <- -2  #detection
K <- 10 # nb capture occasions
J <- nrow(traps) # nb of traps

N.sp <- 40
#Determine sex
sex.sp <- rbinom(N.sp, 1, sexratio) + 1  #assign sex to simulated individuals with equal sex ratio (.5)
#Determine activity centers

b0.sp<- -5.2
b1.sp<- 3
mu.sp=exp(b0.sp + b1.sp*habcov[,3])*pixelarea  #set up habitat to ac relationship
############### Return here and updated if needed
EN.sp <- sum(mu.sp)
EN.sp
#A <- 1
#N <- rpois(1, EN*A)  #come back to this to add stochaticity
pi.sp=mu.sp/sum(mu.sp)   #get probabilities for each cell based on lidar data
ac.sp<-sample(1:SS_size,N.sp, replace=TRUE, pi.sp)


##simulate capture histories
yy.sp <- array(NA, c(N.sp, J, K))#Blank data array #Individuals x Traps as one matrix, then 5 matrices for each occasion (k)
for(j in 1:J) {
  dist <- sqrt((traps$xscale[j] - habcov$xscale[ac.sp])^2 + (traps$yscale[j] - habcov$yscale[ac.sp])^2) #For each trap and individual combo, calculates distance between trap and activity center
  lam.b <- exp(lambda0.sp) 
  lambda <- lam.b * exp(-dist^2 / (2 * sigma.sp[sex.sp]^2)) #Determines the capture probability at each trap as a function of distance
  for(k in 1:K) {
    yy.sp[,j,k] <- rpois(N.sp, lambda) #For each occasion, draws from a poisson dist to see if that animal was caught at that trap
  }
}

tot <- apply(yy.sp, 2, sum) #Total number of detections at each trap across occasions
dat <- data.frame(traps, tot = tot) #Makes it a dataframe
viz_traps <- traps %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(pch = 3)

viz_traps +
  geom_point(data = dat, aes(x = x, y = y, size = tot), alpha = 0.3) +
  scale_size(range = c(0, 20)) +
  labs(x = "",
       y = "",
       size = "# detections")





#Making dataset for the summer
#how many we want to simulate
sexratio <- .5
sigma.su <- c(.2, .3)#sigma for females and males
lambda0.su <- -2  #detection
K <- 10 # nb capture occasions
J <- nrow(traps) # nb of traps

N.su <- 40
#Determine sex
sex.su <- rbinom(N.su, 1, sexratio) + 1  #assign sex to simulated individuals with equal sex ratio (.5)
#Determine activity centers
b0.su<- -3.3
b1.su<- -3
mu.su=exp(b0.su + b1.su*habcov[,3])*pixelarea  #set up habitat to ac relationship
############### Return here and updated if needed
EN.su <- sum(mu.su)
EN.su
#A <- 1
#N <- rpois(1, EN*A)  #come back to this to add stochaticity
pi.su=mu.su/sum(mu.su)   #get probabilities for each cell based on lidar data
ac.su<-sample(1:SS_size,N.su, replace=TRUE, pi.su)


##simulate capture histories
yy.su <- array(NA, c(N.su, J, K))#Blank data array #Individuals x Traps as one matrix, then 5 matrices for each occasion (k)
for(j in 1:J) {
  dist <- sqrt((traps$xscale[j] - habcov$xscale[ac.su])^2 + (traps$yscale[j] - habcov$yscale[ac.su])^2) #For each trap and individual combo, calculates distance between trap and activity center
  lam.b <- exp(lambda0.su) 
  lambda <- lam.b * exp(-dist^2 / (2 * sigma.su[sex.su]^2)) #Determines the capture probability at each trap as a function of distance
  for(k in 1:K) {
    yy.su[,j,k] <- rpois(N.su, lambda) #For each occasion, draws from a poisson dist to see if that animal was caught at that trap
  }
}

tot <- apply(yy.su, 2, sum) #Total number of detections at each trap across occasions
dat <- data.frame(traps, tot = tot) #Makes it a dataframe
viz_traps <- traps %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(pch = 3)

viz_traps +
  geom_point(data = dat, aes(x = x, y = y, size = tot), alpha = 0.3) +
  scale_size(range = c(0, 20)) +
  labs(x = "",
       y = "",
       size = "# detections")





#Making dataset for the fall
#how many we want to simulate
sexratio <- .5
sigma.fa <- c(.2, .3)#sigma for females and males
lambda0.fa <- -2  #detection
K <- 10 # nb capture occasions
J <- nrow(traps) # nb of traps

N.fa <- 40
#Determine sex
sex.fa <- rbinom(N.fa, 1, sexratio) + 1  #assign sex to simulated individuals with equal sex ratio (.5)
#Determine activity centers

b0.fa<- -5.2
b1.fa<- 3
mu.fa=exp(b0.fa + b1.fa*habcov[,3])*pixelarea  #set up habitat to ac relationship
############### Return here and updated if needed
EN.fa <- sum(mu.fa)
EN.fa
#A <- 1
#N <- rpois(1, EN*A)  #come back to this to add stochaticity
pi.fa=mu.fa/sum(mu.fa)   #get probabilities for each cell based on lidar data
ac.fa<-sample(1:SS_size,N.fa, replace=TRUE, pi.fa)


##simulate capture histories
yy.fa <- array(NA, c(N.fa, J, K))#Blank data array #Individuals x Traps as one matrix, then 5 matrices for each occasion (k)
for(j in 1:J) {
  dist <- sqrt((traps$xscale[j] - habcov$xscale[ac.fa])^2 + (traps$yscale[j] - habcov$yscale[ac.fa])^2) #For each trap and individual combo, calculates distance between trap and activity center
  lam.b <- exp(lambda0.fa) 
  lambda <- lam.b * exp(-dist^2 / (2 * sigma.fa[sex.fa]^2)) #Determines the capture probability at each trap as a function of distance
  for(k in 1:K) {
    yy.fa[,j,k] <- rpois(N.fa, lambda) #For each occasion, draws from a poisson dist to see if that animal was caught at that trap
  }
}

tot <- apply(yy.fa, 2, sum) #Total number of detections at each trap across occasions
dat <- data.frame(traps, tot = tot) #Makes it a dataframe
viz_traps <- traps %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(pch = 3)

viz_traps +
  geom_point(data = dat, aes(x = x, y = y, size = tot), alpha = 0.3) +
  scale_size(range = c(0, 20)) +
  labs(x = "",
       y = "",
       size = "# detections")







########################## Model Code ##########################

#Define the nimble model
code <- nimbleCode({
  sigma[1] ~ dunif(0, 3)  #Uniform prior for sigma
  sigma[2] ~ dunif(0, 3)
  lam0 ~ dunif(-5, 5)     #Uniform prior for baseline detection
  sexratio ~ dbeta(1, 1)
  b0.sp ~ dnorm(0, sd=2)
  b0.su ~ dnorm(0, sd=2)
  b0.fa ~ dnorm(0, sd=2)
  b1.sp ~ dnorm(0, sd=2)
  b1.su ~ dnorm(0, sd=2)
  b1.fa ~ dnorm(0, sd=2)
  b2 ~ dnorm(0, sd = 2)
  
  #Spring
  mu.sp[1:SS_size] <- exp(b0.sp + b1.sp*habcov[1:SS_size])*pixelarea
  EN.sp <- sum(mu.sp[1:SS_size])   # expected number of individuals in state-space
  pi.sp[1:SS_size] <- mu.sp[1:SS_size]/EN.sp
  psi.sp <- EN.sp / M # Prob ind is in state space
  
  #Summer
  mu.su[1:SS_size] <- exp(b0.su + b1.su*habcov[1:SS_size])*pixelarea
  EN.su <- sum(mu.su[1:SS_size])   # expected number of individuals in state-space
  pi.su[1:SS_size] <- mu.su[1:SS_size]/EN.su
  psi.su <- EN.su / M # Prob ind is in state space

  
  #Fall
  mu.fa[1:SS_size] <- exp(b0.fa + b1.fa*habcov[1:SS_size])*pixelarea
  EN.fa <- sum(mu.fa[1:SS_size])   # expected number of individuals in state-space
  pi.fa[1:SS_size] <- mu.fa[1:SS_size]/EN.fa
  psi.fa <- EN.fa / M # Prob ind is in state space
  
  

  #Spring
  for(i in 1:M) { #For each individual in our data augmented dataset
    z.sp[i] ~ dbern(psi.sp) #bernouili trial if the individual is in state space
    s.sp[i] ~ dcat(pi.sp[i:SS_size])  #assignment of activity center
    sex.sp[i] ~ dbern(sexratio)
    sex2.sp[i] <- sex.sp[i] + 1
    log(b.lam.sp[i, 1:J]) <- lam0 + b2*size.sp[i]
    lam.sp[i,1:J] <- K*b.lam.sp[i, 1:J] * exp(-d2.sp[s.sp[i], 1:J] / (2 * sigma[sex2.sp[i]]^2)) * z.sp[i]  #Calculates actual detection rate
    
    for(j in 1:J){
      yy.sp[i,j] ~ dpois(lam.sp[i, j])
    }
  }
  
  #Summer
  for(i in 1:M) { #For each individual in our data augmented dataset
    z.su[i] ~ dbern(psi.su) #bernouili trial if the individual is in state space
    s.su[i] ~ dcat(pi.su[i:SS_size])  #assignment of activity center
    sex.su[i] ~ dbern(sexratio)
    sex2.su[i] <- sex.su[i] + 1
    log(b.lam.su[i, 1:J]) <- lam0 + b2*size.su[i]
    lam.su[i,1:J] <- K*b.lam.su[i, 1:J] * exp(-d2.su[s.su[i], 1:J] / (2 * sigma[sex2.su[i]]^2)) * z.su[i]  #Calculates actual detection rate

    for(j in 1:J){
      yy.su[i,j] ~ dpois(lam.su[i, j])
    }
  }

  #Fall
  for(i in 1:M) { #For each individual in our data augmented dataset
    z.fa[i] ~ dbern(psi.fa) #bernouili trial if the individual is in state space
    s.fa[i] ~ dcat(pi.fa[i:SS_size])  #assignment of activity center
    sex.fa[i] ~ dbern(sexratio)
    sex2.fa[i] <- sex.fa[i] + 1
    log(b.lam.fa[i, 1:J]) <- lam0 + b2*size.fa[i]
    lam.fa[i,1:J] <- K*b.lam.fa[i, 1:J] * exp(-d2.fa[s.fa[i], 1:J] / (2 * sigma[sex2.fa[i]]^2)) * z.fa[i]  #Calculates actual detection rate
    
    for(j in 1:J){
      yy.fa[i,j] ~ dpois(lam.fa[i, j])
    }
  }
  
  
  
  
  
  N.sp <- sum(z.sp[1:M]) #This is just producing total abundance of what is left from M
  N.su <- sum(z.su[1:M]) #This is just producing total abundance of what is left from M
  N.fa <- sum(z.fa[1:M]) #This is just producing total abundance of what is left from M
})


J <- nrow(traps) # nb of traps
K <- 10 # nb capture occasions


M <- 100 #Number of data augmented individuals

#Create distance matrix to input as data
d2=matrix(NA, nrow(habcov), nrow(traps))
for(i in 1:nrow(habcov)){
  d2[i,]<-(habcov[i,4]-traps[,3])^2 + (habcov[i,5]-traps[,4])^2
}


#since we don't have any time covariates, sum across occasions
yy.sp.k <- apply(yy.sp, c(1,2), sum)
yy.su.k <- apply(yy.su, c(1,2), sum)
yy.fa.k <- apply(yy.fa, c(1,2), sum)

#Creating data augmented datasets
yaug.sp <- matrix(0, M, n.traps)#Creating augmented dataset of individuals
yaug.sp[1:N.sp, 1:n.traps] <- yy.sp.k
sex.aug.sp <- c(sex.sp-1, rep(NA, M-N.sp))
size.aug.sp <- c(rnorm(N.sp, .5, .2), rep(NA, M-N.sp))

yaug.su <- matrix(0, M, n.traps)#Creating augmented dataset of individuals
yaug.su[1:N.su, 1:n.traps] <- yy.su.k
sex.aug.su <- c(sex.su-1, rep(NA, M-N.su))
size.aug.su <- c(rnorm(N.su, .5, .2), rep(NA, M-N.su))

yaug.fa <- matrix(0, M, n.traps)#Creating augmented dataset of individuals
yaug.fa[1:N.fa, 1:n.traps] <- yy.fa.k
sex.aug.fa <- c(sex.fa-1, rep(NA, M-N.fa))
size.aug.fa <- c(rnorm(N.fa, .5, .2), rep(NA, M-N.fa))

sex.init.sp <- c(rep(NA, N.sp), rbinom(M-N.sp, 1, 0.5))
sex.init.su <- c(rep(NA, N.su), rbinom(M-N.su, 1, 0.5))
sex.init.fa <- c(rep(NA, N.fa), rbinom(M-N.fa, 1, 0.5))

size.init.sp <- c(rep(NA, N.sp), rnorm(M-N.sp, .5, 0.2))
size.init.su <- c(rep(NA, N.su), rnorm(M-N.su, .5, 0.2))
size.init.fa <- c(rep(NA, N.fa), rnorm(M-N.fa, .5, 0.2))




sinit.sp <- round(runif(M, 1, SS_size))
sinit.su <- round(runif(M, 1, SS_size))
sinit.fa <- round(runif(M, 1, SS_size))
sinit.sp[1:40] <- ac.sp
sinit.su[1:40] <- ac.su
sinit.fa[1:40] <- ac.fa



zinit.sp <- c(rep(1, N.sp), rbinom(M-N.sp, 1, 0.5)) #individuals alive state
zinit.su <- c(rep(1, N.su), rbinom(M-N.su, 1, 0.5)) #individuals alive state
zinit.fa <- c(rep(1, N.fa), rbinom(M-N.fa, 1, 0.5)) #individuals alive state


#Defining constanst
constants <- list(M = M, 
                  K = K, #Occasions
                  J = J,
                  SS_size = SS_size,
                  pixelarea = pixelarea) #put back in data #Traps

#Define data
data <- list(yy.sp = yaug.sp, #Data input
             yy.su = yaug.su,
             yy.fa = yaug.fa,
             sex.sp=sex.aug.sp,
             sex.su = sex.aug.su,
             sex.fa = sex.aug.fa,
             size.sp=size.aug.sp,
             size.su = size.aug.su,
             size.fa = size.aug.fa,
             d2.sp = d2,
             d2.su = d2,
             d2.fa = d2,
             habcov = habcov[,3])  #maybe standardize habcov 

#Set inits
inits <- list(sigma = c(1, 1),  #Inputs all these as initial values
              lam0 = -2, 
              s.sp = sinit.sp,
              s.su = sinit.su, 
              s.fa = sinit.fa,
              z.sp = zinit.sp,
              z.su = zinit.su,
              z.fa = zinit.fa,
              sexratio = 0.5,
              sex.sp = sex.init.sp,
              sex.su = sex.init.su,
              sex.fa = sex.init.fa,
              size.sp = size.init.sp,
              size.su = size.init.su,
              size.fa = size.init.fa,
              b0.sp = b0.sp,
              b0.su = b0.su,
              b0.fa = b0.fa,
              b1.sp = b1.sp,
              b1.su = b1.su,
              b1.fa = b1.fa)

#Fit the model
Rmodel <- nimbleModel(code = code,  #This builds the model in R, but is not compiled yet.
                      constants = constants, 
                      data = data, 
                      inits = inits)

Cmodel <- compileNimble(Rmodel)#Compiles the model in c++
conf <- configureMCMC(Rmodel, #Configures MCMC with monitors
                      monitors = c("N.sp", "N.su", "N.fa", "lam0", "psi.sp", "psi.su", "psi.fa", "sigma", "b0.sp", "b0.su", "b0.fa", "b1.sp", "b1.su", "b1.fa", "b2", "sexratio"))
Rmcmc <- buildMCMC(conf)#Building the chains
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)#Compile in c++
samplesList <- runMCMC(Cmcmc, #Run actual compiled model with 2 chains and 10000 iterations
                       niter = 1000,
                       nburnin = 500,
                       nchains = 1)

samples <- rbind(samplesList[[1]], #Combining the two chains
                 samplesList[[2]])

library(basicMCMCplots) #Trace plots
chainsPlot(samplesList,
           var = c("N.ocel", "N.bob", "sigma", "lam0", "sexratio", "b0.ocel", "b1.ocel"))

summary(samplesList)



#Cool! glad that worked, now let's make a rectangular stat space for simulation, very loosely based on El Saux





