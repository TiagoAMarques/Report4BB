# gaussian LMM
# using the data

#for comparizon with a lmer model
glm1<-lmer(formula=crate~1+(1|location)+(1|fyear),data=tags)
summary(glm1)


## define the model
LMMcode <- nimbleCode({
  # the overall intercept
  beta0 ~ dnorm(0,sd=10)
  # random effect standard deviation associated with location, a uniform, might change this to be a gamma latter
  sigmal_RE ~ dunif(0, 2)
  # random effect standard deviation associated with year, a uniform, might change this to be a gamma latter
  sigmay_RE ~ dunif(0, 2)
  # the residual error standard deviation
  errorsd ~ dunif(0, 100)
  #get year random effects
  for(yy in 1:nyears){
    #REy[yy] ~ dnorm(0, sd = sigmay_RE)
    # Using decentered parametrization, a suggestion by Ben Augustine
    REy[yy] ~ dnorm(0, sd = 1)
  }  
  #get location random effects
  for(ll in 1:nlocs){
    #REl[ll] ~ dnorm(0, sd = sigmal_RE)
    # Using decentered parametrization, a suggestion by Ben Augustine
    REl[ll] ~ dnorm(0, sd = 1)
  }  
  for (i in 1:N){
    # Using decentered parametrization, a suggestion by Ben Augustine
    mean[i] <- beta0 + REy[year[i]]*sigmay_RE + REl[loc[i]]*sigmal_RE
    # In bugs the parameterization  of the gamma is done using the shape and the rate
    # via mean and variance, defined as below:
    crate[i] ~ dnorm(mean=mean[i],sd=errorsd)
  }
})

## constants, data, and initial values

#constant, sample size and number of levels for each of the random effects
#number of rows in tags
N<-nrow(tags)
#number of different years
nyears <- length(unique(tags$year))
#number of different locations
nlocs <- length(unique(tags$location))
#the year covariate is passed as a constant
year <- as.numeric(tags$fyear)
#the location covaraite is passed as a constant
loc <- as.numeric(tags$location)
#bundle all in a suitable object
constants <- list(N = N,nyears=nyears,nlocs=nlocs,year = year,loc = loc)

#data
data <- list(crate = tags$crate)

#initial values
#inits <- list(beta0 = 0, sigmal_RE = 1, sigmay_RE = 1, dispersion = 0.8,REy = rep(0,nyears),REl = rep(0,nlocs))
inits <- list(beta0 = 0, sigmal_RE = 1, sigmay_RE = 1, errorsd = 1,REy = rep(0,nyears),REl = rep(0,nlocs))


## create the model object
myLMMModel <- nimbleModel(code = LMMcode, constants = constants, data = data, 
                       inits = inits, check = FALSE)

#things to monitor
#tomon<-c("beta0","dispersion","sigmay_RE","sigmal_RE","crate","REy","REl")
tomon<-c("beta0","errorsd","sigmay_RE","sigmal_RE","crate","REy","REl")

testLMM<-nimbleMCMC(myLMMModel,monitors=tomon,niter=50000,nburnin=10000,progressBar=TRUE,summary=TRUE)

#trace plots to cehck convergence
par(mfrow=c(2,2))
#trace plot intercept
plot(test$samples[,20],pch=".",ylab="intecept")
#trace plot dispersion
plot(test$samples[,135],pch=".",ylab="dispersion")
#trace plot year random effect standard deviation
plot(test$samples[,133],pch=".",ylab="location random effect sigma")
#trace plot location random effect standard deviation
plot(test$samples[,134],pch=".",ylab="year random effect sigma")

#look at main results
testLMM$summary[c(20,133,134,135),]


#check vs lmer
par(mfrow=c(2,2))
# posterior plot intercept
hist(testLMM$samples[,20],pch=".",xlab="intercept",main="")
#the intercept from glmer
abline(v=summary(lmer1)$coeff[1,1],col="red")
abline(v=quantile(testLMM$samples[,20],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testLMM$samples[,20]),col="blue",lty=2)
# posterior plot residual error
hist(testLMM$samples[,133],pch=".",xlab="residual error",main="")
#the variance in glmer (I suspect that is the dispersion, waiting for BB to confirm)
abline(v=summary(lmer1)$sigma,col="red")
abline(v=quantile(testLMM$samples[,133],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testLMM$samples[,133]),col="blue",lty=2)
# posterior plot year random effect standard deviation
hist(testLMM$samples[,134],pch=".",xlab="location random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(lmer1)$varcor))[2],col="red")
abline(v=quantile(testLMM$samples[,134],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testLMM$samples[,134]),col="blue",lty=2)
# posterior plot location random effect standard deviation
hist(testLMM$samples[,135],pch=".",xlab="year random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(lmer1)$varcor))[1],col="red")
abline(v=quantile(testLMM$samples[,135],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testLMM$samples[,135]),col="blue",lty=2)

# the comparison for the lmer results looks better 
# (in the sense of things look much more alike) than that
# of the gamma glmm but there are still some inconsistencies