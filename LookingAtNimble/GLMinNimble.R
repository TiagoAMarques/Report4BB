# GLM
# using the data

#for comparizon with a lmer model
glm1<-glm(crate~1,data=tags,family=Gamma(link="log"))
summary(glm1)


## define the model
GLMcode <- nimbleCode({
  # the overall intercept
  beta0 ~ dnorm(0,sd=10)
  # the gamma variance parameter
  disp ~ dunif(0, 100)
  for (i in 1:N){
    #get the linear predictor, consider a log link function
    log(mean[i]) <- beta0
    # In bugs the parameterization  of the gamma is done using the shape and the rate
    # via mean and variance, defined as below:
    crate[i] ~ dgamma(shape=1/disp,scale=mean[i]*disp)
  }
})

## constants, data, and initial values

#constant, sample size and number of levels for each of the random effects
#number of rows in tags
N<-nrow(tags)
constants <- list(N = N)

#data
data <- list(crate = tags$crate)

#initial values
#inits <- list(beta0 = 0, sigmal_RE = 1, sigmay_RE = 1, dispersion = 0.8,REy = rep(0,nyears),REl = rep(0,nlocs))
inits <- list(beta0 = 0, disp = 1)

## create the model object
myGLMModel <- nimbleModel(code = GLMcode, constants = constants, data = data, 
                       inits = inits, check = FALSE)

#things to monitor
#tomon<-c("beta0","dispersion","sigmay_RE","sigmal_RE","crate","REy","REl")
tomon<-c("beta0","disp","crate")

testGLM<-nimbleMCMC(myGLMModel,monitors=tomon,niter=50000,nburnin=10000,progressBar=TRUE,summary=TRUE)

#trace plots to cehck convergence
par(mfrow=c(1,2))
#trace plot intercept
plot(testGLM$samples[,1],pch=".",ylab="intercept")
#trace plot dispersion
plot(test$samples[,2],pch=".",ylab="dispersion")

#look at main results
testGLM$summary[c(1,114),]


#check vs glm
par(mfrow=c(1,2))
# posterior plot intercept
hist(testGLM$samples[,1],pch=".",xlab="intercept",main="")
#the intercept from glm
abline(v=summary(glm1)$coeff[1,1],col="red")
abline(v=quantile(testGLM$samples[,1],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testGLM$samples[,1]),col="blue",lty=2)
# posterior plot residual error
hist(testGLM$samples[,114],pch=".",xlab="dispersion",main="")
# the glm dispersion
#interesting to see how the GLM estimate is off
abline(v=summary(glm1)$dispersion,col="red")
# but that of gamma dispersion is OK-ish
abline(v=gamma.dispersion(glm1),col="red")
# standard deviation of the gamma GLM given by (see my html on Gamma regression )
#abline(v=sqrt(sum(residuals(glm1, type = "response")^2)/glm1$df.residual),col="red")
abline(v=quantile(testGLM$samples[,114],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(testGLM$samples[,114]),col="blue",lty=2)

# the comparison for the Gamma glm results looks much better 
# (in the sense of things look more alike to the glm results) 
# than that of the Gamma glmm 