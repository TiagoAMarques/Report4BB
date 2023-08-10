# Code edited by Wei Zhang - sent to TAM via gmail on the 4th August 2023


#for GLMMs
library(lme4)
#trying NIMBLE for Bayesian implementation
library(nimble)
# for pipe operator
library(tidyverse)

# set up everything regarding the data to be as in SpermWhaleCueRatesSuppPub1.Rmd
# see also SpermWhaleCueRatesSuppPub1.Rmd for why all this
# which is essentially all the data processing code up to where the GLMM model (crglmer) is fitted in glmer 

# Reading the data that contain the information per deep dive cycle - object ddata1
# file created in Cue_Rates_For_Sperm_Whales.Rmd
# and currently distributed in the parent folder of the github repos at 
# https://github.com/TiagoAMarques/Report4BB
load("../../data_4_article_clickrates_deep_dive.rda")
#removing the tags for animals we know were exposed to sonar
DDCs<-ddata1[ddata1$sonar!="sonar",]
# Creating the data per tag
tags<-DDCs%>%
  group_by(tag)%>%
  summarise(location=unique(location), year=unique(year), sex=unique(sex),
            duration= sum(durations,na.rm=T),nclicks=sum(nclick,na.rm=T),
            crate=sum(nclick,na.rm=T)/sum(durations,na.rm=T),ddc=max(absdives+1,na.rm = T))
#removing the location with a single tag
tags<-tags[tags$location!="Norway Andenes",]
#grouping single tag per year into adjacent years
tags$year[tags$location=="Norway" & tags$year==2009]<- 2010
tags$year[tags$location=="DOMINICA" & tags$year==2017]<- 2016
tags$year[tags$location=="Mediterranean" & tags$year==2001]<- 2003

#making a new variable, year as factor
tags$fyear<-as.factor(tags$year)
CRglm3<-glm(crate~location+fyear,data=tags,family=Gamma(link="log"))

#define DOMINICA as baseline
tags$location<-factor(tags$location,levels=c("DOMINICA","Gulf of Mexico","Azores","Kaikoura","Mediterranean","North Atlantic Delaware","Norway"))
#run model
crglmerDom<-glmer(crate~location+(1|fyear),data=tags,family=Gamma(link="log"))
#define GoM as baseline
tags$location<-factor(tags$location,levels=c("Gulf of Mexico","DOMINICA","Azores","Kaikoura","Mediterranean","North Atlantic Delaware","Norway"))
# run model
crglmerGoM<-glmer(crate~location+(1|fyear),data=tags,family=Gamma(link="log"))

# the lme4 equivalent of what we will do below in NIMBLE
crglmer<-glmer(crate~(1|location)+(1|fyear),data=tags,family=Gamma(link="log"))

# Use glmmTMB, parametrizes the model in the same way as the nimble model code below
library(glmmTMB)
crglmmTMB<-glmmTMB(crate~(1|location)+(1|fyear),data=tags,family=Gamma(link="log"))


## define the model
GLMMcode <- nimbleCode({
  # the overall intercept
  beta0 ~ dnorm(0, sd = 10)
  # random effect standard deviation associated with location, a uniform, might change this to be something else latter
  sigmal_RE ~ dunif(0,10)
  # random effect standard deviation associated with year, a uniform, might change this to be something else latter
  sigmay_RE ~ dunif(0, 10)
  # the gamma dispersion (or variance - see commented parametrization 1) parameter, a uniform, might change this to be something else latter
  #dispersion ~ dunif(0, 10)
  disp ~ dunif(0, 10)
  ## sd ~ dhalfflat()
  #get year random effects
  for(yy in 1:nyears){
    #REy[yy] ~ dnorm(0, sd = sigmay_RE)
    REy[yy] ~ dnorm(0, sd = sigmay_RE)
  }  
  #get location random effects
  for(ll in 1:nlocs){
    #REl[ll] ~ dnorm(0, sd = sigmal_RE)
    REl[ll] ~ dnorm(0, sd = sigmal_RE)
  }  
  for (i in 1:N){
    #get the linear predictor, consider a log link function
    log(mean[i]) <- beta0 + REy[year[i]] + REl[loc[i]]
    # now Using decentered parametrization, a suggestion by Ben Augustine
    # log(mean[i]) <- beta0 + REy[year[i]]*sigmay_RE + REl[loc[i]]*sigmal_RE
    # parametrization 1 - now I know that to not be what is to be used
    #left here for future reference
    # crate[i] ~ dgamma(shape=(mean[i]^2)/disp,scale=disp/mean[i])
    #parametrization 2
    crate[i] ~ dgamma(shape=1/disp, scale=mean[i]*disp)
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
inits <- list(beta0 = 0, sigmal_RE = 1, sigmay_RE = 1, disp = 1,REy = rep(0,nyears),REl = rep(0,nlocs))


## create the model object
myGLMMModel <- nimbleModel(code = GLMMcode, constants = constants, data = data, 
                       inits = inits, check = FALSE, buildDerivs = TRUE) ## Add buildDerivs = TRUE for AD
cmyGLMMModel <- compileNimble(myGLMMModel)

## Build Laplace approximation in nimble
glmmNimLaplace <- buildLaplace(myGLMMModel)
cglmmNimLaplace <- compileNimble(glmmNimLaplace, project = myGLMMModel)
MLEres <- cglmmNimLaplace$findMLE()
summ_MLEres <- cglmmNimLaplace$summary(MLEres)

## nimble and glmmTMB give close answers based on Laplace approximation
summ_MLEres$params$estimates
#output with proper names - via code in https://r-nimble.org/html_manual/cha-AD.html
summaryLaplace(cglmmNimLaplace, MLEres)$params
crglmmTMB

#Consider full MCMC
#things to monitor
#tomon<-c("beta0","dispersion","sigmay_RE","sigmal_RE","crate","REy","REl")
tomon<-c("beta0","disp","sigmay_RE","sigmal_RE","crate","REy","REl")

test<-nimbleMCMC(myGLMMModel,monitors=tomon,niter=50000,nburnin=10000,progressBar=TRUE,summary=TRUE)

#look at main results
test$summary[c(20,133,134,135),]

#the beta0 matches OK-ish with the intercept from glmer
test$summary[20,1]
summary(crglmer)$coeff[1,1]
summary(crglmmTMB)$coefficients$cond[1]
summ_MLEres$params$estimates[1]

# Note that parametrized as 2 the dispersion parameter does not seem to bear a direct relation to the variance reported by glmer
#or... does it? see comment in glmmTMB output below
test$summary[133,1];
summary(crglmer)$sigma^2
#in glmmTMB output it is noted
#Dispersion estimate for Gamma family (sigma^2):
summary(crglmmTMB)$sigma^2 #Closer
#another way
#sd(resid(crglmer,type='pear'))^2
#sd(resid(crglmer,type='response'))^2
summ_MLEres$params$estimates[4]

#but the variances of the random effects are still off
test$summary[134:135,1]
summary(crglmer)$varcor
summary(crglmmTMB)$varcor ## Closer
summ_MLEres$params$estimates[2:3]

par(mfrow=c(2,2))
#trace plot intercept
plot(test$samples[,20],pch=".",ylab="intercept")
#trace plot dispersion
plot(test$samples[,133],pch=".",ylab="dispersion")
#trace plot year random effect standard deviation
plot(test$samples[,134],pch=".",ylab="location random effect sigma")
#trace plot location random effect standard deviation
plot(test$samples[,135],pch=".",ylab="year random effect sigma")

par(mfrow=c(2,2))
#posterior plot intercept
hist(test$samples[,20],pch=".",xlab="intecept",main="")
#the intercept from glmer
abline(v=summary(crglmer)$coeff[1,1],col="red")
abline(v=quantile(test$samples[,20],probs=c(0.025,0.5,0.975)),col=c("green","orange","green"),lty=2,lwd=2)
abline(v=mean(test$samples[,20]),col="blue",lty=2)
#posterior plot intercept response scale
hist(exp(test$samples[,20]),pch=".",xlab="overall mean (i.e. intercept in response scale)",main="")
#the intercept from glmer
abline(v=exp(summary(crglmer)$coeff[1,1]),col="red")
abline(v=quantile(exp(test$samples[,20]),probs=c(0.025,0.5,0.975)),col=c("green","orange","green"),lty=2,lwd=2)
abline(v=mean(exp(test$samples[,20])),col="blue",lty=2)
# #posterior plot dispersion
# hist(test$samples[,133],pch=".",xlab="dispersion",main="please ignore this plot (only used in another parametrization)")
# #the variance in glmer (I suspect that is the dispersion, waiting for BB to confirm)
# abline(v=summary(crglmer)$sigma^2,col="red")
# abline(v=quantile(test$samples[,133],probs=c(0.025,0.5,0.975)),col=c("green","orange","green"),lty=2,lwd=2)
# abline(v=mean(test$samples[,133]),col="blue",lty=2)
#posterior plot location random effect standard deviation
hist(test$samples[,134],pch=".",xlab="location random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(crglmer)$varcor))[1],col="red")
abline(v=quantile(test$samples[,134],probs=c(0.025,0.5,0.975)),col=c("green","orange","green"),lty=2,lwd=2)
abline(v=mean(test$samples[,134]),col="blue",lty=2)
#posterior plot year random effect standard deviation
hist(test$samples[,135],pch=".",xlab="year random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(crglmer)$varcor))[2],col="red")
abline(v=quantile(test$samples[,135],probs=c(0.025,0.5,0.975)),col=c("green","orange","green"),lty=2,lwd=2)
abline(v=mean(test$samples[,135]),col="blue",lty=2)
legend("topright",legend=c("0.025 posterior quantile","posterior median","posterior mean","0.975 posterior quantile","glmer estimate"),
inset=0.05,lwd=2,col=c("green","orange","blue","green","red"),lty=c(2,2,2,2,1))

