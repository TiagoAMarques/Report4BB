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
load("data_4_article_clickrates_deep_dive.rda")
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


## define the model
GLMMcode <- nimbleCode({
  # the overall intercept
  beta0 ~ dnorm(0,sd=10)
  # random effect standard deviation associated with location, a uniform, might change this to be something else latter
  sigmal_RE ~ dunif(0, 1)
  # random effect standard deviation associated with year, a uniform, might change this to be something else latter
  sigmay_RE ~ dunif(0, 1)
  # the gamma dispersion/variance parameter, a uniform, might change this to be a gamma latter
  #dispersion ~ dunif(0, 10)
  disp ~ dunif(0, 100)
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
    #get the linear predictor, consider a log link function
    #log(mean[i]) <- beta0 + REy[year[i]] + REl[loc[i]]
    # Using decentered parametrization, a suggestion by Ben Augustine
    log(mean[i]) <- beta0 + REy[year[i]]*sigmay_RE + REl[loc[i]]*sigmal_RE
    # In bugs the parameterization  of the gamma is done using the shape and the rate
    # via mean and variance, defined as below:
    # note parametrized this way dispersion=variance! 
    #parametrization 1
    # crate[i] ~ dgamma(shape=(mean[i]^2)/disp,scale=disp/mean[i])
    #parametrization 2
    crate[i] ~ dgamma(shape=1/disp,scale=mean[i]*disp)
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
                       inits = inits, check = FALSE)

#things to monitor
#tomon<-c("beta0","dispersion","sigmay_RE","sigmal_RE","crate","REy","REl")
tomon<-c("beta0","disp","sigmay_RE","sigmal_RE","crate","REy","REl")

test<-nimbleMCMC(myGLMMModel,monitors=tomon,niter=50000,nburnin=10000,progressBar=TRUE,summary=TRUE)

#look at main results
test$summary[c(20,133,134,135),]

#the beta0 matches OK-ish with the intercept from glmer
summary(crglmer)$coeff[1,1];test$summary[20,1]

# Note that parametrized as 2 the dispersion parameter does not seem to bear a direct relation to the variance reported by glmer
test$summary[133,1];
summary(crglmer)$sigma^2
summary(crglmer)$sigma
sd(resid(crglmer,type='pear'))

# sqrt(sum(residuals(crglmer, type = "response")^2)/crglmer$df.residual)
# #Dispersion parameter estimation 
# #1.Deviance method
# (crglmer$deviance/crglmer$df.residual)  
# #2.Pearson method
# #resid(model,type='pear')
# sum(resid(crglmer,type='pear')^2)/crglmer$df.residual # estimation used in summary(model)


#but the variances of the random effects are still off
summary(crglmer)$varcor
(test$summary[134:135,1])

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
abline(v=quantile(test$samples[,20],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(test$samples[,20]),col="blue",lty=2)
#posterior plot dispersion
hist(test$samples[,133],pch=".",xlab="dispersion",main="")
#the variance in glmer (I suspect that is the dispersion, waiting for BB to confirm)
abline(v=summary(crglmer)$sigma^2,col="red")
abline(v=quantile(test$samples[,133],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(test$samples[,133]),col="blue",lty=2)
#posterior plot location random effect standard deviation
hist(test$samples[,134],pch=".",xlab="location random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(crglmer)$varcor))[1],col="red")
abline(v=quantile(test$samples[,134],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(test$samples[,134]),col="blue",lty=2)
#posterior plot year random effect standard deviation
hist(test$samples[,135],pch=".",xlab="year random effect sigma",main="")
abline(v=sqrt(as.numeric(summary(crglmer)$varcor))[2],col="red")
abline(v=quantile(test$samples[,135],probs=c(0.025,0.5,0.975)),col="green",lty=2)
abline(v=mean(test$samples[,135]),col="blue",lty=2)

