# sourcing this code will produce a quick simulation that shows, for a specific scenario
# that lme4 presents some issues at the very least estimating the standard deviation 
# of one of the random effects (associated with location) in a Gamma GLMM

#for GLMMs
library(lme4)
# for pipe operator
library(tidyverse)

# This header is a bit nully but it is eventually necessary to set up the data structure over which the issues 
# were originally detected. I could not be bothered identifying what, if something, in the data stucture
# was responsible for th eobserved bias
# (e.g. number of observations per level of each factor random effect, sample size, etc)
# --------------------------------------------------------------------------------------------------
# Nully code begins - 
# --------------------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------------------
# Nully code ends 
# --------------------------------------------------------------------------------------------------

set.seed(1221)
#define number of tags
ntags<-nrow(tags)
# define the 4 key parameters: the overall mean, the dispersion and the two random effect standar deviations
disp<-0.12
sigmaL<-0.2
sigmaY<-0.1
mean<-1
#a stupid way to get the right structure object
yearnames<-rownames(coef(crglmer)$fyear)
locnames <- rownames(coef(crglmer)$location)
B<-999
models<-vector("list",B)
tags.new<-tags
#object to store results
results<-data.frame(ID=1:B,emean=NA,edisp=NA,esdL=NA,esdY=NA)
for (i in 1:B){
  #get random effect values
  byyear <- rnorm(length(yearnames),mean=0,sd=sigmaY)
  byloc <-  rnorm(length(locnames),mean=0,sd=sigmaL)
  #get simulated data
  for(m in 1:ntags){
    tags.new$crate[m]<-exp(mean+byyear[which(yearnames==tags.new$fyear[m])]+byloc[which(locnames==tags.new$location[m])])
  }
  #notice bad coding practice as I am using the same object to hold first mean on the link and then on the response scale
  tags.new$crate<-rgamma(ntags,shape=1/disp,scale=tags.new$crate*disp)
  #fit model
  crglmerB<-glmer(crate~(1|location)+(1|fyear),data=tags.new,family=Gamma(link="log"))
  models[[i]]<-crglmerB
  results[i,2:5]<-c(summary(crglmerB)$coeff[1,1],
  summary(crglmerB)$sigma^2,
  sqrt(as.numeric(summary(crglmerB)$varcor))[1],
  sqrt(as.numeric(summary(crglmerB)$varcor))[2])
}

#inspect results
par(mfrow=c(2,2))
hist(results[,2],xlab="intecept",main="")
abline(v=mean,lty=2,col="green",lwd=2)
abline(v=mean(results[,2]),lty=2,col="red",lwd=2)
hist(results[,3],xlab="dispersion",main="")
abline(v=disp,lty=2,col="green",lwd=2)
abline(v=mean(results[,3]),lty=2,col="red",lwd=2)
hist(results[,4],xlab="location random effect sigma",main="")
abline(v=sigmaL,lty=2,col="green",lwd=2)
abline(v=mean(results[,4]),lty=2,col="red",lwd=2)
hist(results[,5],xlab="year random effect sigma",main="")
abline(v=sigmaY,lty=2,col="green",lwd=2)
abline(v=mean(results[,5]),lty=2,col="red",lwd=2)




