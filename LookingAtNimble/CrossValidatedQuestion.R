#My question is here
#https://stats.stackexchange.com/questions/623051/how-to-simulate-data-for-a-gamma-glm

# Related stack exchange posts
#https://stats.stackexchange.com/questions/484519/how-to-specify-gamma-parameterizations-in-a-generalized-linear-model-setting
#https://stats.stackexchange.com/questions/474326/deviance-for-gamma-glm
#https://stats.stackexchange.com/questions/136909/prior-for-gamma-distribution-in-mean-form/

Possible tweet if no feedback arrives in the next 48 hours:
  
In case you know about #rstast glm's Gamma regression parametrization &/or nimble code to implement gamma regression, just posted this question on stackexchange about how to simulate data for a gamma regression https://stats.stackexchange.com/q/623051/180421?stw=2 



Simulating data for a Gamma regression

I am wondering about whether there might actually be different ways to simulate data for say a Gamma GLM, which in turn relates to what might be the parametrization that the ```glm``` function uses for the ```Gamma``` family.
If one wants to test a Poisson GLM with say a log-link, the way to simulate the data is relatively easy as we only have a mean parameter. Assuming one has an independent variable x one can just simulate mean values of the response and the generate Poisson 
draws with that mean, say

```
#example code: does not work as is
meany<-exp(a+b*x)
ys<-rpois(n,meany)
```

However, in a Gamma (assume same log-link for easier comparison) one could come up with at least two ways in which the mean depends on the regressor x, but then one induces relations between the other parameters, say e.g.

```
#example code: does not work as is
meany<-exp(a+b*x)
# fix variance: dispersion varies with the mean
ys1 <- rgamma(n,shape=(meany^2)/var,scale=var/meany)
# fix dispersion: the variance varies with the mean
ys2<-rgamma(n,shape=1/disp,scale=disp*meany)
```

When one fits a Gamma ```glm``` the output produces a statement about the "dispersion" parameter 

```(Dispersion parameter for Gamma family taken to be x)```

which would hint for the second parametrization. However, on the other hand, when ones fits a Gamma GLMM in ```glmer``` from ```lme4``` the output refers only to ```variance``` components 
which could hint for the first parametrization.

Are there strong reasons to prefer one or the other and what is the one that ```glm``` (and ```glmer```) considers?
  
(in case this induces some additional feedback, the reason I am asking was induced by getting some strange results associatedd with the random effects in a NIMBLE implementation of a gamma GLMM which
 do not match what I got from ```glmer```, want to diagnose what might be the problem which led me to all sorts of rabbit holes, this Gamma (regression) parametrization being one of them. On that tangent, 
 I could not really find examples of NIMBLE/Bugs code for Gamma regression, so if anyone has or knows about them I'd appreciate a pointer to them)
 
A working code example generating data both ways and fitting a Gamma ```glm``` to them follows for completeness
 
```
# GLM Gamma
# parametrization 1 
set.seed(12345)
n<- 1000
minx<--3.2
maxx<-3.2
xs <- runif(n,minx,maxx)
meang <- exp(1.2+0.6*xs)
#define a variance
var<-2.3
#define a dispersion
disp<-0.8
#generate data defining the mean and the (constant) variance
ys1 <- rgamma(n,shape=(meang^2)/var,scale=var/meang)

glmA1<-glm(ys1~xs,family=Gamma(link="log"))
summary(glmA1)
newxs<-seq(minx,maxx,length=200)    
prdys1<-predict(glmA1,newdata=data.frame(xs=newxs),type="response")

# parametrization 2 
# next parametrization is motivated by the statement of
# https://stats.stackexchange.com/questions/247624/dispersion-parameter-for-gamma-family
# "In R GLM assumes shape to be a constant"
# generate data defining the mean and fixing the shape (= fixing dispersion)
# here induces an overall increase variance with mean
ys2 <- rgamma(n,shape=1/disp,scale=meang*disp)

glmA2<-glm(ys2~xs,family=Gamma(link="log"))
summary(glmA2)
prdys2<-predict(glmA2,newdata=data.frame(xs=newxs),type="response")

#plot data and fits
par(mfrow=c(1,2))
plot(xs,ys1)
lines(newxs,prdys1,col="blue",lty=2,lwd=3)
plot(xs,ys2)
lines(newxs,prdys2,col="blue",lty=2,lwd=3)
```
