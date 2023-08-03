# GLM Gamma
# par(mfrow=c(1,3))
par(mfrow=c(1,2))

# parametrization 1 
set.seed(12345)
n<- 1000
minx<--3
maxx<-3
xs <- runif(n,minx,maxx)
meang <- exp(1.2+0.6*xs)
var<-0.3
disp<-0.8
#generate data defining the mean and the (constant) variance
ys1 <- rgamma(n,shape=(meang^2)/var,scale=var/meang)

plot(xs,ys1)
glmA1<-glm(ys1~xs,family=Gamma(link="log"))
summary(glmA1)
# note that the estimated dispersion parameter from the GLM 
# seems to correspond
# to the mean dispersion parameter across all estimated means
mean(var/meang^2);summary(glmA1)$dispersion
#and the variance is well estimated
var;sum(residuals(glmA1, type = "response")^2)/glmA1$df.residual

newxs<-seq(minx,maxx,length=200)    
prdys1<-predict(glmA1,newdata=data.frame(xs=newxs),type="response")
lines(newxs,prdys1,col="blue",lty=2,lwd=3)

# next parametrization is motivated by the statement of
# https://stats.stackexchange.com/questions/247624/dispersion-parameter-for-gamma-family
# "In R GLM assumes shape to be a constant"
# generate data defining the mean and fixing the shape 
# in this example induces increasing variance with mean
ys2 <- rgamma(n,shape=1/disp,scale=meang*disp)
plot(xs,ys2)
glmA2<-glm(ys2~xs,family=Gamma(link="log"))
summary(glmA2)
newxs<-seq(minx,maxx,length=200)    
prdys2<-predict(glmA2,newdata=data.frame(xs=newxs),type="response")
lines(newxs,prdys2,col="blue",lty=2,lwd=3)

# note that the estimated dispersion parameter from the GLM 
# corresponds to dispersion used
disp;summary(glmA2)$dispersion;gamma.dispersion(glmA2)
# but the variance is now no longer constant and depends
# on the mean, so I try the mean of the variances
# to compare with the GLM output - seems to add up
sum(residuals(glmA2, type = "response")^2)/glmA2$df.residual
mean((1/disp)*(meang*disp)^2)
# 
# # generate data defining the mean while inducing 
# # a relation with the variance which in fact means
# # that the paramter var above is not the variance LOL
# ys3 <- rgamma(n,shape=meang*var,scale=1/var)
# plot(xs,ys3)
# glmA3<-glm(ys3~xs,family=Gamma(link="log"))
# summary(glmA3)
# newxs<-seq(minx,maxx,length=200)    
# prdys3<-predict(glmA3,newdata=data.frame(xs=newxs),type="response")
# lines(newxs,prdys3,col="blue",lty=2,lwd=3)
# 
# # note that the estimated dispersion parameter from the GLM 
# # seems to correspond
# # to the mean dispersion parameter across all estimated means
# mean(1/(meang*var));summary(glmA3)$dispersion;gamma.dispersion(glmA3)
# #and the variance seems to be well estimated
# mean(meang/var);sum(residuals(glmA3, type = "response")^2)/glmA3$df.residual
