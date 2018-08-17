#Clear memory
rm(list=ls())

library(spdep)
library(maptools)



#Polygon Shape File (.shp)
aladdin=readShapePoly("huc12_counts.shp")

df=data.frame(aladdin)

class(df)
names(df)
attach(df)

spplot(aladdin,"counts",main="counts")

col.nb1=poly2nb(aladdin)
col.nb2=poly2nb(aladdin,queen=FALSE)


coords=coordinates(aladdin)
summary(coords)
plot(aladdin)
plot(col.nb1,coords,add=TRUE,col=4)
###############
#permutation Test for Moran's I
colw <- nb2listw(col.nb1, style="W")
############
#nb2listw gives a certain neighborhood structure a particular weighting scheme.
nsim <- 99
set.seed(1234)
sim1 <- moran.mc(aladdin$counts, listw=colw, nsim=nsim)
sim1

shapiro.test(aladdin$Shape_Area)
moran.test(aladdin$Shape_Area, listw=colw, randomisation=FALSE)

############
#Construct the spatial connectivity matrix
n=2382
W=matrix(rep(0,len=n*n),nrow=n)
for(i in 1:n){
	slots=rep(0,len=n)
	slots[colw$neighbours[[i]]]=colw$weights[[i]]
	W[i,]=slots
}

S0=sum(W)
S1=0
for(i in 1:n){
	for(j in 1:n){
		S1=S1+0.5*(W[i,j]+W[j,i])^2}}

wi.dot=apply(W,2,sum)
wj.dot=apply(W,1,sum)
S2=sum((wi.dot+wj.dot)^2)

varI= (n^2*S1-n*S2+3*S0^2)/((n-1)*(n+1)*S0^2)
EI=-1/(n-1)

TS=(0.485770914-EI)/sqrt(varI)
TS
1-pnorm(TS)
####################
####################
n=2382
#Local Moran's I
loc.sim=localmoran(aladdin$counts, listw=colw)
loc.sim
#The first column is the local test statistic, then its expected value under the null, then its variance under normality, then its standardized Z-score, and then its p-value.

loc.sim[,5]
pvals=rep(1,n)
alpha=0.05
pvals[which(loc.sim[,5]<alpha)]=rep(0, length(which(loc.sim[,5]<alpha)))

####################

aladdin2=aladdin
aladdin2@data=cbind(aladdin@data,loc.sim[,1],loc.sim[,5],pvals)
spplot(aladdin2,"loc.sim[, 1]",main="Local Moran's I")
spplot(aladdin2,"loc.sim[, 5]",main="Local Moran's I P-Values")
spplot(aladdin2,"pvals",main="Channel Belt A1 Significant Local Moran's I P-Values",scales=list(draw=TRUE))

plot(aladdin$x,loc.sim[,1])
##############
##############
#joint count for categorical data
##############
BT<-cut(aladdin$sorting,breaks=c(-2,1,2,3,4,8),labels=c("1","2","3","4","5"))

joincount.test(BT,nb2listw(col.nb1,style="B"))
set.seed(1234)
joincount.mc(BT,nb2listw(col.nb1,style="B"),nsim=99)






###############
#Getting different adjacency structures
colw <- nb2listw(col.nb1, style="W")
colb <- nb2listw(col.nb1, style="B")
colu <- nb2listw(col.nb1, style="U")
####################
#Fitting models
####################
set.seed(223)

modelresults<-spautolm(aladdin$uvel~aladdin$sorting+aladdin$VERT_POS,listw=colu,family="SAR")
summary(modelresults)
shapiro.test(modelresults$fit$residuals)
moran.test(modelresults$fit$residuals,listw=colw)
set.seed(223)
moran.mc(modelresults$fit$residuals,listw=colw ,nsim=1000)


####################
#Some Extra Plots
####################
#Putting new columns in the object
class(df)
slotNames(df)
dim(aladdin@data)
aladdin@data=cbind(aladdin@data, modelresults$fit$fitted.values, modelresults$fit$residuals, modelresults$fit$signal_trend, modelresults$fit$signal_stochastic)

fitted.vals=modelresults$fit$signal_trend+modelresults$fit$signal_stochastic
fitted.vals-modelresults$fit$fitted.values

spplot(aladdin,"uvel",main="A1 Maximum Grain Size",scales=list(draw=TRUE),xlab="Lateral Distance (m)",ylab="Vertical Distance (m)")

spplot(aladdin,"modelresults$fit$fitted.values",main="A1 Model Fitted = Trend + Spatial",scales=list(draw=TRUE),xlab="Lateral Distance (m)",ylab="Vertical Distance (m)")

spplot(aladdin,"modelresults$fit$residuals",main="A1 Model Residuals",scales=list(draw=TRUE),xlab="Lateral Distance (m)",ylab="Vertical Distance (m)") #Should be iid and normal

spplot(aladdin,"modelresults$fit$signal_trend",main="A1 Model Trend",scales=list(draw=TRUE),xlab="Lateral Distance (m)",ylab="Vertical Distance (m)")

spplot(aladdin,"modelresults$fit$signal_stochastic",main="Model Spatial Component",scales=list(draw=TRUE),xlab="Lateral Distance (m)",ylab="Vertical Distance (m)")

moran.mc(modelresults$fit$signal_stochastic,listw=colw ,nsim=1000) #Should be significant and it is




###########
#spatial lags correlogram
cor8<-sp.correlogram(col.nb1,aladdin$uvel,order=7,method="I",style="W",randomisation=TRUE,zero.policy=TRUE)
plot(cor8,main="Spatial Autocorrelation")

cor8









