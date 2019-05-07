load("../Writing/DensEst.Rdata")

library(lme4)
library(dplyr)
library(ggeffects)
library(readxl)
library(sp)

source('ZebraFuncs.R')

#######################Encounter time and setup time
#time.df <- rbind(LB.distance.est$df, LB.quadrat.est$df, LBL.distance.est$df, LBL.quadrat.est.subset$df, LF.distance.est$df, LF.quadrat.est$df)
time.df <- rbind(LB.distance.est$df, LB.double.est$df, LB.quadrat.est$df, LBL.distance.est$df ,LBL.double.est$df, LBL.quadrat.est.subset$df, LF.distance.est$df, LF.double.est$df, LF.quadrat.est$df)

#contrasts(lm$Lake) <- contr.treatment(3)
#tset.length <- lm(t.set ~ Length + Type + (1|Lake), data=time.df)

tset.length <- lm(t.set ~ Length + Type + Lake, data=time.df,  contrasts=list(Lake="contr.helmert"))
tset.length <- lm(t.set ~ Length + Type, data=time.df)

#tset.length <- lm(t.set ~ Length + Type, data=time.df)
#lm(y ~ f, contrasts = list(f = "contr.helmert"))
detection.rate <- time.df$Detections/time.df$Length
#tenc.length    <- lmer(t.enc ~ Length + detection.rate*Type + (1|Lake), data=time.df)

contrasts(time.df$Lake) <- "contr.helmert"
#tenc.length2   <- lm(t.enc ~ Length + detection.rate*Type + Lake, data=time.df)
#tenc.length    <- lm(t.enc ~ Length + detection.rate*Type + Lake, data=time.df)
tenc.length    <- lm(t.enc ~ Length + detection.rate*Type + Lake, data=time.df)
tenc.length    <- lm(t.enc ~ Length + detection.rate*Type, data=time.df)

#tenc.length <- lm(t.set ~ Length + detection.rate + Type + Lake, data=time.df, contrasts=list(Lake="contr.helmert"))

set.pred  <- ggpredict(tset.length, terms = c("Length", "Type"))
set.pred  <- mutate(set.pred, Type="Setup time")
enc.pred  <- ggpredict(tenc.length, terms = c("Length", "Type"))
enc.pred  <- mutate(enc.pred, Type="Encounter time")
length.pred <- rbind(set.pred, enc.pred)
rate.pred   <- ggpredict(tenc.length, terms = c("detection.rate", "Type"))



############Travel time model

phase2.dat      <- read_xlsx(path="../Data/Season2/Zebra mussel survey_ Day 1 (Responses).xlsx", sheet=1)
phase2.dat <- phase2.dat %>% subset(Observer=="Aislyn") %>% mutate(Lake = as.factor(Lake), TravelTime=NA, TravelDist=NA)

for(i in 1:nlevels(phase2.dat$Lake)) {
  
  phase2.subset <- subset(phase2.dat, Lake == levels(Lake)[i])
  indices       <- which(phase2.dat$Lake == levels(phase2.dat$Lake)[i])
  
  dist.vec    <- time.vec <- vector('numeric', dim(phase2.subset)[1])
  dist.vec[1] <- time.vec[1] <- NA
  
  for(j in 2:(dim(phase2.subset)[1])) {
    dist.vec[j] <- spDists(cbind(phase2.subset$`GPS start latitude`[(j-1):j], phase2.subset$`GPS start longitude`[(j-1):j]), longlat=TRUE)[1,2]*1000
    time.vec[j] <- difftime(phase2.subset$`Start time`[j], phase2.subset$`End time`[j-1], units="secs")
    if(!is.na(phase2.subset$`Break time (duration in minutes)`[j])) {
      time.vec[j] - phase2.subset$`Break time (duration in minutes)`[j]*60
    }
  }
  
  phase2.dat$TravelDist[indices] <- dist.vec
  phase2.dat$TravelTime[indices] <- time.vec
  
}

#travel.lm  <- lmer(TravelTime ~ TravelDist + (1|Lake), data=phase2.dat)
travel.lm  <- lm(TravelTime ~ TravelDist * Lake, data=phase2.dat, contrasts = list(Lake = "contr.helmert"))

###predict travel time at each lake based on distnaces between sites


LB.dist       <- spDists(cbind(subset(phase2.dat, Lake=="Burgen")$`GPS start latitude`, subset(phase2.dat, Lake=="Burgen")$`GPS start longitude`), longlat=TRUE)*1000
LBL.dist      <- spDists(cbind(subset(phase2.dat, Lake=="Little Birch")$`GPS start latitude`, subset(phase2.dat, Lake=="Little Birch")$`GPS start longitude`), longlat=TRUE)*1000
LF.dist       <- spDists(cbind(subset(phase2.dat, Lake=="Florida")$`GPS start latitude`, subset(phase2.dat, Lake=="Florida")$`GPS start longitude`), longlat=TRUE)*1000

#Predict.traveltime2  <- predict(travel.lm, newdata=data.frame(TravelDist=0:100, Lake=rep(phase2.dat$Lake[1], 101)), re.form=NA)
##newdat <- data.frame(Intercept=rep(1, 101), TravelDist=0:100)
#coef(travel.lm)[1:2]
##Predict.traveltime  <- as.matrix(newdat)%*%(coef(travel.lm)[1:2])

#Predict.traveltime  <- predict(travel.lm, newdata=data.frame(TravelDist=0:100, Lake=rep("Test", 101)), re.form=NA)

LF.transectDist <- subSample.dist(LF.dist)
LF.distMean <- NULL
for(i in unique(LF.transectDist$SampleSize)) {
  x <- i
  LF.distMean <- c(LF.distMean, mean(LF.transectDist$Distance[LF.transectDist$SampleSize==i]))
}
names(LF.distMean) <- unique(LF.transectDist$SampleSize)
newdat <- as.matrix(data.frame(Intercept=1, TravelDist=LF.transectDist$Distance))
LF.transectTime <- newdat%*%(coef(travel.lm)[1:2])
  
LB.transectDist <- subSample.dist(LB.dist)
LB.distMean <- NULL
for(i in unique(LB.transectDist$SampleSize)) {
  x <- i
  LB.distMean <- c(LB.distMean, mean(LB.transectDist$Distance[LB.transectDist$SampleSize==i]))
}
names(LB.distMean) <- unique(LB.transectDist$SampleSize)
#LB.transectTime <- predict(travel.lm, newdata=data.frame(TravelDist=LF.transectDist$Distance), re.form=NA)
newdat <- as.matrix(data.frame(Intercept=1, TravelDist=LB.transectDist$Distance))
LB.transectTime <- newdat%*%(coef(travel.lm)[1:2])

LBL.transectDist <- subSample.dist(LBL.dist)
LBL.distMean <- NULL
for(i in unique(LBL.transectDist$SampleSize)) {
  x <- i
  LBL.distMean <- c(LBL.distMean, mean(LBL.transectDist$Distance[LBL.transectDist$SampleSize==i]))
}
names(LBL.distMean) <- unique(LBL.transectDist$SampleSize)

#LBL.transectTime <- predict(travel.lm, newdata=data.frame(TravelDist=LBL.transectDist$Distance), re.form=NA)
newdat <- as.matrix(data.frame(Intercept=1, TravelDist=LBL.transectDist$Distance))
LBL.transectTime <- newdat%*%(coef(travel.lm)[1:2])

TimeDist.frame   <- data.frame(Time=c(LF.transectTime, LB.transectTime, LBL.transectTime), Distance=c(LF.transectDist$Distance, LB.transectDist$Distance, LBL.transectDist$Distance), N=c(LF.transectDist$SampleSize, LB.transectDist$SampleSize, LBL.transectDist$SampleSize), Lake=c(rep("Lake Florida", 77), rep("Lake Burgan", 77), rep("Little Birch Lake", 77))) 

DistN.fit <- lm(log(Distance) ~ log(N) + Lake, data=TimeDist.frame)

nvec.cont    <- seq(2,16,length.out=100)
LF.predict   <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake="Lake Florida")))))%*%(coef(travel.lm)[1:2])
LB.predict   <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake="Lake Burgan")))))%*%(coef(travel.lm)[1:2])
LBL.predict  <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake="Little Birch Lake")))))%*%(coef(travel.lm)[1:2])

#LF.predict  <- predict(travel.lm, newdata=data.frame(TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake=rep(TimeDist.frame$Lake[1], length(nvec.cont)))))), re.form=NA)
#LB.predict  <- predict(travel.lm, newdata=data.frame(TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake=rep(TimeDist.frame$Lake[100], length(nvec.cont)))))), re.form=NA)
#LBL.predict <- predict(travel.lm, newdata=data.frame(TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec.cont, Lake=rep(TimeDist.frame$Lake[200], length(nvec.cont)))))), re.form=NA)

#######Combine the parts to get the total amount of time spent per lake.

nvec       <- 2:16
travel.Time <- matrix(NA, 3, length(nvec))

newDat <- as.matrix(data.frame(Intercept=1, N=nvec))

travel.Time[1,] <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec, Lake="Lake Florida")))))%*%(coef(travel.lm)[1:2])

travel.Time[2,] <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec, Lake="Lake Burgan")))))%*%(coef(travel.lm)[1:2])

travel.Time[3,] <- as.matrix(data.frame(Intercept=1, TravelDist=exp(predict(DistN.fit, newdata=data.frame(N=nvec, Lake="Little Birch Lake")))))%*%(coef(travel.lm)[1:2])


newDat <- data.frame(Intercept=1, Length=30, Type=rep(c("Distance", "Quadrat"), 3))
newDat$Type <- as.numeric(newDat$Type)-1
setup.Time  <- as.matrix(newDat)%*%(coef(tset.length)[1:3])


dr.dist.vec <- c(sum(LF.distance.est$Detections), sum(LB.distance.est$Detections), sum(LBL.distance.est$Detections))/c(sum(LB.distance.est$df$Length), sum(LB.distance.est$df$Length), sum(LBL.distance.est$df$Length))
dr.double.vec <- c(sum(LF.double.est$Detections), sum(LB.double.est$Detections), sum(LBL.double.est$Detections))/c(sum(LB.double.est$df$Length), sum(LB.double.est$df$Length), sum(LBL.double.est$df$Length))
dr.quad.vec <- c(sum(LF.quadrat.est$Mussels), sum(LB.quadrat.est$Mussels), sum(LBL.quadrat.est.subset$Mussels))/c(sum(LF.quadrat.est$df$Length), sum(LB.quadrat.est$df$Length), sum(LBL.quadrat.est.subset$df$Length))



encDat <- data.frame(Intercept=1, Length=30, detection.rate=c(dr.double.vec[1], dr.quad.vec[1], dr.double.vec[2], dr.quad.vec[2], dr.double.vec[3], dr.quad.vec[3]), Type=rep(c("Double", "Quadrat"), 3))
encDat$Type <- as.numeric(encDat$Type)-1
encDat      <- encDat %>% mutate(Interaction=detection.rate*encDat$Type)
#search.Time <- as.matrix(encDat)%*%(coef(tenc.length)[c(1,2,3,4,7)])

#search.Time <- as.matrix(encDat)%*%coef(tenc.length)[c(1,2,3,4,7)]
#search.Time <- as.matrix(encDat)[,1:4]%*%coef(tenc.length)[c(1,2,3,4)]


nbind <- rbind(Distance=nvec, Quadrat=nvec)
total.Time.LF  <- rbind(travel.Time[1,]*nvec, travel.Time[1,]*nvec) + nbind*search.Time[c(1,2)] + nbind*setup.Time[c(1,2)]
total.Time.LB  <- rbind(travel.Time[2,]*nvec, travel.Time[2,]*nvec) + nbind*search.Time[c(3,4)] + nbind*setup.Time[c(3,4)]
total.Time.LBL <- rbind(travel.Time[3,]*nvec, travel.Time[3,]*nvec) + nbind*search.Time[c(5,6)] + nbind*setup.Time[c(5,6)]

totalTime.df <- rbind(data.frame(Time=c(total.Time.LF[1,], total.Time.LF[2,]), N=c(nvec, nvec), Design=c(rep("Distance survey", length(nvec)), rep("Quadrat survey", length(nvec))), Lake=rep("Lake Florida", 2*length(nvec))), 
                      data.frame(Time=c(total.Time.LB[1,], total.Time.LB[2,]), N=c(nvec, nvec), Design=c(rep("Distance survey", length(nvec)), rep("Quadrat survey", length(nvec))), Lake=rep("Lake Burgan", 2*length(nvec))), 
                      data.frame(Time=c(total.Time.LBL[1,], total.Time.LBL[2,]), N=c(nvec, nvec), Design=c(rep("Distance survey", length(nvec)), rep("Quadrat survey", length(nvec))), Lake=rep("Little Birch Lake", 2*length(nvec))))

save.image("../Writing/TimeBudgetEst.Rdata")
