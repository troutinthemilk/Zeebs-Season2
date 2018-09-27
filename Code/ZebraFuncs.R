#set of functions that are useful for analyzing the zebra dataset

twopiece.lik <- function(par.vec, y) {
  
  mu    <- exp(par.vec[1])
  sig1  <- exp(par.vec[2])
  sig2  <- exp(par.vec[3])
  
  ind1 <- y < mu
  ind2 <- y >= mu
  
  integrand <- function(x, mu, sig) { exp(-(x - mu)^2/(2*sig^2)) }
  int1  <- integrate(f=integrand, lower=0, upper=mu, mu=mu, sig=sig1)$value
  int2  <- integrate(f=integrand, lower=mu, upper=1, mu=mu, sig=sig2)$value
  
  ll <- sum(-(y - mu)^2/(sqrt(2)*(ind1*sig1 + ind2*sig2))^2) - length(y)*log(int1+int2) #- sum(ind1*log(int1)) - sum(ind2*log(int2))
  
  return(-ll)
}

twopiece.lik.same <- function(par.vec, y) {
  
  mu    <- exp(par.vec[1])
  sig1  <- exp(par.vec[2])
  sig2  <- sig1
  
  ind1 <- y < mu
  ind2 <- y >= mu

  integrand <- function(x, mu, sig) { exp(-(x - mu)^2/(2*sig^2)) }
  
  int1  <- integrate(f=integrand, lower=0, upper=mu, mu=mu, sig=sig1)$value
  int2  <- integrate(f=integrand, lower=mu, upper=1, mu=mu, sig=sig2)$value
  
  ll <- sum(-(y - mu)^2/(sqrt(2)*(ind1*sig1 + ind2*sig2))^2) - length(y)*log(int1+int2) #sum(ind1*log(int1+int2)) - sum(ind2*log(int2+int1))
  
  return(-ll)
}



onepiece.lik <- function(x, y) {
  
  sig1  <- x
  
  integrand <- function(x, mu, sig) { exp(-x^2/(2*sig^2)) }
  int1      <- integrate(f=integrand, lower=0, upper=1, mu=0, sig=sig1)$value
 
  ll <- sum(-y^2/(sqrt(2)*sig1)^2) - length(y)*log(int1)
  return(-ll)
}


#this interpolates the missing lcoation data for Lake Sylvia by averaging locations of the surrounding transects
interpolate <- function(Trans) {
  
  Sylvia.trans  <- filter(Trans, Region.Label=="Sylvia")
  fillins       <- which(is.na(Sylvia.trans$gps_transect_n))
  nextthree     <- (length(fillins)+1):(length(fillins)+3)
  nextnextthree <- (length(fillins)+4):(length(fillins)+6)
  
  avgN <- mean(Sylvia.trans[c(nextthree, nextnextthree),]$gps_transect_n)
  avgW <- mean(Sylvia.trans[c(nextthree, nextnextthree),]$gps_transect_w)
  bearing <- mean(Sylvia.trans[c(nextthree, nextnextthree),]$bearing)
  
  xy <- data.frame(ID = 1, X = c(avgW), Y = c(avgN))
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS("+proj=utm +zone=15 ellps=WGS84"))
  
  newx <- coordinates(res)[1] + cos(-pi/180*(bearing))*c(-1,-1, -1, -1, 1, 1, 1, 1)*c(1.5, 4.5, 7.5, 10.5, 1.5, 4.5, 7.5, 10.5) 
  newy <- coordinates(res)[2] + sin(-pi/180*(bearing))*c(-1,-1, -1, -1, 1, 1, 1, 1)*c(1.5, 4.5, 7.5, 10.5, 1.5, 4.5, 7.5, 10.5) #just guessin on the order
  
  xynew               <- data.frame(ID = 1:8, X=newx, Y=newy)
  coordinates(xynew)  <- c("X", "Y")
  proj4string(xynew)  <- CRS("+proj=utm +zone=15 ellps=WGS84")  ## for example
  xytrans             <- spTransform(xynew, CRS("+proj=longlat +datum=WGS84"))
  
  Sylvia.trans[fillins,]$gps_transect_n       <- coordinates(xytrans)[,2]
  Sylvia.trans[fillins,]$gps_transect_w       <- coordinates(xytrans)[,1]
  Sylvia.trans[fillins,]$bearing              <- bearing
  Trans[which(Trans$Region.Label=="Sylvia"),] <- Sylvia.trans
  
  return(Trans)
  
}

logit.rand <- function(pars, varcov, x) {
  parDraw <- rmvn(1, pars, varcov)
  plogis(sum(parDraw*x)) 
}


#based on the location of the transect and the distance on the transect that a mussel was discovered, this function calculated the gps coordinates of the detection event.
detection.location <- function(HabZM_entry) {
  gps_detect_w <- gps_detect_n <- NULL
  
  if(is.na(HabZM_entry$gps_transect_w)) {return(c(NA,NA))}
  xy <- data.frame(ID = 1, X = HabZM_entry$gps_transect_w, Y = HabZM_entry$gps_transect_n)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS("+proj=utm +zone=15 ellps=WGS84"))
  
  newx <- coordinates(res)[1] + sin(HabZM_entry$bearing*pi/180)*HabZM_entry$mussel_distance_along + 2*(HabZM_entry$side-1.5)*HabZM_entry$distance*sin((HabZM_entry$bearing-90)*pi/180)
  newy <- coordinates(res)[2] + cos(HabZM_entry$bearing*pi/180)*HabZM_entry$mussel_distance_along + 2*(HabZM_entry$side-1.5)*HabZM_entry$distance*cos((HabZM_entry$bearing-90)*pi/180)
  
  xynew <- data.frame(ID = 1, X=newx, Y=newy)
  coordinates(xynew) <- c("X", "Y")
  proj4string(xynew) <- CRS("+proj=utm +zone=15 ellps=WGS84")  ## for example
  xytrans <- spTransform(xynew, CRS("+proj=longlat +datum=WGS84"))
  
  gps_detect_w <- coordinates(xytrans)[1]
  gps_detect_n <- coordinates(xytrans)[2]
  
  return(c(gps_detect_w, gps_detect_n))
}


#links together habitat and zebra mussel count data so that it is clear what the habitat is at each detection event 
Link.HabCounts <- function(Hab, ZM) {
  
  #clean hab datasets
  Hab <- Hab %>%rename(hab_id=record_id) %>% mutate(record_id = as.factor(str_split(hab_id, pattern="_D", n=2, simplify = TRUE)[,1])) %>% mutate(Observer=str_detect(hab_id, "_2")+1) 
  Hab <- Hab[c(1:61, 63, 62, 64:dim(Hab)[1]), ] 
  
  Hab1  <- Hab %>% filter(Observer==1)
  HabZM <- ZM
  
  for(i in 1:dim(ZM)[1]) {
    
    currRec <- as.character(ZM[i,]$record_id)
    currDis <- ZM[i,]$mussel_distance_along
    currSet <- which(as.character(Hab1$record_id) == currRec)
    
    #uses findinterval to make sure that the current observation is 
    HabZM[i, 56:62]       <- Hab1[findInterval(x=currDis, vec=Hab1[currSet,]$distance_along) + 1, 4:10]
    HabZM[i,]$plantcover  <- Hab1$plant_cover[findInterval(x=currDis, vec=Hab1[currSet,]$distance_along) + 1]
    HabZM$Effort[i]       <- Hab1[currSet[length(currSet)],]$distance_along #multiply the transect length by 2 to get the sampled area. This is the cumulative transect length.
    
  }
  
  #insert missing effort
  #HabZM <- effort.missing(ZM, HabZM)
  
  #get survey areas
  temp <- vector("numeric", 24)
  for(i in 1:24) {
    temp[i] <- filter(HabZM, Region.Label == "Sylvia" & transect_no == i)$Effort[1]
  }
  HabZM$Area[which(HabZM$Region.Label == "Sylvia")] <- 2*sum(temp)
  
  temp <- vector("numeric", 18)
  for(i in 1:18) {
    temp[i] <- filter(HabZM, Region.Label == "Burgen"& transect_no == i)$Effort[1]
  }
  HabZM$Area[which(HabZM$Region.Label == "Burgen")] <- 2*sum(temp)
  
  #Turn cover into factor and rename mussel substrate columns.
  HabZM$plantcover <- as.factor(HabZM$plantcover)
  HabZM$plantcover <- recode(HabZM$plantcover, '0'="Absent", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  HabZM <- HabZM %>% rename(Mud=substrate1, Sand=substrate2, Gravel=substrate3, Pebble=substrate4, Rock=substrate5, Silt=substrate6, Other=substrate7)
  
  #treat these variables as factors
  HabZM$Mud     <- as.factor(HabZM$Mud)
  HabZM$Sand    <- as.factor(HabZM$Sand)
  HabZM$Gravel  <- as.factor(HabZM$Gravel)
  HabZM$Pebble  <- as.factor(HabZM$Pebble)
  HabZM$Rock    <- as.factor(HabZM$Rock)
  HabZM$Silt    <- as.factor(HabZM$Silt)
  HabZM$Other   <- as.factor(HabZM$Other)
  
  HabZM <- HabZM %>% rename(Sample.Label = 'record_id')
  
  returnList <- list()
  returnList[[1]] <- HabZM
  returnList[[2]] <- Hab
  
  return(returnList)

}


#insert missing effort into Sylvia Lake.
#this data was not properly entered into the database for some reason
effort.missing <- function(ZM, HabZM) {
  
  sylvia.missing <- setdiff(1:24, filter(ZM, Region.Label=="Sylvia")$transect_no)
  burgen.missing <- setdiff(1:18, filter(ZM, Region.Label=="Burgen")$transect_no)
  sylvia.missing.distance <- c(15, 10, 15, 15, 15, 15, 15, 27, 15, 18.5)
  burgen.missing.distance <- c(30)
  
  HabZM <- HabZM %>% add_row(record_id="21004900_T10", Region.Label="Burgen", observer___1=1, observer___2=1, observer___6=0, observer___7=0, transect_no=burgen.missing, Effort=burgen.missing.distance)
  
  
  
  for(i in 1:length(sylvia.missing)) {
    HabZM <- HabZM %>% add_row(record_id=paste("73024900_T", sylvia.missing[i], sep=''), Region.Label="Sylvia", observer___1=1, observer___2=1, observer___6=0, observer___7=0, transect_no=sylvia.missing[i], Effort=sylvia.missing.distance[i])
  }
  
  return(HabZM)
  
}

#This function classifies two detections from different observers as either being the same group or not. 
#horDis is the distance away the transect that we will use as criteria for grouping observations
#vertDis is the distance along the transect that we will use as criteria for grouping observations 
#both the actual horizont distance must be less than horDis and the actual vertical distance must be less than vertDis
match.function <- function(BurgenDat, horDis=0.1, vertDis=0.25) {
  BurgenDat <- BurgenDat %>% mutate(nearest=rep(NA,dim(BurgenDat)[1]), horDis=rep(NA,dim(BurgenDat)[1]), vertDis=rep(NA,dim(BurgenDat)[1]))
  
  for(i in 1:dim(BurgenDat)[1]) {
    
    currTrans <- BurgenDat$transect_no[i]
    currObs   <- which(BurgenDat$transect_no == BurgenDat$transect_no[i])
    currLabel <- BurgenDat[i,]$observer
    
    if(length(currObs) == 1) {
      BurgenDat$nearest[i] <- NA
      BurgenDat$horDis[i] <- NA
      BurgenDat$vertDis[i] <- NA
      
      next()
    }
    
    obsOtherDat <- filter(BurgenDat, transect_no==currTrans, observer==switch(currLabel, 2, 1))
    if(dim(obsOtherDat)[1] == 0) { next() }
    
    tmp <- crossdist(X=BurgenDat$distance[i]*2*(BurgenDat$side[i]-1.5), Y=BurgenDat$mussel_distance_along[i], x2=obsOtherDat$distance*2*(obsOtherDat$side-1.5), y2=obsOtherDat$mussel_distance_along)
    min.index             <- which.min(tmp[1,])
    BurgenDat$nearest[i]  <- min(tmp[1,])
    
    #get minimimum distances
    BurgenDat$horDis[i] <- abs(obsOtherDat[min.index,]$distance*2*(obsOtherDat[min.index,]$side-1.5) - BurgenDat$distance[i]*2*(BurgenDat$side[i]-1.5))
    BurgenDat$vertDis[i] <- abs(BurgenDat$mussel_distance_along[i] - obsOtherDat[min.index,]$mussel_distance_along)
    
  }
  
  #I am just defining a value that determines how close two sightings by different observers are that need to be clustered. Working on a new method for this.
  BurgenDat <- BurgenDat %>% mutate(object=rep(NA, dim(BurgenDat)[1]), detected=rep(NA, dim(BurgenDat)[1]))
  if(any(is.na(BurgenDat$distance))) {
    BurgenDat <- BurgenDat[-which(is.na(BurgenDat$distance)),]
  }
  BurgenNew <- NULL
  
  #This section looks at each pair of points seen by the different observers and determines whether these detection events are the same mussel or not.
 
  object.count <- 1
  for(i in unique(BurgenDat$transect_no)) {
    
    tempDat   <- filter(BurgenDat, transect_no == i)
    
    if(dim(tempDat)[1] == 1) {
      
      tempDat$object                            <- object.count
      BurgenNew                                 <- rbind(BurgenNew, tempDat[1,])
      BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
      
      BurgenNew                               <- rbind(BurgenNew, tempDat[1,])
      BurgenNew[dim(BurgenNew)[1],]$object    <- object.count
      BurgenNew[dim(BurgenNew)[1],]$observer  <- switch(tempDat$observer, 2, 1)
      BurgenNew[dim(BurgenNew)[1],]$detected  <- 0
      
      object.count    <- object.count + 1 
    } else {
      for(k in 1:dim(tempDat)[1]) {
        if(!is.na(tempDat[k,]$object)) { next() }
        
        #this is the last observation
        if(k == dim(tempDat)[1]) {
          #print('last one')
          BurgenNew                                 <- rbind(BurgenNew, tempDat[k,])
          BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
          BurgenNew[dim(BurgenNew)[1],]$object      <- object.count
          
          BurgenNew                               <- rbind(BurgenNew, tempDat[k,])
          BurgenNew[dim(BurgenNew)[1],]$object    <- object.count
          BurgenNew[dim(BurgenNew)[1],]$observer  <- switch(tempDat[k,]$observer, 2, 1)
          BurgenNew[dim(BurgenNew)[1],]$detected  <- 0
          
          tempDat[k,]$object  <- object.count
          object.count        <- object.count + 1
          next()
        }
        
        tmp <- crossdist(X=tempDat$distance[k]*2*(tempDat$side[k]-1.5), Y=tempDat$mussel_distance_along[k], x2=tempDat$distance*2*(tempDat$side-1.5), y2=tempDat$mussel_distance_along)
        
        tmpHor <- crossdist(X=tempDat$distance[k]*2*(tempDat$side[k]-1.5), Y=0, x2=tempDat$distance*2*(tempDat$side-1.5), y2=rep(0, length(tempDat$mussel_distance_along)))
        
        tmpVer <- crossdist(X=0, Y=tempDat$mussel_distance_along[k], x2=rep(0, length(tempDat$distance*2*(tempDat$side-1.5))), y2=tempDat$mussel_distance_along)
        
        tmp[k]        <- NA
        tmpHor[k]     <- NA
        tmpVer[k]     <- NA
        
        nearest           <- min(tmp[1, (k+1):dim(tmp)[2]], na.rm=T)
        which.nearest     <- which.min(tmp[1, (k+1):dim(tmp)[2]])
        nearestHor        <- min(tmpHor[1,(k+1):dim(tmp)[2]], na.rm=T)
        which.nearestHor  <- which.min(tmpHor[1, (k+1):dim(tmp)[2]])
        nearestVer        <- min(tmpVer[1, (k+1):dim(tmp)[2]], na.rm=T)
        which.nearestVer  <- which.min(tmpVer[1, (k+1):dim(tmp)[2]])
        
        if(is.na(tempDat[k,]$nearest) & nearestHor > horDis | nearestVer > vertDis) { #there is no nearest neighbor, and this detection has not been already counted
          BurgenNew                                 <- rbind(BurgenNew, tempDat[k,])
          BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
          BurgenNew[dim(BurgenNew)[1],]$object      <- object.count
          
          BurgenNew                               <- rbind(BurgenNew, tempDat[k,])
          BurgenNew[dim(BurgenNew)[1],]$object    <- object.count
          BurgenNew[dim(BurgenNew)[1],]$observer  <- switch(tempDat[k,]$observer, 2, 1)
          BurgenNew[dim(BurgenNew)[1],]$detected  <- 0
          
          tempDat[k,]$object  <- object.count
          object.count        <- object.count + 1
          
        } else { #there is a nearest neighbor and it doesn't yet have a partner
          
          #check that all observations within the criteria are observed by different observers
          valid.match <- which(tmpHor[1, (k+1):dim(tmpHor)[2]] <= horDis & tmpVer[1, (k+1):dim(tmpVer)[2]] <= vertDis)
          valid.match <- valid.match[which(tempDat$observer[k+valid.match] != tempDat$observer[k])]
          if(length(valid.match) == 0) { #there are no valid nearest neighbors
            BurgenNew                                 <- rbind(BurgenNew, tempDat[k,])
            BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
            BurgenNew[dim(BurgenNew)[1],]$object      <- object.count
            
            BurgenNew                               <- rbind(BurgenNew, tempDat[k,])
            BurgenNew[dim(BurgenNew)[1],]$object    <- object.count
            BurgenNew[dim(BurgenNew)[1],]$observer  <- switch(tempDat[k,]$observer, 2, 1)
            BurgenNew[dim(BurgenNew)[1],]$detected  <- 0
            
            tempDat[k,]$object  <- object.count
            object.count        <- object.count + 1
            next()
          }
          #print('match made')
          #assign the nearest neighbor of the above set as the partner of the observation
          best.match <- valid.match[which.min(tempDat$nearest[k + valid.match])]
          BurgenNew <- rbind(BurgenNew, tempDat[k,])
          BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
          BurgenNew[dim(BurgenNew)[1],]$object      <- object.count
          
          BurgenNew                                 <- rbind(BurgenNew, tempDat[k,])
          
          BurgenNew[dim(BurgenNew)[1],]$detected    <- 1
          BurgenNew[dim(BurgenNew)[1],]$object      <- object.count
          BurgenNew[dim(BurgenNew)[1],]$observer    <- switch(tempDat[k,]$observer, 2, 1) 
          tempDat[k,]$object            <- object.count
          tempDat[k+best.match,]$object <- object.count
          
          object.count <- object.count + 1
          
        }#end else 
      }#for k 
    }#else 
  }
  
  BurgenNew$observer <- as.factor(BurgenNew$observer)
  
  return(BurgenNew)
}


#create habitat segments, uded in density/count analysis for mrds. Assumes that habitat does not vary within the segment. See segment.data in ?dsm
create.Seg <- function(Hab, Trans, lakename) {
  
  SegDat <- merge(Hab, Trans, by="record_id") %>% rename(Sample.Label=hab_id) %>% subset(Observer == 1) #assumes habitat measurments from observer 1 are more accurate
  SegDat <- SegDat %>% mutate(Effort = rep(NA, dim(SegDat)[1]), Transect=as.numeric(str_split(Sample.Label, pattern="_T", n=2, simplify = TRUE)[,2])) %>% subset(Region.Label==lakename)
  SegDat$Region.Label <- factor(SegDat$Region.Label)
   
  SegDat$Sample.Label <- as.factor(as.character(SegDat$Sample.Label))
  SegDat <- SegDat[order(SegDat$record_id, SegDat$distance_along),]
  
  for(curr in SegDat$record_id) {
    ind     <- which(SegDat$record_id == curr)
    dis.vec <- c(0, SegDat$distance_along[ind])
    
    SegDat$Effort[ind] <- diff(dis.vec)
    
  }
  
  #get locations of the segments, use utm coordinates for smoothing later
  SegDat <- SegDat %>% mutate(utm_easting=rep(NA, dim(SegDat)[1]), utm_northing=rep(NA, dim(SegDat)[1]))
  
  for(i in 1:dim(SegDat)[1]) {
    
    bearing <- SegDat$bearing[i]
    gpsN    <- SegDat$gps_transect_n[i]
    gpsW    <- SegDat$gps_transect_w[i]
    
    xy              <- data.frame(ID = 1, X = c(gpsW), Y = c(gpsN))
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res             <- spTransform(xy, CRS("+proj=utm +zone=15 ellps=WGS84"))
    
    SegDat$utm_easting[i] <- coordinates(res)[1] + sin(-pi/180*(bearing))*SegDat$distance_along[i]
    SegDat$utm_northing[i] <- coordinates(res)[2] + cos(-pi/180*(bearing))*SegDat$distance_along[i]
    
  }
  #cull missing data
  if(any(is.na(SegDat$clarity))) {
    SegDat <- SegDat[-which(is.na(SegDat$clarity)),]
  }
  
  return(SegDat)
  
}

#this creates the ObsDat format needed by the dsm package. See observation.data in ?dsm
create.Obs <- function(ObsDat, SegDat, lakename) {

  ObsDat      <- ObsDat %>% subset(!is.na(size)) %>% rename(Transect.Label = Sample.Label, Lake=Region.Label) 
  ObsDat      <- ObsDat %>% subset(Lake==lakename)
  ObsDat$Lake <- factor(ObsDat$Lake)
  ObsDat      <- ObsDat %>% mutate(Sample.Label=rep(NA, dim(ObsDat)[1]), object=1:dim(ObsDat)[1])
  
  #get the proper Sample.Label for each observation to link with segdata
  for(i in 1:dim(ObsDat)[1]) {
    currTransect <- ObsDat$Transect.Label[i] #this label
    possHabs     <- which(currTransect == as.character(SegDat$record_id)) #all possible habitats with this transect label
    currInd     <- which(ObsDat$mussel_distance_along[i] <= SegDat$distance_along[possHabs])[1]
    ObsDat$Sample.Label[i] <- as.character(SegDat$Sample.Label[possHabs[currInd]])
    
  }
  
  ObsDat$Sample.Label <- as.factor(ObsDat$Sample.Label)
  
  return(ObsDat)

}

#this creates a total counts in each segment.
#for single observer
create.CountSumm.single <- function(SegDat, ObsDat, region="Sylvia", ds.model) {
  
  segment.data <- subset(SegDat, Region.Label==region)
  segment.data <- segment.data %>% mutate(Area = abs(Effort*2))
  segment.data$plant_cover <- as.factor(segment.data$plant_cover)
  segment.data$plant_cover <- recode(segment.data$plant_cover, '0'="None", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  #for each segment get the total counts
  levels(ObsDat$Sample.Label) <- levels(segment.data$Sample.Label)
  
  count.vec <- vector('numeric', dim(segment.data)[1])
  for(i in 1:dim(segment.data)[1]) {
    currSeg <- segment.data$Sample.Label[i]
    count.vec[i] <- sum(ObsDat$size[(which(ObsDat$Sample.Label == currSeg))])
  }
  phat.pred <- vector('numeric', dim(segment.data)[1])
  
  segment.data <- mutate(segment.data, phat=phat.pred, Counts=count.vec, Effective.area=phat.pred*segment.data$Area)
  
  return(segment.data)
}


#Draws estimated detection probabilities 
#could also redraw segments but then do I need to control for length....
resample.detection <- function(countDat, ddf.dsModel, SegDat, ObsDat, det=F) {
  
  segment.data <- subset(SegDat, Region.Label=="Burgen")
  segment.data <- segment.data %>% mutate(Area = Effort*2)
  segment.data$plant_cover <- as.factor(segment.data$plant_cover)
  segment.data$plant_cover <- recode(segment.data$plant_cover, '0'="None", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  #for each segment get the total counts
  levels(ObsDat$Sample.Label) <- levels(segment.data$Sample.Label)
  count.vec <- vector('numeric', dim(segment.data)[1])
  
  for(i in 1:dim(segment.data)[1]) {
    currSeg <- segment.data$Sample.Label[i]
    count.vec[i] <- sum(ObsDat$size[(which(ObsDat$Sample.Label == currSeg))])
  }
  
  phat.pred <- vector('numeric', dim(segment.data)[1])
  
  #generate p-hats
  if(det == F) {
    varcov      <- solve(ddf.dsModel$mr$hessian)
    pars2.draw  <- mvrnorm(n=1, mu=coef(ddf.dsModel$mr)[,1], Sigma=varcov)
    dispar      <- exp(rnorm(1, coef(ddf.dsModel$ds)$scale[1,1], coef(ddf.dsModel$ds)$scale[1,2]))
  } else {
    pars2.draw <- coef(ddf.dsModel$mr)[,1]
    dispar     <- exp(coef(ddf.dsModel$ds)$scale[1,1])
  }

  for(i in 1:dim(segment.data)[1]) {

    pars2 <- pars2.draw*c(1, 1, segment.data$plant_cover[i]=="Present", segment.data$clarity[i])
    
    pars1 <- pars2[-2]
    p1 <- logit(sum(pars1))
    p2 <- logit(sum(pars2))
    p0 <- p1 + p2 - p1*p2
    
    phat.pred[i] <- sqrt(pi/2)*dispar*erf(1/(sqrt(2)*dispar), 1)*p0
  }
  
  segment.data <- mutate(segment.data, phat=phat.pred, Counts=count.vec, Effective.area=phat.pred*segment.data$Area)
  
  return(segment.data)

}

nbin.regression <- function(pars, X, y, offset) {
  beta <- pars[-length(pars)]
  theta <- exp(pars[length(pars)])
  
  mu <- exp(X%*%beta + log(offset))
  return(-sum(dnbinom(y, size=theta, mu=mu, log=T)))
}

##randomly draw phats from the two-piece normal model
resample.detection.twopiece <- function(countDat, pars2.draw, MCDS.draw, SegDat, ObsDat, det=F) {
  
  segment.data <- subset(SegDat, Region.Label=="Burgen")
  segment.data <- segment.data %>% mutate(Area = Effort*2)
  
  segment.data$plant_cover <- as.factor(segment.data$plant_cover)
  segment.data$plant_cover <- recode(segment.data$plant_cover, '0'="None", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  #for each segment get the total counts
  levels(ObsDat$Sample.Label) <- levels(segment.data$Sample.Label)
  count.vec <- vector('numeric', dim(segment.data)[1])
  
  for(i in 1:dim(segment.data)[1]) {
    currSeg       <- segment.data$Sample.Label[i]
    count.vec[i]  <- sum(ObsDat$size[(which(ObsDat$Sample.Label == currSeg))])
  }
  
  phat.pred <- vector('numeric', dim(segment.data)[1])
  
  #generate p-hats
  mew   <- exp(MCDS.draw[1])
  sig1  <- exp(MCDS.draw[2])
  sig2  <- exp(MCDS.draw[3])
  
  dis <- seq(0, 1, length=1e3)
  dis <- dis[-length(dis)]
  dx  <- dis[2] - dis[1]
  
  p.int <- sum(dx*Two.piece.normal.func(sig1=sig1, sig2=sig2, mu=mew, dis.vec=dis))
  
  cov.mat <- cbind(rep(1, dim(segment.data)[1]), rep(1, dim(segment.data)[1]), segment.data$plant_cover=="Present")
  
  for(i in 1:dim(segment.data)[1]) {
    pars2 <- pars2.draw*cov.vec[1:length(pars2.draw)]
    pars1 <- pars2[-2]
    
    p1 <- logit(sum(pars1))
    p2 <- logit(sum(pars2))
    p0 <- p1 + p2 - p1*p2
    
    phat.pred[i] <- p0*p.int
    
  }

  segment.data <- mutate(segment.data, phat=phat.pred, Counts=count.vec, Effective.area=phat.pred*segment.data$Area)
  
  return(segment.data)
  
}


#error function (integral of normal distribution)
erf <- function(x, sigma) {
  2*pnorm(x*sqrt(2), 0, sigma) - 1
}


#Count up the observations in each segment, then returns dataframe with segment counts
#used to model segment counts.
create.CountSumm <- function(SegDat, ObsDat, BurgenDDF, det=F) {
  
  segment.data <- subset(SegDat, Region.Label=="Burgen")
  segment.data <- segment.data %>% mutate(Area = Effort*2)
  segment.data$plant_cover <- as.factor(segment.data$plant_cover)
  segment.data$plant_cover <- recode(segment.data$plant_cover, '0'="None", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  #for each segment get the total counts
  count.vec <- vector('numeric', dim(segment.data)[1])
  
  for(i in 1:dim(segment.data)[1]) {
    currSeg <- as.character(segment.data$Sample.Label[i])
    count.vec[i] <- sum(ObsDat$size[(which(as.character(ObsDat$Sample.Label) == currSeg))])
  }
  
  ddf.dsModel <- ddf(method="io", dsmodel=~cds(key="hr", formula=~ 1), mrmodel=~glm(link="logit", formula=~ observer + plantcover + clarity), data=BurgenDDF, meta.data=list(width=1))
  
  phat.fitted <- predict(ddf.dsModel)$fitted
  
  phat.pred <- vector('numeric', dim(segment.data)[1])
  
  segment.data <- mutate(segment.data, phat=NA, Counts=count.vec, Effective.area=phat.pred*segment.data$Area)
  
  return(segment.data)
  
}


#Count up the observations in each segment, then returns dataframe with segment counts
#used to model segment counts.
#uses the two-piece normal model
create.CountSumm.twopiece <- function(SegDat, ObsDat, BurgenDDF, ddf.model, MCDS.est, det=F) {
  
  segment.data <- subset(SegDat, Region.Label=="Burgen")
  segment.data <- segment.data %>% mutate(Area = Effort*2)
  segment.data$plant_cover <- as.factor(segment.data$plant_cover)
  segment.data$plant_cover <- recode(segment.data$plant_cover, '0'="None", '1'="Present", '2'="Present", '3'="Present", '4'="Present")
  
  #for each segment get the total counts
  count.vec <- vector('numeric', dim(segment.data)[1])
  
  for(i in 1:dim(segment.data)[1]) {
    currSeg       <- as.character(segment.data$Sample.Label[i])
    count.vec[i]  <- sum(ObsDat$size[(which(as.character(ObsDat$Sample.Label) == currSeg))])
  }
  
  phat.pred   <- vector('numeric', dim(segment.data)[1])
  
  #generate p-hats
  #p0.varcov <- solve(ddf.model$mr$hessian)
  p0.mean   <- ddf.model$mr$par
  
  #MCDS.varcov <- solve(MCDS.model$ds.hessian)
  #MCDS.est    <- c(MCDS.model$Scale.parm.est, MCDS.model$Shape.parm.est)
  mew   <- exp(MCDS.est[1])
  sig1  <- exp(MCDS.est[2])
  sig2  <- exp(MCDS.est[3])
  
  for(i in 1:dim(segment.data)[1]) {
    cov.vec <- c(1, 1, segment.data$plant_cover[i]=="Present", segment.data$clarity[i])
    p0.obs2 <- p0.mean*cov.vec[1:length(p0.mean)]
    p0.obs1 <- p0.obs2[-2]
    p1 <- logit(sum(p0.obs1))
    p2 <- logit(sum(p0.obs2))
    p0 <- p1 + p2 - p1*p2
    
    dis <- seq(0, 1, length=1e3)
    dis <- dis[-length(dis)]
    dx  <- dis[2] - dis[1]
    
    phat.pred[i] <- sum(dx*p0*Two.piece.normal.func(sig1=sig1, sig2=sig2, mu=mew, dis.vec=dis))
  }
  
  segment.data <- mutate(segment.data, phat=phat.pred, Counts=count.vec, Effective.area=phat.pred*segment.data$Area)
  #remove missing data entries.
  if(any(is.na(segment.data$clarity))) {
    segment.data <- segment.data[-which(is.na(segment.data$clarity)),]
  }
  return(segment.data)
  
}


#function to calculate logistic of value x
logit <- function(x) {exp(x)/(1 + exp(x))}



#Gavin Simpsons simulate function
simulate.gam <- function(object, nsim=1, seed=NULL, newdata, freq=FALSE, unconditional=FALSE) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  if (missing(newdata)) {
    newdata <- object$model
  }
  
  ## random multivariate Gaussian
  ## From ?predict.gam in mgcv package
  ## Copyright Simon N Wood
  ## GPL >= 2
  rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m * n), m, n))
  }
  
  Rbeta <- rmvn(n = nsim, mu = coef(object), sig = vcov(object, freq=freq, unconditional=unconditional))
  Xp <- predict(object, newdata = newdata, type = "lpmatrix")
  sims <- Xp %*% t(Rbeta)
  sims
}


#bootstrap count prediction covarainces
create.abundVar <- function(fit.dens4, theta.est, nsim=1e4, offset) {
  
  count.mat    <- exp(simulate.gam(fit.dens4, nsim=nsim))*offset
  count.mean   <- apply(count.mat, 1, mean)
  count.bias   <- count.mean - as.vector(count.est) 
  count.varcov <- matrix(NA, dim(count.mat)[1], dim(count.mat)[1])
  
  for(j in 1:dim(count.mat)[1]) {
    for(k in 1:dim(count.mat)[1]) {
      count.varcov[j,k] <- sum((count.mat[j,] - count.mean[j])*(count.mat[k,] - count.mean[k]))/(dim(count.mat)[2]-1)*theta.est[j]*theta.est[k]
    }
  }
  
  return(count.varcov)
}


##draw new phats
create.phatVar <- function(count.dat, ddf.dsModel, MCDS.fit, SegDat, ObsDat, B=1e2, type="inverse") {
  
  phat.mat  <- matrix(NA, dim(count.dat)[1], B)
  
  p0.varcov <- solve(ddf.dsModel$mr$hessian)
  p0.mean   <- ddf.dsModel$mr$par
  
  MCDS.varcov <- solve(MCDS.fit$hessian)
  MCDS.est    <- c(MCDS.fit$solution)
  
  varcov      <- solve(ddf.dsModel$mr$hessian)
  pars2.draw  <- mvrnorm(n=B, mu=p0.mean, Sigma=p0.varcov)
  MCDS.draw   <- mvrnorm(n=B, mu=MCDS.est, Sigma=MCDS.varcov)
  
  for(i in 1:B) {
    #phat.mat[,i] <- resample.detection.twopiece(countDat=count.dat, pars2.draw=pars2.draw[i,], matrix(MCDS.est, B, 3), SegDat=SegDat, ObsDat=ObsDat, det=F)$phat #
    phat.mat[,i] <- resample.detection.twopiece(countDat=count.dat, pars2.draw=pars2.draw[i,], MCDS.draw[i,], SegDat=SegDat, ObsDat=ObsDat, det=F)$phat
  }
  
  phat.varcov <- diag(count.est/Area)%*%cov(t(1/phat.mat))%*%diag(count.est/Area)
  
  return(list(phat.mat, phat.varcov))
}



Two.piece.normal.func <- function(sig1, sig2, mu, dis.vec) {
  dis1.vec <- dis.vec[dis.vec < mu]
  dis2.vec <- dis.vec[dis.vec >= mu]
  return(c(exp(-(dis1.vec - mu)^2/(2*sig1^2)), exp(-(dis2.vec - mu)^2/(2*sig2^2))))
}

halfnormal.func <- function(sig, dis.vec) {
  return(exp(-(dis.vec)^2/(2*sig^2)))
}


#primary coded as 1
#secondary coded as 2
#secondary sees all that 1 does but then a bit more
#creates 
create.removal.Observer <- function(transect.dat, obs.dat) {
  
  primary   <- transect.dat$`Primary observer (double observer survey)`
  secondary <- transect.dat$`Secondary observer (double observer survey)`
  
  obs.dat <- obs.dat %>% mutate(primary=primary[obs.dat$`Transect #`- min(obs.dat$`Transect #`)+1], secondary=secondary[obs.dat$`Transect #`- min(obs.dat$`Transect #`)+1], observer = rep(1, dim(obs.dat)[1]))
  
  for(i in 1:dim(obs.dat)[1]) {
    curr.tran             <- obs.dat[i,]$`Transect #`
    #obs.dat[i,]$observer  <- as.numeric(obs.dat[i,]$`Observer name` == transect.dat[transect.dat$`Transect number` == curr.tran,]$`Primary observer (double observer survey)`)
    obs.dat[i,]$observer  <- as.numeric(obs.dat[i,]$`Observer name` == transect.dat[transect.dat$`Transect number` == curr.tran,]$`Primary observer (double observer survey)`)
  }

  
  #all observations made by secondary that were not made by primary
  temp01.obs <- obs.dat[obs.dat$observer==0,]
  if(dim(temp01.obs)[1]) {
    temp01.obs$observer <- 1
    temp01.obs$detected <- 0
  }
  #print(dim(temp01.obs)[1])
  temp012.obs <- temp01.obs
  if(dim(temp012.obs)[1]) {
    temp012.obs$observer <- 0
    temp012.obs$detected <- 1
  }
  
  #obs.dat <- rbind(obs.dat, temp.obs)
  
  #copy all observations made by the primary observer to the secondary observer
  temp11.obs <- obs.dat[obs.dat$observer==1,]
  temp11.obs$observer <- 0
  temp11.obs$detected <- 1
  temp112.obs <- temp11.obs
  temp112.obs$observer <- 1
  temp112.obs$detected <- 1
  #obs.dat <- rbind(obs.dat, temp.obs)
  
  obs.dat <- rbind(temp01.obs, temp012.obs, temp11.obs, temp112.obs)
  
  obs.dat$observer[which(obs.dat$observer==0)] <- 2
  obs.dat$primary <- as.factor(obs.dat$primary)
  obs.dat <- obs.dat[order(obs.dat$object),]
  
  return(obs.dat)
}

#what does do?
doubleObs.duplicate <- function(obs.dat) {
  
  num.detect <- which(obs.dat$size > 1)
  if(length(num.detect) == 0) {
    return(obs.dat)
  } 
  
  for(i in 1:length(num.detect)) {
    
    temp.dat  <- obs.dat[num.detect[i],]
    for(j in 2:obs.dat[num.detect[i],]$size) {
      temp.dat <- rbind(temp.dat, temp.dat[1,])
    }
    obs.dat <- rbind(obs.dat, temp.dat[-1,])
  }
  
  return(obs.dat)
  
}



##estimate densities using distnace survey
distance.dens.est <- function(curr.lake, curr.lake2, var.type="design", dist.type="hr", size=FALSE) {
  
  #distance.title    <- suppressMessages(gs_title("Encounters - Double observer - distance survey (Responses)"))
  #distance.dat      <- suppressMessages(gs_read(distance.title, verbose=FALSE))
  distance.dat      <- read_xlsx(path="../Data/Season2/Encounters - Double observer - distance survey (Responses).xlsx", sheet=1)
  distance.dat      <- distance.dat %>% subset(`Lake name` == curr.lake)
  distance.dat      <- distance.dat[order(distance.dat$`Transect #`),]
  distance.dat      <- dplyr::rename(distance.dat, size="Number of mussels in cluster")
  
  #remove 'empty' observations
  if(any(distance.dat$size == 0)){
    rm.vec <- which(distance.dat$size ==0)
    distance.dat <- distance.dat[-rm.vec,]
  }
  if(any(is.na(distance.dat$size))) {
    rm.vec <- which(is.na(distance.dat$size))
    distance.dat <- distance.dat[-rm.vec,]
  }
  varS <- var(distance.dat$size)/length(distance.dat)
  eS   <- mean(distance.dat$size)
  
  distanceorig.dat  <- distance.dat
  
  #transect.title  <- suppressMessages(gs_title("Transect datasheets (Responses)"))
  #transect.dat    <- suppressMessages(gs_read(transect.title, verbose=FALSE))
  transect.dat      <- read_xlsx(path="../Data/Season2/Transect datasheets (Responses).xlsx", sheet=1)
  
  transect.dat <- transect.dat %>% subset(`Lake name:` == curr.lake2) 
  transect.dat <- transect.dat %>% subset(`Survey type` == "Double observer distance")
  transect.dat <- transect.dat[order(transect.dat$"Transect number"),]
  
  setup.time   <- 60^2*hour(transect.dat$`Setup time`) + 60*minute(transect.dat$`Setup time`) + second(transect.dat$`Setup time`)
  hab.time     <- 60^2*hour(transect.dat$`Habitat time`) + 60*minute(transect.dat$`Habitat time`) + second(transect.dat$`Habitat time`)
  enc.time     <- 60^2*hour(transect.dat$`Encounter time`) + 60*minute(transect.dat$`Encounter time`) + second(transect.dat$`Encounter time`)
  
  trans.length  <- sum(transect.dat$`Transect length (if transect survey)`)
  trans.area    <- trans.length*2
  
  dimnames(distance.dat)[[2]][7] <- "distance"
  distance.dat$distance <- distance.dat$distance/100
  
  distance.dat <- distance.dat %>% mutate(detected = rep(1, dim(distance.dat)[1]), object=1:dim(distance.dat)[1])
  distance.dat <- create.removal.Observer(transect.dat=transect.dat, obs.dat=distance.dat)
  
  dis.totdetect <- dim(distance.dat)[1]
  distance.dat$`Observer name` <- as.factor(distance.dat$`Observer name`)
  if(size) {
    dis.detect.model <- ddf(method="rem", dsmodel=~cds(key=dist.type, formula=~1), mrmodel=~glm(formula=~size), data=distance.dat, meta.data=list(width=1))
  } else {
    dis.detect.model <- ddf(method="rem", dsmodel=~cds(key=dist.type), mrmodel=~glm(formula=~1), data=distance.dat, meta.data=list(width=1))
  }
  
  phat    <- summary(dis.detect.model)$average.p
  phat.se <- summary(dis.detect.model)$average.p.se[1,1]
  
  #estimate density
  ##assume poisson model
  count.vec <- detect.vec <- rep(NA, 15)

  for(i in transect.dat$`Transect number`) {
    count.vec[i]  <- sum(distanceorig.dat[which(distanceorig.dat$`Transect #` == i),]$size)
    detect.vec[i] <- sum(distance.dat[which(distance.dat$`Transect #` == i),]$detected)
  }
  names(count.vec) <- names(detect.vec) <- transect.dat$`Transect number`
  if(any(is.na(count.vec))) {
    #count.vec[which(is.na(count.vec))]   <- 0
    #detect.vec[which(is.na(detect.vec))] <- 0
    count.vec  <- count.vec[-which(is.na(count.vec))]  # <- 0
    detect.vec <- detect.vec[-which(is.na(detect.vec))] #<- 0
    
  }
  area.vec <- t(as.matrix(2*transect.dat$`Transect length (if transect survey)`))[1,]
  names(area.vec) <- transect.dat$`Transect number`
  
  lambda <- sum(count.vec, na.rm=T)/sum(area.vec)
  
  dhat <- sum(detect.vec/phat, na.rm=TRUE)*eS/sum(area.vec)
  dhat.detect <- sum(detect.vec/phat, na.rm=TRUE)/sum(area.vec)
  
  df <- data.frame(Detections=detect.vec, Mussels=count.vec, D=detect.vec*eS/area.vec/phat, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, Area=area.vec, Length=area.vec/2, Type=rep("Distance", length(area.vec)), Lake=rep(curr.lake, length(area.vec)), Lat=transect.dat$`GPS latitude`, Long=transect.dat$`GPS longitude`)
  
  if(var.type == "design") {
    
    K  <- length(detect.vec)
    l  <- area.vec/2
    L  <- sum(l)
    n  <- detect.vec
    N  <- sum(detect.vec)
    var.n <- L/(K-1)*sum(l*(n/l - N/L)^2)
    
    dhat.se <- sqrt(dhat^2*(var.n/N^2 + varS/K/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat, Dhat.se=dhat.se, phat=phat, phat.se=phat.se, Transects=length(area.vec), Area=area.vec, Detections=detect.vec, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=dis.detect.model, distance.dat=distance.dat, eS=eS, varS=varS)
    
  } 
  if(var.type == "jack") { #from buckland page 80
    
    K       <- length(detect.vec)
    l       <- area.vec/2
    L       <- sum(l)
    
    dhat            <- sum(count.vec/phat)/sum(area.vec)
    dhat.detect     <- sum(detect.vec/phat)/sum(area.vec)
    
    dhat.K <- pseudo.vec <- vector('numeric', K)
    
    for(k in 1:K) {
      dhat.K[k]     <- sum(count.vec[-k]/phat)/sum(area.vec[-k])
      pseudo.vec[k] <- (L*dhat - (L-l[k])*dhat.K[k])/l[k]
    }
    
    dhat.se.jk  <- sqrt(1/(L*K-L)*sum(l*(pseudo.vec - dhat)^2) )

    b.jk        <- (K-1)*(mean(dhat.K) - dhat)
    dhat.jk     <- K*dhat - (K-1)*mean(dhat.K) 
    
    dhat.se.full <- sqrt(dhat.se.jk^2 + dhat^2*(phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat.detect*eS, Dhat.se=dhat.se.full, phat=phat, phat.se=phat.se, Transects=length(area.vec), Area=area.vec, Detections=detect.vec, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=dis.detect.model, eS=eS, varS=varS)
    
  } 
  if(var.type == "model") {
    dist.glm      <- glm.nb(detect.vec~1 + offset(log(phat*area.vec)))
    dist.predict <- predict(dist.glm, type="response", se.fit=T)
    
    dhat.se <- sqrt(dhat^2*(sum(dist.predict$se.fit^2)/sum(dist.predict$fit)^2 + varS/length(detect.vec)/eS^2 + phat.se^2/phat^2))
   # dhat.se <- sqrt(dhat^2*(var.n/N^2 + varS/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat, Dhat.se=dhat.se, phat=phat, phat.se=phat.se, Transects=length(area.vec), Area=area.vec, Detections=detect.vec, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=dis.detect.model, eS=eS, varS=varS)
    
  }
  
  
  return(return.list)
}


##estimate density from quadrats
quadrat.dens.est <- function(curr.lake, curr.lake2, var.type, quad.side=0.5) {
  
  #quadrat.title    <- suppressMessages(gs_title("Encounters - Quadrats (Responses)"))
  #quadrat.dat      <- suppressMessages(gs_read(quadrat.title, verbose=FALSE))
  quadrat.dat      <- read_xlsx(path="../Data/Season2/Encounters - Quadrats (Responses).xlsx", sheet=1)
  quadrat.dat      <- quadrat.dat %>% subset(`Lake name` == curr.lake)
  if(curr.lake == "Little Birch Lake") {
    quadrat1.dat <- subset(quadrat.dat, (`Transect #` == 1 | `Transect #` == 2 | `Transect #` == 3 | `Transect #` == 4 | `Transect #` == 5 | `Transect #` == 6 | `Transect #` == 7 | `Transect #` == 8 | `Transect #` == 9 | `Transect #` == 10 | `Transect #` == 11) &  `Observer name` =="Austen" )
    #for quadrats 1-11 only use Aislyn
    quadrat2.dat <- subset(quadrat.dat, `Transect #` == 15 | `Transect #` == 14 | `Transect #` == 13 | `Transect #` == 12)
    quadrat.dat <- rbind(quadrat1.dat, quadrat2.dat)
  }

  quadrat.dat <- quadrat.dat[order(quadrat.dat$`Transect #`),]
  
  if(any(is.na(quadrat.dat$`Number of mussels in quadrat`))) {
    quadrat.dat$`Number of mussels in quadrat`[which(is.na(quadrat.dat$`Number of mussels in quadrat`))] <- 0
  }
  #transect.title  <- suppressMessages(gs_title("Transect datasheets (Responses)"))
  #transect.dat    <- suppressMessages(gs_read(transect.title, verbose=FALSE))
  transect.dat      <- read_xlsx(path="../Data/Season2/Transect datasheets (Responses).xlsx", sheet=1)
  
  transect.dat <- transect.dat %>% subset(`Lake name:` == curr.lake2) 
  transect.dat <- transect.dat %>% subset(`Survey type` == "Quadrat survey")
  transect.dat <- transect.dat[order(transect.dat$"Transect number"),]
  
  #if(curr.lake == "Little Birch Lake") {
  #  transect1.dat <- subset(transect.dat, (`Transect number` == 1 | `Transect number` == 2 | `Transect number` == 3 | `Transect number` == 4 | `Transect number` == 5 | `Transect number` == 6 | `Transect number` == 7 | `Transect number` == 8 | `Transect number` == 9 | `Transect number` == 10 | `Transect number` == 11))
    #for transects 1-11 only use Aislyn because sruvyes were done on the same transect
  #  transect2.dat <- subset(transect.dat, `Transect number` == 15 | `Transect number` == 14 | `Transect number` == 13 | `Transect number` == 12)
  #  transect.dat <- rbind(transect1.dat, transect2.dat)
    #transect.dat <- subset(transect.dat, `Transect number` == 15 | `Transect number` == 14 | `Transect number` == 13 | `Transect number` == 12)
  #}
  
  setup.time   <- 60^2*hour(transect.dat$`Setup time`) + 60*minute(transect.dat$`Setup time`) + second(transect.dat$`Setup time`)
  hab.time    <- 60^2*hour(transect.dat$`Habitat time`) + 60*minute(transect.dat$`Habitat time`) + second(transect.dat$`Habitat time`)
  enc.time    <- 60^2*hour(transect.dat$`Encounter time`) + 60*minute(transect.dat$`Encounter time`) + second(transect.dat$`Encounter time`)
  
  #get number of quadrats per transect
  trans.quads <- trans.counts <- rep(NA, max(quadrat.dat$`Transect #`))
  for(i in unique(quadrat.dat$`Transect #`)) {
    trans.quads[i]  <- sum(quadrat.dat$`Transect #` == i)
    trans.counts[i] <- sum(quadrat.dat[quadrat.dat$`Transect #` == i,]$`Number of mussels in quadrat`)
  }
  
  if(any(is.na(trans.quads))) {
    rm.vec         <- which(is.na(trans.quads))
    trans.quads    <- trans.quads[-rm.vec]
    trans.counts   <- trans.counts[-rm.vec]
  }
  trans.num   <- dim(transect.dat)[1]
  quad.counts <- (quadrat.dat$`Number of mussels in quadrat`)
  quad.area   <- quad.side^2
  
  #trans.counts <- quad.counts
  names(trans.counts) <- transect.dat$`Transect number`
  area.vec <- trans.quads*quad.area
  names(area.vec) <- transect.dat$`Transect number`
  
  dhat <- sum(trans.counts)/sum(area.vec)
     
  df <- data.frame(Detections=trans.counts, Mussels=trans.counts, D=trans.counts/area.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, Area=area.vec, Length=transect.dat$`Transect length (if transect survey)`, Type=rep("Quadrat", length(trans.quads)), Lake=rep(curr.lake, length(trans.quads)), Lat=transect.dat$`GPS latitude`, Long=transect.dat$`GPS longitude`)
  
  if(var.type == "design") {
    
    K  <- length(trans.counts)
    l  <- df$Length #sqrt(trans.quads)
    L  <- sum(df$Length) #sum(sqrt(trans.quads))
    nl <- trans.counts/l
    N  <- sum(trans.counts)
    NL <- N/L
    var.n <- L/(K-1)* sum(l*(nl - NL)^2)
    dhat.se <- sqrt(dhat^2*(var.n/N^2))
    
    quad.list <- list(Quadrats=trans.quads, Area=area.vec, Mussels=trans.counts, Dhat=dhat, Dhat.se=dhat.se, Time=setup.time + hab.time + enc.time, df=df)
  } 
  if(var.type == "jack") {
      
    K       <- length(trans.quads)
    l       <- trans.quads
    L       <- sum(l)
    dhat    <- sum(trans.counts)/sum(area.vec)
    
    dhat.K <- pseudo.vec <- vector('numeric', K)
    
    for(k in 1:K) {
      dhat.K[k]     <- sum(trans.counts[-k])/sum(area.vec[-k])
      pseudo.vec[k] <- (L*dhat - (L-l[k])*dhat.K[k])/l[k] #K*dhat - (K-1)*dhat.K[k]
    }
    
    dhat.se.jk  <- sqrt(1/(L*K-L)*sum(l*(pseudo.vec - dhat)^2) )
      
    quad.list <- list(Dhat=dhat, Dhat.se=dhat.se.jk, Transects=length(area.vec), Area=area.vec, Mussels=trans.counts, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df)
      
  }
  if(var.type == "model") {
    
    #quad.pois <- glm(trans.counts[1,] ~ 1 + offset(log(quad.area*trans.quads)), family=poisson)
    quad.nb       <- glm.nb(trans.counts ~ 1 + offset(log(quad.area*trans.quads)))
    quad.predict  <- predict(quad.nb, type="response", se.fit=T)
  
    quad.list <- list(Quadrats=trans.quads, Area=area.vec, Mussels=trans.counts, Dhat=dhat, Dhat.se=sum(quad.predict$se.fit)/sum(quad.area*trans.quads), Time=setup.time + hab.time + enc.time, df=df)
  }
  
  return(quad.list)
}



double.dens.est <- function(curr.lake, curr.lake2, var.type, size=FALSE) {
  
  #double.title    <- suppressMessages(gs_title("Encounters - Double observer - no distance (Responses)"))
  #double.dat      <- suppressMessages(gs_read(double.title, verbose=FALSE))
  double.dat      <- read_xlsx(path="../Data/Season2/Encounters - Double observer - no distance (Responses).xlsx", sheet=1)
  double.dat      <- double.dat %>% subset(`Lake name` == curr.lake)
  double.dat      <- double.dat[order(double.dat$`Transect #`),]
  double.dat      <- dplyr::rename(double.dat, size="Number of mussels in cluster")
  doubleorig.dat <- double.dat
  if(curr.lake == "Little Birch Lake") {
    double.dat <- subset(double.dat, `Transect #` == 15 | `Transect #` == 14 | `Transect #` == 13 | `Transect #` == 12)
  }
  
  if(any(double.dat$size==0)) {
    double.dat <- double.dat[-which(double.dat$size==0),]
  }
  if(any(is.na(double.dat$size))) {
    double.dat <- double.dat[-which(is.na(double.dat$size)),]
  }
  
  eS <- mean(double.dat$size, na.rm=T)
  varS <- var(double.dat$size, na.rm=T)#/length(double.dat$size) 

  #transect.title  <- suppressMessages(gs_title("Transect datasheets (Responses)"))
  #transect.dat    <- suppressMessages(gs_read(transect.title, verbose=FALSE))
  transect.dat      <- read_xlsx(path="../Data/Season2/Transect datasheets (Responses).xlsx", sheet=1)
  transect.dat <- transect.dat %>% subset(`Lake name:` == curr.lake2)
  transect.dat <- transect.dat %>% subset(`Survey type` == "Double observer no distance")
  transect.dat <- transect.dat[order(transect.dat$"Transect number"),]
  if(curr.lake == "Little Birch Lake") {
    transect.dat <- subset(transect.dat, `Transect number` == 15 | `Transect number` == 14 | `Transect number` == 13 | `Transect number` == 12)
  }
  
  
  setup.time   <- 60^2*hour(transect.dat$`Setup time`) + 60*minute(transect.dat$`Setup time`) + second(transect.dat$`Setup time`)
  hab.time     <- 60^2*hour(transect.dat$`Habitat time`) + 60*minute(transect.dat$`Habitat time`) + second(transect.dat$`Habitat time`)
  enc.time     <- 60^2*hour(transect.dat$`Encounter time`) + 60*minute(transect.dat$`Encounter time`) + second(transect.dat$`Encounter time`)
  
  trans.length  <- transect.dat$`Transect length (if transect survey)`
  trans.area    <- trans.length
  
  
  #get number of obsevations made by primary and secondary observers
  count.primary <- count.secondary <- detect.primary <- detect.secondary <- rep(NA, max(transect.dat$`Transect number`))
  for(i in unique(transect.dat$`Transect number`)) {
    curr.prim <- transect.dat[transect.dat$`Transect number` == i,]$`Primary observer (double observer survey)`
    curr.dat  <- double.dat[double.dat$`Transect #`== i,] 
    
    count.primary[i]    <- max(sum(curr.dat[curr.dat$`Observer name` == curr.prim,]$size), 0, na.rm=T)
    count.secondary[i]  <- max(sum(curr.dat[curr.dat$`Observer name` != curr.prim,]$size), 0, na.rm=T)
    
    detect.primary[i]    <- max(dim(curr.dat[curr.dat$`Observer name` == curr.prim,])[1], 0, na.rm=T)
    detect.secondary[i]  <- max(dim(curr.dat[curr.dat$`Observer name` != curr.prim,])[1], 0, na.rm=T)
    
  }
  
  if(any(is.na(count.primary))) {
    rm.vec          <- which(is.na(count.primary))
    count.primary   <- count.primary[-rm.vec]
    count.secondary <- count.secondary[-rm.vec]
    detect.primary  <- detect.primary[-rm.vec]
    detect.secondary <- detect.secondary[-rm.vec]
  }
    
  double.dat <- double.dat %>% mutate(detected = rep(1, dim(double.dat)[1]), distance=rep(0.99, dim(double.dat)[1]), object=1:dim(double.dat)[1])
  
  #double.dat <- doubleObs.duplicate(obs.dat=double.dat)
  double.dat$object <- 1:dim(double.dat)[1]
  double.dat <- create.removal.Observer(transect.dat=transect.dat, obs.dat=double.dat)
  
  if(size) {
    doubleDetect.model <- ddf(method="rem.fi", mrmodel=~glm(formula=~size), dsmodel=~cds(key="unif"), data=double.dat, meta.data=list(width=1))
  } else {
    doubleDetect.model <- ddf(method="rem.fi", mrmodel=~glm(formula=~1), dsmodel=~cds(key="unif"), data=double.dat, meta.data=list(width=1))
  }
  phat    <- summary(doubleDetect.model)$average.p0.1
  phat.se <- summary(doubleDetect.model)$average.p0.1.se[1,1]
#print(doubleDetect.model)

  count.vec <- count.primary + count.secondary
  detect.vec <- detect.primary + detect.secondary
  
  count.vec <- t(as.matrix(count.vec))
  dimnames(count.vec)[[2]] <- transect.dat$`Transect number`
  
  area.vec                <- trans.length
  names(area.vec) <- transect.dat$`Transect number`
  
  trans.area <- sum(transect.dat$`Transect length (if transect survey)`)
  
  dhat.double   <- sum(detect.vec, na.rm=T)*eS/phat/trans.area
  
  df <- data.frame(Detections=detect.vec, Mussels=count.vec[1,], D=detect.vec*eS/area.vec/phat, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, Area=area.vec, Length=area.vec, Type=rep("Double", length(trans.area)), Lake=rep(curr.lake, length(trans.area)), Lat=transect.dat$`GPS latitude`, Long=transect.dat$`GPS longitude`)
  
  if(var.type == "design") {
    K  <- length(detect.vec)
    l  <- area.vec
    L  <- sum(l)
    n  <- detect.vec
    N  <- sum(detect.vec)
    #var.n <- L^2*(K/(L^2*(K-1))*sum(l^2*(n/l - NL)^2))/N^2
    var.n <- L/(K-1)*sum(l*(n/l - N/L)^2)
    
    dhat.se <- sqrt(dhat.double^2*(var.n/N^2 + varS/K/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat.double, Dhat.se=dhat.se, phat=phat, phat.se=phat.se, Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model, eS=eS, varS=varS)
  } 

  if(var.type == "jack") {
    
    K       <- length(detect.vec)
    l       <- area.vec/2
    L       <- sum(l)
    dhat    <- sum(detect.vec*eS/phat)/sum(area.vec)
    dhat.detect    <- sum(detect.vec/phat)/sum(area.vec)
    
    dhat.se.K <- dhat.K <- n.K <- pseudo.vec <- pseudo.nhat <- vector('numeric', K)
    
    for(k in 1:K) {
      n.K[k]        <- sum(detect.vec[-k])
      dhat.K[k]     <- sum(detect.vec[-k]/phat)/sum(area.vec[-k])
      pseudo.vec[k] <- (L*dhat.detect - (L-l[k])*dhat.K[k])/l[k] #K*dhat.detect - (K-1)*dhat.K[k]
    }
    
    dhat.se.jk  <- sqrt(1/(L*K-L)*sum(l*(pseudo.vec - dhat.detect)^2) )
    
    b.jk        <- (K-1)*(mean(dhat.K) - dhat.detect)
    dhat.jk     <- K*dhat.double - (K-1)*mean(dhat.K) 
    
    #dhat.se.full <- sqrt(dhat.se.jk^2*sum(detect.vec)^2/dhat^2/sum(detect.vec)^2 + dhat.jk^2*(varS/eS^2 + phat.se^2/phat^2))
    dhat.se.full <- sqrt(dhat.se.jk^2 + dhat.double^2*(varS/K/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat.double, Dhat.se=dhat.se.full, phat=phat, phat.se=phat.se, Transects=length(area.vec), Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model, eS=eS, varS=varS)
  }
  
  if(var.type == "model") {
    double.glm      <- glm.nb(detect.vec~1 + offset(log(phat*area.vec)))
    double.predict  <- predict(double.glm, type="response", se.fit=T)
    
    se.double     <- sqrt(dhat.double^2*(sum(double.predict$se.fit^2)/sum(double.predict$fit)^2 + varS/length(detect.vec)/eS^2 + phat.se^2/phat^2))#+ varS/eS^2
    
    return.list <- list(Dhat=dhat.double, Dhat.se=se.double, phat=phat, phat.se=phat.se, Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model, eS=eS, varS=varS)
  }
  
  return(return.list)
}



double.dens.sim <- function(curr.lake, curr.lake2, var.type, sd.width=0.1) {
  
  double.dat      <- read_xlsx(path="../Data/Season2/Encounters - Double observer - no distance (Responses).xlsx", sheet=1)
  double.dat      <- double.dat %>% subset(`Lake name` == curr.lake)
  double.dat      <- double.dat[order(double.dat$`Transect #`),]
  double.dat      <- dplyr::rename(double.dat, size="Number of mussels in cluster")
  doubleorig.dat <- double.dat
  
  if(any(double.dat$size==0)) {
    double.dat <- double.dat[-which(double.dat$size==0),]
  }
  if(any(is.na(double.dat$size))) {
    double.dat <- double.dat[-which(is.na(double.dat$size)),]
  }
  
  eS   <- mean(double.dat$size, na.rm=T)
  varS <- var(double.dat$size, na.rm=T) #/length(double.dat$size) #adding term on for numeric purposes
  
  #transect.title  <- suppressMessages(gs_title("Transect datasheets (Responses)"))
  #transect.dat    <- suppressMessages(gs_read(transect.title, verbose=FALSE))
  transect.dat      <- read_xlsx(path="../Data/Season2/Transect datasheets (Responses).xlsx", sheet=1)
  transect.dat <- transect.dat %>% subset(`Lake name:` == curr.lake2)
  transect.dat <- transect.dat %>% subset(`Survey type` == "Double observer no distance")
  transect.dat <- transect.dat[order(transect.dat$"Transect number"),]
  
  setup.time   <- 60^2*hour(transect.dat$`Setup time`) + 60*minute(transect.dat$`Setup time`) + second(transect.dat$`Setup time`)
  hab.time     <- 60^2*hour(transect.dat$`Habitat time`) + 60*minute(transect.dat$`Habitat time`) + second(transect.dat$`Habitat time`)
  enc.time     <- 60^2*hour(transect.dat$`Encounter time`) + 60*minute(transect.dat$`Encounter time`) + second(transect.dat$`Encounter time`)
  
  trans.length  <- transect.dat$`Transect length (if transect survey)`
  trans.area    <- trans.length
  
  
  #get number of obsevations made by primary and secondary observers
  count.primary <- count.secondary <- detect.primary <- detect.secondary <- rep(NA, max(transect.dat$`Transect number`))
  for(i in unique(transect.dat$`Transect number`)) {
    curr.prim <- transect.dat[transect.dat$`Transect number` == i,]$`Primary observer (double observer survey)`
    curr.dat  <- double.dat[double.dat$`Transect #`== i,] 
    
    count.primary[i]    <- max(sum(curr.dat[curr.dat$`Observer name` == curr.prim,]$size), 0, na.rm=T)
    count.secondary[i]  <- max(sum(curr.dat[curr.dat$`Observer name` != curr.prim,]$size), 0, na.rm=T)
    
    detect.primary[i]    <- max(dim(curr.dat[curr.dat$`Observer name` == curr.prim,])[1], 0, na.rm=T)
    detect.secondary[i]  <- max(dim(curr.dat[curr.dat$`Observer name` != curr.prim,])[1], 0, na.rm=T)
    
  }
  
  if(any(is.na(count.primary))) {
    rm.vec          <- which(is.na(count.primary))
    count.primary   <- count.primary[-rm.vec]
    count.secondary <- count.secondary[-rm.vec]
    detect.primary  <- detect.primary[-rm.vec]
    detect.secondary <- detect.secondary[-rm.vec]
  }
  
  double.dat <- double.dat %>% mutate(detected = rep(1, dim(double.dat)[1]), distance=rep(0, dim(double.dat)[1]), object=1:dim(double.dat)[1])
  
  #double.dat <- doubleObs.duplicate(obs.dat=double.dat)
  double.dat$object <- 1:dim(double.dat)[1]
  double.dat <- create.removal.Observer(transect.dat=transect.dat, obs.dat=double.dat)
  
  #doubleDetect.model <- ddf(method="rem", mrmodel=~glm(formula=~1), dsmodel=~cds(key="unif"), data=double.dat, meta.data=list(width=1))
  
  #phat    <- summary(doubleDetect.model)$average.p
  #phat.se <- summary(doubleDetect.model)$average.p.se[1,1]
  doubleDetect.model <- ddf(method="rem.fi", mrmodel=~glm(formula=~1), dsmodel=~cds(key="unif"), data=double.dat, meta.data=list(width=1))
  phat    <- summary(doubleDetect.model)$average.p0.1
  phat.se <- summary(doubleDetect.model)$average.p0.1.se[1,1]
  #print(doubleDetect.model)
  
  count.vec <- count.primary + count.secondary
  detect.vec <- detect.primary + detect.secondary
  
  count.vec <- t(as.matrix(count.vec))
  dimnames(count.vec)[[2]] <- transect.dat$`Transect number`
  
  area.vec        <- trans.length * rnorm(length(trans.length), 1, sd.width)
  names(area.vec) <- transect.dat$`Transect number`
  
  trans.area <- sum(area.vec)
  
  dhat.double   <- sum(detect.vec, na.rm=T)*eS/phat/trans.area
  
  df <- data.frame(Detections=detect.vec, Mussels=count.vec[1,], D=detect.vec*eS/area.vec/phat, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, Area=area.vec, Length=trans.length, Type=rep("Double", length(trans.area)), Lake=rep(curr.lake, length(trans.area)))
  
  if(var.type == "design") {
    K  <- length(detect.vec)
    l  <- area.vec
    L  <- sum(l)
    n  <- detect.vec
    N  <- sum(detect.vec)
    #var.n <- L^2*(K/(L^2*(K-1))*sum(l^2*(n/l - NL)^2))/N^2
    var.n <- L/(K-1)*sum(l*(n/l - N/L)^2)
    
    dhat.se <- sqrt(dhat.double^2*(var.n/N^2 + varS/K/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat.double, Dhat.se=dhat.se, phat=phat, phat.se=phat.se, Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model, eS=eS, varS=varS)
  } 
  
  if(var.type == "jack") {
    
    K       <- length(detect.vec)
    l       <- area.vec/2
    L       <- sum(l)
    dhat    <- sum(detect.vec*eS/phat)/sum(area.vec)
    dhat.detect    <- sum(detect.vec/phat)/sum(area.vec)
    
    dhat.se.K <- dhat.K <- n.K <- pseudo.vec <- pseudo.nhat <- vector('numeric', K)
    
    for(k in 1:K) {
      n.K[k]        <- sum(detect.vec[-k])
      dhat.K[k]     <- sum(detect.vec[-k]/phat)/sum(area.vec[-k])
      pseudo.vec[k] <- (L*dhat.detect - (L-l[k])*dhat.K[k])/l[k] #K*dhat.detect - (K-1)*dhat.K[k]
    }
    
    dhat.se.jk  <- sqrt(1/(L*K-L)*sum(l*(pseudo.vec - dhat.detect)^2) )
    
    b.jk        <- (K-1)*(mean(dhat.K) - dhat.detect)
    dhat.jk     <- K*dhat.double - (K-1)*mean(dhat.K) 
    
    #dhat.se.full <- sqrt(dhat.se.jk^2*sum(detect.vec)^2/dhat^2/sum(detect.vec)^2 + dhat.jk^2*(varS/eS^2 + phat.se^2/phat^2))
    dhat.se.full <- sqrt(dhat.se.jk^2 + dhat.double^2*(varS/K/eS^2 + phat.se^2/phat^2))
    
    return.list <- list(Dhat=dhat.double, Dhat.se=dhat.se.full, phat=phat, phat.se=phat.se, Transects=length(area.vec), Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model, eS=eS, varS=varS)
  }
  
  if(var.type == "model") {
    double.glm      <- glm.nb(detect.vec~1 + offset(log(phat*area.vec)))
    double.predict  <- predict(double.glm, type="response", se.fit=T)
    
    se.double     <- sqrt(dhat.double^2*(sum(double.predict$se.fit^2)/sum(double.predict$fit)^2  + phat.se^2/phat^2))#+ varS/eS^2
    
    return.list <- list(Dhat=dhat.double, Dhat.se=se.double, phat=phat, phat.se=phat.se, Area=area.vec, Detections=detect.primary+detect.secondary, Mussels=count.vec, Time=setup.time+hab.time+enc.time, t.set=setup.time, t.hab=hab.time, t.enc=enc.time, df=df, ddf=doubleDetect.model)
  }
  
  return(return.list)
}



phase2.countrate.pop <- function(phaseTwo.dat, lake.name) {
  
  currLake.phase2 <- subset(phaseTwo.dat, phaseTwo.dat$`Lake name` == lake.name)
  diff1           <- as.numeric(difftime(currLake.phase2$`End time`, currLake.phase2$`Start time`, units="secs"))
  
  currLake.countrate <- sum(currLake.phase2$`Total number of zebra mussels detected`)/sum(diff1 - currLake.phase2$`Break time (duration in minutes)`*60)
  
  return(currLake.countrate)
}

phase2.countrate.transect <- function(phaseTwo.dat, lake.name, distance.est, double.est, quadrat.est) {
  
  currLake.phase2 <- subset(phaseTwo.dat, phaseTwo.dat$Lake == lake.name)

  currLake.phase2.obs1 <- currLake.phase2 %>% subset(., Observer == "Aislyn") 
  currLake.phase2.obs1 <- currLake.phase2.obs1[order(currLake.phase2.obs1$`Traverse number`),]
  currLake.phase2.obs2 <- currLake.phase2 %>% subset(., Observer == "Austin")
  currLake.phase2.obs2 <- currLake.phase2.obs2[order(currLake.phase2.obs2$`Traverse number`),]
  
  diff1           <- as.numeric(difftime(currLake.phase2.obs1$`End time`, currLake.phase2.obs1$`Start time`, units="secs"))
  diff2           <- as.numeric(difftime(currLake.phase2.obs2$`End time`, currLake.phase2.obs2$`Start time`, units="secs"))
  
  newDF <- data.frame(N1=currLake.phase2.obs1$`Total number of zebra mussels detected`, N2=currLake.phase2.obs2$`Total number of zebra mussels detected`, t1=diff1, t2=diff2, Transect=currLake.phase2.obs1$`Traverse number`, N=currLake.phase2.obs1$`Total number of zebra mussels detected`+ currLake.phase2.obs2$`Total number of zebra mussels detected`, T=diff1+diff2)

  newDF <- newDF %>% mutate(D.distance = NA, D.double=NA, D.quadrat=NA)
  
  dist.trans <- as.numeric(names(distance.est$Mussels))
  doub.trans <- as.numeric(names(double.est$Mussels[1,]))
  quad.trans <- as.numeric(names(quadrat.est$Mussels))
  
  z.dist <- intersect(dist.trans, newDF$Transect)
  z.doub <- intersect(doub.trans, newDF$Transect)
  z.quad <- intersect(quad.trans, newDF$Transect)
  
  dist.ind <- match(z.dist, newDF$Transect)
  doub.ind <- match(z.doub, newDF$Transect)
  quad.ind <- match(z.quad, newDF$Transect)
  
  newDF$D.distance[dist.ind]  <- distance.est$Mussels/distance.est$Area/distance.est$phat
  newDF$D.double[doub.ind]    <- double.est$Mussels/double.est$Area/double.est$phat
  newDF$D.quadrat[quad.ind]   <- quadrat.est$Mussels/quadrat.est$Area
  
  if(any(newDF$T < 0)) {
    newDF[-which(newDF$T < 0),]  
  }
  
  return(newDF)
}

chisq <- function(fm) { 
  observed <- getY(fm@data) 
  expected <- fitted(fm) 
  return(sum((observed - expected)^2/expected))
}

EstimatorComparisonPlot <- function(distance.est.design, distance.est.jack, distance.est.model, double.est.design, double.est.jack, double.est.model, quadrat.est.design, quadrat.est.jack, quadrat.est.model, lakename="Burgan") {
  
  library(RColorBrewer)
  col.vec <- brewer.pal(3, "Set1")
  dhat.distance.vec    <- c(distance.est.design$Dhat, distance.est.jack$Dhat, distance.est.model$Dhat)
  dhat.distance.se.vec <- c(distance.est.design$Dhat.se, distance.est.jack$Dhat.se, distance.est.model$Dhat.se)
  
  dhat.double.vec    <- c(double.est.design$Dhat, double.est.jack$Dhat, double.est.model$Dhat)
  dhat.double.se.vec <- c(double.est.design$Dhat.se, double.est.jack$Dhat.se, double.est.model$Dhat.se)
  
  dhat.quadrat.vec    <- c(quadrat.est.design$Dhat, quadrat.est.jack$Dhat, quadrat.est.model$Dhat)
  dhat.quadrat.se.vec <- c(quadrat.est.design$Dhat.se, quadrat.est.jack$Dhat.se, quadrat.est.model$Dhat.se)
  ylims <- c(dhat.double.vec - dhat.double.se.vec, dhat.quadrat.vec - dhat.quadrat.se.vec, dhat.distance.vec - dhat.distance.se.vec, dhat.double.vec + dhat.double.se.vec, dhat.quadrat.vec + dhat.quadrat.se.vec, dhat.distance.vec + dhat.distance.se.vec)
  
  plotCI(x=rep(1,3) + c(-0.1, 0, 0.1), y=dhat.distance.vec, uiw=2*dhat.distance.se.vec, xaxt = 'n', ylab="Estimated density", xlab="Estimator", col=col.vec, pch=19, lwd=2, cex=1.3, sfrac=0, cex.lab=1.3, main=lakename, xlim=c(0.75, 3.25), ylim=range(ylims))
  axis(1, at=c(1,2,3), labels=c("Distance", "Double observer", "Quadrat"))
  abline(h=mean(c(dhat.double.vec, dhat.distance.vec, dhat.quadrat.vec)), lty=3, col="gray")
  plotCI(x=rep(2,3)+c(-0.1,0,0.1), y=dhat.double.vec, uiw=2*dhat.double.se.vec, col=col.vec, pch=19, lwd=2, cex=1.3, sfrac=0, cex.lab=1.3, add=T)
  
  plotCI(x=rep(3,3)+c(-0.1,0,0.1), y=dhat.quadrat.vec, uiw=2*dhat.quadrat.se.vec, col=col.vec, pch=19, lwd=2, cex=1.3, sfrac=0, cex.lab=1.3, add=T)
  
  legend('topleft', legend=c('Design-based estimate', 'Jack-knife estimate', 'Model-based estimate'), pch=19, col=col.vec)
}

rem.logLik <- function(pars, counts, fi=FALSE) {
  p <- plogis(pars[1])
  
  if(fi == TRUE) {
    delta <- 1
  } else { delta <- tanh(pars[2]) }
  
  pdot <- 2*p - delta*p^2
  
  #w11 <- p^2*delta
  #w01 <- p*(1-p*delta)
  
  w11 <- p/(2*p - p^2)
  w01 <- p*(1-p)/(2*p - p^2)
  #w01 <- p*(1-p*delta)
  
  logLik <- counts[1]*log(w11) + counts[2]*log(w01) - sum(counts)*log(delta)
  
  return(-logLik)
} 


rem.optim <- function(counts) {
  pars <- 1 - sum(counts[2]/sum(counts[1]))
  pars[1] <- log(pars[1]/(1 - pars[1]))
  opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8)
  res  <- nloptr(x0=pars, eval_f=rem.logLik, opts=opts, counts=counts, fi=TRUE)
  res$hessian <- hessian(rem.logLik, x=res$solution, counts=counts, fi=TRUE)
  phat <- plogis(res$solution)
}

##
time.predict <- function(time.df, curr.lake="Lake Burgan") {
  
  #first look at which times depend on transect length
  #tset.length <- lm(t.set ~ Length + Lake + Type + D.dens, data=time.df)
  #thab.length <- lm(t.hab ~ Length + Lake + Type + D.dens, data=time.df)
  #tenc.length <- lm(t.enc ~ Length + Lake + Type + D.dens, data=time.df)
  tset.length <- lm(t.set ~ Length + Type + (1|Lake), data=time.df)
  thab.length <- lm(t.hab ~ Length + Type + (1|Lake), data=time.df)
  tenc.length <- lm(t.enc ~ Length + Type + (1|Lake), data=time.df)
  
  #We'll look at Lake Burgan first
  lake.index  <- which(time.df$Lake == curr.lake)
  newdat      <- data.frame(Length=rep(30, 3), Lake=rep(time.df$Lake[lake.index[1]], 3), Type=time.df$Type[c(1,19,37)], D.dens=sum(time.df$Detections[lake.index])/sum(time.df$Length[lake.index]))
  
  tset.pred <- predict(tset.length, newdata=newdat)
  thab.pred <- predict(thab.length, newdata=newdat)
  tenc.pred <- predict(tenc.length, newdata=newdat)
  
  names(tset.pred) <- names(tenc.pred)  <- names(thab.pred) <- time.df$Type[c(1,19,37)]
  
  return(list(setup=tset.pred, habitat=thab.pred, encounters=tenc.pred))
}


hazardDetect.predict <- function(ddf.obj) {
  
  xval      <- seq(0,1, length.out=100)
  p0        <- mean(predict(ddf.obj$mr)$fitted)
  ds.scale  <- as.numeric(exp(coef(ddf.obj$ds)$scale[1]))
  ds.exp  <- as.numeric(exp(coef(ddf.obj$ds)$exponent[1]))
  y <- p0*(1 - exp(-(xval/ds.scale)^(-ds.exp)))

  return(list(x=xval, y=y))
}


subSample.dist <- function(dist.mat) {
 
  dist.vec <- ss.vec <- NULL
  
  for(step in 1:7) {
    
    dist.temp <- NULL
    
    for(i in 1:(15-step)) {
      dist.temp <- c(dist.temp, dist.mat[i + step, i])
    }
    
    dist.vec <- c(dist.vec, dist.temp) #mean(dist.temp)
    ss.vec   <- c(ss.vec, rep(15 %/% step, length(dist.temp)))
    
  }
  return(list(Distance=dist.vec, SampleSize=ss.vec))
}