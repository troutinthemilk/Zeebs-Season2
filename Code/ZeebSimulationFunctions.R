# Simulate zeeb samples
# depends on packages: dplyr

############################################
# generate a zeeb population

zeebPop <- function(len, wid, lam1, lam2, disp = NULL, r.pois = NULL, r.fixed = NULL, type = "none")
{
  # type = "none", "pp" (poisson cluster process), "nbp" (poisson cluster process with negbinom size distribution), "matern" circular matern process
  
  # lake florida density ~ 0.05  use lam2 = 0, lam1 = 4000
  # lake burgan density ~ .5  use lam2 = .4, lam1 = 4000*6
  # little birch density ~ 25 use lam2 = 5, lam1 = 4000*80
  # len <- 20   # transect length
  # wid <- 4000  # width of lake (circumference)
  # lam1 <- 4000*80  # mean number of clusters
  # lam2 <- 5     # mean cluster size = lam2 + 1
  # disp = dispersion (size) param for neg binom (var = mu + mu^2/disp)
  # lam1*(lam2+1)  # expected number of mussels
  # lam1*(lam2+1)/(len*wid)  # expected density
  # r <- .02*1  # Poisson process radius: larger r means less clustering and more randomness (mean dist)
  # r.fixed <- 2*r # fixed max distance for Matern process
  
  if (type == "none")   # no spatial clustering, all clusters size 1
  {
    N.clus <- rpois(1,lambda = lam1) # number of clusters (stage 1)
    N.muss <- rpois(N.clus, lambda = lam2) + 1   # cluster sizes
    N <- sum(N.muss)
    sim.pop <- data.frame(cluster = 1:N, 
                         x.clus = runif(N, 0, wid), 
                         y.clus = runif(N, 0, len), 
                         size = 1)
    sim.pop$x <- sim.pop$x.clus
    sim.pop$y <- sim.pop$y.clus
  }
  
  if (type == "pp")   # poisson cluster process
  {
    N.clus <- rpois(1,lambda = lam1) # number of clusters (stage 1)
    N.muss <- rpois(N.clus, lambda = lam2) + 1   # number with in cluster (stage 2)
    x <- runif(N.clus, 0, wid)
    y <- runif(N.clus, 0, len)
    # random location of clusters
    sim.pop <- data.frame(cluster = rep(1:N.clus, times = N.muss),
                         x.clus = rep(x, times = N.muss), 
                         y.clus = rep(y, times = N.muss),  
                         size = rep(N.muss, times = N.muss))
    # r = radius, theta = angle
    N <- sum(N.muss)
    r.sim <- rexp(N,rate = 1/r.pois)  # random radius with mean r
    theta <- runif(N,0,2*pi)   # random angle
    sim.pop$x <- sim.pop$x.clus + r.sim*cos(theta)
    sim.pop$y <- sim.pop$y.clus + r.sim*sin(theta)
  }

  if (type == "nbp")   # neg binom cluster process
  {
    N.clus <- rpois(1,lambda = lam1) # number of clusters (stage 1)
    N.muss <- rnbinom(N.clus, mu = lam2, size = disp) + 1   # number with in cluster (stage 2)
    x <- runif(N.clus, 0, wid)
    y <- runif(N.clus, 0, len)
    # random location of clusters
    sim.pop <- data.frame(cluster = rep(1:N.clus, times = N.muss),
                          x.clus = rep(x, times = N.muss), 
                          y.clus = rep(y, times = N.muss),  
                          size = rep(N.muss, times = N.muss))
    # r = radius, theta = angle
    N <- sum(N.muss)
    r.sim <- rexp(N,rate = 1/r.pois)  # random radius with mean r
    theta <- runif(N,0,2*pi)   # random angle
    sim.pop$x <- sim.pop$x.clus + r.sim*cos(theta)
    sim.pop$y <- sim.pop$y.clus + r.sim*sin(theta)
  }
  
  if (type == "matern")   # matern process with fixed circular neighborhoods
  {
    N.clus <- rpois(1,lambda = lam1) # number of clusters (stage 1)
    N.muss <- rpois(N.clus, lambda = lam2) + 1   # number with in (stage 2) or negative binomial
    x <- runif(N.clus, 0, wid)
    y <- runif(N.clus, 0, len)
    # random location of clusters
    sim.pop <- data.frame(cluster = rep(1:N.clus, times = N.muss),
                         x.clus = rep(x, times = N.muss), 
                         y.clus = rep(y, times = N.muss),  
                         size = rep(N.muss, times = N.muss))
    # r = radius, theta = angle
    N <- sum(N.muss)
    r.sim <- runif(N, 0 , r.fixed)  # random radius 
    theta <- runif(N,0,2*pi)   # random angle
    sim.pop$x <- sim.pop$x.clus + r.sim*cos(theta)
    sim.pop$y <- sim.pop$y.clus + r.sim*sin(theta)
  }
  
  # restrict individual coords to wid/len bounds:
  sim.pop <- dplyr::filter(sim.pop, x > 0, x < wid, y > 0, y < len)
  sim.pop$type <- type
  return(dplyr::as.tbl(sim.pop))

}

##################################################

# function to take a samples from a pop

# Distance sampling (single and double observer)
# in: 
# zeebPop: population to sample from
# x.samp: x coordinates of sampled transects
# len: transect length
# sigma: half normal detection function parameter
# w: distance max width
# dtype: "double" (removal) or "single"
# cluster: TRUE means sample clusters (with size measured), FALSE means sampling individual mussels 

zeebDistSamp <- function(zeebPop, x.samp, len, sigma, w, dtype = "double", cluster = TRUE)
{
  if (!cluster)  # sample individuals
  { dist.samp <- zeebPop %>% 
    mutate(samp = apply(sapply(x.samp, function(x) between(zeebPop$x,x - w, x + w)), 1, sum)) %>%
    filter(samp ==1) 
  dist.samp <- dist.samp %>%
    mutate(transect = apply(sapply(x.samp, function(x) between(dist.samp$x,x - w, x + w)), 1, function(x) which(x == TRUE))) 
  # get transect ID
  dist.samp$transect.x <- x.samp[dist.samp$transect]
  dist.samp <- dist.samp %>% select(-samp) %>%
    mutate(dist = x - transect.x, 
           prob = exp(-dist^2/(2*sigma^2)))

  if (dtype == "single") # detection indicator
  { dist.samp <- dist.samp %>% 
    filter(rbinom(nrow(dist.samp),1,dist.samp$prob) == 1)  }
  if (dtype == "double")  
  {
    dist.samp$DoubleObserver <- sapply(dist.samp$prob, function(x) sample(c("firstObs","secondObs","neverObs"), size =1, replace = TRUE, prob = c(x, (1-x)*x, (1-x)^2) ))
    dist.samp <- dist.samp %>%
      filter(DoubleObserver != "neverObs")
  }
  }  # end of individual sample
  
  if (cluster)  # sample clusters
  { 
    # just get one response per cluster
    dist.samp <- zeebPop %>% group_by(cluster) %>% slice(1) 
    dist.samp <- dist.samp %>% ungroup() %>%
      mutate(samp = apply(sapply(x.samp, function(x) between(dist.samp$x.clus,x - w, x + w)), 1, sum)) %>%
    filter(samp ==1) 
  dist.samp <- dist.samp %>%
    mutate(transect = apply(sapply(x.samp, function(x) between(dist.samp$x.clus,x - w, x + w)), 1, function(x) which(x == TRUE))) 
  # get transect ID
  dist.samp$transect.x <- x.samp[dist.samp$transect]
  # probs based on cluster distance
  dist.samp <- dist.samp %>% select(-samp) %>%
    mutate(dist = x.clus - transect.x, 
           prob = exp(-dist^2/(2*sigma^2)))

  if (dtype == "single") # detection indicator
  { dist.samp <- dist.samp %>% 
    filter(rbinom(nrow(dist.samp),1,dist.samp$prob) == 1)  }
  if (dtype == "double")  
  {
    dist.samp$DoubleObserver <- sapply(dist.samp$prob, function(x) sample(c("firstObs","secondObs","neverObs"), size =1, replace = TRUE, prob = c(x, (1-x)*x, (1-x)^2) ))
    dist.samp <- dist.samp %>%
      filter(DoubleObserver != "neverObs")
  }
  }  # end of cluster sample

  dist.samp$dtype <- dtype 
  

  
  return(as.tbl(dist.samp))
  
}


############################################

# Double observer removal transect sampling
# in: 
# zeebPop: population to sample from
# x.samp: x coordinates of sampled transects
# len: transect length
# p: individual mussel detection prob
# beta1: cluster probability parameters where p(cluster detect) = logistic(beta0 + beta1*size) with size of cluster. Set so that p(cluster|size =1) = p
# w.dt: distance of strip sampled from transect
# cluster: TRUE means sample clusters (with size measured), FALSE means sampling individual mussels 

zeebDoubleRemSamp <- function(zeebPop, x.samp, len, p, beta1, w.dt, cluster = TRUE)
{
  if (!cluster)  # sample individuals
  { 
  dt.samp <- zeebPop %>% 
    mutate(samp = apply(sapply(x.samp, function(x) between(zeebPop$x,x - w.dt, x + w.dt)), 1, sum)) %>%
    filter(samp ==1) 
  dt.samp <- dt.samp %>%
    mutate(transect = apply(sapply(x.samp, function(x) between(dt.samp$x,x - w.dt, x + w.dt)), 1, function(x) which(x == TRUE))) 
  # get transect ID
  dt.samp$transect.x <- x.samp[dt.samp$transect]
  dt.samp <- dt.samp %>% select(-samp) %>%
    mutate(dist = x - transect.x)
  # get sample
  dt.samp$DoubleObserver <-  sample(c("firstObs","secondObs","neverObs"), size =nrow(dt.samp), replace = TRUE, prob = c(p, (1-p)*p, (1-p)^2))
  dt.samp <- dt.samp %>%
    filter(DoubleObserver != "neverObs")
  }
  
  if (cluster)  # sample clusters
  { 
    dt.samp <- zeebPop %>% group_by(cluster) %>% slice(1) 
    dt.samp <- dt.samp %>% ungroup() %>%
      mutate(samp = apply(sapply(x.samp, function(x) between(dt.samp$x.clus,x - w.dt, x + w.dt)), 1, sum)) %>%
      filter(samp ==1) 
    dt.samp <- dt.samp %>%
      mutate(transect = apply(sapply(x.samp, function(x) between(dt.samp$x.clus,x - w.dt, x + w.dt)), 1, function(x) which(x == TRUE))) 
    # get transect ID
    dt.samp$transect.x <- x.samp[dt.samp$transect]
    dt.samp <- dt.samp %>% select(-samp) %>%
      mutate(dist = x.clus - transect.x)
    # get probs
    beta0 <- log(p/(1-p)) - beta1   # intercept so p.clus = p when size = 1.
    dt.samp$p.clus <- plogis(beta0 + beta1*dt.samp$size)
    # get sample
    dt.samp$DoubleObserver <- sapply(dt.samp$p.clus, function(x) sample(c("firstObs","secondObs","neverObs"), size =1, replace = TRUE, prob = c(x, (1-x)*x, (1-x)^2) ))
    dt.samp <- dt.samp %>%
      filter(DoubleObserver != "neverObs")
     
  }
  
  # # any 0 transects to fill in?
  # no.detect <- !(1:length(x.samp) %in% unique(dt.samp$transect))
  # if (sum(no.detect) > 0)
  # {
  #   missing.trans <- (1:length(x.samp))[no.detect]
  #   n.data <- nrow(dt.samp)
  #   dt.samp[(n.data+1):(n.data+length(missing.trans)),] <- NA
  #   dt.samp[(n.data+1):(n.data+length(missing.trans)),"transect"] <- missing.trans
  #   dt.samp[(n.data+1):(n.data+length(missing.trans)),"size"] <- 0
  #   
  # }
  
  return(as.tbl(dt.samp))
  
}


############################################

# Quadrat removal transect sampling
# in: 
# zeebPop: population to sample from
# x.samp: x coordinates of sampled transects
# len: transect length
# p: individual mussel detection prob (set to 1 as default)
# beta1: cluster probability parameters where p(cluster detect) = logistic(beta0 + beta1*size) with size of cluster. Set so that p(cluster|size =1) = p (NULL as default)
# quad.dim:  square quad side size (0.5)
# x.dist.btw: horizontal distance between quads, centered on transect
# y.dist.btw: vertical distance between quad lower sides (distance unsampled is y.dist.btw - quad.dim)
# cluster: TRUE means sample clusters (with size measured), FALSE means sampling individual mussels 

zeebQuadSamp <- function(zeebPop, x.samp, len, p = 1, beta1 = 0, quad.dim, x.dist.btw, y.dist.btw, cluster = TRUE)
{
  quad.along <- (0:(len/y.dist.btw - 1) * y.dist.btw)  # starting value of quad along transect
  
  if (!cluster)  # sample individuals
  { 
    # get mussels within x (horizontal) bounds along transect
    quad.samp <- zeebPop %>% 
      mutate(samp = apply(sapply(x.samp, function(x) between(zeebPop$x,x - quad.dim - x.dist.btw/2, x + quad.dim + x.dist.btw/2)), 1, sum)) %>%
      filter(samp ==1) 
    quad.samp <- quad.samp %>%
      mutate(transect = apply(sapply(x.samp, function(x) between(quad.samp$x,x - quad.dim - x.dist.btw/2, x + quad.dim + x.dist.btw/2)), 1, function(x) which(x == TRUE))) 
    # get transect ID
    quad.samp$transect.x <- x.samp[quad.samp$transect]
    quad.samp <- quad.samp %>% select(-samp) %>%
      mutate(dist = x - transect.x)
    # get mussels within quads along transect  (keeps all mussels and indicates quad mussels)
    quad.samp <- quad.samp %>% 
      mutate(inquad = apply(sapply(quad.along, function(x) between(quad.samp$y,x , x + quad.dim)), 1, sum) + # get y-indicator
               between(quad.samp$dist, -quad.dim - x.dist.btw/2, -  x.dist.btw/2) +  # get left x indivcator
               between(quad.samp$dist,  x.dist.btw/2,   quad.dim + x.dist.btw/2)) # get right x indicator
    quad.samp <- quad.samp %>%
      mutate(quad = quad.along[1+floor(quad.samp$y) %/% 2], 
             inquad = ifelse(quad.samp$inquad == 2, "in quad", "not in quad"))
    quad.samp <- quad.samp %>% 
      filter(inquad == "in quad") %>%
      select(-inquad)
    
    if (p < 1)  # imperfect detection within quads
    {
      quad.samp <- quad.samp %>% group_by(transect, quad) %>%
        filter(rbinom( n(), 1, p) == 1)
    }
  }
  
  if (cluster)  # sample clusters
  { 
    quad.samp <- zeebPop %>% group_by(cluster) %>% slice(1) 
    # get mussels within x (horizontal) bounds along transect
    quad.samp <- quad.samp %>% ungroup() %>%
      mutate(samp = apply(sapply(x.samp, function(x) between(quad.samp$x.clus,x - quad.dim - x.dist.btw/2, x + quad.dim + x.dist.btw/2)), 1, sum)) %>%
      filter(samp ==1) 
    quad.samp <- quad.samp %>%
      mutate(transect = apply(sapply(x.samp, function(x) between(quad.samp$x.clus,x - quad.dim - x.dist.btw/2, x + quad.dim + x.dist.btw/2)), 1, function(x) which(x == TRUE))) 
    # get transect ID
    quad.samp$transect.x <- x.samp[quad.samp$transect]
    quad.samp <- quad.samp %>% select(-samp) %>%
      mutate(dist = x.clus - transect.x)
    # get mussels within quads along transect  (keeps all mussels and indicates quad mussels)
    quad.samp <- quad.samp %>% 
      mutate(inquad = apply(sapply(quad.along, function(x) between(quad.samp$y.clus,x , x + quad.dim)), 1, sum) + # get y-indicator
               between(quad.samp$dist, -quad.dim - x.dist.btw/2, -  x.dist.btw/2) +  # get left x indivcator
               between(quad.samp$dist,  x.dist.btw/2,   quad.dim + x.dist.btw/2)) # get right x indicator
    quad.samp <- quad.samp %>%
      mutate(quad = quad.along[1+floor(quad.samp$y.clus) %/% 2], 
             inquad = ifelse(quad.samp$inquad == 2, "in quad", "not in quad"))
    quad.samp <- quad.samp %>% 
      filter(inquad == "in quad") %>%
      select(-inquad)
    
    if (p < 1)  # imperfect detection within quads
    {
      # get probs
      beta0 <- log(p/(1-p)) - beta1   # intercept so p.clus=p when size=1.
      quad.samp$p.clus <- plogis(beta0 + beta1*quad.samp$size)
      quad.samp <- quad.samp %>% group_by(transect, quad) %>%
        filter(rbinom( n(), 1, p) == 1)
    }
    
  }
  
  return(as.tbl(quad.samp))
}


########################
## Clean up data: makes the sample data generated by sampling functions look like the actual data collected during surveys:

## Distance data
# omit: x.clus/y.clus or x/y, prob
## Distance data
# omit: x.clus/y.clus or x/y, prob

# in:
# design: "distance", "removal", "quad" 
cleanSample <- function(zsamp, design, cluster = TRUE)
{
  #### Distance
  if (design == "distance")
  {
   if (cluster)
   {
     zsamp.clean <- zsamp %>% 
       select(-prob, -x, -y, -type)
   }
   if (!cluster)  # individuals
    {
      zsamp.clean <- zsamp %>% 
        select(-prob, -x.clus, -y.clus, -cluster, - size, -type)
    }
  }
  
  #### double transect removal
  if (design == "removal")
  {
    if (cluster)
    {
      zsamp.clean <- zsamp %>% 
        select( -x, -y, -x.clus, -dist, -type, -p.clus)
    }
    if (!cluster)  # individuals
    {
      zsamp.clean <- zsamp %>% 
        select( -x.clus, -y.clus, -cluster, - size, -x, -dist, -type)
    }
  }
  
  #### quadrat removal
  if (design == "quad")
  {
    zsamp <- zsamp %>%
      mutate(side = ifelse(dist < 0, "left", "right"))
    if (cluster)
    {
      zsamp.clean <- zsamp %>% 
        select( -x, -y, -x.clus, -dist, -type)
      if ("p.clus" %in% names(zsamp.clean))
        zsamp.clean <- zsamp.clean %>% select(-p.clus)
    }
    if (!cluster)  # individuals
    {
      zsamp.clean <- zsamp %>% 
        select( -x.clus, -y.clus, -cluster, - size, -x, -dist, -type)
    }
  }
  
  
  return(as.tbl(zsamp.clean))
}


###########
#### Fill in NA/0 values for transects/quads with no detections
fillMissing <- function(zsamp, x.samp, design, cluster = TRUE, len = NULL, y.dist.btw = NULL )
{
  if (design != "quad")
  {
  # any 0 transects to fill in?
  no.detect <- !(1:length(x.samp) %in% unique(zsamp$transect))
  if (sum(no.detect) > 0)
  {
    missing.trans <- (1:length(x.samp))[no.detect]
    n.data <- nrow(zsamp)
    zsamp[(n.data+1):(n.data+length(missing.trans)),] <- NA
    zsamp[(n.data+1):(n.data+length(missing.trans)),"transect"] <- missing.trans
    zsamp[(n.data+1):(n.data+length(missing.trans)),"transect.x"] <-x.samp[missing.trans]
    if (cluster)
        zsamp[(n.data+1):(n.data+length(missing.trans)),"size"] <- 0
  }
  }
  
  if (design == "quad")
  {
    quad.along <- (0:(len/y.dist.btw - 1) * y.dist.btw)  # starting value of quad along transect
     frame <-  expand(data_frame(), quad = quad.along, transect = 1:length(x.samp), side  = c("left","right"))
     frame$transect.x <- x.samp[frame$transect]
     if (!cluster)
        zsamp <- full_join(zsamp, frame) 
     if (cluster)
       {
          frame$size <- 0
          zsamp <- full_join(zsamp, frame) 
     }
      
  }
  
  return(as.tbl(zsamp))
  
}


###################   Estimates ##############################
##  Double Distance with removal
distEst <- function(zsamp, len, w, n.tran, wid, type = "model.mrds", cluster = TRUE)
{
  area <- n.tran * len * w*2
  
  if (type == "model.mrds")
  {
    zsamp <- zsamp %>%
      mutate(object = row_number(), 
             observer = ifelse(DoubleObserver == "firstObs", 1, 2), 
             detected = 1, 
             distance = dist, 
             Sample.Label = transect, Effort = len, Region.Label = 1)
    frame <- tidyr::expand(data_frame(), object = 1:nrow(zsamp), observer = c(1,2))

    if (!cluster)  # individual 
      zsamp.rem <- 
        full_join(select(zsamp, object, distance, Sample.Label, Effort, Region.Label), frame, by = c("object")) %>%
        full_join(select(zsamp, detected, observer, object),by = c("object", "observer")) 
    if (cluster)  # cluster
      zsamp.rem <- 
        full_join(select(zsamp, object, distance, size, Sample.Label, Effort, Region.Label), frame, by = c("object")) %>%
        full_join(select(zsamp, detected, observer, object),by = c("object", "observer")) 

    zsamp.rem$detected[is.na(zsamp.rem$detected)] <- 0
    dist.out <- ddf(method = "rem", dsmodel = ~cds(key = "hn"), mrmodel = ~glm(formula = ~1), data = zsamp.rem, meta.data = list(width = w))
    # comment of ddf: detection theta = "scale" is such that sigma = exp(theta)
    #Dest <- dist.out$Nhat/area  # Nhat is number of clusters
    
    # alt method
    dht.out <- dht(dist.out, 
        region.table = data.frame(Region.Label = 1, Area = len*wid), 
        sample.table = data.frame(Region.Label =1, Sample.Label = 1:n.tran, Effort = len))
    Dest <- dht.out$individuals$D$Estimate

  }
  # use unmarked?? I don't think so since I can't find a way to map individual level covariates to detection probs. Only uses site or design based covariates. And distance function assumes J distance intervals instead of a continuous distance covariate.

  
  return(data_frame(Dest, design = "doubleDistance", type, cluster))
  
}

## Double transect removal (no distance)
# designbased depends on FSA package
dtEst <- function(zsamp, len, w.dt, n.tran, type = "designEst.FSA", cluster = TRUE)
{
  area <- n.tran * len * w.dt*2
  
  if (!cluster) 
  {
    zsamp.wide <- zsamp %>% 
      filter(complete.cases(zsamp)) %>% 
      group_by(transect, DoubleObserver) %>%
      summarize(n = n()) %>% 
      spread(key = DoubleObserver, value = n) 
  }
  if (cluster) 
  {
    zsamp.wide <- zsamp %>% 
      filter(complete.cases(zsamp)) %>% 
      group_by(transect, DoubleObserver) %>%
      summarize(n = sum(size)) %>% 
      spread(key = DoubleObserver, value = n) 
  }    
  zsamp.wide[is.na(zsamp.wide)] <- 0
  if (! "secondObs" %in% names(zsamp.wide))
    zsamp.wide$secondObs <- 0
  if (! "firstObs" %in% names(zsamp.wide))
    {  
        zsamp.wide$firstObs <- 0
        zsamp.wide <- zsamp.wide %>% select(1,3,2) 
  }
  
  #add 0 counts for no detection transects
  miss <- which(! 1:n.tran %in% unique(zsamp.wide$transect))
  zsamp.wide <- zsamp.wide %>% bind_rows(data_frame(transect = miss, firstObs = 0, secondObs = 0))

  
   if (type == "designEst.FSA")
  {
    n.list <- split(as.matrix(zsamp.wide[,-1]), f = zsamp.wide$transect)
    remove.out <- lapply(n.list, removal, just.ests = TRUE)
    Dest <- sum(data.frame(remove.out)[1,])/area
  }
  
  if (type == "MultPois.unmarked")
  {
    count.mat <- as.matrix(zsamp.wide[,-1])
    remFrame <- unmarkedFrameMPois(y = count.mat, type = "removal")
    if (mean(zsamp.wide$firstObs - zsamp.wide$secondObs) < 5)
      remove.out <- multinomPois(~ 1 ~ 1, data = remFrame, starts = c(3,-5), control = list(maxit = 1000))
    if (mean(zsamp.wide$firstObs - zsamp.wide$secondObs) >= 5)
      remove.out <- multinomPois(~ 1 ~ 1, data = remFrame, starts = c(3,0), control = list(maxit = 1000))
    abund <- backTransform(linearComb(remove.out, type = "state", coefficients = matrix(c(rep(1,nrow(count.mat))), byrow=FALSE, nrow=nrow(count.mat)) ))  
    # alt: abund <- predict(remove.out, type = "state")$Predicted
    Dest <- sum(abund@estimate)/area
  }
  
  return(data_frame(Dest, design = "doubleRemoval", type, cluster))
  
}

## Quad
quadEst <- function(zsamp, quad.dim, len, y.dist.btw, n.tran, type = "designEst", cluster = TRUE)
{
  if (type == "designEst")
  {
    # area = 2 x quads/side x transects x quad area
    area <- length(0:(len/y.dist.btw - 1) * y.dist.btw)*2*n.tran*quad.dim^2
    if (!cluster)
    {
      tot <- nrow(filter(zsamp , complete.cases(zsamp)))
      Dest <- tot/area
    }
    if (cluster)
    {
      tot <- sum(zsamp$size, na.rm = TRUE)
      Dest <- tot/area
    }
  }
  return(data_frame(Dest, design = "quad", type, cluster))
}