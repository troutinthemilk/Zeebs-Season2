#this function creates a removal survey object that can be analyzed in mrds.
#primary coded as 1
#secondary coded as 2
#secondary sees all that 1 does but then a bit more

suppressPackageStartupMessages(library(tidyverse)) #loads the tidyverse quietly
library(readxl)
library(FSA)
library(mrds)


#' Create data object for removal analysis with distance
#'
#' A helper function to format the data appropriately for a ddf analysis using the R package `MRDS`.
#' @param transect.dat Data table of each transect.
#' @param obs.dat Data table observations in each transect.
#' @return Data table with counts summarized on each transect
#' @export
#' @examples
#' create.removal.Observer(transect.dat=transect.dat, obs.dat=distance.dat)
create.removal.Observer <- function(transect.dat, obs.dat) {

  primary   <- transect.dat$`Primary observer`
  secondary <- transect.dat$`Secondary observer`

  obs.dat <- obs.dat %>% mutate(primary=primary[obs.dat$`Transect #`- min(obs.dat$`Transect #`)+1], secondary=secondary[obs.dat$`Transect #`- min(obs.dat$`Transect #`)+1], observer = rep(1, dim(obs.dat)[1]))

  for(i in 1:dim(obs.dat)[1]) {
    curr.tran             <- obs.dat[i,]$`Transect #`
    obs.dat[i,]$observer  <- as.numeric(obs.dat[i,]$`Observer name` == transect.dat[transect.dat$`Transect number` == curr.tran,]$`Primary observer`)
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

  #remove observations without data
  if(any(obs.dat$size == 0)) {
    obs.dat <- obs.dat[-which(obs.dat$size == 0),]
  }
  if(any(is.na(obs.dat$size))) {
    obs.dat <- obs.dat[-which(is.na(obs.dat$size)),]
  }

  obs.dat <- obs.dat %>% mutate(Region.Label=rep(1, dim(obs.dat)[1])) %>%  rename(Sample.Label = `Transect #`)

  return(obs.dat)

}



#' Create data object for removal analysis without distance
#'
#' A helper function to format the data appropriately for a removal analysis the R package `FSA`.
#' @param transect.dat Data table of each transect.
#' @param obs.dat ata table observations in each transect.
#' @param pool.transects Analyze transects seperately (FALSE) or together (TRUE)
#' @return A list object with counts summarized on each transect or over all transects
#' @export
#' @examples
#' create.removal(transect.dat, obs.dat, pool.transects=TRUE)
create.removal <- function(double.dat, transect.dat, pool.transects=TRUE) {

  double.dat <- left_join(double.dat, transect.dat, by = c("Transect number", "Observer name"= "name"))

  counts <- double.dat %>% group_by(`Transect number`, observer) %>%  summarize(y = sum(`size`))
  counts.wide <- counts %>% spread(key=observer, value=y)

  counts.wide$secondary[which(is.na(counts.wide$secondary))] <- 0

  if(any(is.na(counts.wide))) {
    counts.wide[is.na(counts.wide)] <- 0
  }
  counts.wide$primary
  if(pool.transects) {
    return(c(sum(counts.wide$primary), sum(counts.wide$secondary)))
  }
  n.list      <- split(as.matrix(counts.wide[,2:3]), f = counts.wide$`Transect number`)

  return(n.list)
}
