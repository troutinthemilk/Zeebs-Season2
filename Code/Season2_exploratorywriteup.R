#+ echo=FALSE, warning=FALSE
library(knitr)
opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE, cache=TRUE, include=FALSE, results='hide', error=TRUE)

#+ echo=FALSE, warning=FALSE
library(Distance)
library(mrds)
library(dplyr)
library(dsm)
require(pscl)
library(lubridate)
library(googlesheets)
library(plotrix)
library(ggplot2)
library(gridExtra)
library(readxl)


source('ZebraFuncs.R')
#httr::config(http_version = 1)

curr.lake   <- c("Lake Burgan", "Little Birch Lake", "Lake Florida")
curr.lake2  <- c("Burgan", "Little Birch Lake", "Florida") #"Burgan"

index             <- 1
if(!exists("LB.quadrat.est")) {
LB.distance.est   <- distance.dens.est(curr.lake[index], curr.lake2[index])
LB.double.est     <- double.dens.est(curr.lake[index], curr.lake2[index])
LB.quadrat.est    <- quadrat.dens.est(curr.lake[index], curr.lake2[index])
}

index            <- 2
if(!exists("LS.quadrat.est")) {
LS.distance.est  <- distance.dens.est(curr.lake[index], curr.lake2[index])
LS.double.est    <- double.dens.est(curr.lake[index], curr.lake2[index])
LS.quadrat.est   <- quadrat.dens.est(curr.lake[index], curr.lake2[index])
}

index            <- 3
if(!exists("LF.double.est")) {
LF.distance.est  <- distance.dens.est(curr.lake[index], curr.lake2[index])
LF.double.est    <- double.dens.est(curr.lake[index], curr.lake2[index])
#LF.quadrat.est   <- quadrat.dens.est(curr.lake[index], curr.lake2[index])
}

time.df <- rbind(LB.distance.est$df, LB.double.est$df, LB.quadrat.est$df, LS.distance.est$df, LS.double.est$df, LS.quadrat.est$df, LF.distance.est$df, LF.double.est$df)#, LF.quadrat.est$df)


#' ## Time budgets
#' We are interested in how long it takes to do each task associated with the zebra mussel surveys. Below is the overall amount of time spent doing a given type of survey versus the amount of area that was actually covered.

boxplot(Area/Time ~ Type, data=time.df, xlab="Survey type", ylab="Amount of area surveyed per second")

#' We see that the distance sampling is more efficient at covering that the other types of survey. Given the high detection probalities in the double observer survey

p.enc <- ggplot(time.df, aes(x=Type, y=t.enc, fill=Type))  + coord_trans(y="log10") +  geom_boxplot() + labs(title="Time spent on encounters", x="Survey type", y = "Time") + theme_classic() + guides(fill=FALSE)

p.hab <- ggplot(time.df, aes(x=Type, y=t.hab, fill=Type)) +  geom_boxplot() + labs(title="Time spent on habitat", x="Survey type", y = "Time") + theme_classic() + guides(fill=FALSE)

p.set <- ggplot(time.df, aes(x=Type, y=t.set, fill=Type)) +  geom_boxplot() + labs(title="Time spent on setup", x="Survey type", y = "Time") + theme_classic() + guides(fill=FALSE)

#+ fig.width=5, fig.height=5             
grid.arrange(p.set, p.hab, p.enc, nrow = 1)                
time.dat <- subset(time.df, Type != "Quadrat")


