library(wesanderson)
library(tidyverse)
load('../Writing/DensEst.Rdata')


#dataf <- data.frame(Lake=rep(c("Lake Florida", "Lake Burgan", "Little Birch Lake"),3), Design=c(rep("Quadrat", 3), rep("Double",3), rep("Distance", 3)), Estimate=c(LF.quadrat.est$Dhat, LB.quadrat.est$Dhat, LBL.quadrat.est.subset$Dhat, LF.double.est$Dhat, LB.double.est$Dhat, LBL.double.est$Dhat, LF.distance.est$Dhat, LB.distance.est$Dhat, LBL.distance.est$Dhat), SE=c(LF.quadrat.est$Dhat.se, LB.quadrat.est$Dhat.se, LBL.quadrat.est.subset$Dhat.se, LF.double.est$Dhat.se, LB.double.est$Dhat.se, LBL.double.est$Dhat.se, LF.distance.est$Dhat.se, LB.distance.est$Dhat.se, LBL.distance.est$Dhat.se))

dataf <- data.frame(Lake=rep(c("Lake Florida \n (low density)", "Lake Burgan \n (medium density)", "Little Birch Lake \n (high density)"),2), Design=c(rep("Quadrat", 3), rep("Distance", 3)), Estimate=c(LF.quadrat.est$Dhat, LB.quadrat.est$Dhat, LBL.quadrat.est.subset$Dhat, LF.distance.est$Dhat, LB.distance.est$Dhat, LBL.distance.est$Dhat), SE=c(LF.quadrat.est$Dhat.se, LB.quadrat.est$Dhat.se, LBL.quadrat.est.subset$Dhat.se, LF.distance.est$Dhat.se, LB.distance.est$Dhat.se, LBL.distance.est$Dhat.se))

dataf$Lake <- relevel(dataf$Lake, "Lake Florida \n (low density)")
surveyTypeCols <- wes_palette("Darjeeling1")

ggplot(dataf, aes(x=Design, y=Estimate, colour=Design)) + geom_point(size=4) + facet_wrap(~Lake, ncol=3) + geom_errorbar(aes(ymin=Estimate-2*SE, ymax=Estimate+2*SE, width=0)) + labs(x="", y="Density estimate") + theme_classic() + scale_fill_manual(values=surveyTypeCols) + scale_colour_manual(values=surveyTypeCols) + guides(fill=FALSE)

ggplot(dataf, aes(x=Design, y=Estimate, colour=Design)) + geom_point(size=4) + facet_wrap(~Lake, ncol=3, scale="free") + geom_errorbar(aes(ymin=Estimate-2*SE, ymax=Estimate+2*SE, width=0)) + labs(x="", y="Density estimate") + theme_classic() + scale_fill_manual(values=surveyTypeCols) + scale_colour_manual(values=surveyTypeCols) + guides(fill=FALSE) #+ scale_y_log10()




