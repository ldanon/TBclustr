setwd("~/GitHub/TBclustr")
source('R/all_functions.R')
#UK data up to 2015
rc=read.csv("../realclustersizes2015.csv")

plotclusterdist(rc) # plot the histogram of clustersizes for your data

ABC_output<-fitLogNormalmodel(rc)

mod1<-bootstrapmodel(ABC_output,rc)

examineposteriors(ABC_output)

