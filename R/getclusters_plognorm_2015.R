setwd("~/GitHub/TBclustr")

rc=read.csv("realclustersizes2015.csv")




#system.time(clusterdist(p,r,lambda,clustersizes))

clusterdist_lognorm <- function(noclusters,lnmean, lnsd, overlap, clustersout) {
		rzero = exp(lnmean + 0.5*lnsd^2)
      if (rzero>=1)
      {
          print(paste("logmean=",lnmean,"logsd=",lnsd,"R0=",rzero))
          stop("R0 is greater than 1")
      }
                  
      clustersout=rep(0,noclusters)
#           out <- .C("rcpp_clusterfunction",
#                   numc=as.integer(length(clustersout)), lnmean = as.double(lnmean), nlnsd=as.double(lnsd), oprob=as.double(overlap), clustersizes=as.double(clustersout))
#           return(out)
      rcpp_clusterfunction(length(clustersout),lnmean,lnsd,overlap,clustersout)
}

#test
clusterdist_lognorm(10,-1,0.3,0.3)

mylseq=exp(seq(log(1),log(max(rc$x)),length.out=10))

myseq=seq(1,max(rc$x),1)

h0=hist(rc$x,plot=F,breaks=myseq)
plot(h0$breaks[2:(length(h0$counts)+1)],h0$counts/sum(h0$counts),log="xy",xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,pch=19,col="purple",cex=1.5)

lines(1:200,dlnorm(1:200,-1,2)/sum(dlnorm(1:200,-1,2)),log="xy",type="l")
numberclusters=length(rc$x)
clustersizes=rep(0,numberclusters)
c1=clusterdist_lognorm(-1,0.3,0.3,clustersizes)
system.time(clusterdist(0.2,2,0.3,clustersizes))

# lines(1:200,dlnorm(1:200,10,3)/sum(dlnorm(1:200,10,3)),log="xy")
# lines(1:200,dlnorm(1:200,5,10)/sum(dlnorm(1:200,5,10)),log="xy")


measure = function(pars,data)
{
	prob=pars[1]
	size=pars[2]
	overlap=pars[3]
	rzero = prob*size/(1-prob)
    if (rzero>=1){return(0)}
   	if(prob<0 | prob>1){return(0)}
   	if(size<0){return(0)}
   	if(overlap<0 | overlap>1){return(0)}
	numberclusters=length(data)
	clustersizes=rep(0,numberclusters)
	c1=clusterdist(prob,size,overlap,clustersizes)
	options(warn=-1)
	if(sum(is.na(c1$clustersizes))>0)
		{return(0)}
	k1=ks.test(x=data,y=c1$clustersizes)
	m1=-log(as.numeric(k1$statistic))
	print(paste("pars=",pars[1],pars[2],pars[3],"measure=",round(m1,2)))
	return(m1)
	#return(as.numeric(k1$statistic))
}

modelstats = function(pars)
{
	mu=pars[1]
	sigma=pars[2]
	overlap=pars[3]
	rzero = exp(mu + 0.5*sigma^2)
    
  numberclusters=10968 #this is length(rc$x)
  maxclustersize=300 #max(rc$x=204)
   	mylseq=unique(floor(exp(seq(log(1),log(maxclustersize),length.out=50))))
   	myseq=seq(1,242,1)
	counts=rep(0,length(mylseq))
	
   	if (rzero>=1){return(counts)}
   	
	c1=clusterdist_lognorm(numberclusters,mu,sigma,overlap)
	
	if(sum(is.na(c1$clustersizes))>0){return(counts)}
	if(sum(is.infinite(c1$clustersizes))>0){return(counts)}
	
	modelclusters = c1$clustersizes

	for(i in 1:(length(mylseq)-1))
	{
		counts[i] = sum(modelclusters>=mylseq[i] & modelclusters<mylseq[i+1])
	}				

	return(counts)
}

require("EasyABC")

set.seed(5)

maxclustersize=300 #max(rc$x=204)
mylseq=unique(floor(exp(seq(log(1),log(maxclustersize),length.out=50))))


datacounts=rep(0,length(mylseq))
for(i in 1:(length(mylseq)-1))
{
		datacounts[i] = sum(rc$x>=mylseq[i] & rc$x<mylseq[i+1])
}  	
plot(datacounts,log="xy")


modelstats(c(-2,0.1,0.1))

modelstats(c(-0.1,0.1,0.5))

n=1000
p=0.1
priors=list(c("unif",-10,0),c("unif",0,10),c("unif",0,0.5))
prior_test=c("X1<=0","X2>0","X3>0","X3<1","(X2^2)<(-2*X1)")
ABC_rej<-ABC_rejection(model=modelstats, prior=priors, nb_simul=n, summary_stat_target=datacounts, tol=p,prior_test=prior_test)



ABC_plognorm15<-ABC_mcmc(method="Marjoram", model=modelstats, prior=priors, summary_stat_target=datacounts,  progress_bar=T, prior_test=prior_test,n_rec=1e5)
save(ABC_plognorm15,file="MCMCplognormb15.RData")
load(file="MCMCplognormb15.RData")
N1=length(ABC_plognorm15$param[,1])

#N1=length(ABC_rej$param[,1])
myboot = 1000
mod1=matrix(0,ncol=myboot,nrow=length(mylseq))
for(i in 1:myboot)
{
	myrand=floor(runif(1,1,N1+1))
	
	c1=clusterdist_lognorm(length(rc$x),ABC_plognorm15$param[myrand,1],ABC_plognorm15$param[myrand,2],ABC_plognorm15$param[myrand,3])
	#c1=clusterdist_lognorm(length(rc$x),ABC_rej$param[myrand,1],ABC_rej$param[myrand,2],ABC_rej$param[myrand,3])
	
	modelcounts=rep(0,length(mylseq))
	modelclusters = c1$clustersizes
	#for(j in 1:(length(myseq)-1)){modelcounts[j] = sum(modelclusters>=myseq[j] & modelclusters<myseq[j+1])}
	for(j in 1:(length(mylseq)-1)){modelcounts[j] = sum(modelclusters>=mylseq[j] & modelclusters<mylseq[j+1])}	

	mod1[,i] = modelcounts
}

low1=apply(mod1,MARGIN=1,FUN=function(x){sort(x)[floor(length(x)*0.025)]})
hi1=apply(mod1,MARGIN=1,FUN=function(x){sort(x)[floor(length(x)*0.975)]})


#without importation!
mod0=matrix(0,ncol=myboot,nrow=length(mylseq))
imported=matrix(0,ncol=myboot,nrow=1)
for(i in 1:myboot)
{
	myrand=floor(runif(1,1,N1+1))
	
	c1=clusterdist_lognorm(length(rc$x),ABC_plognorm15$param[myrand,1],ABC_plognorm15$param[myrand,2],0)
	c2=clusterdist_lognorm(length(rc$x),ABC_plognorm15$param[myrand,1],ABC_plognorm15$param[myrand,2],ABC_plognorm15$param[myrand,3])
	#c1=clusterdist_lognorm(length(rc$x),ABC_rej$param[myrand,1],ABC_rej$param[myrand,2],ABC_rej$param[myrand,3])
	
	modelcounts=rep(0,length(mylseq))
	modelclusters = c1$clustersizes
	#for(j in 1:(length(myseq)-1)){modelcounts[j] = sum(modelclusters>=myseq[j] & modelclusters<myseq[j+1])}
	for(j in 1:(length(mylseq)-1)){modelcounts[j] = sum(modelclusters>=mylseq[j] & modelclusters<mylseq[j+1])}	

	mod0[,i] = modelcounts
	
	modelclusters = c2$clustersizes
	#for(j in 1:(length(myseq)-1)){modelcounts[j] = sum(modelclusters>=myseq[j] & modelclusters<myseq[j+1])}
	for(j in 1:(length(mylseq)-1)){modelcounts[j] = sum(modelclusters>=mylseq[j] & modelclusters<mylseq[j+1])}	

	mod1[,i] = modelcounts
	
	imported[i] = (sum(mod1[,i]) - sum(mylseq*mod0[,i] - mod0[,i]))/sum(mod1[,i])
}
mean(imported)
sort(imported)[floor(length(imported)*c(0.025,0.975))]

matplot(x2,mod0,pch=19,log="xy",col=rgb(1,1,0,0.1),xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,cex=1)
matplot(x2,mod1,pch=19,log="xy",col=rgb(1,0,1,0.1),xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,cex=1,add=T)

pdf("figs/plognormmodelfit_15.pdf")
par(mar=c(5,5,1,1))
x1=myseq
x2=mylseq
plot(x2,datacounts,log="xy",xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,pch=19,col="purple",cex=1.5)
#lines(x1,rowMeans(mod1),col=rgb(0,0,0,0.5),cex=1.5,lwd=3)

lines(x2,low1,col=rgb(0,0,0,0.5),cex=1.5,lwd=3)
lines(x2,hi1,col=rgb(0,0,0,0.5),cex=1.5,lwd=3)

polygon(c(x2[low1>0],rev(x2[hi1>0])),c(low1[low1>0],rev(hi1[hi1>0])),col=rgb(0.3,0.5,0.5,0.5),border=NA)
points(x2,datacounts,pch=19,col="purple",cex=1.5)

legend('topright',c("Data","Model"),pch=c(19,22),col=c("purple","black"),cex=1.5,inset=0.01,pt.bg=rgb(0,0,0,0.5),pt.cex=c(1.5,2.5))
dev.off()

pdf("figs/plognormmodelfit2_15.pdf")
par(mar=c(5,5,1,1))
matplot(x2,mod1,pch=19,log="xy",col=rgb(1,1,0,0.1),xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,cex=1)
points(x2,datacounts,pch=22,col="red",cex=1.5,bg="red")
legend('topright',c("Data","Model"),pch=c(22,19),col=c("red",rgb(1,1,0)),cex=1.5,inset=0.01,pt.bg="red",pt.cex=c(2.5,2.5))
dev.off()

p1=ABC_plognorm15$param[,1]
p2=ABC_plognorm15$param[,2]
p3=ABC_plognorm15$param[,3]
> sort(p1)[floor(length(p1)*(0.025))]
[1] -4.681931
> sort(p1)[floor(length(p1)*(0.975))]
[1] -1.518265
> sort(p2)[floor(length(p1)*(0.025))]
[1] 1.180043
> sort(p2)[floor(length(p1)*(0.975))]
[1] 2.836614


rzero = exp(p1 + 0.5*p2^2)
sort(rzero)[floor(length(p1)*(0.5))]
sort(rzero)[floor(length(p1)*(0.025))]
sort(rzero)[floor(length(p1)*(0.975))]
mean(rzero)

sort(p3)[floor(length(p1)*(0.5))]
sort(p3)[floor(length(p1)*(0.025))]
sort(p3)[floor(length(p1)*(0.975))]


uktransmitted = rzero/(rzero+p3) ###????????
sort(uktransmitted)[floor(length(p1)*(0.5))]
sort(uktransmitted)[floor(length(p1)*(0.025))]
sort(uktransmitted)[floor(length(p1)*(0.975))]
mean(uktransmitted)

plot(p3,type="l")

plot(p1,p2)
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(ABC_plognorm15$param,diag.panel=panel.hist)
# x1=10002
# c1=clusterdist(ABC_MCMC3$param[x1,1],ABC_MCMC3$param[x1,2],ABC_MCMC3$param[x1,3],clustersizes)

c1=clusterdist_lognorm(length(rc$x),mean(ABC_rej$param[,1]),mean(ABC_rej$param[,2]),mean(ABC_rej$param[,3]))

modelcounts=rep(0,length(mylseq))
modelclusters = c1$clustersizes
for(i in 1:(length(mylseq)-1)){modelcounts[i] = sum(modelclusters>=mylseq[i] & modelclusters<mylseq[i+1])}	

plot(mylseq,datacounts,log="xy")	
points(mylseq,modelcounts,col="red",pch=19)

myboot = length(ABC_rej$param[,1])
mod1=matrix(0,ncol=myboot,nrow=length(mylseq))
for(i in 1:myboot)
{
	myrand=floor(runif(1,1,myboot+1))
	
	c1=clusterdist_lognorm(length(rc$x),ABC_rej$param[myrand,1],ABC_rej$param[myrand,2],ABC_rej$param[myrand,3])
	
	modelcounts=rep(0,length(mylseq))
	modelclusters = c1$clustersizes
	for(j in 1:(length(mylseq)-1)){modelcounts[j] = sum(modelclusters>=mylseq[j] & modelclusters<mylseq[j+1])}	

	mod1[,i] = modelcounts
}
low1=apply(mod1,MARGIN=1,FUN=function(x){sort(x)[floor(length(x)*0.025)]})
hi1=apply(mod1,MARGIN=1,FUN=function(x){sort(x)[floor(length(x)*0.975)]})


#pdf("figs/lognormmodelfit.pdf")
par(mar=c(5,5,1,1))
x1=mylseq
plot(x1,datacounts,log="xy",xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,pch=19,col="purple",cex=1.5)
#lines(x1,rowMeans(mod1),col=rgb(0,0,0,0.5),cex=1.5,lwd=3)

lines(x1,low1,col=rgb(0,0,0,0.5),cex=1.5,lwd=3)
lines(x1,hi1,col=rgb(0,0,0,0.5),cex=1.5,lwd=3)

polygon(c(x1[low1>0],rev(x1[hi1>0])),c(low1[low1>0],rev(hi1[hi1>0])),col=rgb(0.3,0.5,00.5,0.5),border=NA)
points(x1,datacounts,pch=19,col="purple",cex=1.5)

legend('topright',c("Data","Model"),pch=c(19,22),col=c("purple","black"),cex=1.5,inset=0.01,pt.bg=rgb(0,0,0,0.5),pt.cex=c(1.5,2.5))
#dev.off()

p1=ABC_rej$param[,1]
# hist(p1)

p2=ABC_rej$param[,2]
# hist(p2)
rzero = exp(p1 + 0.5*p2^2)

var2=exp(2*p1+p2^2)*(exp(p2^2)-1)
pdf("figs/plognormrzero15.pdf")
par(mar=c(5,5,1,1))
hist(rzero,col="purple",xlab="Basic reproduction number in UK",main="")
box()
dev.off()
mean(rzero)

p3=ABC_rej$param[,3]

uktransmitted = rzero/(rzero+p3) ###????????

pdf("figs/plognorm_imported15.pdf")
par(mar=c(5,5,1,1))
hist(1-uktransmitted,col=rgb(0,0.5,0.7),xlab="Fraction of cases imported",main="")
box()
dev.off()
mean(1-uktransmitted)

# c1=clusterdist(mean(ABC_MCMC5$param[,1]),mean(ABC_MCMC5$param[,2]),mean(ABC_MCMC5$param[,3]),clustersizes)


# h0=hist(rc$x,plot=F,breaks=seq(0,max(c1$clustersizes,rc$x),1))
# h1=hist(c1$clustersizes,plot=F,breaks=seq(0,max(c1$clustersizes,rc$x),1))
# par(mar=c(5,5,1,1))
# plot(h0$breaks[2:(length(h0$counts)+1)],h0$counts,log="xy",xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,pch=19,col="purple",cex=1.5)
# points(h1$breaks[2:(length(h0$counts)+1)],h1$counts,log="xy",bg=rgb(0,0,0,0.5),pch=22,col=rgb(0,0,0,0.5),cex=1.5)
# legend('topright',c("Data"),pch=19,col="purple",cex=1.5,inset=0.01)


# save(ABC_MCMC3,file="MCMC3output.RData")
# p1=ABC_MCMC3$param[,1]
# hist(p1)

# p2=ABC_MCMC3$param[,2]
# hist(p2)

# p3=ABC_MCMC3$param[,3]
# hist(p3)
# summary(p1*p2/(1-p1))
# hist(p1*p2/(1-p1))
# summary(p3)
# #uk transmitted cases:?
# 1-0.3781152/(0.3781152+0.1644672)

# priors3=list(c("unif",0.0,0.5),c("unif",0,1),c("unif",0,0.5))
# prior_test3=c("X1>0","X2>0","X3>0","X1<0.5","X3<0.75","X2<0.5")
# ABC_MCMC3<-ABC_mcmc(method="Marjoram_original", model=modelstats, prior=priors3, summary_stat_target=data_summary_stats,  progress_bar=T, prior_test=prior_test3,n_rec=1e2)

# #ABC_MCMC<-ABC_mcmc(method="Marjoram", model=modelstats, prior=priors, summary_stat_target=data_summary_stats,  progress_bar=T, prior_test=prior_test,n_rec=100)

