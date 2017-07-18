setwd("~/GitHub/TBclustr")

#UK data up to 2015
rc=read.csv("../../realclustersizes2015.csv")


modelstats = function(pars)
{
	mu=pars[1]
	sigma=pars[2]
	overlap=pars[3]
	numberclusters=pars[4]
	maxclustersize=pars[5] #max(rc$x=204)
	
	rzero = exp(mu + 0.5*sigma^2)
	if (rzero>=1){return(counts)}
	
	clustersout=rep(0,numberclusters)
	rcpp_clusterfunction(numberclusters,mu,sigma,overlap,clustersout)
	
	#no longer using c1$clustersizes
	if(sum(is.na(clustersout))>0){return(counts)}
	if(sum(is.infinite(clustersout))>0){return(counts)}
	
	#log bin cluster sizes 
	mylseq=unique(floor(exp(seq(log(1),log(maxclustersize),length.out=50))))
	counts=rep(0,length(mylseq))
	
	modelclusters = clustersout

	for(i in 1:(length(mylseq)-1))
	{
		counts[i] = sum(modelclusters>=mylseq[i] & modelclusters<mylseq[i+1])
	}				

	return(counts)
}

#test:
modelstats(c(-2,0.1,0.1,length(rc$x),300))

#fitting using ABC
require("EasyABC")

set.seed(5)

maxclustersize=floor(max(rc$x)*1.5)#this is somewhat arbitrary

#log bin data 
mylseq=unique(floor(exp(seq(log(1),log(maxclustersize),length.out=50))))
datacounts=rep(0,length(mylseq))
for(i in 1:(length(mylseq)-1))
{
		datacounts[i] = sum(rc$x>=mylseq[i] & rc$x<mylseq[i+1])
}  	
plot(datacounts,log="xy",xlab="Cluster size",ylab="Frequency")

#define priors (last two variables are fixed)
priors=list(c("unif",-10,0),c("unif",0,10),c("unif",0,0.5),c("unif",length(c1),length(c1)),c("unif",300,300))

#set criteria for parameters (rzero<1)
prior_test=c("X1<=0","X2>0","X3>0","X3<1","(X2^2)<(-2*X1)")

#test ABC code with rejection algorithm:
p=0.1
n=1000
ABC_rej<-ABC_rejection(model=modelstats, prior=priors, nb_simul=n, summary_stat_target=datacounts, tol=p,prior_test=prior_test)

#now run for real with Marjoram MCMC algorithm (may take several hours)
#ABC_plognorm15<-ABC_mcmc(method="Marjoram", model=modelstats, prior=priors, summary_stat_target=datacounts,  progress_bar=T, prior_test=prior_test,n_rec=1e5)

ABC_plognorm15 = ABC_rej
#save(ABC_plognorm15,file="MCMCplognormb15.RData")
#load(file="MCMCplognormb15.RData")


N1=length(ABC_plognorm15$param[,1])
myboot = 1000
mod1=matrix(0,ncol=myboot,nrow=length(mylseq))
clustersout = rep(0,length(rc$x))
for(i in 1:myboot)
{
	myrand=floor(runif(1,1,N1+1))
	rcpp_clusterfunction(length(rc$x),ABC_plognorm15$param[myrand,1],ABC_plognorm15$param[myrand,2],ABC_plognorm15$param[myrand,3],clustersout)
	
	#log bin model output as data
	modelcounts=rep(0,length(mylseq))
	for(j in 1:(length(mylseq)-1)){modelcounts[j] = sum(clustersout>=mylseq[j] & clustersout<mylseq[j+1])}	

	mod1[,i] = modelcounts
}

#draw results figure
#pdf("figs/plognormmodelfit2_15.pdf")
par(mar=c(5,5,1,1))
matplot(mylseq,mod1,pch=19,log="xy",col=rgb(1,1,0,0.1),xlab="Cluster size",ylab="Frequency",cex.lab=1.7,cex.axis=1.5,cex=1)
points(mylseq,datacounts,pch=22,col="red",cex=1.5,bg="red")
legend('topright',c("Data","Model"),pch=c(22,19),col=c("red",rgb(1,1,0)),cex=1.5,inset=0.01,pt.bg="red",pt.cex=c(2.5,2.5))
#dev.off()



#now examine posteriors
p1=ABC_plognorm15$param[,1]
p2=ABC_plognorm15$param[,2]
p3=ABC_plognorm15$param[,3]



rzero = exp(p1 + 0.5*p2^2)
rmed=sort(rzero)[floor(length(p1)*(0.5))]
rlo=sort(rzero)[floor(length(p1)*(0.025))]
rhi=sort(rzero)[floor(length(p1)*(0.975))]
mean(rzero)
print(paste("Median R0 =",round(rmed,3),"95%CI = (",round(rlo,3),",",round(rhi,3),")"))
#pdf("figs/plognormrzero15.pdf")
par(mar=c(5,5,1,1))
hist(rzero,col="purple",xlab="Basic reproduction number",main="")
box()
#dev.off()

uktransmitted = rzero/(1+rzero+p3)
ukmed=sort(uktransmitted)[floor(length(p1)*(0.5))]
uklo=sort(uktransmitted)[floor(length(p1)*(0.025))]
ukhi=sort(uktransmitted)[floor(length(p1)*(0.975))]
mean(uktransmitted)
print(paste("Median % recent transmission =",round(100*ukmed,1),"% 95%CI = (",round(100*uklo,1),",",round(100*ukhi,1),")"))


#pdf("figs/plognorm_imported15.pdf")
par(mar=c(5,5,1,1))
hist(1-uktransmitted,col=rgb(0,0.5,0.7),xlab="Fraction of cases imported",main="")
box()
#dev.off()

