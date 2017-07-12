//#include <time.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector rcpp_clusterfunction(NumericVector numc, 
                                   NumericVector lnmean, 
                                   NumericVector lnsd, 
                                   NumericVector oprob, 
                                   NumericVector clustersizes)
{
  int numclusters=numc[0],i;
  double total_cases, incident_cases, myrand;
  double ln1=lnmean[0], ln2=lnsd[0], overlap=oprob[0];
  double secondarycases,lambda;
  
  RNGScope rngScope;
  
  //GetRNGstate(); //Necessary Random number generator no longer neccesary in Rcpp and generates an error
  
  srand48((long)time(NULL));
  
  //check that R0<1!!!
    for(i=0; i<numclusters; i++){ 
       clustersizes[i]=0; // set cluster sizes to 0 
    } 
    for(i=0; i<numclusters; i++){ // go through each cluster.
   
      total_cases=1; 
      incident_cases=1; 
      do{ 
   	    incident_cases-=1; // carry on doing this until you run out of cases
   	    //secondarycases=rnbinom(nb2,nb1); 
   	    //secondarycases=rlnorm(ln1,ln2); 

   	    lambda=R::rlnorm(ln1,ln2); // pick lambda from a lognormal distribution
	      secondarycases=R::rpois(lambda); // pick the number of secondary cases from a poisson distribution
   	    //secondarycases=gsl_ran_lognormal(gslrand,ln1,ln2); 
   	    total_cases += secondarycases;  // increment the total cases
   	    incident_cases += secondarycases;  // increment the incident cases. 

   	    if(total_cases > 1000){break;}  // carry on until you run out of cases
   	  }while(incident_cases>0);         // or the total number of cases is 1000 (presumably because there are no more than 1000 cases??)
      clustersizes[i]+=total_cases; // save the size of the cluster so you can return it.

//      myrand=R::runif(0,1); 
      myrand=drand48(); 
      //      Rcout << "ok "<<myrand<<" "<<overlap<<std::endl; 
      if(myrand<overlap){i--;} 
      
    } 
 
  //PutRNGstate(); //obsolete RNG stat function call. 
  return 1; 
}

/*** R
clustersout=rep(0,10)
system.time(replicate(1000000,rcpp_clusterfunction(10,-1,0.3,0.3,clustersout)))
clustersout
*/


