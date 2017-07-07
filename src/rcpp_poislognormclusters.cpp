//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
#include <time.h>
//#include <string.h>
#define NRANSI
#include <R.h>
//#include <Rmath.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

int rcpp_clusterfunction(int *numc, double *lnmean, double *lnsd, double *oprob, double *clustersizes)
{
  int numclusters=numc[0],i;
  double total_cases, incident_cases, myrand;
  double ln1=lnmean[0], ln2=lnsd[0], overlap=oprob[0];
  double secondarycases,lambda;

  GetRNGstate(); //Necessary R call
  
  srand48((long)time(NULL));

  //check that R0<1!!!
  for(i=0; i<numclusters; i++) 
     { 
       clustersizes[i]=0; 
    } 
   for(i=0; i<numclusters; i++) 
   { 
      total_cases=1; 
      incident_cases=1; 

      //secondarycases=gsl_ran_lognormal(gslrand,ln1,ln2);
      //total_cases += secondarycases;
      //clustersizes[i]+=total_cases;

      do 
   	{ 
   	  incident_cases-=1; 
   	  //secondarycases=rnbinom(nb2,nb1); 
   	  //secondarycases=rlnorm(ln1,ln2); 
	  lambda=rlnorm(ln1,ln2);
	  secondarycases=rpois(lambda);
   	  //secondarycases=gsl_ran_lognormal(gslrand,ln1,ln2); 
   	  total_cases += secondarycases; 
   	  incident_cases += secondarycases; 

   	  if(total_cases > 1000){break;} 

   	}while(incident_cases>0); 

      clustersizes[i]+=total_cases; 

       myrand=drand48(); 
      if(myrand<overlap){i--;} 
      
    } 
 
  PutRNGstate(); //Necessary R call
  return 1; 
}



