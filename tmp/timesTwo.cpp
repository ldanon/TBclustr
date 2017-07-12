#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x, NumericVector lnmean, NumericVector lnsd) {

  double ln1=lnmean[0], ln2=lnsd[0];
  double lambda;
  lambda=R::rlnorm(-1,0.2); 
  Rcout << lambda <<std::endl;
  R::rpois(lambda);
  
    srand48((long)time(NULL));

  return ln2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42,1,0.3)
*/
