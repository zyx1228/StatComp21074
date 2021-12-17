#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param b the number of burn-in length
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' XC <- myGibbsC(5000,1000)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix myGibbsC(int N, int b) {
  NumericMatrix XC(N, 2);
  NumericMatrix XCN(N-b, 2);
  XC(0,0) = 0.5, XC(0,1) = 0.5;
  for(int i = 1; i < N; i++) {
    XC(i,0) = rbinom(1, 16, XC(i-1,1))[0];
    XC(i,1) = rbeta(1, XC(i,0) + 2, 20 - XC(i,0))[0];
  }
  for(int j = 0; j < N-b;j++){
    XCN(j,0) = XC(j+b,0);
    XCN(j,1) = XC(j+b,1);
  }
  return(XCN);
}