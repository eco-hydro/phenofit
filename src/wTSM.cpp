#include <Rcpp.h>
using namespace Rcpp;

//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/

// Weight updating method in TIMESAT.
// [[Rcpp::export]]
NumericVector rcpp_wTSM(NumericVector y, NumericVector yfit, NumericVector w,
                   int iter, int nptperyear, double wfact) {
  int n  = y.size();
  int m  = sum(w > 0.5); //valid length

  NumericVector w_ceil = ceil(w);
  NumericVector wnew(clone(w));

  double yfitmean = sum(yfit * w_ceil / m);
  double yfitstd  = sqrt((double) sum(pow( (yfit - yfitmean) * w_ceil, 2)/(m - 1)));

  int deltaT = floor(nptperyear/7.0); // fix Solaris error; 20190524
  // Rprintf("deltaT:%d\n", deltaT);

  for (int i = 0; i < n; i++){
      int m1 = std::max(0    , i - deltaT);
      int m2 = std::min(n - 1, i + deltaT);

      IntegerVector idx = Range(m1, m2);
      // Rcout << "idx is: " << idx << std::endl;

      NumericVector yi = yfit[idx];
      double yi_min = min(yi); // local minimum
      double yi_max = max(yi); // local maximum

      // Adjust the weights dependent on if the values are above or below the
      // fitted values
      if (y[i] < yfit[i] - 1e-8){
          // if (yi_min > yfitmean){
          if (yi_min > yfitmean || iter < 2){
            /**
             * If there is a low variation in an interval, i.e. if the interval
             * is at a peak or at a minima compute the normalized distance
             * between the data point and the fitted point.
             */
              double ydiff;
              if (yi_max - yi_min < 0.8 * yfitstd) {
                  ydiff = 2 *(yfit[i] - y[i])/yfitstd;
              } else{
                  ydiff = 0;
              }
              // Use the computed distance to modify the weight. Large distance
              // will give a small weight
              wnew[i] = wfact * w[i] * exp(-ydiff*ydiff);
          }
      }
  }
  return wrap(wnew);
}

/*** R
set.seed(100)
y <- yfit <- rnorm(100)
n <- length(y)
w <- rep(1, n)
nptperyear <- 25

wnew <- rcpp_wTSM(y, yfit, w, 2, nptperyear, 2)
*/
