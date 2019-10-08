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
bool all_finite(NumericVector x) {
    bool ans = is_true(all(is_finite(x)));
    return ans;
}

// [[Rcpp::export]]
double f_goal_cpp(NumericVector par, Function fun, 
    NumericVector y, NumericVector t, NumericVector ypred,
    Nullable<NumericVector> w   = R_NilValue,
    Nullable<NumericVector> ylu = R_NilValue)
{
    // Function FUN(fun);
    if ( !all_finite(par) ) return(9999.0);

    // NumericVector ypred = 
    fun(par, t, ypred);
    if ( !all_finite(ypred) )  return(9999.0);

    double SSE;
    int len = y.size();
    if (w.isNotNull()) {
        NumericVector w_(w);

        if ( ylu.isNotNull() ) {
            NumericVector ylu_(ylu);
            double lower = ylu_[0];
            double upper = ylu_[1];

            for (int i = 0; i < len; i++) {
                if (ypred[i] < lower || ypred[i] > upper) {
                    w_[i] = 0.0;
                }
            }
            // w[ypred < ylu[0] | ypred > ylu[1]] = 0.0;
        }
        SSE = sum( pow(y - ypred, 2.0) * w_ );
    } else {
        SSE = sum( pow(y - ypred, 2.0) );
    }

    double RMSE = sqrt(SSE/len);
    // Rcout << as<bool>(all(is_finite(x))) << std::endl;
    return RMSE;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
funname = "doubleLog.Beck"
FUN = doubleLog.Beck
par = c(
    mn  = 0.1,
    mx  = 0.7,
    sos = 50,
    rsp = 0.1,
    eos = 250,
    rau = 0.1)
t    <- seq(1, 365, 8)
tout <- seq(1, 365, 1)
y <- FUN(par, t)

# initial parameter
par0 <- c(
    mn  = 0.15,
    mx  = 0.65,
    sos = 100,
    rsp = 0.12,
    eos = 200,
    rau = 0.12)

objective <- f_goal # goal function
prior <- as.matrix(par0) %>% t() %>% rbind(., .)
f_goal(par0, y, t, FUN)
f_goal(par, y, t, FUN)
*/
