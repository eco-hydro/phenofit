#include <Rcpp.h>
using namespace Rcpp;


//' Double logistics in Rcpp
//' 
//' @inheritParams doubleLog.Beck
//' 
//' @seealso [doubleLog.Beck()]
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
NumericVector logistic( NumericVector par, NumericVector t){
    double mn  = par[0];
    double mx  = par[1];
    double sos = par[2];
    double rsp = par[3];

    NumericVector pred = mn + (mx - mn)/(1 + exp(-rsp*(t - sos)));
    // # pred <- c/(1 + exp(a + b * t)) + d
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_Zhang( NumericVector par, NumericVector t) {
    double t0  = par[0];
    double mn  = par[1];
    double mx  = par[2];
    double sos = par[3];
    double rsp = par[4];
    double eos = par[5];
    double rau = par[6];

    int n = t.size();
    NumericVector pred(n, 99.0);
    if (t0 - sos <= 1 || t0 - eos >= -1) return(pred);

    pred = ifelse(t <= t0,
        -rsp*(t - sos),
         rau*(t - eos) );
    pred = mn + (mx - mn) / (1 + exp(pred));
    // NumericVector pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_AG( NumericVector par, NumericVector t) {
    double t0  = par[0];
    double mn  = par[1];
    double mx  = par[2];
    double rsp = par[3];
    double a3  = par[4];
    double rau = par[5];
    double a5  = par[6];

    NumericVector pred = ifelse(t <= t0,
        pow( (t0 - t)*rsp, a3),
        pow( (t - t0)*rau, a5));
    pred = mn + (mx - mn) / exp(-pred);
    // NumericVector pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_Beck( NumericVector par, NumericVector t) {
    double mn  = par[0];
    double mx  = par[1];
    double sos = par[2];
    double rsp = par[3];
    double eos = par[4];
    double rau = par[5];

    int n = t.size();
    
    if (eos < sos) return( NumericVector(n, 99.0) );

    NumericVector pred = mn + (mx - mn) * 
        ( 1/(1 + exp(-rsp*(t - sos))) + 
          1/(1 + exp( rau*(t - eos))) - 1);
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_Elmore( NumericVector par, NumericVector t) {
    double mn   = par[0];
    double mx   = par[1];
    double sos  = par[2]; // SOS
    double rsp  = par[3]; // 1/rsp
    double eos  = par[4]; // EOS
    double rau  = par[5]; // 1/rau
    double m7   = par[6];

    NumericVector pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) );
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_Gu( NumericVector par, NumericVector t) {
    double y0  = par[0];
    double a1  = par[1];
    double a2  = par[2];
    double sos = par[3];
    double rsp = par[4];
    double eos = par[5];
    double rau = par[6];
    double c1  = par[7];
    double c2  = par[8];

    NumericVector pred = y0 + 
        (a1/ pow(1 + exp(-rsp*(t - sos)), c1)) - 
        (a2/ pow(1 + exp(-rau*(t - eos)), c2));
    return(pred);
}

//' @rdname logistic
//' @export
// [[Rcpp::export]]
NumericVector doubleLog_Klos( NumericVector par, NumericVector t) {
    double a1 = par[0];
    double a2 = par[1];
    double b1 = par[2];
    double b2 = par[3];
    double c  = par[4];
    double B1 = par[5];
    double B2 = par[6];
    double m1 = par[7];
    double m2 = par[8];
    double q1 = par[9];
    double q2 = par[10];
    double v1 = par[11];
    double v2 = par[12];

    NumericVector pred = (a1*t + b1) + (a2*t*t + b2*t + c) * 
        (1/ pow(1 + q1 * exp(-B1 * (t - m1)), v1) - 
         1/ pow(1 + q2 * exp(-B2 * (t - m2)), v2));
    return(pred);
}

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

# zhang
par = c(
    t0  = 100,
    mn  = 0.1,
    mx  = 0.7,
    sos = 50,
    rsp = 0.1,
    eos = 250,
    rau = 0.1)
doubleLog.Zhang(par, t)
doubleLog_Zhang(par, t)
#
# rbenchmark::benchmark(
#     doubleLog.Zhang(par, t),
#     doubleLog_Zhang(par, t),
#     doubleLog_Zhang2(par, t),
#     replications = 1000000
# )
*/
