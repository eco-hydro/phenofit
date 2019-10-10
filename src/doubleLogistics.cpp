#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
void clogistic( NumericVector par, NumericVector t, NumericVector pred){
    double mn  = par[0];
    double mx  = par[1];
    double sos = par[2];
    double rsp = par[3];

    pred = mn + (mx - mn)/(1 + exp(-rsp*(t - sos)));
    // # pred <- c/(1 + exp(a + b * t)) + d
    // return(pred);
}

// [[Rcpp::export]]
void cdoubleLog_Zhang( NumericVector par, NumericVector t, NumericVector pred) {
    double t0  = par[0];
    double mn  = par[1];
    double mx  = par[2];
    double sos = par[3];
    double rsp = par[4];
    double eos = par[5];
    double rau = par[6];

    if (t0 - sos <= 1 || t0 - eos >= -1) pred = pred*0 + 99.0;

    pred = ifelse(t <= t0,
        -rsp*(t - sos),
         rau*(t - eos) );
    pred = mn + (mx - mn) / (1 + exp(pred));
    // NumericVector pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    // return(pred);
}

// [[Rcpp::export]]
void cdoubleLog_AG( NumericVector par, NumericVector t, NumericVector pred) {
    double t0  = par[0];
    double mn  = par[1];
    double mx  = par[2];
    double rsp = par[3];
    double a3  = par[4];
    double rau = par[5];
    double a5  = par[6];

    pred = ifelse(t <= t0,
        pow( (t0 - t)*rsp, a3),
        pow( (t - t0)*rau, a5));
    pred = mn + (mx - mn) * exp(-pred);
    // NumericVector pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) )
    // return(pred);
}


// [[Rcpp::export]]
void cdoubleLog_Beck( NumericVector par, NumericVector t, NumericVector pred) {
    double mn  = par[0];
    double mx  = par[1];
    double sos = par[2];
    double rsp = par[3];
    double eos = par[4];
    double rau = par[5];

    if (eos < sos) pred = pred*0 + 99.0;

    pred = mn + (mx - mn) *
        ( 1/(1 + exp(-rsp*(t - sos))) +
          1/(1 + exp( rau*(t - eos))) - 1);
    // return(pred);
}

// [[Rcpp::export]]
void cdoubleLog_Elmore( NumericVector par, NumericVector t, NumericVector pred) {
    double mn   = par[0];
    double mx   = par[1];
    double sos  = par[2]; // SOS
    double rsp  = par[3]; // 1/rsp
    double eos  = par[4]; // EOS
    double rau  = par[5]; // 1/rau
    double m7   = par[6];

    // for(int i = 0; i < pred.size(); i++) {
    //     pred[i] = mn + (mx - m7*t[i])*( 1/(1 + exp(-rsp*(t[i]-sos))) - 1/(1 + exp(-rau*(t[i]-eos))) );
    // }
    pred = mn + (mx - m7*t)*( 1/(1 + exp(-rsp*(t-sos))) - 1/(1 + exp(-rau*(t-eos))) );
    // return(pred);
}

// [[Rcpp::export]]
void cdoubleLog_Gu( NumericVector par, NumericVector t, NumericVector pred) {
    double y0  = par[0];
    double a1  = par[1];
    double a2  = par[2];
    double sos = par[3];
    double rsp = par[4];
    double eos = par[5];
    double rau = par[6];
    double c1  = par[7];
    double c2  = par[8];

    pred = y0 +
        (a1/ pow(1 + exp(-rsp*(t - sos)), c1)) -
        (a2/ pow(1 + exp(-rau*(t - eos)), c2));
    // return(pred);
}

// [[Rcpp::export]]
void cdoubleLog_Klos( NumericVector par, NumericVector t, NumericVector pred) {
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

    pred = (a1*t + b1) + (a2*t*t + b2*t + c) *
        (1/ pow(1 + q1 * exp(-B1 * (t - m1)), v1) -
         1/ pow(1 + q2 * exp(-B2 * (t - m2)), v2));
    // return(pred);
}

/*** R
funname = "cdoubleLog.Beck"
FUN = cdoubleLog.Beck
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
cdoubleLog.Zhang(par, t)
ypred <- y*0
cdoubleLog_Zhang(par, t, ypred)
#
# rbenchmark::benchmark(
#     cdoubleLog.Zhang(par, t),
#     cdoubleLog_Zhang(par, t),
#     cdoubleLog_Zhang2(par, t),
#     replications = 1000000
# )
*/
