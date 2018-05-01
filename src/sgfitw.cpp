// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// using namespace Rcpp;

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
arma::mat sgolayB(const arma::mat S, const arma::colvec w) {
    arma::mat Q, R, B, T;
    // Rcout << repmat(sqrt(w), 1, S.n_cols) % S << std::endl;
    arma::qr_econ(Q, R, repmat(sqrt(w), 1, S.n_cols)%S);
    // if (nargout==2) {
    //     mat G = S/(arma::trans(R)*R) ;
    //     B = G*arma::trans(S) ;
    // } else {
        // T = S/R;
        T = solve(R.t(), S.t());
        // Rcout << 'R' <<  R << ;
        B = T.t() * T ;
    // }
    // Rcout << B << std::endl;
    B = repmat(w, 1, B.n_cols).t() % B ;
    return B;
}

// [[Rcpp::export]]
arma::colvec sgfitw_rcpp(const arma::colvec y, const arma::colvec w, const arma::mat S){
    arma::mat B;
    int n       = y.n_rows;
    int frame   = S.n_rows;
    int halfwin = (frame-1)/2;

    B = sgolayB(S, w.subvec(0, frame-1)) ;
    arma::colvec y_head = B.rows(0, halfwin) * y.subvec(0, frame-1);
    arma::colvec y_mid = arma::Col<double>(n-frame) ;
    for (int i=0; i<=n-frame-1; i++) {
        B = sgolayB(S, w.subvec(i, i+frame-1)) ;
        y_mid(i) = as_scalar( B.row(halfwin) * y.subvec(i, i+frame-1) );
    }

    B = sgolayB(S, w.subvec(n-frame, n-1)) ;
    arma::colvec y_tail = B.rows(halfwin+1, frame-1) * y.subvec(n-frame, n-1);
    arma::colvec yfit  = join_vert(y_head, y_mid);
    yfit = join_vert(yfit, y_tail);

    // NumericVector y2 = as<NumericVector>(wrap(yfit));
    // NumericVector(a.begin(),a.end())
    return yfit;
}

/*** R
#' @param frame Odd integer

S <- sgolayS(7, 2)
y <- 1:15
w <- 1:15
# sgolay(S, w)
sgfitw_rcpp(y, w, S)
sgfitw_rcpp(INPUT$y, INPUT$w, S)
*/
