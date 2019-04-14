// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// S matrix of Savitzky Golay
// [[Rcpp::export]]
arma::mat sgmat_S(int halfwin = 5, int d = 2) {
    // int res = halfwin + d;
    int frame = 2*halfwin + 1;
    // IntegerVector vec_f = seq(-halfwin, halfwin);
    // IntegerVector vec_d = seq(0, d);

    arma::mat mat = arma::zeros(2*halfwin+1, d+1);

    for (int i = 0; i < frame; i++) {
        for (int j =0; j <= d; j++) {
            mat(i, j) = pow(i - halfwin, j);
        }
    }
    // Rcout << mat << std::endl;
    return mat;
}

// B matrix of Savitzky Golay
// [[Rcpp::export]]
arma::mat sgmat_B(const arma::mat S) {
    arma::mat Q, R, B, T;
    // Rcout << repmat(sqrt(w), 1, S.n_cols) % S << std::endl;
    arma::qr_econ(Q, R, S);

    T = solve(R.t(), S.t());
    // Rcout << 'R' <<  R << ;
    B = T.t() * T ;
    return B;
}

// B matrix of weighted Savitzky Golay
// [[Rcpp::export]]
arma::mat sgmat_wB(const arma::mat S, const arma::colvec w) {
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

//' Weighted Savitzky-Golay
//'
//' @param y colvec
//' @param w colvec of weight
//' @param halfwin halfwin of Savitzky-Golay
//' @param d polynomial of degree. When d = 0, it becomes moving average.
//'
//' @examples
//' y <- 1:15
//' w <- seq_along(y)/length(y)
//'
//' frame = 5
//' d = 2
//' s1 <- smooth_wSG(y, w, frame, d)
//' s2 <- smooth_SG(y, frame, d)
//' @export
// [[Rcpp::export]]
NumericVector smooth_wSG(const arma::colvec y, const arma::colvec w, const int halfwin, const int d=2){
    int n       = y.n_rows;
    int frame   = halfwin*2 + 1;

    arma::mat S = sgmat_S(halfwin, d);
    arma::mat B = sgmat_wB(S, w.subvec(0, frame-1));
    // Rcpp::Rcout << B << y << std::endl;
    arma::colvec y_head = B.rows(0, halfwin) * y.subvec(0, frame-1);

    arma::colvec y_mid = arma::Col<double>(n-frame-1) ;
    for (int i=1; i<=n-frame-1; i++) {
        B = sgmat_wB(S, w.subvec(i, i+frame-1));
        y_mid(i-1) = as_scalar( B.row(halfwin) * y.subvec(i, i+frame-1) );
    }

    B = sgmat_wB(S, w.subvec(n-frame, n-1));
    arma::colvec y_tail = B.rows(halfwin, frame-1) * y.subvec(n-frame, n-1);

    arma::colvec yfit  = join_vert(y_head, y_mid);
    yfit = join_vert(yfit, y_tail);

    // Rcpp::Rcout << yfit << std::endl;
    // NumericVector y2 = as<NumericVector>(wrap(yfit));
    // NumericVector(a.begin(),a.end())
    return Rcpp::NumericVector(yfit.begin(), yfit.end());
}

//' @rdname smooth_wSG
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector smooth_SG(const arma::colvec y, const int halfwin, const int d=2){
    int n       = y.n_rows;
    int frame   = halfwin*2 + 1;

    arma::mat S = sgmat_S(halfwin, d);
    arma::mat B = sgmat_B(S);
    // Rcpp::Rcout << B << y << std::endl;
    arma::colvec y_head = B.rows(0, halfwin) * y.subvec(0, frame-1);

    arma::colvec y_mid = arma::Col<double>(n-frame-1) ;
    for (int i=1; i<=n-frame-1; i++) {
        y_mid(i-1) = as_scalar( B.row(halfwin) * y.subvec(i, i+frame-1) );
    }

    arma::colvec y_tail = B.rows(halfwin, frame-1) * y.subvec(n-frame, n-1);

    arma::colvec yfit  = join_vert(y_head, y_mid);
    yfit = join_vert(yfit, y_tail);

    return Rcpp::NumericVector(yfit.begin(), yfit.end());
}

/*** R
#' @param frame Odd integer

S <- sgmat_S(7, 2)
y <- 1:15
w <- seq_along(y)/length(y)

frame = 5
d = 2
smooth_wSG(y, w, frame, d)
smooth_SG(y, frame, d)
*/
