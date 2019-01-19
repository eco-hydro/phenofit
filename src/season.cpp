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

/**
 * Fix troughs in \code{season_mov} whose date is later than previous year
 * ending trough, and merge two too close troughs (less than 35 days).
 */
// [[Rcpp::export]]
void fix_dt(DataFrame d) {

    DateVector date_beg = d["beg"];
    DateVector date_end = d["end"];
    NumericVector val_beg = d["y_beg"];
    NumericVector val_end = d["y_end"];

    int n = d.nrow();
    for (int i = 0; i < n - 1; i++){
        // Rcout << i << std::endl;
        int delta_days = date_end[i] - date_beg[i+1];
        bool con = (val_end[i] <= val_beg[i+1]);     // previous end smaller ?

        if (delta_days > 0 || abs(delta_days) <= 50){
            Date  newdate = con? date_end[i] : date_beg[i+1];
            double newval = con? val_end[i] : val_beg[i+1];

            date_end[i]   = newdate; // + (-1);
            date_beg[i+1] = newdate;
            val_end[i]    = newval;
            val_beg[i+1]  = newval;

            // second solution, switch values
            // Date temp_date = date_end[i];
            // date_end[i]    = date_beg[i+1];
            // date_beg[i+1]  =  temp_date;
            //
            // double temp_val = val_end[i];
            // val_end[i]      = val_beg[i+1];
            // val_beg[i+1]    =  temp_val;
        }
    }

    d["beg"] = date_beg;
    d["end"] = date_end;
    d["y_beg"] = val_beg;
    d["y_end"] = val_end;
    d["len"]   = date_end - date_beg + 1;
    // CharacterVector names = d.attr("names");
    // Rcout << names << std::endl;
    // return max_dte;
}


/*** R
# df <- data.frame(
#     id=1:10,
#     date = seq(as.Date("2015-01-01"),
#                as.Date("2015-01-10"),
#                by="day")
# )
# MaxDate(df)
*/
