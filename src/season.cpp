#include <Rcpp.h>
using namespace Rcpp;

// 符合条件的两个相邻(t_diff <= 5) seasons, 合并endtime
void RightCombine_season(NumericVector y_peak, NumericVector y_end, NumericVector len,
                         DateVector date_beg, DateVector date_peak, DateVector date_end, int i)
{
    y_end[i] = y_end[i + 1];
    date_end[i] = date_end[i + 1];
    len[i] = date_end[i] - date_beg[i] + 1;

    if (y_peak[i] < y_peak[i + 1])
    {
        date_peak[i] = date_peak[i + 1];
        y_peak[i] = y_peak[i + 1];
    }
    y_peak[i + 1] = -9999.0; // flag
}

void LeftCombine_season(NumericVector y_peak, NumericVector y_end, NumericVector len,
                        DateVector date_beg, DateVector date_peak, DateVector date_end, int i)
{
    y_end[i + 1] = y_end[i];
    date_end[i + 1] = date_end[i];
    len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;

    if (y_peak[i] > y_peak[i + 1])
    {
        date_peak[i + 1] = date_peak[i];
        y_peak[i + 1] = y_peak[i];
    }
    y_peak[i] = -9999.0; // flag
}

// TODO: ADD SOME TEST FOR `check_season`

/**
 * Fix troughs in \code{season_mov} whose date is later than previous year
 * ending trough, and merge two too close troughs (less than 35 days).
 *  // DateVector t, NumericVector ypred,
 */
//' check_season
//' @param d Data.frame of growing season dividing info
//'
//' @inheritParams season
//' @keywords internal
//' @export
// [[Rcpp::export]]
void check_season(DataFrame d,
                  bool rm_closed = true, double rtrough_max = 0.7, double r_min = 0.02)
{
    DateVector date_beg = d["beg"];
    DateVector date_end = d["end"];
    DateVector date_peak = d["peak"];
    NumericVector len = d["len"];
    NumericVector y_beg = d["y_beg"];
    NumericVector y_end = d["y_end"];
    NumericVector y_peak = d["y_peak"];

    int n = d.nrow();
    Date newdate;
    double newval;

    for (int i = 0; i < n - 1; i++)
    {
        // Rcout << i << std::endl;
        int t_diff = date_beg[i + 1] - date_end[i];
        int T2s    = date_end[i + 1] - date_beg[i];

        double T1_minVal = std::min(y_beg[i], y_end[i]);
        double T2_minVal = std::min(y_beg[i + 1], y_end[i + 1]);
        double maxY = std::max(y_peak[i], y_peak[i + 1]);
        double minY = std::min(T1_minVal, T2_minVal);
        double A = maxY - minY;

        double T2_h_left = y_peak[i + 1] - y_beg[i + 1];
        double T1_h_right = y_peak[i] - y_end[i];

        double trs = A * r_min;
        double trs2 = A * (r_min + 0.1);
        // 1. y_end[i] >= trs || y_beg[i+1] >= trs 定义为高trough
        // double diff = std::abs(y_end[i] - y_beg[i+1]);
        // bool is_HighTrough = y_end[i] >= trs || y_beg[i+1] >= trs;
        bool is_PreEndSmaller = (y_end[i] <= y_beg[i + 1]); // previous end smaller ?
        // Rprintf("trs = %.2f\n", trs);
        // Rcout << i << ": is_HighTrough = " << is_HighTrough << std::endl;

        // 1. growing season日期交差                                  : t_diff < 0
        // 2. 两相邻val_troughs挨得太近（且中间不存在较大的peaks）    : t_diff <= 50
        // 3(temp). 两相邻val_troughs值相差过大; 且不重叠，则进行融合
        if (t_diff < 0) //|| (is_HighTrough && delta_days < 0)
        {
            newdate = is_PreEndSmaller ? date_end[i] : date_beg[i + 1];
            newval = is_PreEndSmaller ? y_end[i] : y_beg[i + 1];

            date_end[i] = newdate; // + (-1);
            date_beg[i + 1] = newdate;
            y_end[i] = newval;
            y_beg[i + 1] = newval;
            // len[i] = date_end[i] - date_beg[i] + 1;
            // len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;
        }

        // 坡脚修复大法，IT-BCi, ending of 2005; 必须是相邻troughs
        // Rcout << peak_temp << std::endl;

        if (t_diff > 0 && t_diff <= 150 &&
            T2_h_left <= trs && y_end[i] < y_beg[i + 1])
        {
            // 向左拱进
            y_beg[i + 1] = y_beg[i];
            date_beg[i + 1] = date_beg[i];
            // len[i + 1] = date_end[i + 1] - date_beg[i + 1] + 1;
        }

        // 处理相邻的peaks，进行融合操作
        // 1. 第二周期，左侧trough过高，超过`0.7*A` (BUG found)
        // 2. 第二周期，左侧r_min过小, （其ypeak > r_min + 0.1）
        // 3. 第一周期， 右侧r_min过小,
        // 3. (not used) 如果相邻生长季返青日期过短，则执行合并操作

        // 这里用max，而非min，意在保护A较小的生长季
        if (rm_closed)
        {
            double is_closed = t_diff >= 0 && t_diff <= 150; // 相差在5天内认为相邻

            // double diff_right = y_peak[i + 1] - y_end[i + 1];
            // Rprintf("T1_h_right = %#4.2f, T2_h_left = % 4.2f, trs = %4.2f\n", T1_h_right, T2_h_left, trs);
            bool con_left = (T1_h_right <= trs && y_beg[i] < y_beg[i + 1] && (y_peak[i] > trs2 + T1_minVal));
            bool con_right = (T2_h_left <= trs && y_end[i] > y_end[i + 1] && (y_peak[i + 1] > trs2 + T2_minVal));

            if (is_closed && T2s <= 650 && 
                ((y_end[i] >= A * rtrough_max + T1_minVal) || con_right || con_left))
            {
                // Rcout << i+1 << "正在进行融合" << std::endl;
                RightCombine_season(y_peak, y_end, len, date_beg, date_peak, date_end, i);
                i++; // 如果进行了融合，则跳过下一生长周期
                continue;
            }
            // if (is_closed && con_left)
            // {
            //     LeftCombine_season(y_peak, y_end, len, date_beg, date_peak, date_end, i);
            //     i++; // 如果进行了融合，则跳过下一生长周期
            //     continue;
            // }
        }
    }
    d["len"] = date_end - date_beg + 1;
    // CharacterVector names = d.attr("names");
    // return max_dte;
}

/*** R
*/
