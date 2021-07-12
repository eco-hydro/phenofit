# list2env(options_season, envir = environment())
options_season <- list(
    rFUN              = "smooth_wWHIT",
    wFUN              = "wTSM",
    iters             = 2,    # rough fitting
    wmin              = 0.1,  # weights updating wmin

    lambda            = NULL,  # only used when rFUN = `smooth_wWHIT`
    .lambda_vcurve    = TRUE, # if lambda not provided, will be optimized by v-curve theory
    frame             = 11, # only used when rFUN = `smooth_wSG`
    nf                = 4,  # only used when rFUN = `smooth_wHANTs`

    minpeakdistance   = NULL, # if `null`, default is `nptperyear/6`
    ypeak_min         = 0.1,

    r_max             = 0.2,
    r_min             = 0.05,
    rtrough_max       = 0.6,

    MaxPeaksPerYear   = 2,
    MaxTroughsPerYear = 3,

    calendarYear      = FALSE,
    adj.param         = TRUE,  # auto adjust rough fitting parameter
    rm.closed         = TRUE,
    is.continuous     = TRUE,
    .check_season     = TRUE,

    # options in `season_mov`
    maxExtendMonth    = 12,   # in month
    nextend           = NULL, # in point, default is `ceiling(maxExtendMonth/12*nptperyear)`

    len_min           = 45,
    len_max           = 650,

    verbose           = FALSE
)

# fine fitting parameters
options_fitting <- list(
    methods            = c("AG", "Beck", "Elmore", "Zhang"),
    wFUN               = "wTSM",
    iters              = 2,   # iterations for fine fitting
    wmin               = 0.1,  # weights updating wmin

    use.y0             = TRUE,
    use.rough          = FALSE,

    nextend            = 2  , # n good-points extend
    minExtendMonth     = 0.5, # in month
    maxExtendMonth     = 1  , # in month
    minPercValid       = 0,   # in 0-1
    verbose            = TRUE
)

.options <- list2env(list(
    # INPUT DATA
    file_vi            = "",
    file_qc            = "",
    qcFUN              = "qc_summary",
    nptperyear         = 23,
    south              = FALSE,
    ymin               = 0.1,
    ws                 = c(0.2, 0.5, 0.8), # initial weights

    # methods
    FUN_season         = "season_mov",
    methods_pheno      = c("TRS", "DER", "Zhang", "Gu"),

    verbose_season_mov = TRUE,
    verbose_season     = FALSE,
    # parameters
    season             = options_season,
    fitting            = options_fitting
))

#' set and get phenofit option
#'
#' @param ... list of phenofit options
#' FUN_season: character, `season_mov` or `season`
#' rFUN: character, rough fitting function. `smooth_wWHIT`, `smooth_wSG` or `smooth_wHANTs`.
#'
#' @examples
#' set_options(verbose_curvefit = FALSE)
#' get_options("verbose_season")
#' @export
set_options <- function(...) {
    opts = list(...)
    modifyList(.options, opts)
    invisible()
}

#' @param names vector of character, names of options
#'
#' @rdname set_options
#' @export
get_options <- function(names = NULL) {
    if (is.null(names)) return(as.list(.options))
    .options[[names]]
}

modifyList <- function (x, val, keep.null = TRUE)
{
    # stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    # if (keep.null) {
    #     for (v in vnames) {
    #         x[v] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]]))
    #             list(modifyList(x[[v]], val[[v]], keep.null = keep.null))
    #         else val[v]
    #     }
    # } else {
        for (v in vnames) {
            if (keep.null) {
                if (is.null(val[[v]])) next()
            }
            x[[v]] <- if (v %in% xnames && is.list(x[[v]]) &&
                is.list(val[[v]]))
                modifyList(x[[v]], val[[v]], keep.null = keep.null)
            else val[[v]]
        }
    # }
    x
}
