# list2env(options_season, envir = environment())
options_season <- list(
    rFUN              = "smooth_wWHIT",
    iters             = 2,    # rough fitting
    wmin              = NULL,  # weights updating wmin
    wFUN              = NULL,
    verbose           = NULL,

    lambda            = NULL,  # only used when rFUN = `smooth_wWHIT`
    .lambda_vcurve    = TRUE, # if lambda not provided, will be optimized by v-curve theory
    frame             = 11, # only used when rFUN = `smooth_wSG`
    nf                = 4,  # only used when rFUN = `smooth_wHANTs`

    minpeakdistance   = NULL, # if `null`, default is `nptperyear/6`
    ypeak_min         = 0.1,

    r_max             = 0.2,
    r_min             = 0.05, # 0
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
    len_max           = 650
)

# fine fitting parameters
options_fitting <- list(
    methods            = c("AG", "Beck", "Elmore", "Zhang"),
    iters              = 2,   # iterations for fine fitting
    wFUN               = NULL,
    wmin               = NULL,  # weights updating wmin
    verbose            = NULL,

    use.y0             = TRUE,
    use.rough          = FALSE,

    nextend            = 2  , # n good-points extend
    minExtendMonth     = 0.5, # in month
    maxExtendMonth     = 1  , # in month
    minPercValid       = 0   # in 0-1
)

.options <- list2env(list(
    # INPUT DATA
    file_vi            = "",
    file_qc            = "",
    qcFUN              = "qc_summary",
    nptperyear         = 23,
    south              = FALSE,
    ymin               = 0.1,
    wFUN               = "wTSM",
    wmin               = 0.1, 
    
    ws                 = c(0.2, 0.5, 0.8), # initial weights

    # methods
    FUN_season         = "season_mov",
    methods_pheno      = c("TRS", "DER", "Zhang", "Gu"),

    verbose            = FALSE,
    debug              = FALSE,
    # parameters
    season             = options_season,
    fitting            = options_fitting, 
    initialized        = FALSE
))

#' set and get phenofit option
#'
#' @param ... list of phenofit options
#' FUN_season: character, `season_mov` or `season`
#' rFUN: character, rough fitting function. `smooth_wWHIT`, `smooth_wSG` or `smooth_wHANTs`.
#'
#' @examples
#' set_options(verbose = FALSE)
#' get_options("season") %>% str()
#' @export
set_options <- function(...) {
    opts = list(...)
    # rm unrelated parameters
    ind = match(names(opts), names(.options)) %>% which.notna()
    # This step might lead to error, but will improve performance
    if (.options$initialized && length(ind) == 0) return()

    .options %<>% modifyList(opts[ind])
    # `season` and `fitting` share the same parameter
    pars_comm = c("wFUN", "wmin", "verbose") 
    for (par in pars_comm) {
        if (is.null(.options$fitting[[par]])) .options$fitting[[par]] <- .options[[par]]
        if (is.null(.options$season[[par]])) .options$season[[par]] <- .options[[par]]
    }
    
    .options$fitting$wFUN %<>% check_function()
    .options$season$wFUN %<>% check_function()
    # .options$season$rFUN %<>% check_function()
    invisible()
}

check_function <- function(fun) {
    if (is.character(fun)) {
        name = fun
        fun = get(name)
        attr(fun, "name") = name
    }
    return(fun)
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
