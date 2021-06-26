.options <- list2env(list(
    # INPUT DATA
    file_vi            = "", 
    file_qc            = "",
    qcFUN              = "qc_summary",
    nptperyear         = 23, 
    ymin               = 0.1,
    ws                 = c(0.2, 0.5, 0.8), # initial weights
    
    # whether to use original `y0` for plot, which is before the rocess of `check_input`.
    # 
    use.y0             = TRUE, 

    # methods
    FUN_season         = "season_mov",
    rFUN               = "smooth_wWHIT",
    methods_fine       = c("AG", "Beck", "Elmore", "Zhang"),
    methods_pheno      = c("TRS", "DER", "Zhang", "Gu"),
    wFUN_rough         = "wTSM",
    wFUN_fine          = "wTSM",

    frame              = 11, # only used when rFUN = `smooth_wSG`
    lambda             = 2,  # only used when rFUN = `smooth_wWHIT`
    nf                 = 4,  # only used when rFUN = `smooth_wHANTs`

    # parameters
    iters_rough        = 2,   # iterations for rough fitting
    iters_fine         = 2,   # iterations for fine fitting
    
    r_max              = 0.2, # see details in growing season dividing
    r_min              = 0,
    rtrough_max        = 0.8,

    nextend            = 2  , # n good-points extend
    minExtendMonth     = 0.5, # in month
    maxExtendMonth     = 1  , # in month
    minPercValid       = 0,   # in 0-1
  
    verbose_season_mov = TRUE,
    verbose_season     = FALSE,
    verbose_curvefit   = TRUE, 
    calendarYear       = TRUE,
    south              = FALSE
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

modifyList <- function (x, val, keep.null = FALSE) 
{
    # stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    vnames <- vnames[nzchar(vnames)]
    if (keep.null) {
        for (v in vnames) {
            x[v] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
                list(modifyList(x[[v]], val[[v]], keep.null = keep.null))
            else val[v]
        }
    } else {
        for (v in vnames) {
            x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && 
                is.list(val[[v]])) 
                modifyList(x[[v]], val[[v]], keep.null = keep.null)
            else val[[v]]
        }
    }
    x
}
