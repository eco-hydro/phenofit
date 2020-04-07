#' tidyFitPheno
#'
#' Tidy for every method with multiple years phenology data
#'
#' @param pheno Phenology metrics extracted from `get_pheno`
#'
#' @keywords internal
#' @examples
#' library(phenofit)
#' # simulate vegetation time-series
#' fFUN = doubleLog.Beck
#' par  = c( mn  = 0.1 , mx  = 0.7 , sos = 50 , rsp = 0.1 , eos = 250, rau = 0.1)
#'
#' t    <- seq(1, 365, 8)
#' tout <- seq(1, 365, 1)
#' y <- fFUN(par, t)
#'
#' methods <- c("AG", "Beck", "Elmore", "Gu", "Zhang") # "Klos" too slow
#' fFITs <- curvefit(y, t, tout, methods)
#'
#' # multiple years
#' fits <- list(`2001` = fFITs, `2002` = fFITs)
#' pheno <- get_pheno(fits, "AG", IsPlot=FALSE)
#' @export
tidyFitPheno <- function(pheno){
    doy2date <- function(datenum) as.Date(unlist(datenum), origin = date.origin)
    # phenonames <- c('TRS2.sos', 'TRS2.eos', 'TRS5.sos', 'TRS5.eos', 'TRS6.sos', 'TRS6.eos',
    #                 'DER.sos', 'DER.pop', 'DER.eos',
    #                 'GU.UD', 'GU.SD', 'GU.DD', 'GU.RD',
    #                 'ZHANG.Greenup', 'ZHANG.Maturity', 'ZHANG.Senescence', 'ZHANG.Dormancy')

    # fix the error of \code{pheno} without name
    # names <- names(pheno)
    # if (is.null(names)){
    #     years <- seq(year(origin), by = 1, length.out = length(pheno))
    #     names(pheno) <- years
    # }
    names <- unlist(pheno[[1]]) %>% names()
    p_date <- ldply(pheno, doy2date, .id = "flag") %>%
        mutate(origin = ymd(paste0(substr(flag, 1, 4), "-01-01"))) %>%
        reorder_name(c("flag", "origin")) %>%
        set_colnames(c("flag", "origin", names)) %>% data.table()

    colnames(p_date) %<>% gsub("GU\\.|ZHANG\\.", "", .)
    phenonames <- setdiff(colnames(p_date),  c("flag", "origin", "meth"))

    df <- p_date %>% melt(c("flag", "origin"), variable.name = "meth",value.name = "date") %>%
        mutate(doy = as.numeric(date - origin + 1))
    # date <- spread(pheno[, c("flag", "meth", "date")], meth, date)
    p_doy <- df[, c("flag", "meth", "doy", "origin")] %>% dcast(flag + origin ~ meth, value.var = "doy")
    vars  <- c("flag", "origin", phenonames)
    list(doy = p_doy[, ..vars], date = p_date[, ..vars])
}
