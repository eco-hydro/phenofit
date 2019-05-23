#' phenofit_TSM.avhrr
#'
#' @inheritParams phenofit_process
#' @param I_part numeric vector, index of pixels needed to run.
#' @param outdir directory of output
#' @param exportType could be one of `"all", "pheno"`. If "all" used, all
#' also will be exported. Note that exported fitting series is daily scale, which
#' is quite large.
#' @param overwrite boolean, whether overwrite previous output?
#' 
#' @importFrom readr read_rds write_rds
#' @export
phenofit_TS.avhrr <- function(
    options,
    dateRange = c(as.Date('2010-01-01'), as.Date('2014-12-31')),
    I_part = NULL, 
    outdir = ".",
    exportType = "all", overwrite = FALSE, 
    .progress = NULL, .parallel = FALSE,
    ...)
{
    file_y  <- options$file_veg_text
    file_qc <- options$file_qc
    nptperyear <- options$nptperyear
    ymin       <- options$ymin

    d_date <- get_date_AVHRR()
    t <- d_date$date

    if (!is.null(dateRange)) {
        I_date <- which(t >= dateRange[1] & t <= dateRange[2])
    } else {
        I_date <- seq_along(t)
    }

    mat_y  <- fread(file_y, skip = 1) %>% as.matrix()
    mat_qc <- fread(file_qc, skip = 1) %>% as.matrix()
    n <- nrow(mat_y)

    if (is.null(I_part) || is.na(I_part)) {
        I_part <- 1:nrow(mat_qc)
    }
    mat_y  <- mat_y[I_part, ]
    mat_qc <- mat_qc[I_part, ]

    showProgress <- !is.null(.progress) # this for shinyapp progress
    if (showProgress){
        on.exit(.progress$close())
        .progress$set(message = sprintf("phenofit (n=%d) | running ", n), value = 0)
    }

    # rv    <- phenofit_loaddata(options, ...)
    # sites <- rv$sites
    # n     <- length(sites)

    # if (nsite > 0) n <- pmin(n, nsite)
    # sites <- seq_len(n) %>% as.character()

    FUN <- ifelse(.parallel, `%dopar%`, `%do%`)
    res <- FUN(foreach(
        y = iter(mat_y, by='row'),
        qc = iter(mat_qc, by='row'),
        i = I_part,
        .packages = c("phenofit")), {

        outfile <- sprintf("%s/phenofit_%05d.RDS", outdir, i)
        if (!file.exists(outfile) || overwrite) {
            # sitename <- rv$sites[i]
            if (showProgress){
                .progress$set(i, detail = paste("Doing part", i))
            }
            fprintf("phenofit (n = %d) | running %03d ... \n", n, i)

            tryCatch({
                wmin <- 0.4
                l_w  <- qc_summary(qc[I_date], wmin = wmin, wmid = 0.5, wmax = 0.8)
                INPUT <- check_input(t[I_date], y[I_date], w = l_w$w, QC_flag = l_w$QC_flag,
                    nptperyear = nptperyear, south = FALSE, ymin = ymin, wmin = wmin)
                # INPUT <- with(rv, getsite_INPUT(df, st, sitename, nptperyear, dateRange))
                brks  <- phenofit_season(INPUT, options, IsPlot = FALSE, verbose = FALSE)
                fits  <- phenofit_finefit(INPUT, brks, options) # multiple methods
                if (exportType == "pheno") { fits <- fits[-(1:3)] }
                # fits
                write_rds(fits, outfile)
            # }, warning = function(w){
            #     message(sprintf('[w] phenofit_process, i=%d: %s', i, w$message))
            }, error = function(e){
                message(sprintf('[e] phenofit_process, i=%d: %s', i, e$message))
            })

        } else {
            fprintf("[file exist] : %s\n", basename(outfile))
        }
    })
    fprintf('Success finished!')
    ############################# CALCULATION FINISHED #####################
    # set_names(res, sites[1:n])
}
