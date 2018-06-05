#' @title fprintf
#' @description print sprintf result into console just like C style fprintf function
#' @export
fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

#' print the running ID in the console
#'
#' @export
runningId <- function(i, step = 1) {
    if (mod(i, step) == 0) cat(sprintf("running %d ...\n", i))
}

#' retry
#' retry to execute expr until success before reaches the max try times (maxTimes)
#' @export
retry <- function(expr, maxTimes = 3){
    eTimes <- 0
    out    <- NULL
    while (eTimes < maxTimes){
      out <- tryCatch({
          expr
          eTimes <- maxTimes
      }, error = function(e){
          eTimes <<- eTimes + 1
          message(sprintf("[try%d]: %s", eTimes, e))
          NULL #If error return NULL
      })
    }
    return(out)
}

#' melt_list
#'
#' @importFrom reshape2 melt
#' @export
melt_list <- function(list, var.name, na.rm = TRUE, ...){
    if (is.data.table(list[[1]])){
        names <- names(list)
        for (i in seq_along(list)){
            x <- list[[i]]
            eval(parse(text = sprintf("x$%s <- names[i]", var.name)))
            list[[i]] <- x
        }
        res <- do.call(rbind, list)#return
    } else{
        id.vars <- colnames(list[[1]])
        res <- reshape2::melt(list, ..., id.vars = id.vars, na.rm = na.rm)
        colnames(res) <- c(id.vars, var.name)
    }
    return(res)
}
# melt_list <- function(data, var.name, ..., na.rm = FALSE, value.name = "value") {
#   id.vars <- colnames(data[[1]])
#   res <- melt(data, ..., id.vars = id.vars, na.rm = na.rm, value.name = value.name)
#   colnames(res) <- c(id.vars, var.name)
#   return(res)
# }

#' listk
#' @export
listk <- function(...){
  # get variable names from input expressions
  cols <- as.list(substitute(list(...)))[-1]
  vars <- names(cols)
  Id_noname <- if (is.null(vars)) seq_along(cols) else which(vars == "")

  if (length(Id_noname) > 0)
    vars[Id_noname] <- sapply(cols[Id_noname], deparse)
  # ifelse(is.null(vars), Id_noname <- seq_along(cols), Id_noname <- which(vars == ""))
  x <- setNames(list(...), vars)
  return(x)
}

#' list.cbind
#' @export
list.cbind <- function(x) do.call(cbind.data.frame, x) %>% set_colnames(names(x))
#' list.rbind
#' @export
list.rbind <- function(x) do.call(rbind.data.frame, x) %>% set_rownames(names(x))#%>% set_rownames(names(x))

#' reorder_name
#' @export
reorder_name <- function(d,
                         headvars = c("site", "date", "year", "doy", "d8", "d16"),
                         tailvars = ""){
    headvars %<>% intersect(colnames(d))
    tailvars %<>% intersect(colnames(d))
    varnames <- c(headvars, setdiff(colnames(d), union(headvars, tailvars)), tailvars)
    if (is.data.table(d)){
        # d[, ..varnames]
        d[, varnames, with = F] #return
    }else{
        d[, varnames]
    }
}

#' rm_empty
#' @export
rm_empty <- function(x){
    if (is.list(x)){
        x[sapply(x, length) > 0]
    }else {
        x[!is.na(x)]
    }
}

#' contain
#' find assigned pattern variable names
#' @export
contain <- function(d, pattern = "NDVI|EVI") {
    names(d) %>% .[grep(pattern, .)]
}

#' merge_pdf
#'
#' rely on python pdfmerge package, `pip install pdfmerge`
#' @export
merge_pdf <- function(outfile = "RPlot.pdf", indir = 'Figs/', pattern = "*.pdf", del = FALSE){
    # "Y:/R/phenofit/Figs/"
    files <- dir(indir, pattern, full.names = T)
    cmd <- sprintf("pdfmerge -o %s %s", outfile, paste(files, collapse = " "))
    shell(cmd, wait = FALSE)
    # shell(sprintf('pdfmerge -o %s %s', outfile, pattern) )
    if (del) shell(sprintf('del %s', pattern))
}

#' weighted CV
#' @export
cv_coef <- function(x, w){
    if (missing(w)) w <- rep(1, length(x))
    if (length(x) == 0){
        return( c(mean = NA, sd = NA, cv = NA) )
    }

    mean <- sum(x * w) / sum(w)
    sd   <- sqrt(sum((x  - mean)^2 * w) /sum(w))
    cv   <- sd / mean
    c(mean = mean, sd = sd, cv = cv) # quickly return
}

#' GOF
#' 
#' Good of fitting
#'
#' @param Y_obs Numeric vector, observations
#' @param Y_sim Numeric vector, corresponding simulated values
#' @param w Numeric vector, weights of every points
#'
#' @export
GOF <- function(Y_obs, Y_sim, w, include.cv = FALSE){
    if (missing(w)) w <- rep(1, length(Y_obs))

    # remove NA values in Y_sim, Y_obs and w
    I <- which(!(is.na(Y_sim) | is.na(Y_obs) | is.na(w)))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)

    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]

    if (include.cv) CV <- cv_coef(Y_obs, w)
    if (is_empty(Y_obs)){
        out <- c(Bias = NA, MAE = NA,RMSE = NA, NSE = NA, R2 = NA,
             pvalue = NA, n_sim = NA, R = NA)
        if (include.cv) out <- c(out, CV)
        return(out) #R = R,
    }

    # R2: the portion of regression explained variance, also known as
    # coefficient of determination
    #
    # https://en.wikipedia.org/wiki/Coefficient_of_determination
    # https://en.wikipedia.org/wiki/Explained_sum_of_squares
    y_mean <- sum(Y_obs * w) / sum(w)

    SSR    <- sum( (Y_sim - y_mean)^2 * w) 
    SST    <- sum( (Y_obs - y_mean)^2 * w)
    R2     <- SSR / SST

    RE     <- Y_sim - Y_obs
    Bias   <- sum ( w*RE)     /sum(w)                     # bias
    MAE    <- sum ( w*abs(RE))/sum(w)                     # mean absolute error
    RMSE   <- sqrt( sum(w*(RE)^2)/sum(w) )                # root mean sqrt error

    # https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
    NSE  <- 1  - sum( (RE)^2 * w) / SST # NSE coefficient

    # Observations number are not same, so comparing correlation coefficient
    # was meaningless.
    # In the current, I have no idea how to add weights `R`.
    R      <- NA
    pvalue <- NA
    tryCatch({
        cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
        R       <- cor.obj$estimate[[1]]
        pvalue  <- cor.obj$p.value
    }, error = function(e){
        message(e$message)
    })

    out <- c(Bias = Bias, MAE = MAE,RMSE = RMSE, NSE = NSE, R2 = R2,
             pvalue = pvalue, n_sim = n_sim, R = R)
    if (include.cv) out <- c(out, CV)
    return(out)
}


R2_sign <- function(n, NumberOfPredictor = 2, alpha = 0.05){
    freedom_r = NumberOfPredictor - 1 # regression
    freedom_e = n - NumberOfPredictor # error 

    F  = qf(1 - alpha, freedom_r, freedom_e)
    R2 = 1 - 1/(1 + F*freedom_r/freedom_e)
    
    # F = 485.1
    # F = R2/freedom_r/((1-R2)/freedom_e)
    # Rc = sqrt(/(qf(1 - alpha, 1, freedom) + freedom)) %T>% print  # 0.11215    
    return(list(F = F, R2 = R2))
}


