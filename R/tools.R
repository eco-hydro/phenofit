#' @title fprintf
#' @description print sprintf result into console just like C style fprintf function
#' @export
fprintf <- function(fmt, ...) cat(sprintf(fmt, ...))

#' @title runningId
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

#' 
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

#' GOF
#' Good of fitting
#' @export
GOF <- function(Y_obs, Y_sim){
    # remove NA values in Y_sim and Y_obs
    I <- which(!(is.na(Y_sim) | is.na(Y_obs)))
    # n_obs <- length(Y_obs)
    n_sim <- length(I)
    
    Y_sim <- Y_sim[I]
    Y_obs <- Y_obs[I]
    
    if (is_empty(Y_obs)){
        return(c(Bias = NA, MAE = NA,RMSE = NA, NASH = NA, 
             pvalue = NA, n_sim = NA, R = NA)) #R = R, 
    }
    RE    <- Y_sim - Y_obs
    Bias  <- mean(RE)                                        # bias
    MAE   <- mean(abs(RE))                                   # mean absolute error
    RMSE  <- sqrt(sum((RE)^2) / length(Y_obs))               # root mean sqrt error
    NASH  <- 1  - sum((RE)^2) / sum((Y_obs - mean(Y_obs))^2) # NASH coefficient
    
    # Observations number are not same, so comparing correlation coefficient
    # was meaningless.
    R      <- NA
    pvalue <- NA
    tryCatch({
        cor.obj <- cor.test(Y_obs, Y_sim, use = "complete.obs")
        R       <- cor.obj$estimate[[1]]
        pvalue  <- cor.obj$p.value  
    }, error = function(e){
        message(e$message)
    })
    return(c(Bias = Bias, MAE = MAE,RMSE = RMSE, NASH = NASH, 
             pvalue = pvalue, n_sim = n_sim, R = R)) #R = R, 
}