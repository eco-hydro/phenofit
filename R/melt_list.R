#' melt_list
#' 
#' @param list Data set to melt
#' @param var.name list names convert into var.name
#' @param na.rm Should NA values be removed from the data set? 
#' This will convert explicit missings to implicit missings.
#' @param ... other parameters to melt.
#' 
#' @examples
#' # data.frame
#' df <- data.frame(year = 2010, day = 1:3, month = 1, site = "A")
#' l  <- list(a = df, b = df)
#' df_new <- melt_list(l, "id")
#' 
#' # data.table
#' df <- data.table::data.table(year = 2010, day = 1:3, month = 1, site = "A")
#' l  <- list(a = df, b = df)
#' df_new <- melt_list(l, "id")
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
