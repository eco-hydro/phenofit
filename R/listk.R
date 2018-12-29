#' listk
#' @param ... objects, possibly named.
#' 
#' @examples
#' a = 1
#' b = 1:2
#' c = 1:3
#' l1 <- listk(a, b, c)
#' l2 <- listk(a, b, c = 1:3)
#' l3 <- listk(a = 1, b = 1:2, c = 1:3)
#' l4 <- listk(1, 1:2, c)
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

#' @param x A list object, with data.frame of element
#' @rdname listk
#' @export
list.cbind <- function(x) do.call(cbind.data.frame, x) %>% set_colnames(names(x))

#' @rdname listk
#' @export
list.rbind <- function(x) do.call(rbind.data.frame, x) %>% set_rownames(names(x))#%>% set_rownames(names(x))
