#' ifelse2
#' 
#' ternary operator just like java `test ? yes : no`
#' 
#' @param test an object which can be coerced to logical mode.
#' @param yes return values for true elements of test.
#' @param no return values for false elements of test.
#' 
#' @examples
#' x <- ifelse2(TRUE, 1:4, 1:10)
#' 
#' @export
ifelse2 <- function(test, yes, no){
    if (test){
        yes
    } else{
        no
    }
}
