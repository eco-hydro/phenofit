lm_info <- function(df){
    m <- lm(total~year, df)
    s <- summary(m)
    
    pvalue <- s$fstatistic %>% {1 - pf(.[1], .[2], .[3])}
    r2 <- s$r.squared
    coef <- s$coefficients[, 1]
    
    txt1 <- sprintf("y = %.3fyear - %.3f", coef[2], coef[1])
    txt2 <- sprintf("R^2 = %.3f, pvalue = %.2e", r2, pvalue)
    return(list(txt1 = txt1, txt2 = txt2))
}

lm_eqn <- function(m) {
    m <- lm(total~year, df)
    s <- summary(m)
    
    pvalue <- s$fstatistic %>% {1 - pf(.[1], .[2], .[3])}
    r2 <- s$r.squared
    coef <- s$coefficients[, 1]
    
    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              r2 = format(summary(m)$r.squared, digits = 3));
    
    if (coef(m)[2] >= 0)  {
        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    } else {
        eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
    }
    
    as.character(as.expression(eq));                 
}