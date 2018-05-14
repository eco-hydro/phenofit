

#' Get the point that Tn is continuously greater (or less) Tmin
findBrks <- function(x, nups = 3, ndowns = 3, zero = "0"){
    if (!zero %in% c("0", "+", "-"))
        stop("Argument 'zero' can only be '0', '+', or '-'.")

    x[is.na(x)] <- 0
    xc <- paste(as.character(sign(x)), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    if (zero != "0")      xc      <- gsub("0", zero, xc)
    peakpat_1 <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
    peakpat_2 <- sprintf("[-]{%d,}", ndowns)

    I   <- gregexpr(peakpat_1, xc)[[1]]
    len <- attr(I, 'match.length')

    n   <- length(I)
    pos <- numeric(n)
    for (i in 1:n){
        i_beg  <- I[i]
        i_end  <- i_beg + len[i] - 1
        pos[i] <- gregexpr(peakpat_2, substr(xc, i_beg, i_end))[[1]][1] + i_beg - 1
    }
    return(pos)
}

nogrowthPolygon <- function(d, Tmin){
    x <- d$Tn#[100:400]
    I_ends   <- findBrks(x > Tmin, zero = "-", nups = 1)
    I_begins <- findBrks(x < Tmin, zero = "-", nups = 1) - 1

    plot(type = "b", x); grid(); abline(a = Tmin, b = 0, col = "red")
    points(I_ends, x[I_ends], col = "red", pch = 19)
    points(I_begins, x[I_begins], col = "blue", pch = 19)

    nptperyear <- 23
    polygons <- list()
    for (i in seq_along(I_ends)){
        i_beg <- I_ends[i]
        i_end <- I_begins[which(I_begins > i_beg & I_begins < i_beg + nptperyear*2/3)][1]
        cat(i, i_end)

        if (!is_empty(i_end)){
            polygons[[i]] <- c(i_beg, -Inf, i_end, -Inf, i_end, Inf,i_beg, Inf,i_beg, -Inf) %>%
                matrix(ncol = 2, byrow = T) %>%
                set_colnames(c("x", "y")) %>% as.data.frame()
        }
    }
    polygons %<>% set_names(seq_along(.))
    df_polygon <- melt_list(polygons, "id")
    df_polygon$date <- d$date[df_polygon$x]
    return(df_polygon)
}

# df_polygon0 <- nogrowthPolygon(d, 0)
# df_polygon5 <- nogrowthPolygon(d, 5)
#
# p + geom_polygon(data = df_polygon5, aes(date, y, color = NULL, shape = NULL),
#                  fill = "red", alpha = 0.2) +
#     geom_polygon(data = df_polygon0, aes(date, y, color = NULL, shape = NULL),
#                  fill = "red", alpha = 0.4)


# ggplot(d, aes(x = seq_along(y), y = y)) + geom_point() +
#     geom_line() +
#     scale_x_continuous(breaks = seq(0, nrow(d), 23)) +
#     theme(
#         panel.grid.minor = element_blank()
#     ) +
#     geom_polygon(data = polygons, aes(x, y), fill = "red", alpha = 0.2)
