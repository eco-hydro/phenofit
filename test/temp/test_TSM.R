## Preparing INPUT data for TIMESAT

dt <- df[date >= ymd("2001-01-01")]

nptperyear <- 23
nyear <- 16
npt   <- nptperyear*nyear
ngrid <- 10

file <- "y.txt"
header <- sprintf("%d\t%d\t%d", nyear, nptperyear, ngrid)

write_TSM <- function(x, file){
    write_lines(header, file)
    write.table(x, file, append = T, row.names = F, col.names = F, sep = "\t")
}

y <- matrix(dt$y, ncol = npt, byrow = T)
w <- matrix(dt$w, ncol = npt, byrow = T)

write_TSM(y, "phenofit_st10_y.txt")
write_TSM(w, "phenofit_st10_w.txt")


par(mfrow = c(10, 1), oma = c(1, 3, 1, 1), mar = rep(0, 4))
for (i in 1:10){
    titlestr <- sprintf("[%02d] %s", i, sites[i])
    plot(y[i, ], type = "b", at = seq(0, ncol(y), 23), xlim = c(0, 23*10), main = titlestr)
    abline(v = seq(0, ncol(y), 23), col = "grey60", lty = 2)
}


Cairo::CairoPDF("test.pdf", 10, 5)
for (i in seq_along(ps)){
    print(ps[[i]])
}
dev.off()

