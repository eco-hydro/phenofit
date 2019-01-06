par(mfrow = c(1, 5),
    oma = c(1, 2, 3, 1), mar = rep(0, 4), yaxt = "n", xaxt = "n")
pheno <- get_pheno.fFITs(fFITs, "AG", IsPlot = TRUE)
# multiple years
fits <- list(`2001` = fFITs, `2002` = fFITs)
pheno <- get_pheno(fits, "AG", IsPlot=TRUE)
