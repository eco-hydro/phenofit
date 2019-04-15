wFUN = wTSM#wBisquare#wTSM

fit  <- curvefits(INPUT, brks2, lambda =lambda,
                  methods = c("zhang"), #,"klos",, 'Gu'， "AG",, "beck", "elmore"
                  nptperyear = nptperyear, debug = F, wFUN = wFUN,
                  nextend = 2, maxExtendMonth = 3, minExtendMonth = 1,
                  qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
                  print = print)
fit$INPUT   <- INPUT
fit$seasons <- brks2

## check the curve fitting parameters
params <- getparam(fit)
print(str(params, 1))
print(params$AG)

## Get GOF information
stat  <- ldply(fit$fits, function(fits_meth){
    ldply(fits_meth, statistic.fFIT, .id = "flag")
}, .id = "meth")
fit$stat <- stat
print(head(stat))
print(mean(stat$NSE, na.rm = T))
## visualization
# svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# Cairo::CairoPDF(file_pdf, 11, 6) #
# dev.off()
g <- plot_phenofit(fit, d, titlestr)
grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
