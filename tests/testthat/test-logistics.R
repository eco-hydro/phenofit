context("logistics")

# source('helper_MOD13A1.R')
# wFUN = wTSM

# lambda   <- init_lambda(INPUT$y) # lambda for whittaker
# # param = listk(
# #     INPUT, nptperyear,
# #     FUN = whitsmw2, wFUN = wBisquare, iters = 2,
# #     lambda,
# #     IsPlot = IsPlot, plotdat = d,
# #     south = sp$lat[1] < 0,
# #     rymin_less = 0.6, ypeak_min = ypeak_min,
# #     max_MaxPeaksperyear =2.5, max_MinPeaksperyear = 3.5
# # )
# brks2 <- season_3y(INPUT, nptperyear, south = sp$lat[1] < 0, 
#     FUN = whitsmw2, wFUN = wFUN,
#     plotdat = d, IsPlot = IsPlot, print = F, partial = F)

# fit  <- curvefits(INPUT, brks2, lambda =lambda,
#                   methods = c("zhang"), #,"klos",, 'Gu'， "AG",, "beck", "elmore"
#                   nptperyear = nptperyear, debug = F, wFUN = wFUN,
#                   nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
#                   qc = as.numeric(dnew$SummaryQA), minPercValid = 0.2,
#                   print = print)
# fit$INPUT   <- INPUT
# fit$seasons <- brks2

# ## check the curve fitting parameters
# params <- getparam(fit)
# # print(str(params, 1))
# # print(params$AG)

# ## Get GOF information
# stat  <- ldply(fit$fits, function(fits_meth){
#     ldply(fits_meth, statistic.phenofit, .id = "flag")
# }, .id = "meth")
# fit$stat <- stat
# print(head(stat))
# # print(mean(stat$NSE, na.rm = T))
# ## visualization
# # svg("Figure1_phenofit_curve_fitting.svg", 11, 7)
# # Cairo::CairoPDF(file_pdf, 11, 6) #
# # dev.off()
# titlestr = "test logistics"
# g <- plot_phenofit(fit, d, titlestr)
# grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
