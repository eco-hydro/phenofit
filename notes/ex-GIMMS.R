library(phenofit)
library(glue)
library(ggplot2)
library(lubridate)
library(latticeGrob)

d = data.table::fread("data-raw/gimms_1pnt.csv")
# sites = colnames(df)[-1]
INPUT <- check_input(d$date, d$y)

wFUN = "wTSM" # "wBisquare"
set_options(
    rFUN = "smooth_wWHIT",
    wFUN_fine  = wFUN,
    wFUN_rough = wFUN,
    methods_fine  = c("AG", "Beck", "Elmore", "Zhang"),
    methods_pheno = c("TRS", "DER", "Zhang", "Gu")
)

brks <- season_mov(INPUT, calendarYear = TRUE, IsPlot = TRUE)
fit  <- curvefits(INPUT, brks, nextend = 2, maxExtendMonth = 3, minExtendMonth = 1, minPercValid = 0.2)

d_fit <- get_fitting(fit)
# g <- plot_phenofit(d_fit, brks, NULL, title.ylab = "NDVI", "Time",
#                    theme = coord_cartesian(xlim = c(ymd("2000-04-01"), ymd("2004-07-31"))))
# grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
# write_fig(g, "a.pdf", 9, 6)

## part 2
l_pheno <- get_pheno(fit, IsPlot = F)$doy #%>% map(~melt_list(., "meth"))

