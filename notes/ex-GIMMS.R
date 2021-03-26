library(phenofit)
library(glue)
library(ggplot2)
library(lubridate)
library(latticeGrob)

d = read_xlsx("data-raw/gimms_1pnt.xlsx")
# sites = colnames(df)[-1]
doy = d$doy[1:24]
date1 = as.Date(glue("2001-{doy}"), "%Y-%j")[1:24]
date2 = as.Date(glue("2002-{doy}"), "%Y-%j")[1:24]
date3 = as.Date(glue("2003-{doy}"), "%Y-%j")[1:24]
d$date <- c(date1, date2, date3)
d$y %<>% divide_by(1e4)

INPUT <- check_input(d$date, d$y, nptperyear = 24)

wFUN = wTSM # wBisquare #
brks <- season_mov(INPUT,
                    rFUN = smooth_wWHIT, wFUN = wFUN, calendarYear = TRUE,
                    IsPlot = TRUE, IsPlot.OnlyBad = F)

fit  <- curvefits(INPUT, brks,
                  methods = c("AG", "Zhang", "Beck", "Elmore", "Gu"), #,"klos",, 'Gu'
                  wFUN = wFUN,
                  nextend = 2, maxExtendMonth = 3, minExtendMonth = 1, minPercValid = 0.2)

d_fit <- get_fitting(fit)
g <- plot_phenofit(d_fit, brks, NULL, title.ylab = "NDVI", "Time",
                   theme = coord_cartesian(xlim = c(ymd("2000-04-01"), ymd("2004-07-31"))))
# grid::grid.newpage(); grid::grid.draw(g)# plot to check the curve fitting
# write_fig(g, "a.pdf", 9, 6)

## part 2
# pheno: list(p_date, p_doy)
l_pheno <- get_pheno(fit, IsPlot = F) #%>% map(~melt_list(., "meth"))

# ratio = 1.15
# file <- "Figure5_Phenology_Extraction_temp.pdf"
# cairo_pdf(file, 8*ratio, 6*ratio)
# temp <- get_pheno(fit$fits$ELMORE[2:6], IsPlot = T)
# dev.off()
# file.show(file)

## check the extracted phenology
inds = 1:min(6, length(brks))
pheno <- get_pheno(fit[inds], "Zhang", IsPlot = T)
# print(str(pheno, 1))
head(l_pheno$doy$AG)
