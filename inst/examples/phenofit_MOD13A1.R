library(phenofit)
data("MOD13A1")

df <- tidy_MOD13(MOD13A1$dt)
st <- MOD13A1$st

date_start <- as.Date('2010-01-01')
date_end   <- as.Date('2016-12-31')

sitename <- 'CA-NS6' # df$site[1]
d     <- df[site == sitename & (date >= date_start & date <= date_end), ]
sp    <- st[site == sitename, ]
south <- sp$lat < 0
nptperyear <- 23

# global parameter
IsPlot = TRUE
print  = FALSE
ypeak_min  = 0.05
wFUN = wTSM

## 1. check_input
# add one year in head and tail
dnew     <- add_HeadTail(d, south = south, nptperyear = nptperyear)
INPUT    <- check_input(dnew$t, dnew$y, dnew$w, QC_flag = dnew$QC_flag,
     nptperyear = nptperyear, south = south,
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)

## 2. Rough fitting and growing season dividing
# Rough fitting and growing season dividing
brks2 <- season_mov(INPUT,
    rFUN = smooth_wWHIT,
    wFUN = wFUN,
    plotdat = d, IsPlot = IsPlot)

## 3. Fine fitting and growing season dividing
fits <- curvefits(
    INPUT, brks2,
    methods = c("AG", "Beck", "Elmore", "Zhang"), #,"klos", "Gu"
    wFUN = wFUN,
    nextend = 2, maxExtendMonth = 2, minExtendMonth = 1, minPercValid = 0.2)

## 4. Phenological metric extraction
l_param   <- get_param(fits)
d_GOF     <- get_GOF(fits)
d_fitting <- get_fitting(fits)
l_pheno   <- get_pheno(fits, "AG", IsPlot=TRUE)

# cairo_pdf('phenofit_MOD13A1.pdf', 12, 6)
cairo_pdf('phenofit_MOD13A1.pdf', 8, 4.5)
g <- plot_curvefits(d_fitting, brks2);
grid::grid.draw(g)
dev.off()

## comparison with TIMESAT and phenopix
# for shiny app: phenofit
# save(df, st, file = "flux10_MOD13A1.rda")
# shiny::runApp("inst/shiny/phenofit/", 80, TRUE)
is_comp <- TRUE # FALSE

if (is_comp) {
    library(data.table)

    d_fit <- merge(d[, .(date, t)], d_fitting)
    d_comp <- fread("../phenofit_comp.txt")
    d_comp$date %<>% ymd()
    d_comp <- d_comp[(date >= date_start & date <= date_end), ]

    dplot <- merge(d_fit, d_comp[, .(date, TIMESAT, phenopix)], all.x = T, by = "date")

    font.size <- 14
    method <- "Beck"
    data <- dplot[meth == method]

    # order  <- c(2, 1, 3)
    colors <- c("red", "blue", "black")
    lgd <- phenofit:::make_legend(linename = c("phenofit", "TIMESAT", "phenopix"),
            linecolor = colors, cex = 1.2)
    p <- ggplot(data, aes_string("t", "y")) +
        geom_point(size = 3, alpha = 0.75,
            aes_string(shape="QC_flag", color = "QC_flag", fill = "QC_flag")) +
        scale_color_manual(values = c(qc_colors,"iter1" = "blue", "iter2" = "red"), drop = F) +
        scale_fill_manual(values = qc_colors, drop = F) +
        scale_shape_manual(values = qc_shapes, drop = F) +
        scale_x_date(date_labels = "%Y/%m", breaks = seq(ymd('20100101'), ymd('20160101'), 'year')) +
        geom_line(aes_string(y = "phenopix"), size = 0.8, alpha = 0.7, color = colors[3]) +
        geom_line(aes_string(y = "TIMESAT"), size = 0.8, alpha = 0.7, color = colors[2]) +
        geom_line(aes_string(y = "ziter2"), size = 0.8, alpha = 0.7, color = colors[1]) +
        theme_gray(base_size = font.size) +
            theme(legend.position="none",
                axis.title = element_text(size = font.size),
                axis.text = element_text(size = font.size - 2)
                # axis.text.x = element_text(angle = 10, hjust = 1, vjust = 1)
                ) +
            labs(x = 'Time', y = 'Vegetation Index')
    g <- arrangeGrob(p, lgd, nrow = 2, heights = c(10, 1),
                padding = unit(1, "line"))

    cairo_pdf('Figure6_comp_TIMESAT_and_phenopix.pdf', 10, 3.5)
    grid::grid.newpage()
    grid::grid.draw(g)
    dev.off()
}
