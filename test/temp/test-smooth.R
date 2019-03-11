source("test/load_pkgs.R")

load('F:/Github/PML_v2/fluxsites_tidy/data/INPUTS/flux129_Phenofit_multi_INPUTS.rda')
stations <- fread('F:/Github/PML_v2/fluxsites_tidy/data/INPUTS/fluxsites_129.csv')
################################################################################

show_legend <- function(){
    legend('topright', c(paste("iteration", 1:3), "0.5% and 99.5% quantile", "peak", "trough"),
           col = c("red", "blue", "green", "red", "red", "blue"),
           lty = c(rep(1, 3), 2, 0, 0), lwd = c(rep(2, 3), 1, 1, 1),
           pch = c(rep(-1, 4), 16, 16))#, bty='n'
}
################################################################################

lst <- map(lst, ~merge(.x, stations[, .(site, lat, IGBP)], by = "site"))

df <- lst$MODGPP#GPPobs
sites <- unique(df$site)

nptperyear <- floor(365/as.numeric(difftime(df$date[2], df$date[1], units = "day")))

var <- colnames(df) %>% .[grep("NDVI|EVI|LAI|GPP", .)]
if (nptperyear > 90){
    cmd <- sprintf("df_season <- df[, .(%s = median(%s, na.rm = T), date = first(date)),
                    .(site,  year, d8)]", var, var)
    eval(parse(text = cmd))
}else{
    df_season <- df
}

## recalculate FAILED points in MODIS dataset
# DE-RuS FI-Lom FI-Sod US-Prr
# 41     49     50     96
I_cal <- 1:length(sites)
I_rem <- c(49, 50, 96)
# seq_along(sites)

fits_site <- list()
get_pheno <- function(site, var = "MODGPP_0m"){
    x <- df[site == sitename]
    tryCatch({
        fit <- curvefit_site(x$date, x$GPP_NT, lambda =1e4,
                             methods = c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos"
                             nptperyear = 365, debug = F, wFUN = bisquare,
                             south = x$lat[1] < 0)
        # p <- plot_phenofit(fit) + ggtitle(file); print(p)
        fit
    }, error = function(e){
        message(sprintf("%s:%s", file, e$message))
    })
}
RES <- llply(sites, get_pheno, .progress = "text")
save(RES, file = infile)

df <- df[lat >= 0,]
sites <- unique(df$site)
nptperyear <- 365


# curvefit_site <- function(t, y, w, lambda, nptperyear = 46, south = FALSE,
#     iters = 2, IsPlot = FALSE,
#     methods = c('AG', 'zhang', 'beck', 'elmore', 'Gu'), ...)
# mar = c(3, 3, 1, 1),

## Fig.2 dividing growing season
cairo_pdf("Fig2_dividing growing season.pdf", 9, 6.5)
win.metafile("Fig2_dividing growing season.emf", 9, 6.5)

## 0. theory illustration, IT-CA1 as an example
sitename <- "IT-CA1"
x <- df[site == sitename]
try_season(x)

show_legend()
dev.off()

write.xlsx(y, "flux136_detailInfo.xlsx")
