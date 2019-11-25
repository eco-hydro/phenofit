source("R/PhenoBrks.R")
source("test/GEE/test-func.R")
source("test/GEE/pkg_seasonality.R")

file <- "file:///C:/Users/kon055/Google Drive/Github/data/MOD13A1_st_20_0m_buffer.csv"
wmin <- 0.1

dt <- tidy_gee_MOD13A1(file)
dt[, per := sum(!is.na(y) & w >= 0.5)/.N, site]

df <- dt[per > 0.3]
df <- df[order(IGBPcode, site, t), ]
# fwrite(df, "phenofit_seasonality_st20.csv")

## global parameters for phenofit
sites      <- unique(df$site)
nptperyear <- 23
wFUN       <- wTSM
ymax_min   <- 0.1

# info <- rm_nonSIPoints(df, IsPlot = T, file = 'SI_2.pdf')
################################################################################

methods <- c("sgfitw", "whitsmw2", "smooth_wHANTS")
method  <- methods[2] #"sgfitw", "whitsmw2" and "smooth_wHANTS".
FUN  <- get(method)
file <- sprintf("st340_%s.pdf", method)

CairoPDF(file, width = 10, height = 12)
par(mfrow = c(6, 1), mar = c(1, 2, 3, 1), mgp = c(1.5, 0.6, 0))

library()
stats <- list()
for (i in 1:length(sites)){
    runningId(i)
    sitename <- sites[i]
    # sitename <- "AU-How"
    d <- df[site == sitename, ]
    stats[[i]] <- whit_brks(d, nptperyear, FUN, frame = 16)
}
dev.off()
file.show(file)


info   <- do.call(rbind, stats) %>% data.table()

min_n  <- nrow(d)*0.8 #333, min valid obs
sites  <- info[NSE < 0.3 & n_sim >= min_n, site] # sites remove

info$n_sim %>% {.[which(. < )]}
################################################################################
## get statistics data
stats <- dlply(df, .(site), function(d){
    whit_brks(d, nptperyear, FUN, IsPlot = F)
}, .progress = "text")

info_hant <- rm_nonSIPoints(df, IsPlot = T, file = 'SI_st20.pdf')

## visualization
# 1. whit  seasonality
ggplot(info_hant[NSE < 0.2], aes(NSE, cv, color = as.factor(IGBPcode))) +
    geom_point() +
    geom_text_repel(aes(label = site), show.legend = F) +
    lims(x = c(0, 0.4))

# 2. HANTS seasonality
ggplot(info_hant[NSE < 0.2], aes(NSE, cv, color = as.factor(IGBPcode))) +
    geom_point() +
    geom_text_repel(aes(label = site), show.legend = F) +
    geom_hline(yintercept = 0.2, color = "red") +
    geom_vline(xintercept = 0  , color = "red") +
    lims(x = c(-0.3, 0.2), y = c(0.05, 0.4))
# writelist_ToXlsx(list(HANTS = info_hant, WHIT = info_whit), "info.xlsx")
