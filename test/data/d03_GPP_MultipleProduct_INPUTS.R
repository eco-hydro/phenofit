source("F:/Github/PML_v2/fluxsites_tidy/R/mainfunc/load_pkgs.R", encoding = "utf-8")
# Write Tn into LAI, NDVI, EVI, for the purpose of eliminating snow illuminated points
# 16 Jan, 2018

################################################################################

files <- "C:/Users/kon055/Google Drive/Github/data/raw/fluxsites_MOD11A2_Tnight_0m_buffer.csv"
df    <- fread(files)
# LST_Night_1km: scale factor, 0.02
df[, `:=`(date = ymd(date),
          Tn = LST_Night_1km*0.02 - 273.15)] #the real value of temperature

df_Tn <- df[, .(date, site, Tn)] %>% add_dn
df_Tn[, c("date", "doy"):=NULL]
df_Tn_16d <- df_Tn[, .(Tn = mean(Tn, na.rm = T)), .(site, year, d16)]

#
indir <- 'C:/Users/kon055/Google Drive/Github/data/phenology/'
df_NDVIv4  = fread(paste0(indir, 'fluxsites_NDVIv4_0m_buffer_1998-2017.csv')) #16DAY
df_LAI     = fread(paste0(indir, 'flux212_LAI.csv'))     #4DAY
df_MOD13A1 = fread(paste0(indir, 'flux212_MOD13A1.csv')) #16DAY
df_MOD13Q1 = fread(paste0(indir, 'flux212_MOD13Q1.csv')) #16DAY
MODGPP     = fread(paste0(indir, 'MOD17A2H_GPP_flux212.csv')) #16DAY

df_NDVIv4$date %<>% ymd()

df_LAI     %<>% add_dn(d16 = FALSE)
df_NDVIv4  %<>% add_dn(d16 = FALSE)
df_MOD13A1 %<>% add_dn(d8 = FALSE)
df_MOD13Q1 %<>% add_dn(d8 = FALSE)
MODGPP %<>% add_dn(d16 = FALSE)

NDVIv4  <- merge(df_NDVIv4 , df_Tn, by = c("site", "year", "d8")) %>% reorder_name
NDVIv4[, `:=`(NDVI = NDVI /1e4)]

LAI     <- merge(df_LAI    , df_Tn, by = c("site", "year", "d8")) %>% reorder_name
MOD13A1 <- merge(df_MOD13A1, df_Tn_16d, by = c("site", "year", "d16")) %>% reorder_name
MOD13Q1 <- merge(df_MOD13Q1, df_Tn_16d, by = c("site", "year", "d16")) %>% reorder_name

# 05. MODIS GPP simulation
MODGPP <- merge(MODGPP, df_Tn, by = c("site", "year", "d8"), all.x = T) %>% reorder_name

vars <- contain(MODGPP, 'GPP')
MODGPP[, (vars) := lapply(.SD, divide_by, 80), .SDcols = vars]

# get the real values of LAI, NDVI, EVI
vars <- contain(df_MOD13A1)
MOD13A1[, (vars) := lapply(.SD, multiply_by, 1e-4), .SDcols = vars]

vars <- contain(MOD13Q1)
MOD13Q1[, (vars) := lapply(.SD, multiply_by, 1e-4), .SDcols = vars]

vars <- contain(LAI, '^LAI')
LAI[, (vars) := lapply(.SD, multiply_by, 1e-1), .SDcols = vars]

source('test/dat/main_QC.R')
# according to QA, initial weights
NDVIv4 [, w := qc_NDVIv4(QA)]
MOD13A1[, w := qc_summary(SummaryQA)]
MOD13Q1[, w := qc_summary(SummaryQA)]
LAI    [, w := qc_5l(FparLai_QC)]
MODGPP [, w := qc_5l(Psn_QC)]

lst <- listk(MODGPP, LAI, MOD13A1, MOD13Q1, NDVIv4)

# match the coresponding date in flux GPP --------------------------------------
# GPPobs = fread(paste0(indir, 'flux129_GPP_16Jan2018.csv')) %>%
#     set_names(c("date", "GPP", "site"))#16DAY
# GPPobs[, date := dmy(date)]
# GPPobs %<>% add_dn

# # mask timeseries
# LAI     %<>% merge(GPPobs[, .(site, date)], by = c("site", "date"))
# NDVIv4  %<>% merge(GPPobs[, .(site, date)], by = c("site", "date"))
# MOD13A1 %<>% merge(GPPobs[, .(site, date)], by = c("site", "date"))
# MOD13Q1 %<>% merge(GPPobs[, .(site, date)], by = c("site", "date"))
# MODGPP  %<>% merge(GPPobs[, .(site, date)], by = c("site", "date"))
#
# ## save files
# fwrite(LAI, 'flux129_Tn_LAI.csv')
# fwrite(MOD13Q1, 'flux129_Tn_MOD13Q1.csv')
# fwrite(MOD13A1, 'flux129_Tn_MOD13A1.csv')
# fwrite(NDVIv4, 'flux129_Tn_NDVIv4.csv')
# fwrite(GPPobs, 'flux129_Tn_GPPobs.csv')
# fwrite(MODGPP, 'flux129_Tn_MODGPP.csv')
#
# lst <- listk(GPPobs, MODGPP, LAI, MOD13A1, MOD13Q1, NDVIv4)
save(lst, file = "Y:/R/phenofit/data/phenofit_MultipleINPUT_flux212.rda")
