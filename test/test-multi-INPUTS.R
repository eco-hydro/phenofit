files <- dir("F:/Github/PML_v2/fluxsites_tidy/data/GPP_Pheno", "*.csv", full.names = T) %>%
    set_names(gsub("flux129_Tn_|.csv", "", basename(.)))

## load stations, and remove sites in South Hemisphere
sites_rm1 <- c()
sites_rm2 <- c("GF-Guy", "BR-Sa3", "US-Whs")
sites_rm  <- union(sites_rm1, sites_rm2)

indir <- 'F:/Github/PML_v2/fluxsites_tidy/data/GPP_Pheno/stations/'
station112 <- fread(paste0(indir, "fluxsites_112.csv"))[lat > 0 & !(site %in% sites_rm), ]
station129 <- fread(paste0(indir, "fluxsites_129.csv"))[lat > 0 & !(site %in% sites_rm), ]

station129

lst <- llply(files, fread, .progress = "text")
# 1. GPPobs  |
# 2. LAI     | FparLai_QC
# 3. MOD13A1 | SummaryQA
# 4. MOD13Q1 | SummaryQA, same as MOD13A1
# 5. MODGPP  | Psn_QC
# 6. NDVIv4  | QA