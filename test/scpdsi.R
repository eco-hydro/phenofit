d <- fread("file:///F:/Github/R_packages/scPDSI/data/weekly/1/potentials")
d[, .(Year, Week, PE)] %>% dcast(Year~Week, value.var = "PE") %>%
    fwrite("weekly_PET", sep = " ", col.names = F)

files <- dir("F:/Github/R_packages/scPDSI-org/data/weekly/", recursive = T, full.names = T) %>%
    set_names({paste0(basename(dirname(.)), "_", basename(.))})
files_new <- gsub("-org", "", files)

i <- 1
diffs <- llply(seq_along(files), function(i){
    org <- fread(files[i])
    new <- fread(files_new[i])

    diff <- new -org
    diff[diff <= 2e-2] <- 0
    res <- unique(diff)
    if(nrow(res) == 1) return(NULL)
    res
}) %>% set_names(names(files)) #%>% rm_empty()
