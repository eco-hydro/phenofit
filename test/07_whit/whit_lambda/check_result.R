library(Ipaper)
library(magrittr)
source('test/07_whit/whit_lambda/smooth_whit_lambda.R')

lst <- get_sbatch("result/whit_lambda/wSELF_grp1/")
# lst <- get_sbatch("result/whit_lambda/wBisquare_grp1/")

I <- sapply(lst, length) %>% {which(. == 1)} %>% as.numeric() %>%
    sort() %T>% print

res <- optim_lambda_FUN(I[12])
res <- par_sbatch(I[1:8], optim_lambda_FUN, wFUN = wSELF,
                  Save = F, return.res = T,
                  outdir = paste0("result/whit_lambda/wSELF", subfix))

sapply(res, length) %>% {which(. == 1)} %>% as.numeric() %>%
    sort() %T>% print

res <- par_sbatch(sites[134], optim_lambda_FUN, wFUN = wSELF,
                  Save = F, return.res = T,
                  outdir = paste0("result/whit_lambda/wSELF", subfix))
