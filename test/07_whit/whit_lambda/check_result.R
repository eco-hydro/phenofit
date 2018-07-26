library(Ipaper)
library(magrittr)
source('test/07_whit/whit_lambda/smooth_whit_lambda.R')

lst <- get_sbatch("result/whit_lambda/wSELF_grp1/")
# lst <- get_sbatch("result/whit_lambda/wBisquare_grp1/")

sapply(lst, length) %>% {which(. == 1)} %>% sort()


res <- par_sbatch(sites[1:100], optim_lambda_FUN, wFUN = wSELF,
                  Save = F, return.res = T,
                  outdir = paste0("result/whit_lambda/wSELF", subfix))

res <- par_sbatch(sites[134], optim_lambda_FUN, wFUN = wSELF,
                  Save = F, return.res = T,
                  outdir = paste0("result/whit_lambda/wSELF", subfix))
