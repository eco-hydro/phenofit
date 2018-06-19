source('test/stable/load_pkgs.R')
source('R/smooth_whit.R')
source('R/smooth_whit_lambda.R')
# source('test/GEE/V-pack.r')
file = 'data_whit_lambda_01.csv'

if (!file.exists(file)){
    dt = fread('test/GEE/data/MOD13A1_st_1e3.csv')
    dt[, `:=`(y    = EVI/1e4,
              t    = ymd(date),
              w    = qc_summary(dt$SummaryQA))]
    dt[, per := sum(!is.na(EVI))/.N, site]
    df <- dt[per > 0.3, .(site, y, t, w, IGBPcode)]
    fwrite(df, file)
}

df         <- fread(file)
sites      <- unique(df$site)
nptperyear <- 23

optim_lambda <- function(sitename, df, IsPlot = F){
    # sitename <- sites[i]#; grp = 1
    d    <- df[site == sitename]
    cat(sprintf('site: %s ...\n', sitename))

    IGBP <- d$IGBPcode[1]
    INPUT <- check_input(d$t, d$y, d$w, trim = T, maxgap = nptperyear / 4, alpha = 0.02)

    # optim lambda by group
    get_lambda <- function(j = NULL){
        tryCatch({
            if (is.null(j)){
                input <- INPUT
            } else{
                I  <- ((j-1)*3*23+1):(j*3*23)
                input <- lapply(INPUT[1:3], `[`, I) %>% c(INPUT[5])
            }

            vc    <- v_curve(input, nptperyear, llas = seq(-2, 3, by = 0.01), d = 2,
                show = IsPlot, iters = 1)
            listk(site = sitename, IGBP, lambda = vc$lambda) #, vc
        }, error = function(e){
            message(sprintf("[e] %s, %d: %s", as.character(sitename), j, e$message))
        })   
    }

    # group = F # three year group
    if (group) {
        temp <- llply(1:6, get_lambda)  
    }else{
        temp <- get_lambda()
    }
    temp # return
}

group = F # three year group
outdir <- ifelse(group, '_grp', '')
outdir <- paste0("result", outdir)

# cpus_per_node <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
par_sbatch(sites, optim_lambda, df = df, save = T, outdir = "result")
# system.time({ res <- pbmclapply(sites, optim_lambda, df = df, mc.cores = cpus_per_node) })
