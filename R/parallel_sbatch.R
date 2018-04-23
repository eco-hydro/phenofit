#' par_sbatch
#'
#' R parallel function, designed for linux sbatch model.
#'
#' @param X a vector (atomic or list) or an expressions vector. Other objects
#' (including classed objects) will be coerced by as.list.
#' @param FUN the function to be applied to (mclapply) each element of X or
#' (mcmapply) in parallel to ....
#' @param nodes The (maximum) number of cluster nodes to spread the calculation
#' over. slurm_apply automatically divides params in chunks of approximately
#' equal size to send to each node. Less nodes are allocated if the parameter
#' set is too small to use all CPUs on the requested nodes.
#' @param cpus_per_node The number of CPUs per node on the cluster;
#' determines how many processes are run in parallel per node.
#' @param ... other parameters passed to mclapply
#'
#' @export
par_sbatch <- function(X, FUN, ..., nodes = 8, cpus_per_node = 16){
    nparams       <- length(X)

    if (nparams < cpus_per_node * nodes) {
        nchunk <- cpus_per_node
    }else {
        nchunk <- ceiling(nparams/nodes)
    }
    nodes <- ceiling(nparams/nchunk)

    I_node <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) #node ID
    I_beg  <- I_node * nchunk + 1
    I_end  <- min((I_node + 1) * nchunk, nparams)

    # cl  <- parallel::makeCluster(cpus_per_node, type = "FORK") #fork only work in linux
    # print(length(cl))

    # If error, still can close cluster works
    # https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/mclapply.html
    res <- tryCatch({
        ## parallel function for parLapplyLB
        # parLapplyLB(cl, sites[I_beg:I_end], calsite_pheno, df = df) #mcmapply,
        mclapply(X[I_beg:I_end], FUN, mc.cores = cpus_per_node, ...)  #SIMPLIFY
    }, error = function(e){
        message(sprintf('[error]: %sn', e$message))
        # stopCluster(cl)
    })

    saveRDS(res, file = paste0('results_', I_node, '.RDS'))
    # stopCluster(cl)
}

#' get_slurm_out
#' Merge the slurm result.
#' @export
get_slurm_out <- function(indir = '.', pattern = 'result.*.RDS',
                          outfile = "result.rda", IsSave = TRUE){
    files <- dir(indir, pattern, full.names = T)
    cat('OUTPUTs:', "\n")
    print(basename(files))
    RES <- lapply(files, readRDS) %>% do.call(c, .)

    if (IsSave) save(RES, file = outfile)
    return(RES)
}
