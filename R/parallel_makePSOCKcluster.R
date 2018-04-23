makePSOCKcluster <- function(names, ...)
{
    newPSOCKnode <- function(machine = "localhost", ...,
                         options = defaultClusterOptions, rank)
    {
        options <- addClusterOptions(options, list(...))
        if (is.list(machine)) {
            options <- addClusterOptions(options, machine)
            machine <- machine$host
        }
        outfile <- getClusterOption("outfile", options)
        master <- if (machine == "localhost") "localhost"
        else getClusterOption("master", options)
        port <- getClusterOption("port", options)
        manual <- getClusterOption("manual", options)
        timeout <- getClusterOption("timeout", options)
        methods <- getClusterOption("methods", options)
        useXDR <- getClusterOption("useXDR", options)

        ## build the local command for starting the worker
        env <- paste0("MASTER=", master,
                     " PORT=", port,
                     " OUT=", outfile,
                     " TIMEOUT=", timeout,
                     " XDR=", useXDR)
        arg <- "parallel:::.slaveRSOCK()"
        rscript <- if (getClusterOption("homogeneous", options)) {
            shQuote(getClusterOption("rscript", options))
        } else "Rscript"
        rscript_args <- getClusterOption("rscript_args", options)
        if(methods) rscript_args <-c("--default-packages=datasets,utils,stats,methods", rscript_args)

        ## in principle we should quote these,
        ## but the current possible values do not need quoting
        cmd <- if(length(rscript_args))
            paste(rscript, paste(rscript_args, collapse = " "),
                  "-e", shQuote(arg), env)
        else paste(rscript, "-e", shQuote(arg), env)

        ## We do redirection of connections at R level once the process is
        ## running.  We could instead do it at C level here, at least on
        ## a Unix-alike.
        renice <- getClusterOption("renice", options)
        if(!is.na(renice) && renice) ## ignore 0
            cmd <- sprintf("nice +%d %s", as.integer(renice), cmd)

        if (manual) {
            cat("Manually start worker on", machine, "with\n    ", cmd, "\n")
            utils::flush.console()
        } else {
            ## add the remote shell command if needed
            if (machine != "localhost") {
                ## This assumes an ssh-like command
                rshcmd <- getClusterOption("rshcmd", options)
                user <- getClusterOption("user", options)
                ## this assume that rshcmd will use a shell, and that is
                ## the same shell as on the master.
                cmd <- shQuote(cmd)
                cmd <- paste(rshcmd, "-l", user, machine, cmd)
            }

            if (.Platform$OS.type == "windows") {
                ## snow said:
                ## On Windows using input = something seems needed to
                ## disconnect standard input of an ssh process when run
                ## from Rterm (at least using putty's plink).  In
                ## principle this could also be used for supplying a
                ## password, but that is probably a bad idea. So, for now
                ## at least, on Windows password-less authentication is
                ## necessary.
                ##
                ## (Not clear if that is the current behaviour: works for me)
              # system(cmd, intern = T)
              cat(cmd, sep = "\n")
              system(cmd, wait = FALSE, input = "")
            }
            else system(cmd, wait = FALSE)
        }

        con <- socketConnection("localhost", port = port, server = TRUE,
                                blocking = TRUE, open = "a+b", timeout = timeout)
        structure(list(con = con, host = machine, rank = rank),
                  class = if(useXDR) "SOCKnode" else "SOCK0node")
    }#end function: newPSOCKnode

    if (is.numeric(names)) {
        names <- as.integer(names[1L])
        if(is.na(names) || names < 1L) stop("numeric 'names' must be >= 1")
        names <- rep('localhost', names)
    }
    .check_ncores(length(names))
    options <- addClusterOptions(defaultClusterOptions, list(...))
    cl <- vector("list", length(names))
    for (i in seq_along(cl)){
      cat(sprintf("[%2d] --- connecting to %s ---\n", i, names[[i]][1]))
      cl[[i]] <- newPSOCKnode(names[[i]], options = options, rank = i)
    }
        
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}
environment(makePSOCKcluster)<- environment(parallel::makePSOCKcluster)
# environment(newPSOCKnode)<- environment(parallel:::newPSOCKnode)