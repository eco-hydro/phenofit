# Apply V-curve to production data Dongdong Kong

library(readxl)
library(spam)

library(data.table)
library(magrittr)
library(lubridate)

library(plyr)
library(Cairo)
library(purrr)

source('V-pack.R')

file <- "data/fluxsites_LAI_MCD15A3H_0m_buffer.csv"
dt   <- fread(file)
dt$date %<>% ymd()

stations <- fread("data/flux-212.txt")
dt <- merge(dt, stations)

lst <- split(dt, dt$site)
df <- lst$`AR-SLu`

# test stations order
# match(stations$site, names(lst)) %>% diff() %>% unique()


ylab  <- expression(log[10]*' ('*lambda*')')

lambda     <- 1000
nptperyear <- 92*3

CairoPDF(sprintf("optim_lambda=%d, nptperyear=%d.pdf", lambda, nptperyear),
         width = 10, height = 6)

# par(mfrow = c(2, 1), mar = c(4, 3, 2, 1), mgp = c(1.5, 0.6, 0))
op <- par(mfrow = c(2, 2),
          oma = c(1, 2, 2, 1), mar = c(3, 2, 1, 1)) #, yaxt = "n", xaxt = "n"

temp <- llply(lst, fit_site, lambda = lambda, .progress = "text")
dev.off()


temp <- unlist(temp)
median(10^temp, na.rm = T)
boxplot(temp, ylab = ylab, main = 'yearly lambda'); grid()



## 2. test yearly data

# ror]: DE-Akm | 'ylim'值不能是无限的
# |============================                                              |
# 38%[error]: DE-RuS | 'ylim'值不能是无限的
# |==============================================================            |
# 84%[error]: US-ORv | 'ylim'值不能是无限的
# |========================================================================= |
# 99%[error]: US-WPT | 'ylim'值不能是无限的

fit_site <- function(df, lambda){
    site   <- df$site[1]
    nyear  <- floor(nrow(df)/nptperyear)

    df     <- df[1:(nyear*nptperyear), ]
    df$grp <- rep(1:nyear, rep(nptperyear, nyear))

    tryCatch({
        fit <- dlply(df, .(grp), fit_group, lambda = lambda)
        z <- map(fit, "z") %>% unlist
        t <- df$date

        title = sprintf('%s, %s', site, df$IGBP[1])
        plot(t, df$Lai/10, type = "l", main = title); grid();

        lines(t, z, col = "red")
        map_dbl(fit, "opt_lambda")#return
    }, error = function(e){
        message(sprintf("[e] %s|%s", site, e$message))
        return(NA)
    })
}

fit_group <- function(x, lambda){
    grp   <- x$grp[1]
    t     <- x$date
    y_raw <- x$Lai/10
    y <- y_raw
    n <- length(t)

    w <- rep(1, n)
    I_na <- is.na(y)
    y[I_na] <- 0
    w[I_na] <- 0

    lambda <- v_opt(y, w = 0 * y + 1, d = 2, llas = c(0, 4), tol = 0.01)

    list(z = whit(y, lambda, w)$z,
         opt_lambda = lambda) #quickly return
}


optim_lambda <- function(df, IsPlot = T, lims = c(0, 4), ...){
    subfun <- function(df, IsPlot = FALSE){
        tryCatch({
            # Find best lambda with V-curve
            if (IsPlot){
                d = 2
                llas = seq(lims[1], lims[2], by = 0.01)
                vc = v_curve(y, w = w, llas, d = d, show = IsPlot)

                grid()
                # Plot data and smooth
                plot(t, y_raw, type = 'l', col = 'darkgrey', xlab = '', ylab = 'LAI'); grid()
                lines(t, vc$z, col = 'blue', lwd = 1.5)
                str_title <- sprintf("%s group = %d, log10(lambda) = %.3f", site, grp, log10(vc$lambda))
                title(str_title)
            }
            lambda <- v_opt(y, w = 0 * y + 1, d = 2, llas = lims, tol = 0.01)
            return(lambda)
        },
        error = function(e){
            message(sprintf("[error]: %s | %s", site, e$message))
            return(NA);
        })
    }

    # df <- lst$`AR-SLu`
    site   <- df$site[1]
    nyear  <- floor(nrow(df)/nptperyear)

    df     <- df[1:(nyear*nptperyear), ]
    df$grp <- rep(1:nyear, rep(nptperyear, nyear))

    daply(df, .(grp), subfun, IsPlot = IsPlot)
    # subfun(df, IsPlot)
}
# subfun(x, T)

# length | lambda
# -------| ------
# 1y     | 372.4539 , 355.0303
# 2y     | 796.9794 , 762.3611
# 3y     | 928.0381 , 928.0381
# ...
# 15y    | 4679.206
