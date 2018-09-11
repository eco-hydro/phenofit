# library(rgdal)
library(stars)
library(pracma)
library(Ipaper)

file_y      <- "D:/Document/GoogleDrive/phenofit/whit_lambda/img_ymat.tif"
file_w      <- "D:/Document/GoogleDrive/phenofit/whit_lambda/img_w_2.tif"
file_lambda <- "D:/Document/GoogleDrive/phenofit/whit_lambda/lambda.tif"


# info <- GDALinfo(file)
# x <- readGDAL(file, band = 1)

region.dim = c(1200, 720)

ymat <- read_stars(file_y)[[1]]
wmat <- read_stars(file_w)[[1]]
lambdas <- read_stars(file_lambda)[[1]]

whit <- function (y, lambda = 1600, w, d = 2)
{
    m <- length(y)
    E <- eye(m)
    D <- diff(E, lag = 1, differences = d)
    B <- diag(w) + (lambda * t(D) %*% D)
    z <- solve(B, w*y)
    return(z)
}

dims <- dim(ymat)
nx <- dims[1] # col
ny <- dims[2] # row
prefix <- "whit"

for (i in 1:nx){
    runningId(i, 10)
    # cat(sprintf("%s running %4d ...\n", prefix, i))
    for (j in 1:ny){
        y <- ymat[i, j, ]
        w <- wmat[i, j, ]
        lambda <- lambdas[i, j]

        if (!is.na(lambda)){
            tryCatch({
                z <- whit(y, lambda, w = w)
            }, error = function(e){
                message(sprintf("[i=%d, j=%d]: %s", i, j, e$message))
            })
        }
    }
}
# +proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
