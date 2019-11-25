library(data.table)
library(magrittr)
library(foreach)
library(phenofit)
# devtools::load_all()
# Rcpp::sourceCpp("src/season.cpp")

# {
#     x <- readRDS("x.rds")
#     dt = readRDS("x.rds")
#     print(dt)
#     check_seasons(dt, rtrough_max = 0.6, r_min = 0.1)
#     # check_seasons(dt)
#     print(dt)
#     all.equal(x, dt)
# }

{
    # source("../PhenoAsync/R/tidy_seasons.R")
    tidy_seasons <- function(l, rtrough_max = 0.6, r_min = 0.1) {
        INPUT <- l$INPUT
        brks <- l$brks
        titlestr <- l$titlestr

        year_grps <- names(brks)

        lst_dt <- foreach(dt = brks$dt) %do% {
            dt <- data.frame(dt) %>% data.table()
            check_season(dt, rtrough_max = rtrough_max, r_min = r_min)
            dt <- dt[y_peak != -9999.0 & (len > 45 & len < 650), ]
            dt
        }
        dt <- do.call(rbind, lst_dt)
        check_season(dt, rtrough_max = rtrough_max, r_min = r_min)
        dt <- dt[y_peak != -9999.0 & (len > 45 & len < 650), ]

        brks$dt <- dt
        plot_season(INPUT, brks, title = titlestr, show.legend = FALSE)
        dt
    }

    load("lst_brks.rda")
    # names <- map_chr(lst, "titlestr") %>% str_extract("(?<=th ).{6}(?=,)")
    lst2 = lst[c(15:24, 3:14)]

    file <- "phenofit_fluxnet10_multi_seasons_v3.pdf"
    Cairo::CairoPDF(file, 10, 8)
    par(mfrow = c(5, 1))
    # grps <- c(15:24, 3:14) # CN-Sw2
    grps = seq_along(lst2)
    # grps = 14
    for (i in grps) {
        # Rcpp::sourceCpp("src/season.cpp")
        Ipaper::runningId(i)
        l <- lst2[[i]]
        dt <- tidy_seasons(l, rtrough_max = 0.6, r_min = 0.1)
    }
    dev.off()
    SumatraPDF(file)
}
