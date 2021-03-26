#! /usr/bin/Rscript
# library(data.table)
# library(magrittr)
# library(foreach)
# library(phenofit)
devtools::load_all()
# library(data.table)
# library(magrittr)
# library(testthat)

# read_data <- function(file) {
#     dt = data.table::fread(file)
#     dt$beg  %<>% as.Date()
#     dt$peak %<>% as.Date()
#     dt$end  %<>% as.Date()
#     dt
# }

# test_check_season <- function(dt, len_min = 45, len_max = 650, 
#     verbose = FALSE) 
# {
#     check_season(dt, TRUE, rtrough_max = 0.6, r_min = 0.1)
#     if (verbose) print(dt)

#     dt[y_peak != -9999.0 & (len > len_min & len < len_max), ]
# }

# Rcpp::sourceCpp("src/season.cpp")
# {
#     # 1. debug: 错误高值融合，生长季丢失 (LAI_US-KS2)
#     x = read_data("cheak_season-LAI-(US-KS2).csv")
#     dt = read_data("cheak_season-LAI-(US-KS2).csv")
#     dt = test_check_season(dt)
    
#     expect_equal(nrow(x), nrow(dt))
#     d_diff = dt - x
#     # print(d_diff)
# }
