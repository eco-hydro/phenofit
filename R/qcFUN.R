# source('R/GPP_Pheno/main_QC.R')
# 
# 1. GPPobs  |
# 2. LAI     | FparLai_QC
# ------------------------------------------------------------------------------
# Bits 5-7: SCF_QC (five-level confidence score)
# bitwShiftR(255L, 5)
# 000 (0): Main (RT) method used, best result possible (no saturation)
# 001 (1): Main (RT) method used with saturation. Good, very usable
# 010 (2): Main (RT) method failed due to bad geometry, empirical algorithm used
# 011 (3): Main (RT) method failed due to problems other than geometry, empirical algorithm used
# 100 (4): Pixel not produced at all, value couldn't be retrieved (possible reasons: bad L1B data, unusable MOD09GA data)
# 
# 2. MOD13A1 | SummaryQA
# ------------------------------------------------------------------------------
#    SummaryQA      : Pixel reliability summary QA
#    -1 Fill/No data: Not processed
#    0 Good data    : Use with confidence
#    1 Marginal data: Useful but look at detailed QA for more information
#    2 Snow/ice     : Pixel covered with snow/ice
#    3 Cloudy       : Pixel is cloudy
# 
# 3. MOD13Q1 | SummaryQA, same as MOD13A1
# 4. MODGPP  | Psn_QC
# ------------------------------------------------------------------------------
#    Bits 5, 6, 7: 5-level Confidence Quality score.
#    000 (0): Very best possible
#    001 (1): Good,very usable, but not the best
#    010 (2): Substandard due to geometry problems - use with caution
#    011 (3): Substandard due to other than geometry problems - use with caution
#    100 (4): couldn't retrieve pixel (not produced at all - non-terrestrial biome)
#    111 (7): Fill Value
#    
# 5. NDVIv4  | QA:
# QC bit flags, bit #: description (1 = yes, 0 = no)
# 0 : Unused
# 1 : Pixel is cloudy
# 2 : Pixel contains cloud shadow
# 3 : Pixel is over water
# 4 : Pixel is over sunglint
# 5 : Pixel is over dense dark vegetation
# 6 : Pixel is at night (high solar zenith)
# 7 : Channels 1-5 are valid
# 8 : Channel 1 value is invalid
# 9 : Channel 2 value is invalid
# 10: Channel 3 value is invalid
# 11: Channel 4 value is invalid
# 12: Channel 5 value is invalid
# 13: RHO3 value is invalid
# 14: BRDF correction is invalid
# 15: Polar flag, latitude over 60 degrees (land) or 50 degrees (ocean)

## QC control for MOD09A1
# https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD09A1
# Bits 0-1: Cloud state
# 0: Clear
# 1: Cloudy
# 2: Mixed
# 3: Not set, assumed clear
#
# Bit 2: Cloud shadow
# 0: No
# 1: Yes
#
# Bits 6-7: Aerosol quantity
# 0: Climatology
# 1: Low
# 2: Average
# 3: High
#
# Bit 10: Internal cloud algorithm flag
# 0: No cloud
# 1: Cloud
#
# Bit 12: MOD35 snow/ice flag
# 0: No
# 1: Yes
#
# Bit 15: Internal snow mask
# 0: No snow
# 1: Snow

#' Initial weights according to qc
#' 
#' @description
#' \describe{
#'   \item{getBits}{Extract bitcoded QA information from bin value}
#'   \item{qc_summary}{Initial weigths based on Quality reliability of VI pixel, 
#' suit for MOD13A1, MOD13A2 and MOD13Q1 (SummaryQA band).}
#'   \item{qc_5l}{Initial weights based on Quality control of five-level 
#' confidence score, suit for MCD15A3H(LAI, FparLai_QC), MOD17A2H(GPP, Psn_QC) 
#' and MOD16A2(ET, ET_QC).}
#'   \item{qc_NDVIv4}{For NDVIv4}
#'   \item{qc_StateQA}{Initial weights based on `StateQA`, suit for MOD09A1, MYD09A1. }
#' }
#' 
#' @param x Binary value
#' @param start Bit starting position, count from zero
#' @param end Bit ending position
#' @param wmin Double, minimum weigth (i.e. weight of snow, ice and cloud).
#' @return A list object with
#' \item{weigths}{Double vector, initial weights}
#' \item{QC_flag}{Factor vector, with the level of \code{c("snow", "cloud", "shadow", "aerosol", "marginal", "good")}}
#' 
#' @rdname qcFUN
#' @export
getBits <- function(x, start, end = start){
    # Geometric progression Sn = a1*(1 - q^n)/(1-q)
    n   <- end - start + 1
    a1  <- 2^start; q <- 2
    Sn  <- a1*(1 - q^n)/(1-q)
    
    bitwAnd(x, Sn) %>% bitwShiftR(start) #quickly return
}


#' @param QA quality control variable
#' 
#' @rdname qcFUN
#' @export
qc_summary <- function(QA, wmin = 0.2){
    ## 1. initial weights
    w <- rep(NA, length(QA)) # default weight is zero
    
    w[QA == 0] <- 1             # clear, good
    w[QA == 1] <- 0.5           # margin
    
    w[QA >= 2 & QA <=3] <- wmin # Snow/ice, or cloudy

    ## 2. initial QC_flag
    QC_flag <- factor(QA, 0:3, c("good", "marginal", "snow", "cloud"))

    list(w = w, QC_flag = QC_flag) # quickly return
}

#' @rdname qcFUN
#' @export
qc_StateQA <- function(QA, wmin = 0.2){
    qc_cloud   = getBits(QA, 0, 1)
    qc_shadow  = getBits(QA, 2, 2)
    qc_aerosol = getBits(QA, 6, 7)   # climatology treated as good values
    qc_snow    = getBits(QA, 12, 12)

    ## 1. initial weights
    w <- rep(0.5, length(QA)) # default weight is zero

    I_good <- qc_cloud %in% c(0, 3) & qc_aerosol %in% c(0, 1, 2) & !qc_snow
    I_bad  <- qc_snow | qc_cloud %in% c(1, 2) | qc_aerosol == 3

    # others are marginal
    w[I_good] <- 1
    w[I_bad]  <- wmin

    ## 2. initial QC_flag
    QC_flag = rep("marginal", length(QA))

    I_aerosol <- qc_aerosol == 3
    I_shadow  <- qc_shadow == 1
    I_cloud   <- qc_cloud %in% c(1, 2)
    I_snow    <- qc_snow == 1
    I_good    <- qc_cloud %in% c(0, 3) & qc_aerosol %in% c(0, 1, 2) & !qc_snow

    QC_flag[I_aerosol] <- "aerosol"
    QC_flag[I_shadow]  <- "shadow"
    QC_flag[I_cloud]   <- "cloud"
    QC_flag[I_snow]    <- "snow"
    QC_flag[I_good]    <- "good"

    levels  <- c("snow", "cloud", "shadow", "aerosol", "marginal", "good")
    QC_flag <- factor(QC_flag, levels)

    list(w = w, QC_flag = QC_flag) # quickly return
}

#' @export
#' @rdname qcFUN
qc_5l <- function(QA, wmin = 0.2){
    # bit5-7, five-level confidence score
    # QA <- bitwShiftR(bitwAnd(QA, 224), 5) #1110 0000=224L
    QA <- getBits(QA, 5, 7)
    w  <- rep(NA, length(QA)) #default zero
    
    w[QA <= 1] <- 1            #clear, good
    w[QA >= 2 & QA <=3] <- 0.5 #geometry problems or others
    w[QA >  4] <- wmin
    return(w)
}

#' @rdname qcFUN
qc_NDVIv4 <- function(QA){
    # bit1-2: cloudy, cloud shadow
    QA <- bitwShiftR(bitwAnd(QA, 7), 1) 
    
    w  <- rep(NA, length(QA))
    w[QA == 0] <- 1   #clear, good
    w[QA == 2] <- 0.5 #cloud shadow
    return(w)
}
