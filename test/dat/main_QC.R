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
qc_summary <- function(QA){
    w <- numeric(length(QA)) #default zero
    
    w[QA == 0] <- 1    #clear, good
    w[QA == 1] <- 0.5  #margin
    
    w[QA >= 2 & QA <=3] <- 0 #cloud shadow
    return(w)
}
#' qc_5l
#' Quality control of five-level confidence score
#' 
#' This function is suit for MCD15A3H(LAI) and MODGPP.
qc_5l <- function(QA){
    # bit5-7, five-level confidence score
    QA <- bitwShiftR(bitwAnd(QA, 224), 5) #1110 0000=224L
    w  <- numeric(length(QA)) #default zero
    
    w[QA <= 1] <- 1            #clear, good
    w[QA >= 2 & QA <=3] <- 0.5 #geometry problems or others
    return(w)
}
qc_NDVIv4 <- function(QA){
    # bit1-2: cloudy, cloud shadow
    QA <- bitwShiftR(bitwAnd(QA, 7), 1) 
    
    w  <- numeric(length(QA))
    w[QA == 0] <- 1   #clear, good
    w[QA == 2] <- 0.5 #cloud shadow
    return(w)
}

# fix_snow
# Roughly assume growing season is Apr-Oct.
fix_snow <- function(x){

}
