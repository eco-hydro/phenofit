
findpeaks_season_jl <- function(ypred,
                                r_max = 0, r_min = 0,
                                minpeakdistance = 0, minpeakheight = 0L,
                                nups = 1L, ndowns = nups,
                                # other param
                                nyear = 1) {
  A <- max(ypred) - min(ypred)
  ans <- JuliaCall::julia_call("phenofit.findpeaks_season", ypred,
    r_max = r_max, r_min = r_min,
    # r_max = h_max/A, r_min = h_min/A,
    minpeakdistance = as.integer(minpeakdistance), minpeakheight = minpeakheight,
    nups = nups, ndowns = ndowns
  )
  ans$threshold <- data.table(h_max = r_max * A, h_min = r_min * A, r_max, r_min)
  ans
}
