data("CA_NS6")
d = CA_NS6
head(d)

nptperyear <- 23
INPUT <- check_input(d$t, d$y, d$w, QC_flag = d$QC_flag,
     nptperyear = nptperyear, south = FALSE, 
     maxgap = nptperyear/4, alpha = 0.02, wmin = 0.2)
plot_input(INPUT)
