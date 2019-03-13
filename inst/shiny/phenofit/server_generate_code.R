# AUTO generating code by `phenofit`
# Dongdong Kong, 20181213

params <- list(
	# 1.1 load input
	txt_varVI, # variable of vegetation index
	file_type, 
	file_veg, 
	file_site, 
	file_rda,

	# 1.2 check_input
	check_QC2weight,
	txt_varQC, 
	qcFUN, 

	# 2 Rough fitting and gs dividing
	FUN_season     = input$FUN_season, 
    rFUN           = input$rFUN,
    iters          = input$iters, 
    lambda         = input$lambda, 
    nf             = input$nf, 
    frame          = input$frame,
    wFUN           = input$wFUN,
    maxExtendMonth = input$maxExtendMonth, 
    rtrough_max   = input$rtrough_max,
    r_max  = input$r_max, 
    r_min  = input$r_min

    # 3 Fine curve fitting
    methods = input$FUN, #c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
    debug = F,
    wFUN2 = get(input$wFUN2),
    nextent = 2, 
    maxExtendMonth = 3, minExtendMonth = 1,
    QC_flag = NULL, # if QC2weight, it's QC_flag
    minPercValid = 0.2, # If valid points less than 20%, will have no phenological
    	                # metrics output.
    print = TRUE
)

# re-define parameter's files, in the principle of easy-editing
# params <- reactiveValuesToList(input)

df <- st <- sites <- nptperyear <- NULL

# all sites input

# updateINPUT()
# updateY(input)
# convert_QC2weight(input)
