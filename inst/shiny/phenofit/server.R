# shiny::runApp('inst/shiny/phenofit')

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ## define reactiveValues
    rv <- reactiveValues(df = df, st = st, sites = sites)
    # df_rv    <- reactive(rv$df)
    # sites_rv <- reactive(rv$sites)

    # INPUTall   <- reactive({ updateINPUT() })
    ############################################################################
    ################################ observeEvent ##############################
    
    observeEvent(input$txt_varVI , {
        update_VI(rv, input$txt_varVI)
    })

    observeEvent(input$qcFUN     , convert_QC2weight(input))
    observeEvent(input$txt_varQC , convert_QC2weight(input))

    observeEvent(input$nptperyear, {
        # Update wSG parameter
        updateNumericInput(session,
            "frame", "moving window size (frame):",
            floor(nptperyear / 5 * 2 + 1),
            floor(nptperyear / 12),
            floor(nptperyear / 2),
            floor(nptperyear / 12)
        )

        # Update Whittaker parameter
        if (nptperyear >= 300) {
            lambda <- 1e4
        } else if (nptperyear >= 40){
            lambda <- 15
        } else if (nptperyear >= 20){
            lambda <- 5
        } else {
            lambda <- 2
        }
        updateNumericInput(session, "lambda", value = lambda)
    })

    # file_veg changed
    observeEvent(input$file_veg, {
        file_veg  <- input$file_veg$datapath
        file_site <- input$file_site$datapath

        if (check_file(file_veg)) {
            df    <- fread(file_veg) %>% check_datestr()
            sites <- unique(df$site) %>% sort()

            if (!check_file(file_site)){
                rv$st <- data.table(ID = seq_along(sites), site = sites, lat = 30)
            }

            rv$df <- df
            rv$site <- sites

            ## update `nptperyear`
            var_time   <-  intersect(c("t", "date"), colnames(df))[1]
            deltaT     <-  as.numeric(diff(df[[var_time]][c(1, 2)]))
            nptperyear <<- ceiling(365/deltaT)

            updateNumericInput(session, 'nptperyear', value = nptperyear)

            ## update `var_VI``
            updateSelectInput(session, "txt_varVI",
                choices = select_var_VI(rv$df),
                selected = select_var_VI(rv$df)[1])

            # print("here")
            # browser()
            # update_VI(rv, input$txt_varVI) # update_VI right now.

            ## update `var_QC`
            varQCs <- colnames(rv$df) %>% .[grep("QA|QC|qa|qc", .)]
            if (length(varQCs) == 0){
                sel_qc_title <- paste0("vairable of QC: ",
                    "No QC variables!")
                seq_qc <- ""; varQCs <- ""
            } else {
                sel_qc_title <- "vairable of QC:"
                sel_qc <- varQCs[1]
            }

            # Do not update values when varQC changes. Because, it also rely on
            # qcFUN
            updateSelectInput(session, "txt_varQC", sel_qc_title,
                choices = varQCs, varQCs[1])

            ## Update `site` in main panel
            updateSelectInput(session, "site",
                              choices = sites, sites[1])
        }
    })

    observeEvent(input$file_site, {
        file_veg  <- input$file_veg$datapath
        file_site <- input$file_site$datapath

        if (check_file(file_veg)) {
            if (check_file(file_site)){
                rv$st <- fread(file_site)
            }
        }
    })

    ############################################################################
    output$t_input_veg  <- DT::renderDataTable( DT_datatable(rv$df, scrollX = TRUE) )
    output$t_input_site <- DT::renderDataTable( DT_datatable(rv$st))

    d          <- reactive({ getDf.site(rv$df, input$site, input$dateRange) })
    date_range <- reactive({ range(d()$t) })

    INPUT <- reactive({ getINPUT.site(rv$df, rv$st, input$site, input$dateRange) })
    brks  <- reactive({
        cal_season(input, INPUT())
    })

    # Figure height of fine curve fitting
    heightSize <- reactive({
        n <- length(input$FUN)
        height <- 250
        ifelse(n == 1, height+50, height)*n
    }) # number of fine curve fitting

    params_fineFitting <- reactive({
        list(
            INPUT(), brks(),
            methods = input$FUN, #c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
            debug = F,
            wFUN = get(input$wFUN2),
            nextent = 2, maxExtendMonth = 3, minExtendMonth = 1,
            qc = 1,
            minPercValid = 0.2,
            print = TRUE
        )
    })

    # Fine Curve Fitting
    fineFitting <- reactive({
        # params <- params_fineFitting()
        # print(str(params))
        fit  <- do.call(curvefits, params_fineFitting())
        # params_fineFitting <- getparam(fit)
        # print(str(params_fineFitting))

        ## Get GOF information
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.fFIT, .id = "flag")
        }, .id = "meth")

        # get phenology information
        p <- lapply(fit$fits, PhenoExtract)
        pheno <- map(p, tidyFitPheno, origin = INPUT()$t[1]) %>% purrr::transpose()

        c(fit, list(INPUT = INPUT(), seasons = brks(), stat = stat, pheno = pheno))
    })

    lst_metrics <- reactive({
        fineFitting()$pheno %>% map(function(lst){
            d <- melt_list(lst, "meth") %>% reorder_name(c("site", "meth", "flag", "origin"))
            colnames(d) %<>% gsub("GU.|ZHANG.", "", .)
            d
        })
    })
    ############################################################################
    # output$txt_title <- renderText({ INPUT()$titlestr })

    ## Rough fitting and growing season dividing
    output$plot_gs <- renderPlot({
        sitename <- input$site
        par(setting)

        plot_season(INPUT(), brks(), d(), INPUT()$ylu)
        abline(h = 1, col = "red")
        title(INPUT()$titlestr)
        # rv$brks = do.call(check_season, param)
    })

    output$t_gs <- DT::renderDataTable({
        DT_datatable(brks()$dt) %>%
            DT::formatRound(c(4:6), 3) %>%
            DT::formatStyle(columns = c(4:6), 'text-align' = 'center')
    })

    # Fine Curve Fitting Output Figure
    # adjust size of fineFitting output
    output$sized_plot_fineFitting <- renderUI({
        plotOutput("plot_fineFitting", height = heightSize())
    })

    output$plot_fineFitting <- renderPlot({
        g <- plot_phenofit(fineFitting(), d(), INPUT()$titlestr)
        grid::grid.draw(g)
    })

    output$t_fineFitting <- DT::renderDataTable({
        # d_fineParams <- getparam(fineFitting()) %>% melt_list("method")
        d <- fineFitting()$stat
        DT_datatable(d) %>% DT::formatRound(c(3:6), 3) #%>%
            # DT::formatStyle(columns = c(3:6), 'text-align' = 'center')
    })

    output$t_phenoMetrics_date <- DT::renderDataTable({
        d <- lst_metrics()$date
        d <- d[, 3:ncol(d)] %>% lapply(., format.Date, "%Y/%m/%d") %>%
            as.data.frame() %>%
            cbind(d[, 1:2], .)

        DT_datatable(d,
            columnDefs = list(list(width = '20px', targets = 3:ncol(d) )),
            info = FALSE)
    })

    output$t_phenoMetrics_doy <- DT::renderDataTable({
        d <- lst_metrics()$doy
        d[[3]] %<>% format("%Y/%m/%d")
        DT_datatable(d,
            autoWidth = TRUE,
            # scrollX = TRUE, #"900px",
            columnDefs = list(list(width = '20%', targets = 3:ncol(d) ))
        )
    })
    # output$console_phenoMetrics <- renderPrint({ lst_metrics })

    output$download_allsites <- downloadHandler(
        filename = function() {
            paste('data-', Sys.Date(), '.RData', sep='')
        },
        content = function(file) {
            fprintf("debug | n = %d\n", length(sites))
            progress <- shiny::Progress$new(min = 0, max = length(sites))

            # future::future({
                res <- phenofit_all(input, progress)
                # # browser()
                save(res, file = file)
            # })
            # write.csv(mtcars, file)
        }
    )
}

# shinyApp(ui = ui, server = server)
