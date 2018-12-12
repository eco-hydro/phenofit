# runApp("test/phenology_async/check_season/")
# load("data/shiny_flux115.rda")
# sites <- sort(sites)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ## define reactiveValues
    # INPUTall   <- reactive({ updateINPUT() })
    output$t_input_veg <- DT::renderDataTable({
        DT_datatable(df, scrollX = TRUE)
    })

    output$t_input_site <- DT::renderDataTable({
        DT_datatable(st)    
    })

    d          <- reactive({ getDf.site(df, input$site) })
    date_range <- reactive({ range(d()$t) })
    INPUT <- reactive({ getINPUT.site(df, st, input$site) })
    

    brks  <- reactive({
        param <- list(
            FUN_season     = input$FUN_season, 
            rFUN           = input$rFUN,
            iters          = input$iters, 
            lambda         = input$lambda, 
            nf             = input$nf, 
            frame          = input$frame,
            wFUN           = input$wFUN,
            maxExtendMonth = input$maxExtendMonth, 
            rytrough_max   = input$rytrough_max,
            threshold_max  = input$threshold_max, 
            threshold_min  = input$threshold_min
        )
        # param <- lapply(varnames, function(var) input[[var]])
        param <- c(list(INPUT()), param)
        # print(str(param))
        do.call(check_season, param) # brk return
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
        fit$INPUT   <- INPUT()
        fit$seasons <- brks()

        # params_fineFitting <- getparam(fit)
        # print(str(params_fineFitting))

        ## Get GOF information
        stat  <- ldply(fit$fits, function(fits_meth){
            ldply(fits_meth, statistic.phenofit, .id = "flag")
        }, .id = "meth")
        fit$stat <- stat

        # get phenology information
        p <- lapply(fit$fits, ExtractPheno)
        pheno <- map(p, tidyFitPheno, origin = INPUT()$t[1]) %>% purrr::transpose()

        fit$pheno  <- pheno
        fit
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
            columnDefs = list(list(width = '20px', targets = 3:nrow(d) )),
            info = FALSE)
    })

    output$t_phenoMetrics_doy <- DT::renderDataTable({
        d <- lst_metrics()$doy
        d[[3]] %<>% format("%Y/%m/%d")
        DT_datatable(d,
            autoWidth = TRUE,
            # scrollX = TRUE, #"900px",
            columnDefs = list(list(width = '20%', targets = 3:nrow(d) ))
        )
    })
    # output$console_phenoMetrics <- renderPrint({ lst_metrics })
}

# Run the application
# shinyApp(ui = ui, server = server)
