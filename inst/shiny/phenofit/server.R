# shiny::runApp('inst/shiny/phenofit')
roots <- c('wd' = getwd(), Home = "~", getVolumes()())
nrun  <- 0

#' tidy_shinyFiles
tidy_shinyFiles <- function(selection) {
    paths <- parseFilePaths(roots, selection)$datapath
    # update working dir
    wd <- dirname(paths[1])

    if (dir.exists(wd)) {
        roots['wd'] <<- wd
    }
    return(paths)
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ## file select  ------------------------------------------------------------
    shinyFileChoose(input, 'file_veg_text', roots=roots,
        # defaultPath='',
        defaultRoot='wd',
        filetypes=c('', 'txt', 'csv'))
    shinyFileChoose(input, 'file_site', roots=roots,
        # defaultPath='',
        defaultRoot='wd',
        filetypes=c('', 'txt', 'csv'))
    shinyFileChoose(input, 'file_veg_rda', roots=roots,
        # defaultPath='',
        defaultRoot='wd',
        filetypes=c('', 'RData', 'rda'))

    # Figure height of fine curve fitting
    heightSize <- reactive({
        n <- length(input$FUN_fine)
        if (n == 1) n <- 1.2 # increase height for single plot
        fig.height*n + lgd.height
    }) # number of fine curve fitting

    ## define reactiveValues
    # rv    <- reactiveValues(df = df, st = st, sites = sites)
    rv    <- reactiveValues(df = dataIN$df, st = dataIN$st, sites = dataIN$sites,
        nptperyear = dataIN$nptperyear,
        filepaths = options[c("file_veg_text", "file_veg_rda", "file_site", "file_type", "nptperyear")],
        nrun = 0)

    observeEvent(rv$sites, {
    #     browser()
        updateSelectInput(session, "site",
                      choices = rv$sites, selected = rv$sites[1])
    })
    ############################################################################
    ## observeEvent ------------------------------------------------------------

    # observe({
    #     load_data(options, rv, input)
    #     # load_data(filepaths, rv)
    #     convert_QC2weight(input, rv)
    #     updateSelectInput(session, "site",
    #                   choices = rv$sites, rv$sites[1])
    #     # updateInput_phenofit(session, rv, init = TRUE)
    # })

    observeEvent(input$load_data, {
        rv$filepaths <-
            list(file_veg_text = input$file_veg_text %>% tidy_shinyFiles(),
                file_veg_rda   = input$file_veg_rda %>% tidy_shinyFiles(),
                file_site      = input$file_site %>% tidy_shinyFiles(),
                file_type      = input$file_type,
                nptperyear     = input$nptperyear)

        load_data(rv$filepaths, rv)
    })

    observeEvent(input$pre_process, {
        convert_QC2weight(input, rv)
        update_VI(rv, input$var_y)
    })

    ## INPUT -------------------------------------------------------------------
    ############################################################################
    output$t_input_veg  <- DT::renderDataTable( DT_datatable(rv$df, scrollX = TRUE) )
    output$t_input_site <- DT::renderDataTable( DT_datatable(rv$st))

    d          <- reactive({ getDf.site(rv$df, input$site, input$dateRange) })
    date_range <- reactive({ range(d()$t) })

    INPUT <- reactive({ getINPUT.site(rv$df, rv$st, input$site, rv$nptperyear,
        input$dateRange) })
    brks  <- reactive({ cal_season(input, INPUT()) })

    params_fineFitting <- reactive({
        list(
            INPUT(), brks(),
            iters          = input$iters_fine,
            methods        = input$FUN_fine, #c("AG", "zhang", "beck", "elmore", 'Gu'), #,"klos",
            verbose        = FALSE,
            wFUN           = get(input$wFUN_fine),
            nextend        = input$nextend_fine,
            maxExtendMonth = input$max_extend_month_fine,
            minExtendMonth = 1,
            minPercValid   = 0.2,
            print          = TRUE,
            use.rough      = input$use.rough
        )
    })

    ## WRITE options to json ---------------------------------------------------
    # add option to write setting
    options_phenofit <- reactive({
        # browser()
        # nrun <- isolate(rv$nrun)
        # print(filepaths() %>% str()) # who induced?
        # print(rv$filepaths)
        setting.get(input, rv$filepaths)
    })

    observeEvent(options_phenofit(), {
        if (rv$nrun >= 1) {
            timestr <- Sys.time() %>% format("%Y%m%d_%H0000")
            outfile <- sprintf("perference/phenofit_setting.json")
            # outfile <- sprintf("perference/phenofit_setting_%s.json", timestr)

            sprintf('[s] write setting %s ... \n', basename(outfile)) %>% cat()
            setting.write(options_phenofit(), outfile)
        }
    })

    ## Fine Curve Fitting ------------------------------------------------------
    fineFitting <- reactive({
        # params <- params_fineFitting()
        # print(str(params))
        fit  <- do.call(curvefits, params_fineFitting())
        # params_fineFitting <- getparam(fit)

        stat  <- get_GOF(fit)                       # Goodness-Of-Fit
        pheno <- get_pheno(fit, IsPlot=FALSE)   # Phenological metrics

        list(fit = fit, INPUT = INPUT(), seasons = brks(), stat = stat, pheno = pheno)
    })

    lst_metrics <- reactive({
        isolate({
            rv$nrun <- rv$nrun + 1 # increase running times
        })
        nrun <<- nrun + 1

        fineFitting()$pheno %>% map(function(lst){
            d <- melt_list(lst, "meth") %>% reorder_name(c("site", "meth", "flag", "origin"))
            colnames(d) %<>% gsub("GU.|ZHANG.", "", .)
            d
        })
    })

    ###################### END OF REACTIVE VALUES ##############################
    ############################################################################
    # output$txt_title <- renderText({ INPUT()$titlestr })

    ## Rough fitting and growing season dividing
    output$plot_gs <- renderPlot({
        sitename <- input$site
        par(par_setting)

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
        df_fit <- get_fitting(fineFitting()$fit)

        g <- plot_phenofit(df_fit, brks(), INPUT()$titlestr)
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
