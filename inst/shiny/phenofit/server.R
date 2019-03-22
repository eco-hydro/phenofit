# shiny::runApp('inst/shiny/phenofit')
roots <- c('wd' = getwd(), Home = "~", getVolumes()())
nrun  <- 0

if(.Platform$OS.type == "windows") {
    windowsFonts(
      # lishu = windowsFont(family = "LiSu"),           # 隶书
      yahei = windowsFont(family = "Microsoft YaHei") # 微软雅黑
    )
}


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
        updateSelectInput(session, "site",
                      choices = rv$sites, selected = rv$sites[1])
    })
    ############################################################################
    ## observeEvent ------------------------------------------------------------

    observeEvent(input$load_data, {
        filepaths <-
            list(file_veg_text = input$file_veg_text %>% tidy_shinyFiles(),
                file_veg_rda   = input$file_veg_rda %>% tidy_shinyFiles(),
                file_site      = input$file_site %>% tidy_shinyFiles(),
                file_type      = input$file_type,
                nptperyear     = input$nptperyear)

        phenofit_loaddata(filepaths, rv)
        rv$filepaths <- filepaths
    })

    observeEvent(input$pre_process, {
        convert_QC2weight(input, rv)
        update_VI(rv, input$var_y)
    })

    ## INPUT -------------------------------------------------------------------
    ############################################################################
    output$t_input_veg  <- DT::renderDataTable( DT_datatable(rv$df, scrollX = TRUE) )
    output$t_input_site <- DT::renderDataTable( DT_datatable(rv$st))

    d          <- reactive({ getsite_data(rv$df, input$site, input$dateRange) })
    date_range <- reactive({ range(d()$t) })

    INPUT <- reactive({
        getsite_INPUT(rv$df, rv$st, input$site, rv$nptperyear, input$dateRange)
    })
    brks  <- reactive({ phenofit_season(INPUT(), input) })

    ## WRITE options to json ---------------------------------------------------
    # add option to write setting
    options_phenofit <- reactive({
        # nrun <- isolate(rv$nrun)
        # print(filepaths() %>% str()) # who induced?
        setting.get(input, rv$filepaths)
    })

    observeEvent(options_phenofit(), {
        # if (rv$nrun >= 1) {
        fprintf('nrun = %d\n', nrun)
        if (nrun >= 1) {
            timestr <- Sys.time() %>% format("%Y%m%d_%H0000")
            file_json <- sprintf("perference/phenofit_setting.json")
            # outfile <- sprintf("perference/phenofit_setting_%s.json", timestr)
            if (!dir.exists(dirname(file_json)))
                dir.create(dirname(file_json), recursive = TRUE)

            sprintf('[s] write setting %s ... \n', basename(file_json)) %>% cat()
            setting.write(options_phenofit(), file_json)
        }
    })

    ## Fine Curve Fitting ------------------------------------------------------
    fineFitting <- reactive({
        phenofit_finefit(INPUT(), brks(), input)
    })

    lst_metrics <- reactive({
        # isolate({ rv$nrun <- rv$nrun + 1 }) # increase running times
        nrun <<- nrun + 1
        fineFitting()$pheno %>% map(function(lst){
            d <- melt_list(lst, "meth") %>% reorder_name(c("site", "meth", "flag", "origin"))
            colnames(d) %<>% gsub("GU.|ZHANG.", "", .)
            d
        })
    })

    ###################### END OF REACTIVE VALUES ##############################
    ############################################################################

    ## Rough fitting and growing season dividing
    output$plot_gs <- renderPlot({
        sitename <- input$site
        par(par_setting)

        plot_season(INPUT(), brks(), d(), INPUT()$ylu)
        abline(h = 1, col = "red")
        title(INPUT()$titlestr, family = "yahei")
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
        d <- dcast(d, flag~meth, value.var = "NSE")
        # browser()
        DT_datatable(d) %>% DT::formatRound(c(2:ncol(d)), 3) #%>%
            # DT::formatStyle(columns = c(3:6), 'text-align' = 'center')
    })

    # output$t_phenoMetrics_date <- DT::renderDataTable({
    #     d <- lst_metrics()$date
    #     d <- d[, 3:ncol(d)] %>% lapply(., format.Date, "%Y/%m/%d") %>%
    #         as.data.frame() %>%
    #         cbind(d[, 1:2], .)

    #     DT_datatable(d,
    #         columnDefs = list(list(width = '20px', targets = 3:ncol(d) )),
    #         info = FALSE)
    # })

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
        filename = function() { paste('data-', Sys.Date(), '.RData', sep='') },
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
