# runApp("test/phenology_async/check_season/")
# load("data/shiny_flux115.rda")

# sites <- sort(sites)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # define reactiveValues
    d          <- reactive({ getSiteData(df, input$site) })
    date_range <- reactive({ range(d()$t) })
    INPUT <- reactive({
        getINPUT_GPPobs(df, st, input$site)
    })
    brks  <- reactive({
        varnames <- c("FUN_season", "FUN_fit",
            "iters", "lambda", "nf", "frame",
            "wFUN",
            "maxExtendMonth", "rytrough_max",
            "threshold_max", "threshold_min", "threshold_max") %>%
            intersect(names(input)) %>% set_names(., .)
        param <- lapply(varnames, function(var) input[[var]])
        param <- c(list(INPUT()), param)
        do.call(check_season, param) # brk return
    })
    ############################################################################
    # output$txt_title <- renderText({
    #     INPUT()$titlestr
    # })
    output$plot_GPPobs <- renderPlot({
        sitename <- input$site
        # threshold_max <- input$threshold_max
        # threshold_min <- input$threshold_min
        # rytrough_max  <- input$rytrough_max
        # iters  <- input$iters
        # lambda <- input$lambda
        # FUNname <- input$FUN
        par(setting)
        plot_season(INPUT(), brks(), d(), INPUT()$ylu)
        abline(h = 1, col = "red")
        title(INPUT()$titlestr)
        # rv$brks = do.call(check_season, param)
    })

    output$t_gs <- DT::renderDataTable({
        site <- input$site
        DT::datatable(brks()$dt, options = list(
            autoWidth = TRUE,
            # columnDefs = list(list(width = '10px', targets = c(4:10)))
            pageLength = 20
        )) %>%
            formatRound(c(4:6), 3) %>%
            formatStyle(columns = c(4:6), 'text-align' = 'center')
    })
    # output$gppPlot <- renderPlot({
    #     sitename <- input$site
    #     d <- INPUT_lst$MODGPP[site == sitename, .(site, t = date, y = MODGPP, w)]
    #     plot_data(d, "MOD17A2 GPP")
    #     abline(h = 1, col = "red")
    # })
    # output$eviPlot <- renderPlot({
    #     sitename <- input$site
    #     d <- INPUT_lst$MOD13A1[site == sitename, .(site, t = date, y = EVI, w)]
    #     plot_data(d, "MOD13A1 EVI")
    # })
    # output$ndviPlot <- renderPlot({
    #     sitename <- input$site
    #     d <- INPUT_lst$MOD13A1[site == sitename, .(site, t = date, y = NDVI, w)]
    #     plot_data(d, "MOD13A1 NDVI")
    # })
    output$plot_GPP_mod <- renderPlot({
        sitename <- input$site
        d <- lst_sm$GPP_mod[site == sitename & scale == "0m", .(site, t = date, y = GPP_mod, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]

        plot_data(d, "MOD17A2 GPP")
        abline(h = 1, col = "red")
    })
    output$plot_GPP_vpm <- renderPlot({
        sitename <- input$site
        d <- lst_sm$GPP_vpm[site == sitename, .(site, t = date, y = GPP_vpm)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]

        plot_data(d, "GPP_vpm (Yao 2017)")
        abline(h = 1, col = "red")
    })
    output$plot_MOD13A1_EVI <- renderPlot({
        sitename <- input$site
        d <- lst_sm$MOD13A1[site == sitename & scale == "0m", .(site, t = date, y = EVI, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]
        plot_data(d, "MOD13A1 EVI")
    })
    output$plot_MOD13A1_NDVI <- renderPlot({
        sitename <- input$site
        d <- lst_sm$MOD13A1[site == sitename & scale == "0m", .(site, t = date, y = NDVI, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]
        plot_data(d, "MOD13A1 NDVI")
    })
    output$plot_MOD13Q1_EVI <- renderPlot({
        sitename <- input$site
        d <- lst_sm$MOD13Q1[site == sitename & scale == "0m", .(site, t = date, y = EVI, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]
        plot_data(d, "MOD13Q1 EVI")
    })
    output$plot_MOD13Q1_NDVI <- renderPlot({
        sitename <- input$site
        d <- lst_sm$MOD13Q1[site == sitename & scale == "0m", .(site, t = date, y = NDVI, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]
        plot_data(d, "MOD13Q1 NDVI")
    })
    output$plot_MCD15A3H_LAI <- renderPlot({
        sitename <- input$site
        d <- lst_sm$LAI[site == sitename & scale == "0m", .(site, t = date, y = LAI, w)]
        d <- d[t >= date_range()[1] & t <= date_range()[2]]
        plot_data(d, "MCD15A3H LAI")
    })
}

# Run the application
# shinyApp(ui = ui, server = server)
