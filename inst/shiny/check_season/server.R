library(shiny)
library(DT)
# runApp("test/phenology_async/check_season/")
# load("GPPobs_check_season.rda")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    # define reactiveValues
    d     <- reactive({ getSiteData(df, input$site) })
    INPUT <- reactive({
        get_input(df, st, input$site)
    })
    brks  <- reactive({
        varnames <- c("FUN_season", "FUN_fit", "iters", "lambda", "nf", "frame",
            "maxExtendMonth",
            "wFUN", "threshold_max", "threshold_min", "threshold_max") %>%
            intersect(names(input)) %>% set_names(., .)
        param <- lapply(varnames, function(var) input[[var]])
        param <- c(list(INPUT()), param)
        do.call(check_season, param) # brk return
    })
    ############################################################################
    # output$txt_title <- renderText({
    #     INPUT()$titlestr
    # })
    output$obsPlot <- renderPlot({
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
        datatable(brks()$dt, options = list(
            autoWidth = TRUE,
            # columnDefs = list(list(width = '10px', targets = c(4:10)))
            pageLength = 20
        )) %>%
            formatRound(c(4:6), 3) %>%
            formatStyle(columns = c(4:6), 'text-align' = 'center')
    })

    output$gppPlot <- renderPlot({
        sitename <- input$site
        d <- INPUT_lst$MODGPP[site == sitename, .(site, t = date, y = MODGPP, w)]
        plot_data(d, "MOD17A2 GPP")
        abline(h = 1, col = "red")
    })
    output$eviPlot <- renderPlot({
        sitename <- input$site
        d <- INPUT_lst$MOD13A1[site == sitename, .(site, t = date, y = EVI, w)]
        plot_data(d, "MOD13A1 EVI")
    })
    output$ndviPlot <- renderPlot({
        sitename <- input$site
        d <- INPUT_lst$MOD13A1[site == sitename, .(site, t = date, y = NDVI, w)]
        plot_data(d, "MOD13A1 NDVI")
    })
}

# Run the application
# shinyApp(ui = ui, server = server)
