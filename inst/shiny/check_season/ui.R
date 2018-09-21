# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            # tags$head(
            tags$script(src = "selectInput.js"),
            # selectInput(
            #     "source",
            #     "Choose a data source:",
            #     choices = sites,
            #     selected = "MOD17A2"
            # ),
            titlePanel("Phenofit: growing season dividing"),
            numericInput("iters", "iters:", 2, 1, 3),
            selectInput("FUN_season","Choose a season dividing function (FUN_season):",
                        choices = c('season', 'season_3y'), selected = "season_3y"),
            selectInput("FUN_fit","Choose a Rough fitting function (FUN):",
                        choices = c('wWHIT', 'wSG', 'wHANTS'), selected = "wHANTS"),
            # selectInput("sites", "Choose sites group:",
            #             choices = ,
            #             selected = "single season"
            # ),
            conditionalPanel(
                condition = "input.FUN_fit == 'wWHIT'",
                numericInput("lambda", "lambda:", 1e4, 2, 1e4)
            ),
            conditionalPanel(
                condition = "input.FUN_fit == 'wSG'",
                numericInput("frame", "moving window size (frame):", floor(nptperyear/5*2+1),
                             floor(nptperyear/12), floor(nptperyear/2), floor(nptperyear/12))
            ),
            conditionalPanel(
                condition = "input.FUN_fit == 'wHANTS'",
                numericInput("nf", "number of frequencies (nf):", 3, 1, 6)
            ),
            conditionalPanel(
                condition = "input.FUN_season == 'season_3y'",
                numericInput("maxExtendMonth", "Include n previous and subsequent month (maxExtendMonth):",
                    2, 0, 12)
            ),
            selectInput("wFUN","Choose a weights updating function (wFUN):",
                        choices = c("wTSM", "wBisquare", "wChen"), selected = "wBisquare"),
            selectInput("site","Choose a site:", choices = sites, selected = "US-Me2"), #CH-Fru, FR-LBr
            sliderInput("threshold_max", "threshold_max:",
                min = 0, max = 1, value = 0.2, param_step ),
            sliderInput( "threshold_min", "threshold_min:",
                min = 0, max = 0.2, value = 0, 0.02 ),
            sliderInput( "rytrough_max", "rytrough_max:",
                min = 0, max = 1, value = 0.8, param_step)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            # textOutput('txt_title'),
            tabsetPanel(type = "tabs",
                tabPanel("Plot",
                    plotOutput("plot_GPPobs", height = fig.height),
                    plotOutput("plot_GPP_mod", height = fig.height),
                    plotOutput("plot_GPP_vpm", height = fig.height),
                    plotOutput("plot_MOD13A1_EVI" , height = fig.height),
                    plotOutput("plot_MOD13A1_NDVI", height = fig.height),
                    plotOutput("plot_MOD13Q1_EVI" , height = fig.height),
                    plotOutput("plot_MOD13Q1_NDVI", height = fig.height),
                    plotOutput("plot_MCD15A3H_LAI", height = fig.height)
                ),
                tabPanel("Table",
                    h3("Growing season dividing information:"),
                    br(),
                    DT::dataTableOutput("t_gs") # table of growing season), , width = "100%"
                )
            )
        )
    )
)
