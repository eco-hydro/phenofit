# Define UI for application that draws a histogram
source('ui_tab_1_loading.R')
# ui <- fluidPage(

ui <- navbarPage(
    "phenofit",
    selected = "Main", #"Load data",
    # tags$head(tags$style(HTML( "hr {border-top: 1px solid black;}" ))),
    tags$script(src = "selectInput.js"),
    tab_loading,
    tabPanel(
        "Main",
        fluidRow(
            column(
                width_sidebar,
                h4("1. Growing season dividing"),
                dateRangeInput("dateRange", "Date range:",
                    start = date_begin,
                    end   = date_end,
                    startview = "decade"),
                selectInput("site", "Choose a site:",
                    choices = sites, selected = sites[1]
                ),
                numericInput("iters", "iters of Rough fitting:", 2, 1, 6),
                checkboxInput("calendarYear", "calendarYear", FALSE),
                selectInput(
                    "FUN_season", "Choose a season dividing function (FUN_season):",
                    choices = c('season', 'season_mov'),
                    selected = "season_mov"
                ),
                selectInput(
                    "wFUN", "Choose a weights updating function for Rough Fitting (wFUN):",
                    choices = options_wFUN,
                    selected = "wBisquare"
                ),
                selectInput(
                    "rFUN", "Choose a Rough fitting function (rFUN):",
                    choices = c('wWHIT', 'wSG', 'wHANTS'),
                    selected = "wHANTS"
                ),
                conditionalPanel(condition = "input.rFUN == 'wWHIT'",
                                 numericInput("lambda", "lambda:", 1e4, 2, 1e4)),
                conditionalPanel(condition = "input.rFUN == 'wSG'",
                    numericInput(
                        "frame", "moving window size (frame):",
                        floor(nptperyear / 5 * 2 + 1),
                        floor(nptperyear / 12),
                        floor(nptperyear / 2),
                        floor(nptperyear / 12)
                    )
                ),
                conditionalPanel(condition = "input.rFUN == 'wHANTS'",
                    numericInput("nf", "number of frequencies (nf):", 2, 1, 6)
                ),
                conditionalPanel(
                    condition = "input.FUN_season == 'season_mov'",
                    numericInput(
                        "maxExtendMonth",
                        "Include n previous and subsequent month in Rough fitting (maxExtendMonth):",
                        2, 0, 12
                    )
                ),
                sliderInput( "r_max", "r_max:",
                    min = 0,  max = 1, value = 0.2, param_step
                ),
                sliderInput( "r_min", "r_min:",
                    min = 0, max = 0.2, value = 0, 0.02
                ),
                sliderInput( "rtrough_max", "rtrough_max:",
                    min = 0, max = 1, value = 0.8, param_step
                ),
                hr(),
                br(),

                ################################################################
                ## curve fitting, select curve fitting methods
                h3("2. Fine Curve fitting"),
                checkboxInput("use.rough", "use.rough", FALSE),
                numericInput("iters2", "iters of Fine fitting:", 2, 1, 6),
                checkboxGroupInput("FUN", 
                    "Choose Fine fitting functions (fFUN):",
                    choiceNames  = list("Asymmetric Gaussian (AG)", 
                        "Beck logistic (Beck)", 
                        "Elmore logistic (Elmore)", 
                        "Gu logistic (Gu)", 
                        "Piecewise logistic (Zhang)"), 
                    choiceValues = list("AG", "Beck", "Elmore", "Gu", "Zhang"),
                    selected = list("Elmore")
                ),
                selectInput(
                    "wFUN2", "Choose a weights updating function for Fine Fitting (wFUN2):",
                    choices = options_wFUN,
                    selected = "wBisquare"
                ),
                numericInput(
                    "nextent2",
                    "Extend curve fitting window, until n good or marginal elements are found in previous and subsequent growing season (nextent).",
                    1, 0, 10
                ),
                numericInput(
                    "maxExtendMonth2",
                    "Max extend window size (month) in Fine fitting (maxExtendMonth):",
                    2, 0, 12
                ),
                hr(),
                br(),
                
                ################################################################
                h3("3. Phenology extraction"),
                checkboxGroupInput("pheno_methods", 
                    "Choose Phenology Extraction methods:",
                    choiceNames  = list(
                        "Threshold (TRS)", 
                        "Derivate (DER)", 
                        "Inflection (Zhang)", 
                        "Gu (Gu)"),
                    choiceValues = list("TRS", "DER", "Zhang", "Gu"),
                    selected = list("TRS", "DER", "Zhang", "Gu")
                ),
                style = "overflow-x: scroll; overflow-y: scroll"
            ),
            column(
                12 - width_sidebar,
                # textOutput('txt_title'),
                # tabsetPanel(type = "tabs",
                # tabPanel("Plot",
                h4("1.1 Rough Fitting"),
                plotOutput("plot_gs", height = fig.height+lgd.height),
                conditionalPanel(condition = "output.warnstat == TRUE",
                           verbatimTextOutput("warnmsg")),

                h4("1.2 Growing season dividing info:"),
                DT::dataTableOutput("t_gs", width = "100%"), # table of growing season), , width = "100%"
            
                h4("2.1 Fine Fitting:"),     
                uiOutput("sized_plot_fineFitting"),
                h4("2.2 Good of fitting:"),     
                DT::dataTableOutput("t_fineFitting", width = "40%"),
               
                h4("3 Phenological metrics (date & doy):"),
                # selectInput(
                #     "t_pheno_methods",
                #     "chose aphenological methods to inspect:",
                #     choices = c("wTSM", "wBisquare", "wChen"),
                #     selected = "wBisquare"
                # ),
                # verbatimTextOutput("console_phenoMetrics")
                column(12,
                    DT::dataTableOutput("t_phenoMetrics_date"),
                    DT::dataTableOutput("t_phenoMetrics_doy"),
                    style = "overflow-x: scroll")
            )
        )),
    navbarMenu("Export",
        tabPanel("Generate code", 
            downloadButton("download_allsites", 
                "Export all sites phenological metrics")
        )
    ),
    tabPanel("help")
)
