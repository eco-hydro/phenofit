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
                    # choices = ""
                    choices = sites, selected = options$site
                ),
                numericInput("iters_rough", "iters of Rough fitting:", 2, 1, 6),
                checkboxInput("calendarYear", "calendarYear", options$calendarYear),
                selectInput(
                    "FUN_season", "Choose a season dividing function (FUN_season):",
                    choices = c('season', 'season_mov'),
                    selected = options$FUN_season
                ),
                selectInput(
                    "FUN_rough", "Choose a Rough fitting function (FUN_rough):",
                    choices = c('wWHIT', 'wSG', 'wHANTS'),
                    selected = options$FUN_rough
                ),
                selectInput(
                    "wFUN_rough", "Choose a weights updating function for Rough Fitting (wFUN_rough):",
                    choices = options_wFUN,
                    selected = options$wFUN_rough,
                ),
                conditionalPanel(condition = "input.FUN_rough == 'wWHIT'",
                                 numericInput("lambda", "lambda:", options$lambda, 2, 1e4)),
                conditionalPanel(condition = "input.FUN_rough == 'wSG'",
                    numericInput(
                        "frame", "moving window size (frame):",
                        options$frame,
                        floor(nptperyear / 12),
                        floor(nptperyear / 2),
                        floor(nptperyear / 12)
                    )
                ),
                conditionalPanel(condition = "input.FUN_rough == 'wHANTS'",
                    numericInput("nf", "number of frequencies (nf):", options$nf, 1, 6)
                ),
                conditionalPanel(
                    condition = "input.FUN_season == 'season_mov'",
                    numericInput(
                        "max_extend_month_rough",
                        "Include n previous and subsequent month in Rough fitting (max_extend_month_rough):",
                        options$max_extend_month_rough, 0, 12
                    )
                ),
                sliderInput( "r_max", "r_max:",
                    min = 0,  max = 1, value = options$r_max, param_step
                ),
                sliderInput( "r_min", "r_min:",
                    min = 0, max = 0.2, value = options$r_min, 0.02
                ),
                sliderInput( "rtrough_max", "rtrough_max:",
                    min = 0, max = 1, value = options$rtrough_max, param_step
                ),
                hr(),
                br(),

                ################################################################
                ## curve fitting, select curve fitting methods
                h3("2. Fine Curve fitting"),
                checkboxInput("use.rough", "Use rough fitting smoothed series?", options$use.rough),
                numericInput("iters_fine", "iters of Fine fitting:", options$iters_fine, 1, 6),
                checkboxGroupInput("FUN_fine", 
                    "Choose Fine fitting functions (FUN_fine):",
                    choiceNames  = list("Asymmetric Gaussian (AG)", 
                        "Beck logistic (Beck)", 
                        "Elmore logistic (Elmore)", 
                        "Gu logistic (Gu)", 
                        "Piecewise logistic (Zhang)"), 
                    choiceValues = list("AG", "Beck", "Elmore", "Gu", "Zhang"),
                    selected = options$FUN_fine
                ),
                selectInput(
                    "wFUN_fine", "Choose a weights updating function for Fine Fitting (wFUN_fine):",
                    choices = options_wFUN,
                    selected = options$wFUN_fine
                ),
                numericInput(
                    "nextend_fine",
                    "Extend curve fitting window, until n good or marginal elements are found in previous and subsequent growing season (nextend_fine).",
                    options$nextend_fine, 0, 10
                ),
                numericInput(
                    "max_extend_month_fine",
                    "Max extend window size (month) in Fine fitting (max_extend_month_fine):",
                    options$max_extend_month_fine, 
                    0, 12
                ),
                hr(),
                br(),
                
                ################################################################
                h3("3. Phenology extraction"),
                checkboxGroupInput("meths_pheno", 
                    "Choose Phenology Extraction methods:",
                    choiceNames  = list(
                        "Threshold (TRS)", 
                        "Derivate (DER)", 
                        "Inflection (Zhang)", 
                        "Gu (Gu)"),
                    choiceValues = list("TRS", "DER", "Zhang", "Gu"),
                    selected = options$meths_pheno
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
                h4("2.2 Good of fitting (NSE):"),     
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
                    # DT::dataTableOutput("t_phenoMetrics_date"),
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
