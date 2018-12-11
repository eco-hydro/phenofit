# Define UI for application that draws a histogram
width_sidebar <- 3
# ui <- fluidPage(

# sidebar
# Application title

ui <- navbarPage(
    "phenofit",
    selected = "Season dividing",
    tags$head(tags$style(HTML(
        "hr {border-top: 1px solid black;}"
    ))),
    tags$script(src = "selectInput.js"),

    tabPanel("Load data", 
        radioButtons("file_type", "file type:", 
            choices = c("text", ".rda | .RData"), 
            selected = "text"),
        conditionalPanel(condition = "input.file_type == 'text'",
            fileInput("file_veg", "File of vegetation time-series",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
            fileInput("file_site", "File of site information",
                multiple = FALSE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
        ),
        conditionalPanel(condition = "input.file_type == '.rda | .RData'",
            fileInput("file_rda", "INPUT file",
                multiple = FALSE,
                accept = c(".rda", ".RData"))
        ),
        numericInput("nptperyear", "nptperyear:", 1, 366, 23)
        # Input: Select a file ----     
    ),
    tabPanel(
        "Season dividing",
        fluidRow(
            column(
                width_sidebar,
                h4("1. Growing season dividing"),
                selectInput(
                    "site", "Choose a site:",
                    choices = sites, selected = "US-Me2"
                ),
                numericInput("iters", "iters:", 2, 1, 3),
                selectInput(
                    "FUN_season", "Choose a season dividing function (FUN_season):",
                    choices = c('season', 'season_3y'),
                    selected = "season_3y"
                ),
                selectInput(
                    "wFUN", "Choose a weights updating function for Rough Fitting (wFUN):",
                    choices = c("wTSM", "wBisquare", "wChen"),
                    selected = "wBisquare"
                ),
                selectInput(
                    "rFUN", "Choose a Rough fitting function (FUN):",
                    choices = c('wWHIT', 'wSG', 'wHANTS'),
                    selected = "wHANTS"
                ),
                conditionalPanel(condition = "input.FUN_fit == 'wWHIT'",
                                 numericInput("lambda", "lambda:", 1e4, 2, 1e4)),
                conditionalPanel(
                    condition = "input.FUN_fit == 'wSG'",
                    numericInput(
                        "frame",
                        "moving window size (frame):",
                        floor(nptperyear / 5 * 2 + 1),
                        floor(nptperyear / 12),
                        floor(nptperyear / 2),
                        floor(nptperyear / 12)
                    )
                ),
                conditionalPanel(
                    condition = "input.rFUN == 'wHANTS'",
                    numericInput("nf", "number of frequencies (nf):", 3, 1, 6)
                ),
                conditionalPanel(
                    condition = "input.FUN_season == 'season_3y'",
                    numericInput(
                        "maxExtendMonth",
                        "Include n previous and subsequent month (maxExtendMonth):",
                        2, 0, 12
                    )
                ),
                sliderInput(
                    "threshold_max", "r_max:",
                    min = 0,  max = 1, value = 0.2,
                    param_step
                ),
                sliderInput(
                    "threshold_min", "r_min:",
                    min = 0, max = 0.2, value = 0, 0.02
                ),
                sliderInput(
                    "rytrough_max", "rytrough_max:",
                    min = 0, max = 1, value = 0.8, param_step
                ),
                style = "overflow-x: scroll; overflow-y: scroll"
            ),
            column(
                12 - width_sidebar,
                # textOutput('txt_title'),
                # tabsetPanel(type = "tabs",
                # tabPanel("Plot",
                plotOutput("plot_gs", height = fig.height),
                h4("Growing season dividing information:"),
                br(),
                DT::dataTableOutput("t_gs") # table of growing season), , width = "100%"
                # )
                # )
            )
        )),
    ## curve fitting, select curve fitting methods
    tabPanel(
        "Curve fitting & Phenology Extraction",
        fluidRow(
            column(width_sidebar, 
                h3("3. Fine Curve fitting"),
                checkboxGroupInput("FUN", 
                    "Choose Fine fitting functions (fFUN):",
                    choiceNames  = list("Asymmetric Gaussian (AG)", 
                        "Piecewise logistic (zhang)", 
                        "Beck logistic (beck)", 
                        "Elmore logistic (elmore)", 
                        "Gu logistic (Gu)"), 
                    choiceValues = list("AG", "zhang", "beck", "elmore", "Gu"),
                    selected = list("elmore")
                ),
                selectInput(
                    "wFUN2",
                    "Choose a weights updating function for Fine Fitting (wFUN2):",
                    choices = c("wTSM", "wBisquare", "wChen"),
                    selected = "wBisquare"
                ),
                ################################################################
                h3("4. Phenology extraction"),
                checkboxGroupInput("pheno_methods", 
                    "Choose Phenology Extraction methods:",
                    choiceNames  = list(
                        "Threshold (TRS)", 
                        "Derivate (DER)", 
                        "Inflection (Zhang)", 
                        "Gu (Gu)"),
                    choiceValues = list("TRS", "DER", "Zhang", "Gu"),
                    selected = list("TRS", "DER", "Zhang", "Gu")
                )
            ),
            column(12 - width_sidebar,
                verticalLayout(
                    uiOutput("sized_plot_fineFitting"),
                    br(),
                    h3("3.1 Good of fitting:"),
                    DT::dataTableOutput("t_fineFitting", width = "40%"),
                    br(),

                    h3("4.1 Phenological metrics (date & doy):"),
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
            )
        )
    ),
    navbarMenu("Export",
        tabPanel("generate script"),
        tabPanel("export phenological metrics")),
    tabPanel("help")
)
# titlePanel("Curve Fitting")
