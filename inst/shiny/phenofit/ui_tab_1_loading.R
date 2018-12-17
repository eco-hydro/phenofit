## 1. loading data table
# source('ui_tab_1_loading.R')
width_sidebar <- 3

tab_loading <- tabPanel("Load data",
    # 1.1 loading data
    fluidRow(
        column(width_sidebar + 2, 
            h3('1.1 load data'),         
            radioButtons("file_type", "file type:", 
                choices = c("text", ".rda | .RData"), 
                selected = "text"),
            conditionalPanel(condition = "input.file_type == 'text'",
                fileInput("file_veg", "File of vegetation time-series (file_veg):",
                    multiple = FALSE,
                    accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
                fileInput("file_site", "File of site information (file_site):",
                    multiple = FALSE,
                    accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
            ),
            conditionalPanel(condition = "input.file_type == '.rda | .RData'",
                fileInput("file_rda", 
                    "RData of vegetation time-series and site information (file_rda):",
                    multiple = FALSE,
                    accept = c(".rda", ".RData"))
            ),
            numericInput("nptperyear", "nptperyear:", 365, 12, 366, 1),

            ## 1.2 check_input
            h3('1.2 check_input'),
            textInput("txt_varVI", "vairable of vegetation index", "y"),

            checkboxInput("check_QC2weight", "Convert QC to weight?", FALSE),
            conditionalPanel(condition = "input.check_QC2weight", 
                textInput("txt_varQC", "vairable of QC:", ""),
                selectInput(
                    "qcFUN", "function of initializing weights according to QC (qcFUN):",
                    choices = c("qc_summary", "qc_5l", "qc_StateQA", "qc_NDVIv4"),
                    selected = "qc_summary"
                )
            ),
            actionButton("btn_updateInput", strong("Update INPUT"))
        ),
        column(6, 
            # verbatimTextOutput("console_phenoMetrics", "help info")
            br(), br(),
            h3(em("If no input data assigned, the default is Eddy covariance daily GPP data.", 
                style="color:red")),
            br(), br(), br(), 
            conditionalPanel(condition = "input.file_type == 'text'", 
                strong("File of vegetation time-series:"),
                p(code('file_veg'), " should have the column of 'site', 'y', 't', and 'w' (optional).", tags$br(),  
                    "If w is missing, weights of all points are 1.0."),
                br(),
                strong("File of vegetation time-series:"),
                p(code('file_site'), " should have the column of ", 
                    "'ID (numeric)', 'site (string)', 'lat (double)'. ", tags$br(), 
                    "IGBPname (string) is optional.")
            ),
            conditionalPanel(condition = "input.file_type == '.rda | .RData'", 
                strong("RData of vegetation time-series and site information:"),
                p(code('file_rda'), "should have the variable of ", 
                    code("df"), " (data.frame of vegetation time-series) and ",
                    code("st"), " (data.frame of site information)."), 
                br(),
                p(code('df'), " should have the column of 'site', 'y', 't', and 'w' (optional).", tags$br(),  
                    "If w is missing, weights of all points are 1.0."),
                p(code('st'), " should have the column of ", 
                    "'ID (numeric)', 'site (string)', 'lat (double)'. ", tags$br(), 
                    "IGBPname (string) is optional.")
            )
        )
    ), 

    br(), 
    fluidRow(
        # 2.3 preview input data
        hr(),
        h3("1.1 Vegetation time-series:"),
        DT::dataTableOutput("t_input_veg", width = "50%"), 
        
        h3("1.2 Site information:"),
        DT::dataTableOutput("t_input_site")
    )
)
