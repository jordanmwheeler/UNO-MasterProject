packages = c("shiny", "markdown", "shinythemes", "dplyr", "ggplot2", 
             "plotly", "tidyr", "kableExtra", "reshape2", "quantmod",
             "lubridate", "shinyWidgets", "Hmisc", "fGarch", "parallel",
             "shinycssloaders")
invisible(lapply(packages, library, character.only = TRUE))


shinythemes::themeSelector()

tags$style(type="text/css",
           ".shiny-output-error { visibility: hidden; }",
           ".shiny-output-error:before { visibility: hidden; }"
)

navbarPage("Praedicere", theme = shinytheme("darkly"),
           navbarMenu("Stock Overview",
                      tabPanel("Daily",
                               fluidPage(
                                 fluidRow(column(4, align="center", textInput(inputId = "stockSymbol", label = "Stock Symbol:", value = "AAPL")),
                                          column(4, align="center", dateRangeInput(inputId = "overviewDateRan", label = "Start and End Date:", start = Sys.Date() %m-% months(1), end = Sys.Date())),
                                          column(4, align="center", selectInput(inputId = "overviewGraphic", label = "Graph:", choices = c("Candlestick" = "cs", "Time Series" = "ts")))),
                                 hr(style = "border-color: #ffffff;"),
                                 tags$br(),
                                 withSpinner(plotly::plotlyOutput("overviewPlot")),
                                 tags$br(),
                                 tags$br(),
                                 tags$br(),
                                 hr(style = "border-color: #ffffff;"),
                                 tags$br(),
                                 tags$br(),
                                 tags$br(),
                                 tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
                                                 color: #ffffff !important;}")),
                                 withSpinner(DT::dataTableOutput("overviewTable")),
                                 tags$br(),
                                 tags$br(),
                                 tags$br()
                                 )),
                      "----",
                      tabPanel("Intraday",
                               fluidPage(
                                 fluidRow(align = "center", column(4, align="center", textInput(inputId = "intra.ov.stockSymbol", label = "Stock Symbol:", value = "AAPL")),
                                          column(4, align="center", dateRangeInput(inputId = "intra.ov.overviewDateRan", label = "Start and End Date:", start = Sys.Date(), end = Sys.Date()), selectInput("intra.ov.intraday.step", label = "Choose Intraday Step:", choices = c("1 min." = "1min", "5 min." = "5min", "15 min." = "15min", "30 min." = "30min", "60 min." = "60min"))),
                                          column(4, align="center", selectInput(inputId = "intra.ov.overviewGraphic", label = "Graph:", choices = c("Candlestick" = "cs", "Time Series" = "ts")))),
                                 hr(style = "border-color: #ffffff;"),
                                 withSpinner(plotly::plotlyOutput("intra.ov.overviewPlot")),
                                 hr(style = "border-color: #ffffff;"),
                                 tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
                                                 color: #ffffff !important;}")),
                                 DT::dataTableOutput("intra.ov.overviewTable")
                                 ))),
           tabPanel("Daily Forecast",
                    fluidPage(
                      fluidRow(column(4, align="center", textInput(inputId = "stockSymbolForecast", label = "Stock Symbol:", value = "AAPL")),
                               column(4, align="center", dateRangeInput(inputId = "forecastDateRan", label = "Start and End Date:", start = Sys.Date() %m-% months(1), end = Sys.Date())),
                               column(4, align="center", numericInput(inputId = "numForecast", label = "Number of Forecasts:", value = 1))),
                      fluidRow(column(4, align="center", actionButton(inputId = "forecastButton", label = "Auto GARCH", width = '75%')),
                               column(4, align="center", actionButton(inputId = "model.performance", label = "Model Performance", width = '75%')),
                               column(4, align="center", actionButton(inputId = "man.forecast.button", label = "Manual GARCH", width = '75%'))),
                      hr(style = "border-color: #ffffff;"),
                      withSpinner(plotly::plotlyOutput("forecastPlot")),
                      hr(style = "border-color: #ffffff;"),
                      fluidRow(column(8, offset = 2, align="center", selectInput(inputId = "forecastDisplay.button", label = "Display", choices = c("Actuals and Forecast", "Only Forecast")))),
                      fluidRow(column(4, offset = 4, align="center", downloadButton(outputId = "forecast.download", label = "Download Forecast .csv"))),
                      hr(style = "border-color: #ffffff;"),
                      tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
                                      color: #ffffff !important;}")),
                      withSpinner(DT::dataTableOutput("forecastTable")),
                      hr(style = "border-color: #ffffff;")
                      )),
           tabPanel("Intraday Forecast",
                    fluidPage(
                      fluidRow(column(4, align="center", textInput(inputId = "intra.stockSymbolForecast", label = "Stock Symbol:", value = "AAPL")),
                               column(4, align="center", dateRangeInput(inputId = "intra.forecastDateRan", label = "Start and End Date:", start = Sys.Date(), end = Sys.Date()), selectInput("intra.intraday.step", label = "Choose Intraday Step:", choices = c("1 min." = "1min", "5 min." = "5min", "15 min." = "15min", "30 min." = "30min", "60 min." = "60min"), selected = "5 min.")),
                               column(4, align="center", numericInput(inputId = "intra.numForecast", label = "Number of Forecasts:", value = 1))),
                      fluidRow(column(4, align="center", actionButton(inputId = "intra.forecastButton", label = "Auto GARCH", width = '75%')),
                               column(4, align="center", actionButton(inputId = "intra.model.performance", label = "Model Performance", width = '75%')),
                               column(4, align="center", actionButton(inputId = "intra.man.forecast.button", label = "Manual GARCH", width = '75%'))),
                      hr(style = "border-color: #ffffff;"),
                      withSpinner(plotly::plotlyOutput("intra.forecastPlot")),
                      hr(style = "border-color: #ffffff;"),
                      fluidRow(column(8, offset = 2, align="center", selectInput(inputId = "intra.forecastDisplay.button", label = "Display", choices = c("Actuals and Forecast", "Only Forecast")))),
                      fluidRow(column(4, offset = 4, align="center", downloadButton(outputId = "intra.forecast.download", label = "Download Forecast .csv"))),
                      hr(style = "border-color: #ffffff;"),
                      tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
                                      color: #ffffff !important;}")),
                      withSpinner(DT::dataTableOutput("intra.forecastTable")),
                      hr(style = "border-color: #ffffff;")
                    ))
)
