packages = c("shiny", "markdown", "shinythemes", "dplyr", "ggplot2", 
             "plotly", "tidyr", "kableExtra", "reshape2", "quantmod",
             "lubridate", "shinyWidgets")
invisible(lapply(packages, library, character.only = TRUE))

shinythemes::themeSelector()
navbarPage("Stock Market App", theme = shinytheme("darkly"),
           tabPanel("Overview",
                    fluidPage(
                      fluidRow(column(4, align="center", textInput(inputId = "stockSymbol", label = "Stock Symbol:", value = "AAPL")),
                               column(4, align="center", dateRangeInput(inputId = "overviewDateRan", label = "Start and End Date:", start = Sys.Date() %m-% months(1), end = Sys.Date())),
                               column(4, align="center", selectInput(inputId = "overviewGraphic", label = "Graph:", choices = c("Candlestick" = "cs", "Time Series" = "ts")))),
                      hr(style = "border-color: #ffffff;"),
                      plotly::plotlyOutput("overviewPlot"),
                      hr(style = "border-color: #ffffff;"),
                      tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
            color: #ffffff !important;}")),
                      DT::dataTableOutput("overviewTable")
                    )),
           tabPanel("Forecast",
                    fluidPage(
                      fluidRow(column(4, align="center", textInput(inputId = "stockSymbolForecast", label = "Stock Symbol:", value = "AAPL")),
                               column(4, align="center", dateRangeInput(inputId = "forecastDateRan", label = "Start and End Date:", start = Sys.Date() %m-% months(1), end = Sys.Date())),
                               column(4, align="center", numericInput(inputId = "numForecast", label = "Number of Forecasts:", value = 1))),
                      fluidRow(column(4, offset = 4, align="center", actionButton(inputId = "forecastButton", label = "GARCH Forecast"))),
                      hr(style = "border-color: #ffffff;"),
                      plotly::plotlyOutput("forecastPlot"),
                      hr(style = "border-color: #ffffff;"),
                      tags$style(HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
                                      color: #ffffff !important;}")),
                      DT::dataTableOutput("forecastTable")
                      )),
           tabPanel("Graphics"),
           tabPanel("Comparison")
)