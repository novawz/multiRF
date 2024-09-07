library(shiny)
library(shinythemes)
library(shinycssloaders)

ui <- fluidPage(
  theme = shinytheme("yeti"),  # Using a theme close to our desired style
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    # Import Font Awesome for icons
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.3/css/all.min.css")
  ),
  navbarPage(
    title = div(style = "color: white; font-weight: bold;", "multiRF"),
    id = "navbar",
    tabPanel(
      HTML('<i class="fas fa-info-circle"></i> About'),  # Using Font Awesome icon for "About"
      fluidPage(
        h2("About"),
        p("This Shiny app demonstrates the use of the shinydashboard package with three sections: About, Variable Selection, and Clustering Analysis.")
      )
    ),
    tabPanel("Variable Selection",
             fluidPage(
               h2("Variable Selection"),
                sidebarLayout(sidebarPanel(
                  h3("Step 1: Input Data"),
                  fileInput("file1", "Upload your own data", accept = ".csv", multiple = T),
                  actionButton("submit", "Submit", icon = NULL, width = NULL),
                  helpText("Note: Upload the .csv file"),
                  h5("Use Example Data"),
                  actionButton("example", "TCGA-BRCA", icon = NULL, width = NULL, disabled = FALSE),
                  h3("Step 2: Initial Model Fitting"),
                  numericInput("ntree", "Number of trees", 300, min = 0, width = "50%"),
                  numericInput("yprob", "Response proportion", 0.5, min = 0.1, max = 1, width = "50%"),
                  checkboxGroupInput("para", "Other parameters", choices = c("Scale", "Normalized weights", "Direct"),
                                     selected = "Scale", inline = T),
                  actionButton("run", "Run", icon = NULL, width = NULL),
                  helpText("This could take a few minites"),
                  h3("Step 3: Variable Selection"),
                  selectInput("method", "Method", c("test", "filter", "mixture"),  width = "50%"),
                  uiOutput("newSelect1"),
                  uiOutput("newSelect2"),
                  actionButton("run2", "Run", icon = NULL, width = NULL),
                  width = 3),
               mainPanel( width = 9,
                 tabsetPanel(
                   tabPanel("Data Summary",
                            br(),
                            tableOutput("datasummary"),
                            plotOutput("dataplot1"), align = "center",
                            br(),
                            span(textOutput("initialModel") %>% withSpinner(), style="color:navy")  ),
                   tabPanel("Results",
                            br(),
                            h4("Variable weights"),
                            br(),
                            plotOutput("scatterPlot1") %>% withSpinner(),
                            br(),
                            h4("Selection summary"),
                            br(),
                            column(4, tableOutput("table2") ),
                            column(8, plotOutput("scatterPlot2") )
                            ),
                   tabPanel("Visualization",
                            br(),
                            h4("Heatmap"),
                            br(),
                            plotOutput("dataplot2"),
                            br(),
                            h4("Pairwise importance"),
                            br(),
                            plotOutput("dataplot3")

                   )
                 )
                )))),
    tabPanel("Clustering Analysis",
             fluidPage(
               h2("Clustering Analysis"),
               sidebarLayout(
                 sidebarPanel(
                   selectInput("xvar_cluster", "X-axis variable", names(mtcars)),
                   selectInput("yvar_cluster", "Y-axis variable", names(mtcars)),
                   numericInput("clusters", "Number of clusters", 3, min = 1, max = 10)
                 ),
                 mainPanel(
                   plotOutput("clusterPlot")
                 )
               )
             )
    )
  )
)
