library(shiny)
library(Seurat)
library(SeuratData)
library(shinyjs)

# Define UI
ui <- fluidPage(
  useShinyjs(),
  headerPanel('pearsonCorrelateCells'),
  div(
    tags$p("pearsonCorrelateCells does [insert].
         This Shiny app demonstrates the capabilities of the
         pearsonCorrelateCells package. The posterior slide 10X Genomics
         Visium Sagittal Mouse Brain is used here. You may upload a csv file
           with genes and corresponding expression. You will then see
           whether each spot in the Mouse Brain spatial data has a significant
           positive correlation with the genes expression in the file you
           uploaded."),
    style="margin-left: 15px"
  ),
  br(),
  sidebarPanel(
    tags$p("You may upload a csv file
           with two columns, the first column containing gene names and the
           second column containing gene expression values. Here is an example file"),
    downloadButton("downloadData", "Download example csv"),
    br(),
    br(),
    fileInput(inputId="uploadedFile", label="Choose gene expression CSV file",
          accept = c("text/csv",
                     "text/comma-separated-values,text/plain",
                     ".csv"))
  ),
  mainPanel(
    div(
      id = "loadingIndicator",
      h3("Loading brain slice..."),
      style="color:blue"
    ),
    hidden(
      div(
        id = "processingIndicator",
        h3("Processing correlations..."),
        style="color:blue"
      )
    ),
    hidden(
      div(id = "brainPlot",
        imageOutput('image1')
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {

  InstallData("stxBrain")
  brain <- LoadData("stxBrain", type = "posterior1")
  hide("loadingIndicator")
  show("brainPlot")

  # https://stackoverflow.com/questions/43619046/how-to-initialise-default-values-for-rendered-items-in-a-shiny-app
  values <-  reactiveValues(imgPlot = renderImage({
    img <- htmltools::capturePlot({
      plot(1, 1)
    }, height = 500, width = 500)
    list(src = img, width = 500, height = 500)
  }, deleteFile = TRUE))

  output$image1 <- renderImage({
    img <- htmltools::capturePlot({
      plot(GetImage(brain, mode ="raster"))
    }, height = 500, width = 500)
    list(src = img, width = 500, height = 500)
  }, deleteFile = TRUE)

  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = "sampleGeneExpr.csv",
    content = function(file) {
      download.file("https://raw.githubusercontent.com/jessielam415/pearsonCorrelateCells/master/inst/extdata/sampleGeneExpr.csv", file)
    }
  )

  # https://github.com/rstudio/shiny/issues/3170
  observeEvent(input$uploadedFile, {
    show("processingIndicator")
    req(input$uploadedFile)
    tryCatch(
      {
        uploadedFileDf <- read.csv(input$uploadedFile$datapath, header=FALSE)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    brainCorrelated <- findCorrelatedCells(comparator = uploadedFileDf, visiumData = brain,
                                           assay="Spatial", topGenes = nrow(uploadedFileDf))
    hide("processingIndicator")
    values$imgPlot <- renderImage({
      img <- htmltools::capturePlot({
        visualizeCells(brainCorrelated)
      }, height = 500, width = 500)
      list(src = img, width = 500, height = 500)
    }, deleteFile = TRUE)
    output$image1 <- values$imgPlot
  })
}

# Run the app ----
shinyApp(ui, server)
