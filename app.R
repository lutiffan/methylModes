# Shiny-specific
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinycssloaders)

# Tables and plots
library(data.table)
library(dplyr)
library(ggplot2)
library(plotly)

# Parallelization
library(foreach)
library(iterators)
library(parallelly)
library(doParallel)

# Load default parameters
source(file = "standardParams.R", local = T)

# Load functions
source(file = "methylModes.R", local = T)
source(file = "fillPeakSummaryParallel.R", local = T)
source(file = "visualizeBetaPeaks.R", local = T)

# Global variables
KERNEL_TYPE = "gaussian"
BANDWIDTH_TYPE = NA
NUM_BREAKS = 500

# TODO: let user set max file upload size via slider??
ui <- fluidPage( # create a display that automatically adjusts to the dimensions of your users browser window
  useShinyjs(), # Gives us tricks like disabling an element until condition satisfied
  navbarPage("MethylModes",
    tabPanel("Get Started",
             h3("Data Upload"),
             sidebarPanel(
               radioButtons("probeType", "Probe type",
                            selected = character(0),
                            c("450k" = "il450k",
                              "EPIC" = "ilepic")),
               fileInput("betaFile", "Choose Beta Matrix File",
                         multiple = FALSE,
                         accept = c(".RDS", ".RDA")), # Most people are not gonna have an RDS format
               helpText("Accepted file formats: .RDS, .RDA", br(),
                        "Max file size: 20 GB"),
               p("Possible features: upload phenotype , select annotation file to use",
                 "Slider bar for histogram, show kernel density estimate")
            ),
             h4("Data summary (generated upon successful file upload)"),
             mainPanel(
               # Histogram of all probes
               textOutput(outputId = "dimensions"),
               withSpinner(plotOutput(outputId = "betaOverview")),
               withSpinner(plotlyOutput(outputId = "chromosomeBar")),
               withSpinner(plotlyOutput(outputId = "islandBar"))
              )
    ),
    
    tabPanel("Run MethylModes",
               h3("Settings"),
               sidebarPanel(
                 numericInput(inputId = "proportionSample",
                                "Minimum proportion of sample considered to be a peak",
                                value = 0.05),
                   numericInput(inputId = "peakDistance",
                                "Minimum distance between adjacent peaks",
                                value = 0.1),
                   h5("Smoothing parameters"),
                   numericInput("densityAdjust", "density() 'adjust' parameter",
                                value = 1.5),
                   numericInput("pushToZero", 
                                "Threshold for numbers small enough to be set to zero",
                                value = 1e-6),
                   helpText("The density() function fits very small floating-point",
                            "values to the data, creating tiny perturbations that ",
                            "reduce the accuracy of MethylModes. Setting a threshold ",
                            "under which small values are considered equivalent ",
                            "to zero mitigates this issue."),
                   shinyjs::disabled(numericInput("rangeStart", "Start of range of beta matrix rows",
                                     value = 1)),
                   shinyjs::disabled(numericInput("rangeEnd", "End of range of beta matrix rows",
                                     value = 1)),
                 add_busy_spinner(spin = "self-building-square", 
                                  timeout = 500,
                                  position = "top-right"),
                   shinyjs::disabled(actionButton("run", "Run MethylModes"))
               ),
             mainPanel(
               h4("Summary of MethylModes Results"),
               shinyjs::disabled(downloadButton("downloadPeakSummary",
                                                "Download Results")),
               plotlyOutput("peakCountBar"),
               h3("Coming soon: view sorted MethylModes results (e.g. among multimodal probes, show probes with highest number of modes to fewest")
                       )
             ),
    tabPanel("Review Previously Generated Results",
              h3("Upload Previously Calculated Results"),
              sidebarPanel(
                fileInput("peakSummaryFile", "Choose MethylModes Result File",
                          multiple = FALSE,
                           accept = c(".RDS",
                                      ".csv"))
              )
             )
    )
    
)

server <- function(input, output) {
  # Set maximum file upload size to 20 GB
  # Set to max memory?
  options(shiny.maxRequestSize = 20*1024^3)
  
  getBetas <- reactive({
    req(input$betaFile$datapath)
    readRDS(input$betaFile$datapath)
  })
  
  output$dimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    nProbe <- nrow(betas)
    nSample <- ncol(betas)
    paste0("Found ", nProbe, " probes and ", nSample, " samples.")
  })
  
  output$betaOverview <- renderPlot({
    betas <- getBetas()
    if (is.null(betas)) return()
    
    probeMeans <- rowMeans(betas)

    hist(probeMeans, col = "#75AADB", border = "white",
         xlab = "Beta value",
         main = "Histogram of methylation proportions across probes")

  })
  
  getAnnotation <- reactive({
    req(input$probeType)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    if (input$probeType == "il450k") {
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      annotationData <- data.frame(Chromosome = 
                                     sub("chr", "", Locations[rownames(betas), 
                                                              "chr"]),
                                   Island = Islands.UCSC[rownames(betas), 
                                                         "Relation_to_Island"])
    } else if (input$probeType == "ilepic") {
      library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      annotationData <- data.frame(Chromosome = 
                                     sub("chr", "", Locations[rownames(betas), "chr"]))
    }
    annotationData
    
    # if (input$probeType == "il450k") {
    #   library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    #   
    # }
    # else if (input$probeType == "ilepic") {
    #   library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    # }
  })
  
  # TODO: annotations. Which one do I use for EPIC?
  # Make an accompanying table so you know there is data for the Y chromosome?
  output$chromosomeBar <- renderPlotly({
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # relevantLabels <- as.data.frame(
    #   sub(pattern = "chr", replacement = "", x = Locations[rownames(betas), "chr"])
    # )
    # names(relevantLabels) <- "Chromosome"
    
    relevantLabels <- getAnnotation()
    if (is.null(relevantLabels)) return(NULL)
    
    relevantLabels$Chromosome <- factor(relevantLabels$Chromosome,
                                        levels = c(as.character(1:22), "X", "Y"))
    
    p <- ggplot(relevantLabels, aes(x = Chromosome)) +
      geom_bar() +
      labs(title = "Barplot of Chromosome Counts", x = "Chromosome", y = "Count")
    ggplotly(p)
  })
  
  output$islandBar <- renderPlotly({
    betas <- getBetas()
    if (is.null(betas)) return()
    
    relevantLabels <- getAnnotation()
    if (is.null(relevantLabels)) return(NULL)
    
    relevantLabels$Island <- factor(relevantLabels$Island,
                                    levels = c("N_Shelf", "N_Shore", "Island",
                                               "S_Shore", "S_Shelf", "OpenSea"))
    
    p <- ggplot(relevantLabels, aes(x = Island)) +
      geom_bar() +
      labs(title = "Barplot of Probe Relations to Islands", x = "Relation", y = "Count")
    ggplotly(p)
  })
  
  observeEvent(input$betaFile, {
    if (!is.null(input$betaFile$datapath)) {
      
      betas <- getBetas()
      shinyjs::enable("run")  # Enable the button when the file is uploaded
      shinyjs::enable("rangeStart")
      shinyjs::enable("rangeEnd")
      updateNumericInput(inputId = "rangeEnd", value = nrow(betas))
    }
  })
  
  getPeakSummary <- reactive({
    req(input$run)
    
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # betas <- readRDS("/home/lutiffan/betaMatrix/smolBetas.RDS")
    
    # Set MethylModes parameters
    proportionSample <- input$proportionSample
    personalSpace <- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    rangeStart <- input$rangeStart
    rangeEnd <- input$rangeEnd
    
    totalRows = rangeEnd - rangeStart + 1
    
    peakSummary <<- fillPeakSummaryParallel(betas)
  })
  
  observeEvent(input$run, {
    req(input$betaFile$datapath)
    
    getPeakSummary()
    
    shinyjs::enable("downloadPeakSummary")
  }) # End of observeEvent() for "run" button
  
  
  output$peakCountBar <- renderPlotly({
    peakSummary <- getPeakSummary()
    if (is.null(peakSummary)) return()
    
    peakCounts <- data.frame(Peaks = peakSummary$numPeaks)
    
    p <- ggplot(peakCounts, aes(x = Peaks)) +
      geom_bar() +
      labs(title = "Barplot of Peak Counts per Probe", x = "Number of peaks detected", y = "Count")
    ggplotly(p)
  })
  
  # Define a download handler
  output$downloadPeakSummary <- downloadHandler(
    filename = function() {
      paste("peakSummary", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      fwrite(getPeakSummary(), file)
    }
  )
  
}

shinyApp(ui, server)
