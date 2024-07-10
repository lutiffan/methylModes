function(input, output) {
  # Arbitrarily set maximum file upload size to 100 GB
  # The real limit comes from the user's local memory limit
  # Set to max memory?
  options(shiny.maxRequestSize = 100*1024^3)
  
  ##### Get Started #####
  
  # Implicity returns betas as a reactive object
  getBetas <- reactive({
    req(input$betaFile$datapath)
    
    fileType <- file_ext(input$betaFile$datapath)
    
    if (fileType %in% c("RDS", "rds")) {
      readRDS(input$betaFile$datapath) 
    } else if (fileType %in% c("csv", "TXT", "txt", "tsv")) {
      as.matrix(data.table::fread(input$betaFile$datapath), rownames = 1)
    } else if (fileType %in% c("RDA", "rda")) {
      load(file = input$betaFile$datapath)
    } else {
      warning("Invalid file type.")
    }
    
  })
  
  output$dimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    nProbe <- nrow(betas)
    nSample <- ncol(betas)
    paste0("Found ", nProbe, " probes and ", nSample, " samples.")
  })
  
  getAnnotation <- reactive({
    req(input$arrayType)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    if (input$arrayType == "il450k") {
      if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
        BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
      }
      
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

    } else if (input$arrayType == "ilepic1") {
      if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
        BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
      }
      
      # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      # I think this one is more up-to-date
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      
    } else if (input$arrayType == "ilepic2") {
      if (!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38")) {
        BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
      }
      
      library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    }
    annotationData <- data.frame(Chromosome = 
                                   sub("chr", "", Locations[rownames(betas), "chr"]),
                                 Island = Islands.UCSC[rownames(betas), 
                                                       "Relation_to_Island"])
    annotationData
  })
  
  # output$betaOverview <- renderPlot({
  #   betas <- getBetas()
  #   if (is.null(betas)) return()
  #   
  #   probeMeans <- rowMeans(betas)
  #   
  #   hist(probeMeans, col = "#75AADB", border = "white",
  #        xlab = "Beta value",
  #        main = "Average methylation proportions across probes")
  #   
  # })
  
  output$betaOverview <- renderPlotly({
    betas <- getBetas()
    if (is.null(betas)) return()

    probeMeans <- data.frame("Beta" = round(rowMeans(betas), 2))
    
    p <- ggplot(probeMeans, aes(Beta)) +
      # geom_histogram(aes(y = ..density.., text = paste('Histogram Bar Density:', round(..density.., 2))), binwidth = 0.05) +
      # geom_line(aes(y = ..density.., text = paste('Smoothed Density:', round(..density.., 2))), stat = 'density', size = 1) +
      geom_histogram(aes(y = ..density..)) +
      geom_line(aes(y = ..density.., ), stat = 'density') +
      labs(title = "Density of Average Methylation Proportion Per Probe") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    ggplotly(p)

    # hist(probeMeans, col = "#75AADB", border = "white",
    #      xlab = "Beta value",
    #      main = "Average methylation proportions across probes")

  })
  
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
      labs(title = "Probe Locations by Chromosome", x = "Chromosome", y = "Count") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
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
      labs(title = "CpG Island Coverage", x = "", y = "Count") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    ggplotly(p)
  })
  
  observeEvent(input$betaFile, {
    if (!is.null(input$betaFile$datapath) & 
        file_ext(input$betaFile$datapath) %in% 
        c("RDS", "rds", "csv", "TXT", "txt", "tsv", "RDA", "rda")) {
      
      # Probe analysis inputs
      betas <- getBetas()
      numProbes <- nrow(betas)
      firstProbe <- rownames(betas)[1]
      
      # Remove disabling from buttons if the appropriate analysis type is selected
      # if (input$analysisType == "individual") {
        shinyjs::enable("runProbe")
        shinyjs::enable("runProbeRandom")
        shinyjs::enable("probeId")
        updateTextInput(inputId = "probeId", value = firstProbe)
      # } else {
        shinyjs::enable("runBetaMatrix")
        shinyjs::enable("rangeStart")
        shinyjs::enable("rangeEnd")
        updateNumericInput(inputId = "rangeStart", value = 1)
        updateNumericInput(inputId = "rangeEnd", value = numProbes)
      #}
    }
  })
  
  ##### Probe-level analysis #####
  
  observeEvent(input$runProbeRandom, {
    req(input$betaFile$datapath)
    betas <- getBetas()
    
    updateTextInput(inputId = "probeId", value = rownames(betas)[sample(1:nrow(betas), 1)])
  })
  
  # TODO: use peakSummary when available
  getProbeVisual <- reactive({
    req(input$runProbe | input$runProbeRandom)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # Set MethylModes parameters
    proportionSample <- input$proportionSample
    personalSpace <- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    probeId <- input$probeId
    # probeId <- "cg27399079" # for testing
    # histogramBins <- 50
    rowId <- which(rownames(betas) == probeId)
    histogramBins <- input$numHistogramBins
    
    calculate <- methylModes(row.data = betas[rowId,])
    
    detectedPeaks <- calculate$detected
    fittedDensity <- calculate$probeDensityEst
    
    rowId <- which(rownames(betas) == probeId)
    # histData <- betas[rowId,]
    histData <- data.frame("beta" = betas[rowId,])
    histogramBins <- input$numHistogramBins
    
    # Generate the histogram data
    ggHist <- ggplot(histData, aes(x = beta)) +
      geom_histogram(binwidth = 1 / histogramBins, fill = "black", color = "white") +
      labs(title = paste("Beta distribution for probe", probeId),
           x = "Beta",
           y = "Density") +
      xlim(-0.05, 1.05) + ylim(0, max(fittedDensity$y))
    
    fittedDensityDF <- data.frame(x = fittedDensity$x, y = fittedDensity$y)
    
    # Add density lines
    ggHist <- ggHist +
      geom_line(data = fittedDensityDF, aes(x = x, y = y), color = "orchid", size = 1)
    
    dataFrameForMaxima <- data.frame(beta = round(fittedDensity$x[detectedPeaks$maximaIdx], 2),
                                     density = fittedDensity$y[detectedPeaks$maximaIdx])
    dataFrameForLeftMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$leftMinIdx], 2),
                                      density = fittedDensity$y[detectedPeaks$leftMinIdx])
    dataFrameForRightMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$rightMinIdx], 2),
                                       density = fittedDensity$y[detectedPeaks$rightMinIdx])
    
    
    # Add points for the peaks (maxima and minima) using the correctly structured data frames
    ggHist <- ggHist +
      geom_point(data = dataFrameForMaxima, aes(x = beta, y = density), color = "red", size = 3) +
      geom_point(data = dataFrameForLeftMin, aes(x = beta, y = density), color = "blue", size = 3) +
      geom_point(data = dataFrameForRightMin, aes(x = beta, y = density), color = "blue", size = 3)
    
    ggHistPlotly <- ggplotly(ggHist, tooltip = "beta")
    
    return(ggHistPlotly)
    

  })

  getProbeVisualBaseR <- reactive({
    req(input$runProbe)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # Set MethylModes parameters
    proportionSample <- input$proportionSample
    personalSpace <- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    probeId <- input$probeId
    # probeId <- "cg27399079" # for testing
    # histogramBins <- 50
    rowId <- which(rownames(betas) == probeId)
    histogramBins <- input$numHistogramBins
    
    calculate <- methylModes(row.data = betas[rowId,])
    
    detectedPeaks <- calculate$detected
    fittedDensity <- calculate$probeDensityEst
    
    title <- paste("Beta distribution for probe", probeId)
    labelHeight <- 0.05*max(fittedDensity$y)
    hist(betas[rowId,], breaks = histogramBins, xlim = c(0,1),
         probability = TRUE, main = title, col = "black")
    lines(fittedDensity, col = "orchid", lwd = 2, pch = 19)
    points(fittedDensity$x[detectedPeaks$maximaIdx], fittedDensity$y[detectedPeaks$maximaIdx], col = "red", pch = 19, cex = 1.25)
    points(fittedDensity$x[detectedPeaks$leftMinIdx], fittedDensity$y[detectedPeaks$leftMinIdx], col = "blue", pch = 19, cex = 1.25)
    points(fittedDensity$x[detectedPeaks$rightMinIdx], fittedDensity$y[detectedPeaks$rightMinIdx], col = "blue", pch = 19, cex = 1.25)
    text(fittedDensity$x[detectedPeaks$maximaIdx],
         rep(labelHeight, nrow(detectedPeaks)), col = "black",
         labels = round(fittedDensity$x[detectedPeaks$maximaIdx], 2), pos = 3)
  })
  
  output$probeVisual <- renderPlotly({
    req(getProbeVisual())
    getProbeVisual()  # Return the plotly object
  })
  
  output$probeVisualBaseR <- renderPlot({
    req(input$runProbe | input$runProbeRandom)
    # req(getProbeVisual())
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # Set MethylModes parameters
    proportionSample <- input$proportionSample
    personalSpace <- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    probeId <- input$probeId
    # probeId <- "cg27399079" # for testing
    # histogramBins <- 50
    rowId <- which(rownames(betas) == probeId)
    histogramBins <- input$numHistogramBins
    
    calculate <- methylModes(row.data = betas[rowId,])
    
    detectedPeaks <- calculate$detected
    fittedDensity <- calculate$probeDensityEst

    title <- paste("Beta distribution for probe", probeId)
    labelHeight <- 0.05*max(fittedDensity$y)
    hist(betas[rowId,], breaks = histogramBins, xlim = c(0,1),
         probability = TRUE, main = title, xlab = "Beta", ylab = "Density", 
         col = "black")
    lines(fittedDensity, col = "orchid", lwd = 2, pch = 19)
    points(fittedDensity$x[detectedPeaks$maximaIdx], fittedDensity$y[detectedPeaks$maximaIdx], col = "red", pch = 19, cex = 1.25)
    points(fittedDensity$x[detectedPeaks$leftMinIdx], fittedDensity$y[detectedPeaks$leftMinIdx], col = "blue", pch = 19, cex = 1.25)
    points(fittedDensity$x[detectedPeaks$rightMinIdx], fittedDensity$y[detectedPeaks$rightMinIdx], col = "blue", pch = 19, cex = 1.25)
    # text(fittedDensity$x[detectedPeaks$maximaIdx],
    #      rep(labelHeight, nrow(detectedPeaks)), col = "black",
    #      labels = round(fittedDensity$x[detectedPeaks$maximaIdx], 2), pos = 3)
    # req(getProbeVisualBaseR())
    # getProbeVisualBaseR()
    # hist(1:10)
  })
  
  ##### Beta matrix-level analysis #####
  
  getPeakSummary <- reactive({

    req(input$runBetaMatrix)
    betas <- getBetas()
    if (is.null(betas)) return()

    # TODO require start and end range
    
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
    print(paste("Running in parallel on", availableCores(), "cores"))
    # Note that I had a bug here where I forgot to restrict the beta matrix
    # to the user-requested range! Fixed 6/27/24
    peakSummary <- fillPeakSummaryParallel(betas[rangeStart:rangeEnd,])
    
    # Sort results by number of detected peaks, descending
    # peakSummary <- peakSummary[order(peakSummary$numPeaks, decreasing = TRUE),]
    
    return(peakSummary)
  })
  
  observeEvent(input$runBetaMatrix, {
    req(input$betaFile$datapath)
    
    getPeakSummary()
    
    shinyjs::enable("downloadPeakSummary")
  }) 
  
  output$peakCountBar <- renderPlotly({
    peakSummary <- getPeakSummary()
    if (is.null(peakSummary)) return()

    peakCounts <- data.frame(Peaks = peakSummary$numPeaks)

    p <- ggplot(peakCounts, aes(x = Peaks)) +
      geom_bar() +
      labs(title = "Barplot of Peak Counts per Probe", x = "Number of peaks detected", y = "Count")
    ggplotly(p)
  })
  
  output$peakSummaryPreview <- renderPlot({
    peakSummary <- getPeakSummary()
    if (is.null(peakSummary)) return()
    
    betas <- getBetas()
    if (is.null(betas)) return()
    
    whichProbes <- sample(1:nrow(peakSummary), 9)
    par(mfrow = c(3,3))
    for (i in 1:9) {
      peakSummaryPlot(beta.data = betas[whichProbes[i],,drop = FALSE], 
                      peak.summary = peakSummary[whichProbes[i],])
    }
    
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
  
  # Conditional features
  # Render different action buttons in the sidebar based on selection
  # output$conditionalSidebarButtons <- renderUI({
  #   if (input$analysisType == "individual") {
  #     # Individual Probe Analysis buttons
  #     list(
  #       shinyjs::disabled(textInput("probeId", "Probe Id",
  #                                   value = "cg27399079")),
  #       # add_busy_spinner(spin = "self-building-square",
  #       #                  timeout = 500,
  #       #                  position = "top-right"),
  #       shinyjs::disabled(actionButton("runProbe", "Run on Selected Probe")),
  #       div(style = "margin-bottom: 10px;"),
  #       shinyjs::disabled(actionButton("runProbeRandom", "Run on Randomly Selected Probe"))
  #     )
  #   } else if (input$analysisType == "multi_probe") {
  #     # Whole Genome Analysis buttons
  #     list(
  #       actionButton("action3", "Whole Genome 1"),
  #       actionButton("action4", "Whole Genome 2")
  #       # ... Add other UI elements here for Whole Genome Analysis
  #     )
  #   }
  # })
  
  # output$conditionalPlots <- renderUI({
  #   if (input$analysisType == "individual") {
  #     list(
  #       h4("Summary of Probe-Level MethylModes Results"),
  #       plotlyOutput(outputId = "probeVisual"),
  #       plotOutput(outputId = "probeVisualBaseR")
  #     )
  #   } else if (input$analysisType == "multi_probe") {
  #   }
  # })
  
}