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
  
  output$wholeDataDimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    nProbe <- nrow(betas)
    nSample <- ncol(betas)
    paste0("Found ", nProbe, " probes and ", nSample, " samples.")
  })
  
  # Check whether annotation packages are already installed or not
  missing450kAnnotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylation450kanno.ilmn12.hg19"))))
  missingEPICAnnotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))))
  missingEPICV2Annotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38"))))
  
  output$missing450kAnnotation <- reactive({
    missing450kAnnotation()
  })
  
  outputOptions(output, "missing450kAnnotation", suspendWhenHidden = FALSE)
  
  output$missingEPICAnnotation <- reactive({
    missingEPICAnnotation()
  })
  
  outputOptions(output, "missingEPICAnnotation", suspendWhenHidden = FALSE)
  
  output$missingEPICV2Annotation <- reactive({
    missingEPICV2Annotation()
  })
  
  outputOptions(output, "missingEPICV2Annotation", suspendWhenHidden = FALSE)
  
  getAnnotation <- reactive({
    req(input$arrayType)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    if (input$arrayType == "il450k") {
      if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
        BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
        missing450kAnnotation(FALSE)
      }
      
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

    } else if (input$arrayType == "ilepic1") {
      if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
        BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
        missingEPICAnnotation(FALSE)
      }
      
      # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      # I think this one is more up-to-date
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      
    } else if (input$arrayType == "ilepic2") {
      if (!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38")) {
        BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
        missingEPICV2Annotation(FALSE)
      }
      
      library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    }

    annotationData <- data.frame(Chromosome = 
                                   sub("chr", "", Locations[rownames(betas), "chr"]),
                                 Island = Islands.UCSC[rownames(betas), 
                                                       "Relation_to_Island"],
                                 Position = Locations[rownames(betas), "pos"])

    annotationData
  })
  
  # Will set to TRUE when probe-average plot on "get started" page is created
  plotCreatedBetaOverview <- reactiveVal(FALSE)
  
  # Updates according to the plotCreatedBetaOverview reactive value
  output$plotCreatedBetaOverview <- reactive({
    plotCreatedBetaOverview()
  })
  
  outputOptions(output, "plotCreatedBetaOverview", suspendWhenHidden = FALSE)
  
  output$betaOverview <- renderPlotly({
    betas <- getBetas()
    if (is.null(betas)) return()

    probeMeans <- data.frame("Beta" = round(rowMeans(betas), 2))
    p <- ggplot(probeMeans, aes(Beta)) +
      geom_histogram(aes(y = after_stat(density)), 
                     binwidth = 1 / input$numHistogramBinsBetaOverview) +
      geom_line(aes(y = after_stat(density), ), stat = 'density') +
      labs(title = "Density of Average Methylation Proportion Per Probe",
           x = "Beta",
           y = "Density") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    # Enables display of slider bar
    plotCreatedBetaOverview(TRUE)
    ggplotly(p)
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
        shinyjs::enable("runMultiProbe")
        shinyjs::enable("rangeStart")
        shinyjs::enable("rangeEnd")
        updateNumericInput(inputId = "rangeStart", value = 1)
        updateNumericInput(inputId = "rangeEnd", value = numProbes)
      #}
    }
  })
  
  ##### Probe-level analysis #####
  
  # Reactive value to track updates needed for the eventReactive
  probeUpdate <- reactiveVal(NULL)
  
  observeEvent(input$runProbeRandom, {
    req(input$betaFile$datapath)
    betas <- getBetas()
    randomProbeId <- rownames(betas)[sample(1:nrow(betas), 1)]

    updateTextInput(inputId = "probeId", value = randomProbeId)
    
    probeUpdate(randomProbeId)
  })
  
  observeEvent(input$runProbe, {
    probeUpdate(input$probeId)
  })
  
  # Will set to TRUE when individual-probe visual is created
  plotCreatedProbeVisual <- reactiveVal(FALSE)
  
  output$plotCreatedProbeVisual <- reactive({
    plotCreatedProbeVisual()
  })
  
  outputOptions(output, "plotCreatedProbeVisual", suspendWhenHidden = FALSE)
  
  #### Set up graph and table after analyzing a single probe ####
  getSingleProbeSummary <- eventReactive(probeUpdate(), {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # Set MethylModes parameters
    proportionSample <<- input$proportionSample
    peakDistance <<- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    
    probeId <- probeUpdate()
    rowId <- which(rownames(betas) == probeId)

    methylModes(row.data = betas[rowId,])
  })
  
  # Example of trimodal: cg26261358
  output$probeTable <- renderTable({
    req(input$runProbe | input$runProbeRandom)
    betas <- getBetas()
    if (is.null(betas)) return()
    
    mmResults <- getSingleProbeSummary()
    if (is.null(mmResults)) return()
    
    rowId <- which(rownames(betas) == input$probeId)
    data.frame("numPeaks" = nrow(mmResults$detected), 
                              "meanBeta" = mean(betas[rowId,]),
                              "peakLocations" = mmResults$probeDensityEst$x[mmResults$detected$maximaIdx],
                              "leftMin" = mmResults$probeDensityEst$x[mmResults$detected$leftMinIdx],
                              "rightMin" = mmResults$probeDensityEst$x[mmResults$detected$rightMinIdx],
                              "proportionSample" = mmResults$detected$propSample,
                              "peakVariance" = var(betas[rowId,])) 
  })
  
  # problematic trimodal: cg04837091
  getProbeVisualAdHoc <- eventReactive(getSingleProbeSummary(), {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    mmResults <- getSingleProbeSummary()
    if (is.null(mmResults)) return()
    
    detectedPeaks <- mmResults$detected
    fittedDensity <- mmResults$probeDensityEst
    
    fittedDensityDF <- data.frame(x = fittedDensity$x, y = fittedDensity$y)
    dataFrameForMaxima <- data.frame(beta = round(fittedDensity$x[detectedPeaks$maximaIdx], 2),
                                     density = fittedDensity$y[detectedPeaks$maximaIdx])
    dataFrameForLeftMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$leftMinIdx], 2),
                                      density = fittedDensity$y[detectedPeaks$leftMinIdx])
    dataFrameForRightMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$rightMinIdx], 2),
                                       density = fittedDensity$y[detectedPeaks$rightMinIdx])
    
    probeId <- probeUpdate()
    rowId <- which(rownames(betas) == probeId)

    # histData <- betas[rowId,]
    histData <- data.frame("beta" = betas[rowId,])
    # hist(histData$beta, breaks = 50)
    
    # Colors chosen using the following code: 
    # hcl.colors(3, palette = 'viridis')
    # Add points for the peaks (maxima and minima) using the correctly structured data frames
    ggHist <- ggplot(histData, aes(x = beta)) +
      geom_histogram(aes(y = after_stat(density)), 
                     binwidth = 1 / input$numHistogramBinsOneProbe) +
      labs(title = paste("Beta distribution for Probe", probeId),
           x = "Beta",
           y = "Density") +
      xlim(-0.05, 1.05) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      geom_point(data = dataFrameForMaxima, aes(x = beta, y = density, color = "Maxima"), size = 3) +
      scale_color_manual(
        name = "Legend",
        values = c(
          "Density Estimate" = "#00588B", 
          "Maxima" = "#B2DC3C", 
          "Minima" = "#009B95"
        )
      )
    
    if (input$showDensitySingleProbe) {
      ggHist <- ggHist + 
        geom_line(data = fittedDensityDF, aes(x = x, y = y, color = "Density Estimate"))
    }
    if (input$showMinimaSingleProbe) {
      ggHist <- ggHist + 
        geom_point(data = dataFrameForLeftMin, aes(x = beta, y = density, color = "Minima"), size = 3) +
        geom_point(data = dataFrameForRightMin, aes(x = beta, y = density, color = "Minima"), size = 3)
    }
    
    ggHistPlotly <- ggplotly(ggHist, tooltip = "beta")
    
    # Enable the slider bar
    plotCreatedProbeVisual(TRUE)
    return(ggHistPlotly)
  })
  

  
  output$probeVisual <- renderPlotly({
    req(getProbeVisualAdHoc())
    getProbeVisualAdHoc()  # Return the plotly object
  })
  
  ##### Beta matrix-level analysis #####
  getMultiProbeSummary <- eventReactive(input$runMultiProbe, {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # betas <- readRDS("/home/lutiffan/betaMatrix/smolBetas.RDS")

    # Set MethylModes parameters
    proportionSample <<- input$proportionSample
    peakDistance <<- input$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero
    # rangeStart <- input$rangeStart
    # rangeEnd <- input$rangeEnd
    
    # totalRows = rangeEnd - rangeStart + 1
    print(paste("Running in parallel on", availableCores(), "cores"))
    # Note that I had a bug here where I forgot to restrict the beta matrix
    # to the user-requested range! Fixed 6/27/24
    
    if (input$region == "wholeGenome") {
      peakSummary <- fillPeakSummaryParallel(betas)
    } else {
      peakSummary <- fillPeakSummaryParallel(betas[betaFilter(),])
    }
    
    # Sort results by number of detected peaks, descending
    # peakSummary <- peakSummary[order(peakSummary$numPeaks, decreasing = TRUE),]
    
    return(peakSummary)
  })
  
  observeEvent(input$runMultiProbe, {
    req(input$betaFile$datapath)
    
    getMultiProbeSummary()
    reset("runMultiProbe")
    shinyjs::enable("downloadPeakSummary")
  }) 
  
  output$peakCountBar <- renderPlotly({
    req(input$runMultiProbe)
    peakSummary <- getMultiProbeSummary()
    if (is.null(peakSummary)) return()

    peakCounts <- data.frame(Peaks = peakSummary$numPeaks)

    p <- ggplot(peakCounts, aes(x = Peaks)) +
      geom_bar() +
      labs(title = "Number of Peaks Detected Per Probe", x = "Peaks", y = "Count")
    ggplotly(p)
  })
  
  output$peakSummaryPreview <- renderPlot({
    peakSummary <- getMultiProbeSummary()
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

  output$peakSummaryTable <- DT::renderDataTable({
    peakSummary <- getMultiProbeSummary()
    if (is.null(peakSummary)) return()
    
    listToString <- function(listData, decimalPlaces) {
      paste(as.character(round(unlist(listData), decimalPlaces)), collapse = ", ")
    }

    peakLocationsPreviewFriendly <- unlist(lapply(peakSummary$peakLocations, 
                                                  FUN = listToString,
                                                  decimalPlaces = 2))
    proportionSamplePreviewFriendly <- unlist(lapply(peakSummary$proportionSample, 
                                                     FUN = listToString,
                                                     decimalPlaces = 3))
    peakVariancePreviewFriendly <- unlist(lapply(peakSummary$peakVariance, 
                                                 FUN = listToString,
                                                 decimalPlaces = 6))
    
    datatable(data.frame("probeName" = peakSummary$probeName,
                         "numPeaks" = as.factor(peakSummary$numPeaks),
                         "meanBeta" = round(peakSummary$meanBeta, 2),
                         "peakLocations" = peakLocationsPreviewFriendly,
                         "proportionSample" = proportionSamplePreviewFriendly,
                         "peakVariance" = peakVariancePreviewFriendly),
              selection = 
                list(mode = "single",
                     target = 'row+column',
                     selected = list(
                       rows = 1,
                       cols = c())),
              filter = "top",
              options = list(columnDefs = list(
                list(searchable = FALSE, targets = c(4, 5, 6)) 
              ))
    )
  })
  
  output$probeVisualFromPeakSummary <- renderPlot({
    betas <- getBetas()
    if (is.null(betas)) return()
    
    betaFilter <- betaFilter()
    if (is.null(betas)) return()
    
    peakSummary <- getMultiProbeSummary()
    if (is.null(peakSummary)) return()

    req(getMultiProbeSummary())
    
    # output$selectedRow <- renderPrint({
    #   selected <- input$peakSummaryTable_rows_selected
    #   if(length(selected)) {
    #     print(selected)
    #   } else {
    #     print("No row selected")
    #   }
    # })
    
    rowIdx <- input$peakSummaryTable_rows_selected
    peakSummaryPlot(beta.data = betas[betaFilter,][rowIdx, , drop = FALSE], peak.summary = peakSummary[rowIdx,])
    
    # beta <- betas[betaFilter,][rowIdx, , drop = FALSE]
    # probeId <- rownames(beta)
    # 
    # histData <- data.frame("beta" = beta)
    # 
    # bandwidthType = NULL
    # if (is.na(BANDWIDTH_TYPE)) {
    #   bandwidthType <- "nrd0"
    # } else if (bandwidthType == "sheatherJones") {
    #   bandwidthType <- bw.SJ(row.data)
    # }
    # 
    # fittedDensity <- density(beta, from = 0, to = 1, n = NUM_BREAKS, 
    #         adjust = input$densityAdjust, bw = bandwidthType,
    #         kernel = KERNEL_TYPE)
    # 
    # fittedDensityDF <- data.frame(x = fittedDensity$x, y = fittedDensity$y)
    # dataFrameForMaxima <- data.frame(beta = round(fittedDensity$x[detectedPeaks$maximaIdx], 2),
    #                                  density = fittedDensity$y[detectedPeaks$maximaIdx])
    # dataFrameForLeftMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$leftMinIdx], 2),
    #                                   density = fittedDensity$y[detectedPeaks$leftMinIdx])
    # dataFrameForRightMin <- data.frame(beta = round(fittedDensity$x[detectedPeaks$rightMinIdx], 2),
    #                                    density = fittedDensity$y[detectedPeaks$rightMinIdx])
    # 
    # ggHist <- ggplot(histData, aes(x = beta)) +
    #   geom_histogram(aes(y = after_stat(density)), 
    #                  binwidth = 1 / input$numHistogramBinsOneProbe) +
    #   labs(title = paste("Beta distribution for Probe", probeId),
    #        x = "Beta",
    #        y = "Density") +
    #   xlim(-0.05, 1.05) + 
    #   theme(panel.grid.major = element_blank(), 
    #         panel.grid.minor = element_blank(),
    #         panel.background = element_blank(), 
    #         axis.line = element_line(colour = "black")) +
    #   geom_line(data = fittedDensityDF, aes(x = x, y = y, color = "Density Estimate"), 
    #             linewidth = 1) +
    #   geom_point(data = dataFrameForMaxima, aes(x = beta, y = density, color = "Maxima"), size = 3) +
    #   geom_point(data = dataFrameForLeftMin, aes(x = beta, y = density, color = "Minima"), size = 3) +
    #   geom_point(data = dataFrameForRightMin, aes(x = beta, y = density, color = "Minima"), size = 3) +
    #   scale_color_manual(
    #     name = "Legend",
    #     values = c(
    #       "Density Estimate" = "#00588B", 
    #       "Maxima" = "#B2DC3C", 
    #       "Minima" = "#009B95"
    #     )
    #   )
    # 
    # ggHistPlotly <- ggplotly(ggHist, tooltip = "beta")
  })
  
  # output$selectedRow <- renderPrint({
  #   selected <- input$peakSummaryTable_rows_selected
  #   if(length(selected)) {
  #     print(selected)
  #   } else {
  #     print("No row selected")
  #   }
  # })
  
  # output$probeVisualFromPeakSummary <- renderPlot({
  #   req(getProbeVisualFromPeakSummary())
  #   getProbeVisualFromPeakSummary()  # Return the plotly object
  # })
  
  # Define a download handler
  output$downloadPeakSummary <- downloadHandler(
    filename = function() {
      paste("peakSummary", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      fwrite(getMultiProbeSummary(), file)
    }
  )
  
  # Create a list of the unique chromosomes represented in data
  uniqueChromosomes <- reactive({
    annotationData <- getAnnotation()
    if (is.null(annotationData)) return(NULL)
    
    allChromosomes <- c(as.character(1:22), "X", "Y")
    intersect(allChromosomes, unique(annotationData$Chromosome))
  })
  
  output$chromosomeSelect <- renderUI({
    req(uniqueChromosomes()) 
    selectInput("selectedChromosome", "Choose Chromosome:",
                choices = uniqueChromosomes(),
                selected = uniqueChromosomes()[1])
  })
  
  # Get minimum and maximum base pair locations from data
  basePairMinMax <- reactive({
    annotationData <- getAnnotation()
    if (is.null(annotationData)) return(NULL)
    
    list("bpMin" = min(annotationData$Position), 
         "bpMax" = max(annotationData$Position))
  })
  
  output$basePairRangeSelect <- renderUI({
    req(basePairMinMax())

    numericRangeInput("selectedBasePairRange", 
                      "Enter CpG location range",
                      value = c(basePairMinMax()$bpMin, basePairMinMax()$bpMax),
                      min = basePairMinMax()$bpMin,
                      max = basePairMinMax()$bpMax)
  })
  
  # Create a vector used to subset the beta matrix
  betaFilter <- reactive({
    annotationData <- getAnnotation()
    if (is.null(annotationData)) return(NULL)
    
    if (input$region == "chromosome") {
      req(input$selectedChromosome)
      rowsToKeep <- annotationData$Chromosome == input$selectedChromosome
    } else if (input$region == "basePair") {
      req(input$selectedBasePairRange)
      rowsToKeep <- annotationData$Position >= input$selectedBasePairRange[1] & 
        annotationData$Position <= input$selectedBasePairRange[2]
    } else {
      rowsToKeep <- TRUE
    }

    rowsToKeep
  })
  
  output$subsetDataDimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    betaFilter <- betaFilter()
    
    nProbe <- sum(betaFilter)
    paste(nProbe, "probes selected.")
  })
}