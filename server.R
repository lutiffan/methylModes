function(input, output) {
  # Arbitrarily set maximum file upload size to 100 GB
  # The real limit comes from the user's local memory limit
  # Is there a simple way to set this to max memory based on local machine?
  options(shiny.maxRequestSize = 100*1024^3)
  
  ##### Get Started #####
  
  # Implicitly returns betas as a reactive object
  getBetas <- reactive({
    req(input$betaFile$datapath)
    
    fileType <- file_ext(input$betaFile$datapath)
    
    if (fileType %in% c("RDS", "rds")) {
      betas <- readRDS(input$betaFile$datapath) 
    } else if (fileType %in% c("csv", "TXT", "txt", "tsv")) {
      betas <-as.matrix(data.table::fread(input$betaFile$datapath), rownames = 1)
    } else if (fileType %in% c("RDA", "rda")) {
      betas <- load(file = input$betaFile$datapath)
    } else if (fileType == "qs") {
      betas <- qread(file = input$betaFile$datapath)
    } else {
      warning("Invalid file type.")
    }
    # Sort beta matrix (important for matching with annotation data)
    betas[order(rownames(betas)),]
  })
  
  output$wholeDataDimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    nProbe <- nrow(betas)
    nSample <- ncol(betas)
    paste0("Found ", nProbe, " probes and ", nSample, " samples.")
  })
  
  output$numProbesAnalyzed <- renderUI({
    betaFilter <- betaFilter()
    nProbe <- sum(betaFilter)
    h4(paste0("A total of ", nProbe, " probes were analyzed."))
  })
  
  # # Check whether annotation packages are already installed or not
  # missing450kAnnotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylation450kanno.ilmn12.hg19"))))
  # missingEPICAnnotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))))
  # missingEPICV2Annotation <- reactiveVal(suppressMessages(suppressWarnings(!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38"))))
  # 
  # output$missing450kAnnotation <- reactive({
  #   missing450kAnnotation()
  # })
  # 
  # outputOptions(output, "missing450kAnnotation", suspendWhenHidden = FALSE)
  # 
  # output$missingEPICAnnotation <- reactive({
  #   missingEPICAnnotation()
  # })
  # 
  # outputOptions(output, "missingEPICAnnotation", suspendWhenHidden = FALSE)
  # 
  # output$missingEPICV2Annotation <- reactive({
  #   missingEPICV2Annotation()
  # })
  # 
  # outputOptions(output, "missingEPICV2Annotation", suspendWhenHidden = FALSE)
  
  getAnnotationLocal <- reactive({
    if (is.null(input$arrayType)) return(NULL)
    betas <- getBetas()
    if (is.null(betas)) return(NULL)
    
    if (input$arrayType == "il450k") {
      manifestFile <- fread(paste0(getwd(), "/IlluminaManifest450k.csv"))
      
    } else if (input$arrayType == "ilepic1") {
      manifestFile <- fread(paste0(getwd(), "/IlluminaManifestEPIC.csv"))
      
    } else {
      simpleError("Invalid annotation package selected.")
    }
    # TODO: get Epicv2 manifest data
    # else if (input$arrayType == "ilepic2") {
    #   if (!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38")) {
    #     BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
    #     missingEPICV2Annotation(FALSE)
    #   }
    #   
    #   library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    # } 
    
    
    # annotationData <- data.frame(Chromosome = 
    #                                sub("chr", "", Locations[rownames(betas), "chr"]),
    #                              Island = Islands.UCSC[rownames(betas), 
    #                                                    "Relation_to_Island"],
    #                              Position = Locations[rownames(betas), "pos"])
    
    # Keep only rows corresponding to probes that we have
    manifestFile <- manifestFile[manifestFile$IlmnID %in% rownames(betas),]
    # Now sort the probes in order of name
    manifestFile <- manifestFile[order(manifestFile$IlmnID),]
    
    annotationData <- data.frame(Chromosome = manifestFile$CHR,
                                 Island = manifestFile$Relation_to_UCSC_CpG_Island,
                                 Position = manifestFile$MAPINFO)
    
    if (input$arrayType == "il450k") {
      annotationData$SNP_within_10Bp <- manifestFile$Probe_SNPs_10
      annotationData$SNP_10Bp_and_beyond <- manifestFile$Probe_SNPs
    } else if (input$arrayType == "ilepic1") {
      annotationData$SNP_distance <- manifestFile$SNP_DISTANCE
    }
    
    # Replace empty strings in Island column with "OpenSea"
    annotationData$Island[annotationData$Island == ""] <- "OpenSea"
    annotationData
  })
  
  # getAnnotation <- reactive({
  #   if (is.null(input$arrayType)) return(NULL)
  #   betas <- getBetas()
  #   if (is.null(betas)) return(NULL)
  #   
  #   if (input$arrayType == "il450k") {
  #     if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  #       BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
  #       missing450kAnnotation(FALSE)
  #     }
  #     
  #     library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # 
  #   } else if (input$arrayType == "ilepic1") {
  #     if (!require("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")) {
  #       BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  #       missingEPICAnnotation(FALSE)
  #     }
  #     
  #     # library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  #     # I think this one is more up-to-date
  #     library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  #     
  #   } 
  #   # else if (input$arrayType == "ilepic2") {
  #   #   if (!require("IlluminaHumanMethylationEPICv2anno.20a1.hg38")) {
  #   #     BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  #   #     missingEPICV2Annotation(FALSE)
  #   #   }
  #   #   
  #   #   library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  #   # } 
  #   else {
  #     simpleError("Invalid annotation package selected.")
  #   }
  # 
  #   annotationData <- data.frame(Chromosome = 
  #                                  sub("chr", "", Locations[rownames(betas), "chr"]),
  #                                Island = Islands.UCSC[rownames(betas), 
  #                                                      "Relation_to_Island"],
  #                                Position = Locations[rownames(betas), "pos"])
  # 
  #   annotationData
  # })
  
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
      geom_histogram(aes(y = after_stat(density), fill = after_stat(x)), # Use x (Beta) to map to fill
                     binwidth = 1 / input$numHistogramBinsBetaOverview, 
                     color = "black") + # Add a border color for clarity
      geom_line(aes(y = after_stat(density)), stat = 'density', color = "black") + # Add density line
      labs(title = "Mean Methylation Across Probes",
           x = "Beta",
           y = "Density") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      scale_fill_viridis_c() # Apply viridis color scale
    
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
    
    relevantLabels <- getAnnotationLocal()
    if (is.null(relevantLabels)) return(NULL)
    
    relevantLabels$Chromosome <- factor(relevantLabels$Chromosome,
                                        levels = c(as.character(1:22), "X", "Y"))
    
    p <- ggplot(relevantLabels, aes(x = Chromosome)) +
      geom_bar(color = "#155087", fill = "#007996") +  # Set the bar fill color to cerulean
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
    
    relevantLabels <- getAnnotationLocal()
    if (is.null(relevantLabels)) return(NULL)
   
    relevantLabels$Island <- factor(relevantLabels$Island,
                                    levels = c("N_Shelf", "N_Shore", "Island",
                                               "S_Shore", "S_Shelf", "OpenSea"))

    islandColors <- c("Island" ="#FDE333"  ,       
                      "N_Shore" = "#C1DE35",      
                      "S_Shore" = "#00C376",      
                      "N_Shelf" = "#008498",      # Blue
                      "S_Shelf" = "#006892",      # Blue
                      "OpenSea" = "#363D7C")      # Dark Blue (or the closest that matches)
    
    p <- ggplot(relevantLabels, aes(x = Island, fill = Island)) +
      geom_bar(color = "black") +  # Add black border to the bars
      scale_fill_manual(values = islandColors) +  # Assign the custom colors
      labs(title = "CpG Island Coverage", x = "", y = "Count") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
    ggplotly(p)
  })
  
  observeEvent(input$betaFile, {
    if (!is.null(input$betaFile$datapath)) {
      
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
        # shinyjs::enable("rangeStart")
        # shinyjs::enable("rangeEnd")
        # updateNumericInput(inputId = "rangeStart", value = 1)
        # updateNumericInput(inputId = "rangeEnd", value = numProbes)
      #}
    }
  })
  
  # Disable the Analysis Type radio buttons on startup
  shinyjs::disable("analysisType")
  
  # Enable the radio buttons upon file upload
  observeEvent(input$betaFile, {
    shinyjs::enable("analysisType")
  })
  
  # Disable results cross-tabulation until results are available
  shinyjs::disable("crossTabMultimodal")
  shinyjs::disable("crossTabCpgIsland")
  
  observeEvent(c(input$peakSummaryFile, input$runMultiProbe), {
    shinyjs::enable("crossTabMultimodal")
    shinyjs::enable("crossTabCpgIsland")
  })
  
  ##### Probe-level analysis #####
  
  # # Reactive value to track updates needed for the eventReactive
  # probeUpdate <- reactiveVal(NULL)
  
  # Reactive value to track the probe parameters
  probeParams <- reactive({
    list(
      probeId = input$probeId,
      proportionSample = input$proportionSample,
      peakDistance = input$peakDistance,
      triggeredBy = input$runProbe + input$runProbeRandom,  # Summing button clicks to create a new trigger
      showDensity = input$showDensitySingleProbe,
      showMinima = input$showMinimaSingleProbe,
      numHistogramBins = input$numHistogramBinsOneProbe
    )
  })
  
  observeEvent(input$runProbeRandom, {
    req(input$betaFile$datapath)
    betas <- getBetas()
    randomProbeId <- rownames(betas)[sample(1:nrow(betas), 1)]

    updateTextInput(inputId = "probeId", value = randomProbeId)
    
    # probeUpdate(randomProbeId)
  })
  
  # observeEvent(input$runProbe, {
  #   probeUpdate(input$probeId)
  # })
  
  # Will set to TRUE when individual-probe visual is created
  plotCreatedProbeVisual <- reactiveVal(FALSE)
  
  output$plotCreatedProbeVisual <- reactive({
    plotCreatedProbeVisual()
  })
  
  outputOptions(output, "plotCreatedProbeVisual", suspendWhenHidden = FALSE)
  
  #### Set up graph and table after analyzing a single probe ####
  getSingleProbeSummary <- eventReactive(probeParams(), {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    # Set MethylModes parameters
    proportionSample <<- probeParams()$proportionSample
    peakDistance <<- probeParams()$peakDistance
    kernelType <<- KERNEL_TYPE
    bandwidthType <<- BANDWIDTH_TYPE
    numBreaks <- NUM_BREAKS
    densityAdjust <- input$densityAdjust
    pushToZero <- input$pushToZero

    probeId <- probeParams()$probeId
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
    
    probeId <- probeParams()$probeId
    rowId <- which(rownames(betas) == probeId)

    histData <- data.frame("beta" = betas[rowId,])
    
    # Colors chosen using the following code: 
    # hcl.colors(3, palette = 'viridis')
    # Add points for the peaks (maxima and minima) using the correctly structured data frames
    ggHist <- ggplot(histData, aes(x = beta)) +
      geom_histogram(aes(y = after_stat(density)), 
                     binwidth = 1 / probeParams()$numHistogramBins) +
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
          "Minima (peak boundaries)" = "#009B95"
        )
      )
    
    if (probeParams()$showDensity) {
      ggHist <- ggHist + 
        geom_line(data = fittedDensityDF, aes(x = x, y = y, color = "Density Estimate"))
    }
    if (probeParams()$showMinima) {
      ggHist <- ggHist + 
        geom_point(data = dataFrameForLeftMin, aes(x = beta, y = density, color = "Minima (peak boundaries)"), size = 3) +
        geom_point(data = dataFrameForRightMin, aes(x = beta, y = density, color = "Minima (peak boundaries)"), size = 3)
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
  
  # Reactive expression that responds to file upload
  getUploadedPeakSummary <- reactive({
    req(input$peakSummaryFile$datapath) # Ensure the file is uploaded
    peakSummary <- data.table::fread(input$peakSummaryFile$datapath)
    peakSummary <- peakSummary[order(peakSummary$probeName),]
    reset("runMultiProbe")
    
    # Users need to restart the app if they want to run a new analysis after uploading past results
    shinyjs::disable("runMultiProbe")
    shinyjs::disable("region")
    
    peakSummary
  })
  
  # Object that gets peakSummary from either upload or run-button
  selectedPeakSummary <- reactive({
    if (input$runMultiProbe == 0 && is.null(input$peakSummaryFile$datapath)) {
      return(NULL)  # Neither triggered initially
    }
    
    if (input$runMultiProbe > 0) {
      return(getMultiProbeSummary())
    } else {
      return(getUploadedPeakSummary())
    }
  })
  
  ##### Beta matrix-level analysis #####
  getMultiProbeSummary <- eventReactive(input$runMultiProbe, {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    req(betaFilter())

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
    print(paste("Running in parallel on", availableCores() - 1, "cores"))
    # Note that I had a bug here where I forgot to restrict the beta matrix
    # to the user-requested range! Fixed 6/27/24
    
    if (input$region == "wholeGenome") {
      peakSummary <- fillPeakSummaryParallel(betas)
    } else {
      peakSummary <- fillPeakSummaryParallel(betas[betaFilter(),])
    }
    
    # Sort results by number of detected peaks, descending
    # peakSummary <- peakSummary[order(peakSummary$numPeaks, decreasing = TRUE),]

    # Label invariant and hypo/hypermethylated CpG sites
    lowVariance <- logical(nrow(peakSummary))
    for (i in 1:nrow(peakSummary)) {
      if (peakSummary$numPeaks[i] > 1) {
        lowVariance[i] <- NA
      } else {
        lowVariance[i] <- peakSummary$peakVariance[i] < input$varianceThreshold
      }
    }
    
    hypoMethylated <- unlist(peakSummary$meanBeta) < input$hypoThreshold
    hyperMethylated <- unlist(peakSummary$meanBeta) > input$hyperThreshold
    
    peakSummary[, lowVariance := lowVariance]
    peakSummary[, hypoMethylated := hypoMethylated]
    peakSummary[, hyperMethylated := hyperMethylated]

    return(peakSummary)
  })
  
  observeEvent(input$runMultiProbe, {
    req(input$betaFile$datapath)
    
    getMultiProbeSummary()
    reset("runMultiProbe")
    shinyjs::enable("downloadPeakSummary")
  }) 
  
  # output$peakCountBar <- renderPlotly({
  #   req(input$runMultiProbe)
  #   peakSummary <- getMultiProbeSummary()
  #   if (is.null(peakSummary)) return()
  # 
  #   peakCounts <- data.frame(Peaks = peakSummary$numPeaks)
  # 
  #   p <- ggplot(peakCounts, aes(x = Peaks)) +
  #     geom_bar() +
  #     labs(title = "Number of Peaks Detected Per Probe", x = "Peaks", y = "Count")
  #   ggplotly(p)
  # })
  
  # Will set to TRUE when probe-average plot on "get started" page is created
  tableCreatedResultSummary <- reactiveVal(FALSE)
  
  # Updates according to the plotCreatedBetaOverview reactive value
  output$tableCreatedResultSummary <- reactive({
    tableCreatedResultSummary()
  })
  
  outputOptions(output, "tableCreatedResultSummary", suspendWhenHidden = FALSE)
  
  output$modalityTable <- renderTable(align = 'l', {
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()
 
    if (is.null(peakSummary)) return()
    
    counts <- peakSummary[, .N, by = numPeaks]
    counts$numPeaks <- as.integer(counts$numPeaks)
    
    colnames(counts) <- c("Peak Count/Modality", "Probes")

    # Set this to TRUE so title display is triggered
    tableCreatedResultSummary(TRUE)
    
    counts
  })
  
  output$flaggedProbesTableCounts <- renderTable(align = 'l', {
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()

    if (is.null(peakSummary)) return()
    
    varianceString = paste0(sum(peakSummary$lowVariance, na.rm = TRUE),
                            " (",
                            round(mean(peakSummary$lowVariance, na.rm = TRUE), 2) * 100,
                            "%)")
    hypoMethString = paste0(sum(peakSummary$hypoMethylated),
                            " (",
                            round(mean(peakSummary$hypoMethylated), 2) * 100,
                            "%)")
    hyperMethString = paste0(sum(peakSummary$hyperMethylated),
                             " (",
                             round(mean(peakSummary$hyperMethylated), 2) * 100,
                             "%)")
    
    
    results <- data.frame("Low Variance" = varianceString,
               "Hypomethylated" = hypoMethString,
               "Hypermethylated" = hyperMethString)

    colnames(results)[1] <- "Low-Variance"
    
    results
  })
  
  # output$flaggedProbesTableCounts2 <- renderTable(align = 'l', {
  #   req(selectedPeakSummary())
  #   peakSummary <- selectedPeakSummary()
  # 
  #   if (is.null(peakSummary)) return()
  #   
  #   nearCutoffPropSampleString = paste0(round(mean(peakSummary$nearCutoffPropSample, na.rm = TRUE), 2) * 100, "%")
  #   nearCutoffPeakDistanceString = paste0(round(mean(peakSummary$nearCutoffPeakDistance, na.rm = TRUE), 2) * 100, "%")
  #   
  #   results <- data.frame("Near proportionSample Threshold" = nearCutoffPropSampleString,
  #                         "Near peakDistance Threshold" = nearCutoffPeakDistanceString)
  #   colnames(results) <- c("Near proportionSample Threshold", "Near peakDistance Threshold")
  #   
  #   results
  # })
  
  # output$peakSummaryPreview <- renderPlot({
  #   peakSummary <- selectedPeakSummary()
  # 
  #   if (is.null(peakSummary)) return()
  #   
  #   betas <- getBetas()
  #   if (is.null(betas)) return()
  #   
  #   whichProbes <- sample(1:nrow(peakSummary), 9)
  #   par(mfrow = c(3,3))
  #   for (i in 1:9) {
  #     peakSummaryPlot(beta.data = betas[whichProbes[i],,drop = FALSE], 
  #                     peak.summary = peakSummary[whichProbes[i],])
  #   }
  #   
  # })

  output$peakSummaryTable <- DT::renderDataTable({
    peakSummary <- selectedPeakSummary()
    if (is.null(peakSummary)) return()
    
    listToString <- function(listData, decimalPlaces) {
      paste(as.character(round(unlist(listData), decimalPlaces)), collapse = ", ")
    }

    if (is.list(peakSummary$peakLocations)) {

      peakLocationsPreviewFriendly <- unlist(lapply(peakSummary$peakLocations, 
                                                    FUN = listToString,
                                                    decimalPlaces = 2))
      proportionSamplePreviewFriendly <- unlist(lapply(peakSummary$proportionSample, 
                                                       FUN = listToString,
                                                       decimalPlaces = 3))
      peakVariancePreviewFriendly <- unlist(lapply(peakSummary$peakVariance, 
                                                   FUN = listToString,
                                                   decimalPlaces = 6))
    } else {

      peakLocationsPreviewFriendly <- strsplit(peakSummary$peakLocations, split = "[|]")
      
      proportionSamplePreviewFriendly <- strsplit(peakSummary$proportionSample, split = "[|]")
      
      peakVariancePreviewFriendly <- strsplit(peakSummary$peakVariance, split = "[|]")
    }
    
    datatable(data.table("probeName" = peakSummary$probeName,
                         "numPeaks" = as.factor(peakSummary$numPeaks),
                         "meanBeta" = round(peakSummary$meanBeta, 2),
                         "peakLocations" = peakLocationsPreviewFriendly,
                         "proportionSample" = proportionSamplePreviewFriendly,
                         "peakVariance" = peakVariancePreviewFriendly),
              selection = 
                list(mode = "single",
                     target = 'row+column',
                     selected = list(
                       rows = which.max(peakSummary$numPeaks),
                       cols = c())),
              filter = "top",
              options = list(columnDefs = list(
                list(searchable = FALSE, targets = c(4, 5, 6)) 
              ),
              order = list(list(2, 'desc')))
    )
  })
  
  # Reactive expression to track the selected row
  multiProbeParams <- reactive({
    list(
      selectedRow = input$peakSummaryTable_rows_selected,
      showDensity = input$showDensityMultiProbe,
      showMinima = input$showMinimaMultiProbe
    )
  })
  
  getProbeVisualFromPeakSummary <- eventReactive(multiProbeParams(), {

    betas <- getBetas()
    if (is.null(betas)) return()
    
    betaFilter <- betaFilter()
    
    peakSummary <- selectedPeakSummary()

   # peakSummary <- getMultiProbeSummary()
    
    if (is.null(peakSummary)) return()
    
    req(selectedPeakSummary())
    # req(getMultiProbeSummary())
    
    # Here is where multiProbeParams() fails to update when peakSummary is created via file upload
    req(multiProbeParams()$selectedRow)  # Ensures that the reactive expression has a value
    
    beta.data <- betas[betaFilter,][multiProbeParams()$selectedRow, , drop = FALSE]
    peak.summary <- peakSummary[multiProbeParams()$selectedRow,]
    
    if (is.na(BANDWIDTH_TYPE)) {
      bandwidth = "nrd0"
    } else if (bandwidthType == "sheatherJones") {
      bandwidth = bw.SJ(betas[row.id,])
    }

    fittedDensity <- density(beta.data, from = 0, to = 1, n = numBreaks, 
                             adjust = densityAdjust, bw = bandwidth)

    if (is.list(peakSummary$peakLocations)) {
      detectedPeaks <- unlist(peak.summary$peakLocations)
  
      detectedMins <- c(unlist(peak.summary$leftMin), 
                        unlist(peak.summary$rightMin)[peak.summary$numPeaks])
    } else {
      detectedPeaks <- as.numeric(unlist(strsplit(peak.summary$peakLocations, split = "[|]")))
      
      detectedMins <- as.numeric(c(unlist(strsplit(peak.summary$leftMin, split = "[|]")),
                        unlist(strsplit(peak.summary$rightMin, split = "[|]")[peak.summary$numPeaks])))
    }
    
    # Place the detected peaks and minima on the fitted line
    fittedPeaks <- numeric(length(detectedPeaks))
    for (p in 1:length(detectedPeaks)) {
      closestIndex <- which.min(abs(fittedDensity$x - detectedPeaks[p]))
      fittedPeaks[p] <- fittedDensity$y[closestIndex]
    }
    
    fittedMins <- numeric(length(detectedMins)) 
    for (m in 1:length(detectedMins)) {
      closestIndex <- which.min(abs(fittedDensity$x - detectedMins[m]))
      fittedMins[m] <- fittedDensity$y[closestIndex]
    }
    
    fittedDensityDF <- data.frame(x = fittedDensity$x, y = fittedDensity$y)
    dataFrameForMaxima <- data.frame(beta = detectedPeaks,
                                     density = fittedPeaks)
    dataFrameForMinima <- data.frame(beta = detectedMins,
                                      density = fittedMins)
    
    histData <- data.frame("beta" = as.numeric(beta.data))
    
    # Colors chosen using the following code: 
    # hcl.colors(3, palette = 'viridis')
    # Add points for the peaks (maxima and minima) using the correctly structured data frames
    ggHist <- ggplot(histData, aes(x = beta)) +
      geom_histogram(aes(y = after_stat(density)), 
                     binwidth = 1 / probeParams()$numHistogramBins) +
      labs(title = paste("Beta distribution for Probe", peak.summary$probeName),
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
          "Minima (peak boundaries)" = "#009B95"
        )
      )

    if (multiProbeParams()$showDensity) {
      ggHist <- ggHist + 
        geom_line(data = fittedDensityDF, aes(x = x, y = y, color = "Density Estimate"))
    }
    if (multiProbeParams()$showMinima) {
      ggHist <- ggHist + 
        geom_point(data = dataFrameForMinima, aes(x = beta, y = density, color = "Minima (peak boundaries)"), size = 3)
    }
    
    ggHistPlotly <- ggplotly(ggHist, tooltip = "beta")

    # peakSummaryPlot(beta.data = betas[betaFilter,][multiProbeParams()$selectedRow, , drop = FALSE], 
    #                 peak.summary = peakSummary[multiProbeParams()$selectedRow,])
    ggHistPlotly
  })
  
  output$probeVisualFromPeakSummary <- renderPlotly({
    getProbeVisualFromPeakSummary()
  })
  
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
    annotationData <- getAnnotationLocal()
    if (is.null(annotationData)) return(NULL)
    
    allChromosomes <- c(as.character(1:22), "X", "Y")
    intersect(allChromosomes, unique(annotationData$Chromosome))
  })
  
  output$chromosomeSelectWhole <- renderUI({
    req(uniqueChromosomes()) 
    selectInput("selectedChromosomeWhole", "Choose Chromosome:",
                choices = uniqueChromosomes(),
                selected = uniqueChromosomes()[1])
  })
  
  output$chromosomeSelectPartial <- renderUI({
    req(uniqueChromosomes()) 
    selectInput("selectedChromosomePartial", "Choose Chromosome:",
                choices = uniqueChromosomes(),
                selected = uniqueChromosomes()[1])
  })
  
  # chromosomeToSubset <- reactive({
  #   input$selectedChromosomePartial
  # })
  # 
  # # Get minimum and maximum base pair locations from data
  # basePairMinMax <- eventReactive(chromosomeToSubset(), {
  #   annotationData <- getAnnotation()
  #   if (is.null(annotationData)) return(NULL)
  # 
  #   list("bpMin" = min(annotationData$Position[annotationData$Chromosome == chromosomeToSubset()]), 
  #        "bpMax" = max(annotationData$Position[annotationData$Chromosome == chromosomeToSubset()]))
  # })
  
  # Get minimum and maximum base pair locations from data
  basePairMinMax <- eventReactive(selectedRegion()$chromosomeToSubset, {
    annotationData <- getAnnotationLocal()
    if (is.null(annotationData)) return(NULL)
    
    list("bpMin" = min(annotationData$Position[annotationData$Chromosome == selectedRegion()$chromosomeToSubset]), 
         "bpMax" = max(annotationData$Position[annotationData$Chromosome == selectedRegion()$chromosomeToSubset]))
  })
  
  output$basePairRangeSelect <- renderUI({
    req(basePairMinMax())

    numericRangeInput("selectedBasePairRange", 
                      "Enter CpG location range",
                      value = c(basePairMinMax()$bpMin, basePairMinMax()$bpMax),
                      min = basePairMinMax()$bpMin,
                      max = basePairMinMax()$bpMax)
  })
  
  selectedRegion <- reactive({
    list(
      region = input$region,
      chromosomeWhole = input$selectedChromosomeWhole,
      chromosomeToSubset = input$selectedChromosomePartial,
      annotation = input$arrayType,
      existingFile = input$peakSummaryFile$datapath
    )
  })
  
  # Create a vector used to subset the beta matrix
  betaFilter <- eventReactive(selectedRegion(), {
    betas <- getBetas()
    if (is.null(betas)) return(NULL)

    if (!is.null(selectedRegion()$existingFile)) {
      peakSummaryProbes <- selectedPeakSummary()$probeName
      return(rownames(betas) %in% peakSummaryProbes)
    }
    
    annotationData <- getAnnotationLocal()
    
    if (is.null(annotationData)) return(rep(TRUE, nrow(betas)))
    
    if (selectedRegion()$region == "chromosome") {
      req(input$selectedChromosomeWhole)
      rowsToKeep <- annotationData$Chromosome == input$selectedChromosomeWhole
    } else if (selectedRegion()$region == "basePair") {
      req(input$selectedBasePairRange)
      rowsToKeep <- annotationData$Position >= input$selectedBasePairRange[1] & 
        annotationData$Position <= input$selectedBasePairRange[2] & 
        annotationData$Chromosome == input$selectedChromosomePartial
    } else {
      return(rep(TRUE, nrow(betas)))
    }

    rowsToKeep
  })
  
  output$subsetDataDimensions <- renderText({
    betas <- getBetas()
    if (is.null(betas)) return("")
    betaFilter <- betaFilter()

    annotationData <- getAnnotationLocal()
    if (is.null(annotationData)) simpleError("Annotation package required.")

    nProbe <- sum(betaFilter)
    paste(nProbe, "probes selected.")
  })
  
  # Cross-tabulate for multimodal tables (SNP under probe)
  output$tableMultimodalSNP <- DT::renderDataTable({
    req(selectedPeakSummary())
    req(getAnnotationLocal())
    peakSummary <- selectedPeakSummary()
    annotationData <- getAnnotationLocal()

    peakCounts <- sort(unique(peakSummary$numPeaks))
    
    # Apply betaFilter to ensure matching data lengths
    filter <- betaFilter()
    annotationData <- annotationData[filter,]

    if (input$arrayType == "il450k") {
      snpClose <- numeric(length(peakCounts))
      snpFar <- numeric(length(peakCounts))
      
      for (i in seq_along(peakCounts)) {
        snpClose[i] <- sum(peakSummary$numPeaks == peakCounts[i] 
                           & annotationData$SNP_within_10Bp != "")
        snpFar[i] <- sum(peakSummary$numPeaks == peakCounts[i] 
                         & annotationData$SNP_10Bp_and_beyond != "")
      }

      return(DT::datatable(data.frame(`Number of Modes` = peakCounts,
                                      `SNP within 10 bp` = snpClose,
                                      `SNP in 10-50 bp` = snpFar),
                           rownames = FALSE))
      
    } else if (input$arrayType == "ilepic1") {
      annotationData <- annotationData %>%
        mutate(SNP_distance_bin = cut(SNP_distance, breaks = c(0, 1, 10, 50), labels = c("1", "2-10", "11-50")))
      
      subset <- which(!is.na(annotationData$SNP_distance_bin))
      dataSubset <- peakSummary[subset, ]
      
      crossTab <- table(
        `Number of Modes Detected` = dataSubset$numPeaks,
        `SNP Distance Bin` = annotationData$SNP_distance_bin[subset]
      )
    }
    
    DT::datatable(as.data.frame.matrix(crossTab), rownames = TRUE)
  })
  
  # Cross-tabulate for multimodal tables (Relation to CpG island)
  output$tableMultimodalCpG <- DT::renderDataTable({
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()
    req(getAnnotationLocal())
    annotationData <- getAnnotationLocal()
    annotationData <- annotationData[betaFilter(),]
    
    crossTab <- table("Number of Modes Detected" = peakSummary$numPeaks, annotationData$Island)
    DT::datatable(as.data.frame.matrix(crossTab), rownames = TRUE)
  })
}