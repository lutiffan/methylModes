function(input, output) {
  # Arbitrarily set maximum file upload size to 100 GB
  # The real limit comes from the user's local memory limit
  # Is there a simple way to set this to max memory based on local machine?
  options(shiny.maxRequestSize = MAX_FILE_SIZE)
  
  ##### Get Started #####
  
  # Implicitly returns betas as a reactive object
  getBetas <- reactive({
    req(input$betaFile$datapath)
    
    fileType <- tools::file_ext(input$betaFile$datapath)
    
    if (fileType %in% c("RDS", "rds")) {
      betas <- readRDS(input$betaFile$datapath) 
    } else if (fileType %in% c("csv", "TXT", "txt", "tsv")) {
      warning("Converting to matrix format. 
              Interpreting first column as row names.")
      betas <- as.matrix(fread(input$betaFile$datapath), rownames = 1)
    } else if (fileType %in% c("RDA", "rda")) {
      objname <- load(input$betaFile$datapath)
      betas <- get(objname)
    } else if (fileType == "qs") {
      betas <- qread(file = input$betaFile$datapath)
    } else {
      warning("Invalid file type.")
    }
    
    if(!is.matrix(betas)) {
      warning("Converting data to matrix format.")
      betas <- as.matrix(betas, rownames = 1)
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

  
  getAnnotationLocal <- reactive({
    if (is.null(input$arrayType)) return(NULL)
    betas <- getBetas()
    if (is.null(betas)) return(NULL)
    
    if (input$arrayType == "il450k") {
      # Load data from RDA file
      data("IlluminaManifest450k", package = "methylModes", envir = environment())
      manifestFile <- IlluminaManifest450k
      
    } else if (input$arrayType == "ilepic1") {
      # Load data from RDA file
      data("IlluminaManifestEPIC", package = "methylModes", envir = environment())
      manifestFile <- IlluminaManifestEPIC
      
    } else if (input$arrayType == "ilepic2") {
      data("IlluminaManifestEPICv2", package = "methylModes", envir = environment())
      manifestFile <- IlluminaManifestEPICv2
    } else {
      simpleError("Invalid annotation package selected.")
    } 
    
    # Keep only rows corresponding to probes that we have
    manifestFile <- manifestFile[manifestFile$IlmnID %in% rownames(betas),]
    # Now sort the probes in order of name
    manifestFile <- manifestFile[order(manifestFile$IlmnID),]

    annotationData <- data.frame(Chromosome = manifestFile$CHR,
                                 Island = manifestFile$Relation_to_UCSC_CpG_Island,
                                 Position = manifestFile$MAPINFO,
                                 IlmnID = manifestFile$IlmnID)
    
    if (input$arrayType == "il450k") {
      annotationData$SNP_within_10Bp <- manifestFile$Probe_SNPs_10
      annotationData$SNP_10Bp_and_beyond <- manifestFile$Probe_SNPs
    } else { # if (input$arrayType == "ilepic1") {
      annotationData$snpDistance0_1 <- manifestFile$snpDistance0_1
      annotationData$snpDistance2_10 <- manifestFile$snpDistance2_10
      annotationData$snpDistance11_50 <- manifestFile$snpDistance11_50
    } 
    
    # Replace empty strings in Island column with "OpenSea"
    annotationData$Island[annotationData$Island == ""] <- "OpenSea"
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
                     binwidth = 1 / input$numHistogramBinsBetaOverview, 
                     color = "black", 
                     fill = "#2FA4E7") +
      geom_line(aes(y = after_stat(density)), stat = 'density', color = "black") +
      labs(title = "Mean Methylation Across Probes",
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
    
    relevantLabels <- getAnnotationLocal()
    if (is.null(relevantLabels)) return(NULL)
    
    relevantLabels$Chromosome <- factor(relevantLabels$Chromosome,
                                        levels = c(as.character(1:22), "X", "Y"))
    
    p <- ggplot(relevantLabels, aes(x = Chromosome)) +
      geom_bar(color = "black", fill = "#2FA4E7") +  # Set the bar fill color to cerulean
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

    p <- ggplot(relevantLabels, aes(x = Island, fill = Island)) +
      geom_bar(color = "black") +  # Add black border to the bars
      scale_fill_manual(values = ISLAND_COLORS) +  # Assign the custom colors
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

    # Get parameters from UI inputs
    params <- list(
      proportionSample = probeParams()$proportionSample,
      peakDistance = probeParams()$peakDistance,
      kernelType = FIXED_KERNEL_TYPE,
      bandwidthType = FIXED_BANDWIDTH_TYPE,
      numBreaks = FIXED_NUM_BREAKS,
      densityAdjust = DEFAULT_DENSITY_ADJUST,
      pushToZero = DEFAULT_PUSH_TO_ZERO
    )

    probeId <- probeParams()$probeId
    rowId <- which(rownames(betas) == probeId)

    methylModes::methylModes(row.data = betas[rowId,],
                proportionSample = params$proportionSample,
                peakDistance = params$peakDistance,
                kernelType = params$kernelType,
                bandwidthType = params$bandwidthType,
                numBreaks = params$numBreaks,
                densityAdjust = params$densityAdjust,
                pushToZero = params$pushToZero)
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
      geom_point(data = dataFrameForMaxima, aes(x = beta, y = density, color = "Maxima"), size = 3, alpha = 0.6) +
      scale_color_manual(
        name = "Legend",
        values = DENSITY_PLOT_COLORS
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
    
    fileType <- tools::file_ext(input$peakSummaryFile$datapath)

    if (fileType %in% c("RDS", "rds")) {
      peakSummary <- readRDS(input$peakSummaryFile$datapath) 
    } else if (fileType %in% c("csv", "TXT", "txt", "tsv")) {
      peakSummary <- fread(input$peakSummaryFile$datapath)
    } else if (fileType %in% c("RDA", "rda")) {
      objname <- load(input$peakSummaryFile$datapath)
      peakSummary <- get(objname)
    } else if (fileType == "qs") {
      peakSummary <- qread(file = input$peakSummaryFile$datapath)
    } else {
      warning("Invalid file type.")
    }

    peakSummary <- peakSummary[order(peakSummary$probeName),]
    shinyjs::reset("runMultiProbe")
    
    # Users need to restart the app if they want to run a new analysis after uploading past results
    shinyjs::disable("runMultiProbe")
    shinyjs::disable("region")
    
    # Check for threshold columns in uploaded results
    existingLowVarCols <- grep("^lowVariance_", colnames(peakSummary), value = TRUE)
    existingHypoCols <- grep("^hypoMethylated_", colnames(peakSummary), value = TRUE)
    existingHyperCols <- grep("^hyperMethylated_", colnames(peakSummary), value = TRUE)

    if (length(existingLowVarCols) > 0 && length(existingHypoCols) > 0 && length(existingHyperCols) > 0) {
      # Extract thresholds from column names
      uploadedVarianceThreshold <- sub("^lowVariance_", "", existingLowVarCols[1])
      uploadedHypoThreshold <- sub("^hypoMethylated_", "", existingHypoCols[1])
      uploadedHyperThreshold <- sub("^hyperMethylated_", "", existingHyperCols[1])

    } else {
      # If columns are missing, add them for the current thresholds
      peakSummary <- peakSummaryPostProcessing(peakSummary, input$varianceThreshold, input$hypoThreshold, input$hyperThreshold)
    }
    
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
  
peakSummaryPostProcessing <- function(peakSummary, varianceThreshold, hypoThreshold, hyperThreshold) {
  # Create dynamic column names
  lowVarCol <- paste0("lowVariance_", varianceThreshold)
  hypoCol <- paste0("hypoMethylated_", hypoThreshold)
  hyperCol <- paste0("hyperMethylated_", hyperThreshold)

  # Calculate values
  lowVariance <- ifelse(
    peakSummary$numPeaks > 1,
    NA,
    unlist(peakSummary$peakVariance) < varianceThreshold
  )
  hypoMethylated <- unlist(peakSummary$meanBeta) < hypoThreshold
  hyperMethylated <- unlist(peakSummary$meanBeta) > hyperThreshold 

  # Assign with dynamic names
  peakSummary[, (lowVarCol) := lowVariance]
  peakSummary[, (hypoCol) := hypoMethylated]
  peakSummary[, (hyperCol) := hyperMethylated]
  
  peakSummary
}

  ##### Beta matrix-level analysis #####
  getMultiProbeSummary <- eventReactive(input$runMultiProbe, {
    betas <- getBetas()
    if (is.null(betas)) return()
    
    req(betaFilter())

    # Get parameters from UI inputs
    params <- list(
      proportionSample = input$proportionSample,
      peakDistance = input$peakDistance,
      kernelType = FIXED_KERNEL_TYPE,
      bandwidthType = FIXED_BANDWIDTH_TYPE,
      numBreaks = FIXED_NUM_BREAKS,
      densityAdjust = DEFAULT_DENSITY_ADJUST,
      pushToZero = DEFAULT_PUSH_TO_ZERO
    )
    
    # totalRows = rangeEnd - rangeStart + 1
    print(paste("Running in parallel on", availableCores() - 1, "cores"))
    # Note that I had a bug here where I forgot to restrict the beta matrix
    # to the user-requested range! Fixed 6/27/24
    
    if (input$region == "wholeGenome") {
      peakSummary <- methylModesBatch(betas,
                                      proportionSample = params$proportionSample,
                                      peakDistance = params$peakDistance,
                                      kernelType = params$kernelType,
                                      bandwidthType = params$bandwidthType,
                                      numBreaks = params$numBreaks,
                                      densityAdjust = params$densityAdjust,
                                      pushToZero = params$pushToZero)
    } else {
      peakSummary <- methylModesBatch(betas[betaFilter(),],
                                      proportionSample = params$proportionSample,
                                      peakDistance = params$peakDistance,
                                      kernelType = params$kernelType,
                                      bandwidthType = params$bandwidthType,
                                      numBreaks = params$numBreaks,
                                      densityAdjust = params$densityAdjust,
                                      pushToZero = params$pushToZero)
    }

    # Sort results by number of detected peaks, descending
    # peakSummary <- peakSummary[order(peakSummary$numPeaks, decreasing = TRUE),]

    peakSummaryPostProcessing(peakSummary, input$varianceThreshold, input$hypoThreshold, input$hyperThreshold)
  })
  
  observeEvent(input$runMultiProbe, {
    req(input$betaFile$datapath)
    
    getMultiProbeSummary()
    shinyjs::reset("runMultiProbe")
    shinyjs::enable("downloadPeakSummary")
  }) 
  
  # output$peakCountBar <- plotly::renderPlotly({
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

    # Find the actual columns present in the data
    lowVarCol <- grep("^lowVariance_", colnames(peakSummary), value = TRUE)
    hypoCol <- grep("^hypoMethylated_", colnames(peakSummary), value = TRUE)
    hyperCol <- grep("^hyperMethylated_", colnames(peakSummary), value = TRUE)

    # Extract thresholds from column names
    extract_threshold <- function(colname, prefix) sub(paste0("^", prefix, "_"), "", colname)
    varianceThreshold <- if (length(lowVarCol)) extract_threshold(lowVarCol[1], "lowVariance") else NA
    hypoThreshold <- if (length(hypoCol)) extract_threshold(hypoCol[1], "hypoMethylated") else NA
    hyperThreshold <- if (length(hyperCol)) extract_threshold(hyperCol[1], "hyperMethylated") else NA

    # Use the found columns for summary
    varianceString = if (length(lowVarCol)) paste0(sum(peakSummary[[lowVarCol]], na.rm = TRUE),
                                                 " (",
                                                 round(mean(peakSummary[[lowVarCol]], na.rm = TRUE), 2) * 100,
                                                 "%)") else "N/A"
    hypoMethString = if (length(hypoCol)) paste0(sum(peakSummary[[hypoCol]]),
                                               " (",
                                               round(mean(peakSummary[[hypoCol]]), 2) * 100,
                                               "%)") else "N/A"
    hyperMethString = if (length(hyperCol)) paste0(sum(peakSummary[[hyperCol]]),
                                                 " (",
                                                 round(mean(peakSummary[[hyperCol]]), 2) * 100,
                                                 "%)") else "N/A"

    # Set column names with thresholds
    colnames_with_thresholds <- c(
      paste0("Low-Variance (Unimodal and variance < ", varianceThreshold, ")"),
      paste0("Hypomethylated (Mean beta < ", hypoThreshold, ")"),
      paste0("Hypermethylated (Mean beta > ", hyperThreshold, ")")
    )

    results <- data.frame(varianceString, hypoMethString, hyperMethString)
    colnames(results) <- colnames_with_thresholds
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

  output$peakSummaryTable <- renderDataTable({
    peakSummary <- selectedPeakSummary()
    if (is.null(peakSummary)) return()
    
    listToString <- function(listData, decimalPlaces = 2) {
      paste(format(round(unlist(listData), decimalPlaces), nsmall = 2), collapse = ", ")
    }
    
    unpackExcelColumn <- function(listData, decimalPlaces = 2) {
      paste(format(round(as.numeric(unlist(strsplit(listData, 
                                                    split = "[|]"))), 
                         decimalPlaces), 
                   nsmall = 2), 
            collapse = ", ")
    }

    if (is.list(peakSummary$peakLocations)) { # After running methylModesBatch

      peakLocationsPreviewFriendly <- unlist(lapply(peakSummary$peakLocations, 
                                                    FUN = listToString,
                                                    decimalPlaces = 2))
      proportionSamplePreviewFriendly <- unlist(lapply(peakSummary$proportionSample, 
                                                       FUN = listToString,
                                                       decimalPlaces = 2))
      peakVariancePreviewFriendly <- unlist(lapply(peakSummary$peakVariance, 
                                                   FUN = listToString,
                                                   decimalPlaces = 6))
    } else { # When you've read in existing results
      peakLocationsPreviewFriendly <- unlist(lapply(peakSummary$peakLocations,
                                                    FUN = unpackExcelColumn,
                                                    decimalPlaces = 2))
      proportionSamplePreviewFriendly <- unlist(lapply(peakSummary$proportionSample,
                                                       FUN = unpackExcelColumn,
                                                       decimalPlaces = 2))
      peakVariancePreviewFriendly <- unlist(lapply(peakSummary$peakVariance,
                                                   FUN = unpackExcelColumn,
                                                   decimalPlaces = 6))
      # peakLocationsPreviewFriendly <- round(as.numeric(unlist(strsplit(peakSummary$peakLocations, split = "[|]"))), 2)
      # 
      # proportionSamplePreviewFriendly <- round(as.numeric(unlist(strsplit(peakSummary$proportionSample, split = "[|]"))), 2)
      # 
      # peakVariancePreviewFriendly <- round(as.numeric(unlist(strsplit(peakSummary$peakVariance, split = "[|]"))), 6)
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
              order = list(list(2, 'desc')),
              server = TRUE)
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
    
    # Use fixed bandwidth type
    if (is.na(FIXED_BANDWIDTH_TYPE)) {
      bandwidth = "nrd0"
    } else {
      bandwidth = FIXED_BANDWIDTH_TYPE
    }

    fittedDensity <- density(beta.data, from = 0, to = 1, n = FIXED_NUM_BREAKS, 
                             adjust = DEFAULT_DENSITY_ADJUST, bw = bandwidth)

    if (is.list(peak.summary$peakLocations)) {
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
        values = DENSITY_PLOT_COLORS
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
      paste("methylModes_PS_", 
            gsub("\\.", "_", input$proportionSample), "_", 
            "PD_",
            gsub("\\.", "_", input$peakDistance), "_",
            "VT_",
            gsub("\\.", "_", input$varianceThreshold), "_",
            "HT_",
            gsub("\\.", "_", input$hyperThreshold), "_",
            Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Get the peak summary data
      peakSummary <- getMultiProbeSummary()
      
      # Add annotation info columns if available
      annotationData <- getAnnotationLocal()
      if (!is.null(annotationData)) {
        # Convert to data.table for consistent operations
        peakSummary <- as.data.table(peakSummary)
        annotationData <- as.data.table(annotationData)
        
        # Select only the columns we need
        annotationSubset <- annotationData[, .(IlmnID, Island)]
        
        # Use data.table merge instead of dplyr left_join
        peakSummary <- merge(peakSummary, annotationSubset, 
                            by.x = "probeName", by.y = "IlmnID", 
                            all.x = TRUE)
      }
      
      # Write to file - fwrite should handle list columns automatically
      fwrite(peakSummary, file)
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

    shinyWidgets::numericRangeInput("selectedBasePairRange", 
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
  
  # Cross-tabulate multimodality (SNP under probe)
  output$tableMultimodalSNP <- renderDataTable({
    req(selectedPeakSummary())
    req(getAnnotationLocal())
    peakSummary <- selectedPeakSummary()
    annotationData <- getAnnotationLocal()

    peakCounts <- sort(unique(peakSummary$numPeaks))

    # Apply betaFilter to ensure matching data lengths
    filter <- betaFilter()
    annotationData <- annotationData[betaFilter(),]

    if (input$arrayType == "il450k") {
      snpClose <- numeric(length(peakCounts))
      snpFar <- numeric(length(peakCounts))

      for (i in seq_along(peakCounts)) {
        snpClose[i] <- sum(peakSummary$numPeaks == peakCounts[i]
                           & annotationData$SNP_within_10Bp != "")
        snpFar[i] <- sum(peakSummary$numPeaks == peakCounts[i]
                         & annotationData$SNP_10Bp_and_beyond != "")
      }
      
      crossTab <- data.frame(`Number of Modes` = peakCounts,
                             `SNP within 10 bp` = snpClose,
                             `SNP in 10-50 bp` = snpFar)
      colnames(crossTab) <- c("Number of Modes",
                              "SNP within 10 bp",
                              "SNP in 10-50 bp")

    } else { # if (input$arrayType == "ilepic1") {
      snp1_10 <- numeric(length(peakCounts))
      snp2_10 <- numeric(length(peakCounts))
      snp11_50 <- numeric(length(peakCounts))

      for (i in seq_along(peakCounts)) {
        snp1_10[i] <- sum(annotationData$snpDistance0_1[peakSummary$numPeaks == peakCounts[i]], na.rm = T)
        snp2_10[i] <- sum(annotationData$snpDistance2_10[peakSummary$numPeaks == peakCounts[i]], na.rm = T)
        snp11_50[i] <- sum(annotationData$snpDistance11_50[peakSummary$numPeaks == peakCounts[i]], na.rm = T)
      }

      crossTab <- data.frame(`Number of Modes` = peakCounts,
                             `SNP at or neighboring CpG site` = snp1_10,
                             `SNP between 2-10 bp from CpG` = snp2_10,
                             `SNP between 11-50 bp from CpG` = snp11_50)
      colnames(crossTab) <- c("Number of Modes", 
                              "SNP at or neighboring CpG site", 
                              "SNP between 2-10 bp from CpG", 
                              "SNP between 11-50 bp from CpG")
    }
    datatable(crossTab,
              rownames = FALSE,
              escape = FALSE)
  })
  
  # Cross-tabulate for multimodality (Relation to CpG island)
  output$tableMultimodalCpG <- renderDataTable({
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()
    req(getAnnotationLocal())
    annotationData <- getAnnotationLocal()
    annotationData <- annotationData[betaFilter(),]
    
    crossTab <- table("Number of Modes Detected" = peakSummary$numPeaks, annotationData$Island)
    datatable(as.data.frame.matrix(crossTab), rownames = TRUE)
  })

  # Cross-tabulate hypo- and hypermethylation status (SNP under probe)
  output$tableHypoHyperSNP <- renderDataTable({
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()
    req(getAnnotationLocal())
    annotationData <- getAnnotationLocal()
    annotationData <- annotationData[betaFilter(),]

    # Create a data table counting probes with hypo- and hypermethylation status along rows
    # and SNP under probe along columns
    if (input$arrayType == "il450k") {
      snpClose <- numeric(2)
      snpFar <- numeric(2)
      snpClose[1] <- sum(peakSummary$hypoMethylated & annotationData$SNP_within_10Bp != "")
      snpClose[2] <- sum(peakSummary$hyperMethylated & annotationData$SNP_within_10Bp != "")
      snpFar[1] <- sum(peakSummary$hypoMethylated & annotationData$SNP_10Bp_and_beyond != "")
      snpFar[2] <- sum(peakSummary$hyperMethylated & annotationData$SNP_10Bp_and_beyond != "")

      crossTab <- data.frame(`status` = c("Hypo", "Hyper"),
                             `SNP within 10 bp` = snpClose,
                             `SNP in 10-50 bp` = snpFar)
      colnames(crossTab) <- c("Hypo/hypermethylation",
                              "SNP within 10 bp",
                              "SNP in 10-50 bp")
      
    } else {
      snp1_10 <- numeric(2)
      snp2_10 <- numeric(2)
      snp11_50 <- numeric(2)
      snp1_10[1] <- sum(annotationData$snpDistance0_1[peakSummary$hypoMethylated], na.rm = T)
      snp1_10[2] <- sum(annotationData$snpDistance0_1[peakSummary$hyperMethylated], na.rm = T)
      snp2_10[1] <- sum(annotationData$snpDistance2_10[peakSummary$hypoMethylated], na.rm = T)
      snp2_10[2] <- sum(annotationData$snpDistance2_10[peakSummary$hyperMethylated], na.rm = T)
      snp11_50[1] <- sum(annotationData$snpDistance11_50[peakSummary$hypoMethylated], na.rm = T)
      snp11_50[2] <- sum(annotationData$snpDistance11_50[peakSummary$hyperMethylated], na.rm = T)

      crossTab <- data.frame(`status` = c("Hypo", "Hyper"),
                             `SNP at or neighboring CpG site` = snp1_10,
                             `SNP between 2-10 bp from CpG` = snp2_10,
                             `SNP between 11-50 bp from CpG` = snp11_50)
      colnames(crossTab) <- c("Hypo/hypermethylation",
                              "SNP at or neighboring CpG site",
                              "SNP between 2-10 bp from CpG",
                              "SNP between 11-50 bp from CpG")
    }
    datatable(crossTab,
              rownames = FALSE,
              escape = FALSE)
  })

  # Cross-tabulate hypo- and hypermethylation status (Relation to CpG island)
  output$tableHypoHyperCpG <- renderDataTable({
    req(selectedPeakSummary())
    peakSummary <- selectedPeakSummary()
    req(getAnnotationLocal())
    annotationData <- getAnnotationLocal()
    annotationData <- annotationData[betaFilter(),]

  # For every value of annotationData$Island, count the number of probes that are hypo- or hypermethylated
  islandLabels <- unique(annotationData$Island)
  hypo <- numeric(length(islandLabels))
  hyper <- numeric(length(islandLabels))
  for (i in seq_along(islandLabels)) {
    hypo[i] <- sum(peakSummary$hypoMethylated & annotationData$Island == islandLabels[i], na.rm = T)
    hyper[i] <- sum(peakSummary$hyperMethylated & annotationData$Island == islandLabels[i], na.rm = T)
  }

  crossTab <- matrix(c(hypo, hyper), nrow = 2, byrow = TRUE)
  colnames(crossTab) <- islandLabels
  rownames(crossTab) <- c("Hypo", "Hyper")
    datatable(as.data.frame.matrix(crossTab), rownames = TRUE)
  })
}

