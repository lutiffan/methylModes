# TODO: let user set max file upload size via slider??
# create a display that automatically adjusts to the dimensions of your user's browser window
fluidPage(theme = shinythemes::shinytheme("cerulean"), 
  useShinyjs(), # Gives us tricks like disabling an element until condition satisfied
  navbarPage("MethylModes",
             tabPanel("Get Started",
                      sidebarPanel(
                        fileInput("betaFile", "Upload file containing beta values",
                                  multiple = FALSE,
                                  accept = c(".txt", ".csv", ".RDS", ".RDA")), 
                        helpText("Accepted file formats: .RDS, .RDA, .csv, .txt (tab-delimited), .tsv", br(), br(),
                                 "Max file size: 100 GB or local memory limit, 
                                 whichever is smaller"),
                        conditionalPanel(
                          condition = "output.missing450kAnnotation ||
                          output.missingEPICAnnotation ||
                          output.missingEPICV2Annotation",
                          h5("The following array annotation R packages are not 
                             installed:"),
                          conditionalPanel(
                            condition = "output.missing450kAnnotation",
                            h5("450k")
                          ),
                          conditionalPanel(
                            condition = "output.missingEPICAnnotation",
                            h5("EPIC")
                          ),
                          conditionalPanel(
                            condition = "output.missingEPICV2Annotation",
                            h5("EPIC V2")
                          ),
                          h5("Selecting a missing package below will 
                             begin installing it, which can take a long time.")
                        ),
                        radioButtons("arrayType", "Array type (Optional)",
                                     selected = character(0),
                                     c("450k" = "il450k",
                                       "EPIC v1.0" = "ilepic1",
                                       "EPIC v2.0" = "ilepic2")),
                        helpText("Loads array-specific annotations (e.g. chromosome, 
                                 base pair, etc.. Required for running MethylModes on 
                                 a subset of data.)")
                      ),
                      mainPanel(
                        conditionalPanel(
                          condition = "output.plotCreatedBetaOverview === false",
                          h4("Data summary will be displayed upon successful file upload.")
                        ),
                        # Histogram of all probes
                        textOutput(outputId = "wholeDataDimensions"),
                        plotlyOutput(outputId = "betaOverview"),
                        conditionalPanel(
                          condition = "output.plotCreatedBetaOverview",
                          sliderInput(inputId = "numHistogramBinsBetaOverview", 
                                      "Number of histogram bins", 
                                      min = 10, max = 100,
                                      value = 50)
                        ),
                        withSpinner(plotlyOutput(outputId = "chromosomeBar")),
                        plotlyOutput(outputId = "islandBar")
                      )
              ),
             #### Run MethylModes ####
             tabPanel("Run MethylModes",
                      sidebarPanel(
                        radioButtons(
                          "analysisType",
                          "Select Analysis Type:",
                          choices = c("Individual Probe" = "individual", 
                                      "Multiple Probe" = "multiProbe")
                        ),
                        h4("Peak Detection Thresholds"),
                        numericInput(label = "ProportionSample",
                                     inputId = "proportionSample",
                                     value = proportionSample,
                                     min = 0.00001,
                                     max = 0.99999),
                        helpText("Minimum proportion of sample considered to be a peak"),
                        numericInput(label = "PeakDistance",
                                     inputId = "peakDistance",
                                     value = peakDistance,
                                     min = 0.00001,
                                     max = 0.99999),
                        helpText("Minimum distance between adjacent peaks"),
                        # h4("Smoothing parameters"),
                        # numericInput(label = "DensityAdjust",
                        #              inputId = "densityAdjust",
                        #              value = 1.5),
                        # numericInput(label = "Epsilon",
                        #              inputId = "pushToZero",
                        #              value = 1e-6),
                        # helpText("Default values are recommended for smoothing parameters."),
                        # helpText("The density() function fits very small floating-point",
                        #          "values to the data, creating tiny perturbations that ",
                        #          "reduce the accuracy of MethylModes. Setting a threshold ",
                        #          "under which small values are considered equivalent ",
                        #          "to zero mitigates this issue."),
                        #### Individual probe options ####
                        conditionalPanel(
                          condition = "input.analysisType == 'individual'",
                          textInput("probeId", "Enter probe ID", value = ""),
                          # actionButton("runProbe", "Run on Selected Probe", 
                          #              disabled = TRUE),
                          div(style = "margin-bottom: 10px;"),
                          h5("Alternatively, you can..."),
                          actionButton("runProbeRandom", "View a random probe", disabled = TRUE),
                          div(style = "margin-bottom: 10px;")
                        ),
                        #### Multi-probe options ####
                        conditionalPanel(
                          condition = "input.analysisType == 'multiProbe'",
                          radioButtons("region", "Select region", 
                                       choices = c("Whole genome" = "wholeGenome",
                                                   "Chromosome*" = "chromosome",
                                                   "Base pair range*" = "basePair")),
                          conditionalPanel(
                            condition = "input.region == 'chromosome'",
                            uiOutput("chromosomeSelectWhole")
                          ),
                          conditionalPanel(
                            condition = "input.region == 'basePair'",
                            uiOutput("chromosomeSelectPartial"),
                            uiOutput("basePairRangeSelect")
                          ),
                          conditionalPanel(
                            condition = "input.region == 'chromosome' ||
                            input.region == 'basePair'",
                            textOutput(outputId = "subsetDataDimensions"),
                          ),
                          helpText("*Requires that annotations for the array type are loaded (select on 'Get Started') page."),
                          h4("Optional Thresholds"),
                          helpText("Label low-variance, hypomethylated, and hypermethylated CpG sites"),
                          numericInput(inputId = "varianceThreshold",
                                       label = "Variance",
                                       value = 1e-5,
                                       min = 0,
                                       max = 1),
                          numericInput(inputId = "hypoThreshold",
                                       label = "Hypomethylation",
                                       value = 0.3,
                                       min = 0,
                                       max = 0.5),
                          numericInput(input = "hyperThreshold",
                                       label = "Hypermethylation",
                                       value = 0.7,
                                       min = 0.5,
                                       max = 1),
                          actionButton("runMultiProbe", "Run multiple-probe analysis", disabled = TRUE)
                        )
                      ), # sidebarPanel
                      mainPanel(
                        #### Individual probe display ####
                        conditionalPanel(
                          condition = "input.analysisType == 'individual'",
                          h3("Individual Probe Analysis"),
                          checkboxInput("showDensitySingleProbe", 
                                        label = "Display density estimate curve",
                                        value = TRUE),
                          checkboxInput("showMinimaSingleProbe",
                                        label = "Display detected peak boundaries",
                                        value = FALSE),
                          withSpinner(plotlyOutput(outputId = "probeVisual")),
                          conditionalPanel(
                            condition = "output.plotCreatedProbeVisual",
                            sliderInput(inputId = "numHistogramBinsOneProbe", 
                                        "Number of histogram bins", 
                                        min = 10, max = 100,
                                        value = 50)
                          ),
                          tableOutput('probeTable')
                        ),
                        #### Multi-probe display ####
                        conditionalPanel(
                          condition = "input.analysisType == 'multiProbe'",
                          h3("Multiple Probe Analysis"),
                          shinyjs::disabled(downloadButton("downloadPeakSummary",
                                                           "Download results")),
                          checkboxInput("showDensityMultiProbe", 
                                        label = "Display density estimate curve",
                                        value = TRUE),
                          checkboxInput("showMinimaMultiProbe",
                                        label = "Display detected peak boundaries",
                                        value = FALSE),
                          conditionalPanel(
                            condition = "output.tableCreatedResultSummary",
                            h4("Result Summary")
                          ),
                          tableOutput("modalityTable"),
                          tableOutput("flaggedProbesTableCounts"),
                          plotlyOutput("probeVisualFromPeakSummary"),
                          DT::dataTableOutput("peakSummaryTable")
                        )
                        )
                      ), 
             tabPanel("Review and Analyze Results",
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
