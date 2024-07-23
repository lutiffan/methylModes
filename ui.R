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
                        h4("Descriptive Thresholds for Peak Detection"),
                        numericInput(label = "ProportionSample",
                                     inputId = "proportionSample",
                                     value = 0.05),
                        helpText("Minimum proportion of sample considered to be a peak"),
                        numericInput(label = "PeakDistance",
                                     inputId = "peakDistance",
                                     value = 0.1),
                        helpText("Minimum distance between adjacent peaks"),
                        h4("Smoothing parameters"),
                        numericInput(label = "DensityAdjust",
                                     inputId = "densityAdjust",
                                     value = 1.5),
                        numericInput(label = "Epsilon",
                                     inputId = "pushToZero",
                                     value = 1e-6),
                        helpText("Default values are recommended for smoothing parameters."),
                        # helpText("The density() function fits very small floating-point",
                        #          "values to the data, creating tiny perturbations that ",
                        #          "reduce the accuracy of MethylModes. Setting a threshold ",
                        #          "under which small values are considered equivalent ",
                        #          "to zero mitigates this issue."),
                        #### Individual probe options ####
                        conditionalPanel(
                          condition = "input.analysisType == 'individual'",
                          textInput("probeId", "Probe ID", value = ""),
                          actionButton("runProbe", "Run on Selected Probe", 
                                       disabled = TRUE),
                          div(style = "margin-bottom: 10px;"),
                          actionButton("runProbeRandom", "Run on Randomly 
                                       Selected Probe", disabled = TRUE),
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
                            uiOutput("chromosomeSelect")
                          ),
                          conditionalPanel(
                            condition = "input.region == 'basePair'",
                            uiOutput("basePairRangeSelect")
                          ),
                          conditionalPanel(
                            condition = "input.region == 'chromosome' ||
                            input.region == 'basePair'",
                            textOutput(outputId = "subsetDataDimensions"),
                          ),
                          helpText("*Requires that annotations for the array type are loaded (select on 'Get Started') page."),
                          actionButton("runMultiProbe", "Run Multiple Probe Analysis", disabled = TRUE)
                          # numericInput("rangeStart", "Range Start", value = NULL),
                          # numericInput("rangeEnd", "Range End", value = NULL),
                        )
                      ), # sidebarPanel
                      mainPanel(
                        #### Individual probe display ####
                        conditionalPanel(
                          condition = "input.analysisType == 'individual'",
                          h4("Individual Probe Analysis"),
                          withSpinner(plotlyOutput(outputId = "probeVisual")),
                          conditionalPanel(
                            condition = "output.plotCreatedProbeVisual",
                            sliderInput(inputId = "numHistogramBinsOneProbe", 
                                        "Number of histogram bins", 
                                        min = 10, max = 100,
                                        value = 50)
                          ),
                          tableOutput('probeTable')
                          # plotOutput(outputId = "probeVisualBaseR")
                        ),
                        #### Multi-probe display ####
                        conditionalPanel(
                          condition = "input.analysisType == 'multiProbe'",
                          h4("Multiple Probe Analysis"),
                          shinyjs::disabled(downloadButton("downloadPeakSummary",
                                                           "Download Results")),
                          withSpinner(plotlyOutput("peakCountBar")),
                          # plotOutput("peakSummaryPreview"),
                          DT::dataTableOutput("peakSummaryTable")
                        )
                        )
                      ), # Run MethylModes
             # tabPanel("Run MethylModes",
             #          sidebarPanel(
             #            radioButtons(
             #              "analysisType",
             #              "Select Analysis Type:",
             #              choices = c("Individual Probe" = "individual", 
             #                          "Whole Genome" = "multiProbe")
             #            ),
             #            h4("Hyperparameters for Probe-Level Analysis"),
             #            numericInput(label = "proportionSample",
             #                         inputId = "proportionSampleProbe",
             #                         value = 0.05),
             #            helpText("Minimum proportion of sample considered to be a peak"),
             #            numericInput(label = "peakDistance",
             #                         inputId = "peakDistanceProbe",
             #                         value = 0.1),
             #            helpText("Minimum distance between adjacent peaks"),
             #            h4("Smoothing parameters"),
             #            numericInput(inputId = "densityAdjustProbe", "density() 'adjust' parameter",
             #                         value = 1.5),
             #            numericInput(inputId = "pushToZeroProbe",
             #                         "Threshold for numbers small enough to be set to zero",
             #                         value = 1e-6),
             #            helpText("The density() function fits very small floating-point",
             #                     "values to the data, creating tiny perturbations that ",
             #                     "reduce the accuracy of MethylModes. Setting a threshold ",
             #                     "under which small values are considered equivalent ",
             #                     "to zero mitigates this issue."),
             #            uiOutput("conditionalSidebarButtons")
             #          ), # sidebarPanel
             #          mainPanel(
             #            uiOutput("conditionalPlots")
             #          )
             # ), # Run MethylModes
             # tabPanel("Individual Probe Analysis",
             #          sidebarPanel(
             #            h4("Hyperparameters for Probe-Level Analysis"),
             #            numericInput(inputId = "proportionSampleProbe",
             #                         "Minimum proportion of sample considered to be a peak",
             #                         value = 0.05),
             #            numericInput(inputId = "peakDistanceProbe",
             #                         "Minimum distance between adjacent peaks",
             #                         value = 0.1),
             #            h4("Smoothing parameters"),
             #            numericInput(inputId = "densityAdjustProbe", "density() 'adjust' parameter",
             #                         value = 1.5),
             #            numericInput(inputId = "pushToZeroProbe",
             #                         "Threshold for numbers small enough to be set to zero",
             #                         value = 1e-6),
             #            helpText("The density() function fits very small floating-point",
             #                     "values to the data, creating tiny perturbations that ",
             #                     "reduce the accuracy of MethylModes. Setting a threshold ",
             #                     "under which small values are considered equivalent ",
             #                     "to zero mitigates this issue."),
             #            h4("Graph Display Options"),
             #            numericInput(inputId = "numHistogramBins", "Number of histogram bins", value = 50),
             #            shinyjs::disabled(textInput("probeId", "Probe Id",
             #                                           value = "cg27399079")),
             #            # add_busy_spinner(spin = "self-building-square",
             #            #                  timeout = 500,
             #            #                  position = "top-right"),
             #            shinyjs::disabled(actionButton("runProbe", "Run on Selected Probe")),
             #            div(style = "margin-bottom: 10px;"),
             #            shinyjs::disabled(actionButton("runProbeRandom", "Run on Randomly Selected Probe"))
             #          ),
             #          mainPanel(
             #            # h4("Summary of Probe-Level MethylModes Results"),
             #            # withSpinner(plotlyOutput(outputId = "probeVisual")),
             #            # plotOutput(outputId = "probeVisualBaseR")
             #          )
             # ), # End of "Individual probe analysis" tab
             ##### Beta Matrix-level analysis #####
             # tabPanel("Whole Genome Analysis",
             #          sidebarPanel(
             #            h4("Hyperparameters for Beta Matrix-Level Analysis"),
             #            numericInput(inputId = "proportionSampleBetaMatrix",
             #                         "Minimum proportion of sample considered to be a peak",
             #                         value = 0.05),
             #            numericInput(inputId = "peakDistanceBetaMatrix",
             #                         "Minimum distance between adjacent peaks",
             #                         value = 0.1),
             #            h4("Smoothing parameters"),
             #            numericInput(inputId = "densityAdjustBetaMatrix", "density() 'adjust' parameter",
             #                         value = 1.5),
             #            numericInput(inputId = "pushToZeroBetaMatrix",
             #                         "Threshold for numbers small enough to be set to zero",
             #                         value = 1e-6),
             #            helpText("The density() function fits very small floating-point",
             #                     "values to the data, creating tiny perturbations that ",
             #                     "reduce the accuracy of MethylModes. Setting a threshold ",
             #                     "under which small values are considered equivalent ",
             #                     "to zero mitigates this issue."),
             #            shinyjs::disabled(numericInput(inputId = "rangeStartBetaMatrix", "Start of range of beta matrix rows",
             #                                           value = 1)),
             #            shinyjs::disabled(numericInput(inputId = "rangeEndBetaMatrix", "End of range of beta matrix rows",
             #                                           value = 1)),
             #            add_busy_spinner(spin = "self-building-square",
             #                             timeout = 500,
             #                             position = "top-right"),
             #            shinyjs::disabled(actionButton("runMultiProbe", "Run MethylModes"))
             #          ),
             #          mainPanel(
             #            h4("Preview of Beta Matrix-Level MethylModes Results"),
             #            shinyjs::disabled(downloadButton("downloadPeakSummary",
             #                                             "Download Results")),
             #            withSpinner(plotlyOutput("peakCountBar")),
             #            plotOutput("peakSummaryPreview"),
             #            h3("Coming soon: view sorted MethylModes results (e.g. among multimodal probes, show probes with highest number of modes to fewest")
             #          )
             # ), # End of "Whole Genome Analysis" tab
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
