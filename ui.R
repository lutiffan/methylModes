# TODO: let user set max file upload size via slider??
# create a display that automatically adjusts to the dimensions of your user's browser window
fluidPage(theme = shinythemes::shinytheme("cerulean"), 
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
                                  accept = c(".txt", ".csv", ".idat", ".RDS", ".RDA")), # Most people are not gonna have an RDS format
                        helpText("Accepted file formats: .csv, .txt (tab-delimited), .idat, RGSET, .RDS, .RDA", br(),
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
             ##### Probe-level analysis #####
             tabPanel("Probe-Level Analysis",
                      sidebarPanel(
                        h4("Hyperparameters for Probe-Level Analysis"),
                        numericInput(inputId = "proportionSampleProbe",
                                     "Minimum proportion of sample considered to be a peak",
                                     value = 0.05),
                        numericInput(inputId = "peakDistanceProbe",
                                     "Minimum distance between adjacent peaks",
                                     value = 0.1),
                        h4("Smoothing parameters"),
                        numericInput(inputId = "densityAdjustProbe", "density() 'adjust' parameter",
                                     value = 1.5),
                        numericInput(inputId = "pushToZeroProbe",
                                     "Threshold for numbers small enough to be set to zero",
                                     value = 1e-6),
                        helpText("The density() function fits very small floating-point",
                                 "values to the data, creating tiny perturbations that ",
                                 "reduce the accuracy of MethylModes. Setting a threshold ",
                                 "under which small values are considered equivalent ",
                                 "to zero mitigates this issue."),
                        h4("Graph Display Options"),
                        numericInput(inputId = "numHistogramBins", "Number of histogram bins", value = 50),
                        shinyjs::disabled(actionButton("runProbeRandom", "Run on Random Probe")),
                        shinyjs::disabled(textInput("probeId", "Probe Id",
                                                       value = "cg27399079")),
                        # add_busy_spinner(spin = "self-building-square",
                        #                  timeout = 500,
                        #                  position = "top-right"),
                        shinyjs::disabled(actionButton("runProbe", "Run on Selected Probe"))
                      ),
                      mainPanel(
                        h4("Summary of Probe-Level MethylModes Results"),
                        withSpinner(plotlyOutput(outputId = "probeVisual")),
                        plotOutput(outputId = "probeVisualBaseR")
                      )
             ),
             ##### Beta Matrix-level analysis #####
             tabPanel("Beta Matrix-Level Analysis",
                      sidebarPanel(
                        h4("Hyperparameters for Beta Matrix-Level Analysis"),
                        numericInput(inputId = "proportionSampleBetaMatrix",
                                     "Minimum proportion of sample considered to be a peak",
                                     value = 0.05),
                        numericInput(inputId = "peakDistanceBetaMatrix",
                                     "Minimum distance between adjacent peaks",
                                     value = 0.1),
                        h4("Smoothing parameters"),
                        numericInput(inputId = "densityAdjustBetaMatrix", "density() 'adjust' parameter",
                                     value = 1.5),
                        numericInput(inputId = "pushToZeroBetaMatrix",
                                     "Threshold for numbers small enough to be set to zero",
                                     value = 1e-6),
                        helpText("The density() function fits very small floating-point",
                                 "values to the data, creating tiny perturbations that ",
                                 "reduce the accuracy of MethylModes. Setting a threshold ",
                                 "under which small values are considered equivalent ",
                                 "to zero mitigates this issue."),
                        shinyjs::disabled(numericInput(inputId = "rangeStartBetaMatrix", "Start of range of beta matrix rows",
                                                       value = 1)),
                        shinyjs::disabled(numericInput(inputId = "rangeEndBetaMatrix", "End of range of beta matrix rows",
                                                       value = 1)),
                        add_busy_spinner(spin = "self-building-square",
                                         timeout = 500,
                                         position = "top-right"),
                        shinyjs::disabled(actionButton("runBetaMatrix", "Run MethylModes"))
                      ),
                      mainPanel(
                        h4("Preview of Beta Matrix-Level MethylModes Results"),
                        shinyjs::disabled(downloadButton("downloadPeakSummary",
                                                         "Download Results")),
                        withSpinner(plotlyOutput("peakCountBar")),
                        plotOutput("peakSummaryPreview"),
                        h3("Coming soon: view sorted MethylModes results (e.g. among multimodal probes, show probes with highest number of modes to fewest")
                      )
             ),
             tabPanel("Upload Previously Generated Results",
                      h3("Upload Previously Calculated Results"),
                      sidebarPanel(
                        fileInput("peakSummaryFile", "Choose MethylModes Result File",
                                  multiple = FALSE,
                                  accept = c(".RDS",
                                             ".csv"))
                      )
              ),
             tabPanel("Explore Results",
                      sidebarPanel(),
                      mainPanel())
   )
  
)
