./app.R                                                                                             0000644 0001750 0001750 00000065324 13605370703 013124  0                                                                                                    ustar 00abhinandan                      abhinandan                                                                                                                                                                                                             library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(maftools)
library(shinycssloaders)

options(shiny.maxRequestSize=40*1024^2)

ui <- dashboardPagePlus(
  dashboardHeaderPlus(title = "EasyMAF", titleWidth = 315),
  
  dashboardSidebar( width = 315,
                    sidebarMenu(
                      menuItem("Introduction", tabName = "introduction", icon = icon("fas fa-user-o")),
                      menuItem("Input", tabName = "input", icon = icon("fas fa-upload")),
                      menuItem("MAF Files", tabName = "Maf Summary", icon = icon("fas fa-file"),
                               menuItem("Summary", tabName = "summary", icon = icon("fas fa-calendar-check-o"),
                                        menuSubItem("Sample Summary", tabName = "samplesummary", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Gene Summary", tabName = "genesummary", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Clinical Data", tabName = "clinicaldata", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Fields", tabName = "fields", icon = icon("fas fa-caret-right"))),
                               menuItem("Visualization", tabName = "visualization", icon = icon("fas fa-line-chart"),
                                        menuSubItem("Plotting MAF Summary", tabName = "mafsummary", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Oncoplots", tabName = "oncoplots", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Oncostrip", tabName = "oncostrip", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Transition and Transversions", tabName = "transitionandtransversion", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Lollipop Plot", tabName = "lollipopplot", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Rainfall Plot", tabName = "rainfallplot", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Compare mutation load against TCGA cohorts", tabName = "compare", icon = icon("fas fa-caret-right")),
                                        menuSubItem("VAF Plots", tabName = "vafplot", icon = icon("fas fa-caret-right")),
                                        menuSubItem("Genecloud", tabName = "genecloud", icon = icon("fas fa-caret-right")))),
                      menuItem("Copy-Number Data", tabName = "copynumberdata", icon = icon("fas fa-list-ol"),
                               menuItem("GISTIC Input", tabName = "gisticinput", icon = icon("fas fa-caret-right"),
                                        menuItem("Summary", tabName = "gisticsum", icon = icon("fas fa-calendar-check-o"),
                                                 menuSubItem("Sample Summary", tabName = "gissamplesum"),
                                                 menuSubItem("Gene Summary", tabName = "gisgenesum"),
                                                 menuSubItem("Cytoband Summary", tabName = "cytobandsum")),
                                        menuItem("Visualization", tabName = "gisticplot", icon = icon("fas fa-line-chart"),
                                                 menuSubItem("Genome Plot", tabName = "genomeplot", icon = icon("fas fa-caret-right")),
                                                 menuSubItem("Bubble Plot", tabName = "bubbleplot", icon = icon("fas fa-caret-right")),
                                                 menuSubItem(" Oncoplot", tabName = "oncoplot", icon = icon("fas fa-caret-right")),
                                                 menuSubItem("Visualize CBS segments", tabName = "cbssegments", icon = icon("fas fa-caret-right"))))),
                      menuItem("Analysis", tabName = "analysis", icon = icon("fas fa-pencil-square-o"),
                               menuSubItem("Somatic Interactions", tabName = "somaticinteraction", icon = icon("fas fa-caret-right")),
                               menuSubItem("Oncodrive", tabName = "oncodrive", icon = icon("fas fa-caret-right")),
                               menuSubItem("pfam Domains", tabName = "pfamdomains", icon = icon("fas fa-caret-right")),
                               menuSubItem("Pan-Cancer Comparision", tabName = "comparision", icon = icon("fas fa-caret-right")),
                               menuSubItem("Survival Analysis", tabName = "survivalanalysis", icon = icon("fas fa-caret-right")),
                               menuSubItem("Clinical Enrichment Analysis", tabName = "clinicalenrichmentanalysis", icon = icon("fas fa-caret-right")),
                               menuSubItem("Oncogenic Signalling Pathways", tabName = "signallingpathways", icon = icon("fas fa-caret-right")))
                    )
  ),
  dashboardBody(
                 tabItems(
                   tabItem(tabName = "introduction",
                              h2("Introduction"),
                               p("EasyMAF is an interactive web based application developed in R using Shiny, to analyze and visualize Mutation Annotation Format (MAF) files from large scale sequencing studies.
                                 The application aims to provide users with an easy and interactive user interface to perform most commonly used analyses in cancer genomics and to create feature rich customizable visualizations with minimal effort."),
                                
                              h3("1. Data Input", style="padding-left: 1em"),
                               p("The input data for this app requires 2 files in '.maf.gz' and '.tsv' format.", style="padding-left: 2em"),
                               p("Required input files:", style="padding-left: 2em"),
                                tags$ul(style="padding-left: 4em",
                                  tags$li("an MAF file - can be gz compressed."),
                                  tags$li("an optional but recommended clinical data associated with each sample/Tumor_Sample_Barcode in MAF."),
                                  tags$li("an optional copy number data if available. Can be GISTIC output or a custom table containing sample names, gene names and copy-number status."),
                               p("An example of MAF file named as 'sample_maf.tsv' can be accessed", a("here.", href="https://docs.google.com/spreadsheets/d/1m_x-MfieBtw4H0iYqkXX0Z8La_321LDXtF_T07koza8/edit#gid=623536838", style="padding-left: 0em"), 
                                 "Similarly an example of copy-number data named as 'sample_copy-number_data.tsv' can be accessed", a("here.", href="https://docs.google.com/spreadsheets/d/14XSllwr-zcwDDW_q4kILNWHTTg_bTEpvXF6aps-dyXQ/edit#gid=2080309177"))),
                                  
                              h3("2. Summary and Visualization", style="padding-left: 1em"),
                               p(a("maftools", href="https://bioconductor.org/packages/release/bioc/html/maftools.html"), "package provides various functions to perform widely used analyses in cancer genomics and to create feature rich customizable visualizations.", style="padding-left: 2em"),
                               p("Information on sample summary, gene summary, clinical data associated with samples, fields in MAF file can be obtained. It also provides various functions to visualize MAF files. 
                                  These include:", style="padding-left: 2em"),
                                    tags$ul(style="padding-left: 4em",
                                      tags$li("Plotting MAF summary"),
                                      tags$li("Oncoplots"),
                                      tags$li("Oncostrip"),
                                      tags$li("Transitions and Transversions"),
                                      tags$li("Lollipop Plot"),
                                      tags$li("Rainfall Plot"),
                                      tags$li("Compare mutation load against TCGA cohorts"),
                                      tags$li("VAF Plots"),
                                      tags$li("Genecloud")),
                           
                              h3("3. Processing copy-number data", style="padding-left: 1em"),
                               p("Varoius functions are available in the package for processing copy-number data. Summary of output files generated by GISTIC programme can be obtained. As mentioned earlier, we need four files that are generated by GISTIC, i.e, all_lesions.conf_XX.txt, amp_genes.conf_XX.txt, del_genes.conf_XX.txt and scores.gistic, where XX is the confidence level. See ", 
                                  a("GISTIC documentation", href="ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTICDocumentation_standalone.htm"), "for details. 
                                  There are 3 types of plots available to visualize gistic results:", style="padding-left: 2em"),
                                     tags$ul(style="padding-left: 4em",
                                       tags$li("Genome Plot"),
                                       tags$li("Bubble Plot"),
                                       tags$li("Oncoplot")),
                           
                              h3("4. Analysis", style="padding-left: 1em"),
                               p("maftools package also provides numerous functions for analyzing copy-number data. These functions are implemented in the interface for easy understanding of the analyses process.", style="padding-left: 2em"),
                               p("Functions available in the package for analysis include:", style="padding-left: 2em"),
                                     tags$ul(style="padding-left: 4em",
                                       tags$li("Somatic Interactions"),
                                       tags$li("Detecting cancer driver genes based on positional clustering"),
                                       tags$li("Adding and summarizing pfam domains"),
                                       tags$li("Pan-Cancer comparision"),
                                       tags$li("Survival analysis"),
                                       tags$li("Clinical Enrichment Analysis"),
                                       tags$li("Oncogenic signalling pathways")),
                                 ),
                   tabItem(tabName = "input",
                           boxPlus(title = "MAF Files Upload", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   boxPlus(title = "MAF File", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                           fileInput("file1", "Upload MAF File", multiple = FALSE, accept = c(".maf",".maf.gz"))),
                                   boxPlus(title = "Clinical Data", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                           fileInput("file2", "Upload Clinical File (Optional)", accept = c(".csv",".tsv"))))),
                   
                   
                   
                   tabItem(tabName = "samplesummary",
                           withSpinner(dataTableOutput("samplesummary"))),
                   
                   tabItem(tabName = "genesummary",
                           withSpinner(dataTableOutput("genesummary"))),
                   
                   tabItem(tabName = "clinicaldata",
                           withSpinner(dataTableOutput("clinicalsummary"))),
                   
                   tabItem(tabName = "fields",
                           withSpinner(tableOutput("fields")),
                           uiOutput("bugg")),
                   
                   tabItem(tabName = "mafsummary",
                           boxPlus(title = "Plot Maf Summmary", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("mafsummary"))),
                           boxPlus(title = "Add Stat", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto", collapsed = TRUE,
                                   radioButtons(inputId = "addstat", label = "Add Stat:" ,choiceNames = c("Mean","Median"), choiceValues = c("mean","median"))),
                           boxPlus(title = "Remove Outlier", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto", collapsed = TRUE,
                                   radioButtons(inputId = "rmoutlier", label = "Remove Outlier:" ,choiceNames = c("Yes","No"), choiceValues = c("TRUE","FALSE")))),
                   
                   tabItem(tabName = "oncoplots",
                           boxPlus(title = "Plot Oncoplot", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("oncoplots"))),
                           boxPlus(title = "Number of top mutations to plot:", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   sliderInput("oncoslider", "Number of top mutations to display:", 0.0, 100, 10, step = 10))),
                   
                   
                   tabItem(tabName = "oncostrip",
                           boxPlus(title = "Plot Oncostrip", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("oncostrip"))),
                           boxPlus(title = "Enter gene names to plot:", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto", collapsed = TRUE,
                                   textInput("oncotext", "Gene name(s) separated by a comma" , value = "TTN,MUC2,FLT3,NPM1"))),
                   
                   tabItem(tabName = "transitionandtransversion",
                           boxPlus(title = "Plot Transitions and Transvaersions", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("transitionandtransversion")))),
                   
                   tabItem(tabName = "lollipopplot",
                           boxPlus(title = "Lollipop Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("lollipopplot"))),
                           boxPlus(title = "Enter gene name to plot:", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto", collapsed = TRUE,
                                   textInput("lollipopgene" , "Enter gene name" , value = "TTN"))),
                   
                   tabItem(tabName = "rainfallplot",
                           boxPlus(title = "Rainfall Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("rainfallplot"))),
                           boxPlus(title = "Select point size", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, collapsed = TRUE, width = "auto",
                                   sliderInput("rainslider", "", 0.0, 1, 0.5, step = .1)),
                           boxPlus(title = "Detect change points", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, collapsed = TRUE, width = "auto",
                                   radioButtons(inputId = "rainradio", label = "Detect change points:", choiceNames =list("Yes","No"), choiceValues =list("TRUE","FALSE")))),
                   
                   tabItem(tabName = "compare",
                           boxPlus(title = "Compare mutation load against TCGA cohorts", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("compare")))),
                   
                   tabItem(tabName = "vafplot",
                           boxPlus(title = "VAF Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("vafplot")))),
                   
                   tabItem(tabName = "genecloud",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("genecloud")))),
                   tabItem(tabName = "gisticinput",
                           fluidRow(
                             column(12,
                                    boxPlus(title = "GISTIC Files Upload", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, collapsed = TRUE, width = "auto",
                                            boxPlus(title = "", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                                    fileInput("file3", "Upload 'all lesions' file", multiple = FALSE, accept = c(".txt"))),
                                            boxPlus(title = "", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                                    fileInput("file4", "Upload amplified genes file", accept = c(".txt"))),
                                            boxPlus(title = "", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                                    fileInput("file5", "Upload deleted genes file", multiple = FALSE, accept = c(".txt"))),
                                            boxPlus(title = "", status = "warning", solidHeader = TRUE, collapsible = TRUE, closable = FALSE,
                                                    fileInput("file6", "Upload gistic scores file", accept = c(".gis",".gistic"))))))),
                   tabItem(tabName = "gissamplesum",
                           withSpinner(dataTableOutput("gissamplesum"))),
                   
                   tabItem(tabName = "gisgenesum",
                           withSpinner(dataTableOutput("gisgenesum"))),
                   
                   tabItem(tabName = "cytobandsum",
                           withSpinner(dataTableOutput("cytobandsum"))),
                   
                   tabItem(tabName = "genomeplot",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("genomeplot")))),
                   
                   tabItem(tabName = "bubbleplot",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("bubbleplot")))),
                   
                   tabItem(tabName = "gisticoncoplot",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("gisticoncoplot")))),
                   
                   tabItem(tabName = "cbssegments",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("cbssegments")))),
                   
                   tabItem(tabName = "somaticinteraction",
                           boxPlus(title = "Somatic Interaction", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("somaticinteractions")))),
                   
                   tabItem(tabName = "oncodrive",
                           boxPlus(title = "Oncodrive", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("oncodrive")))),
                   
                   tabItem(tabName = "pfamdomains",
                           boxPlus(title = "pfam Domains", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("pfamdomains")))),
                   
                   tabItem(tabName = "comparision",
                           boxPlus(title = "Pan-Cancer Comparision", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("pancancer")))),
                   
                   tabItem(tabName = "survivalanalysis",
                           boxPlus(title = "Survival Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("survivalanalysis")))),
                   
                   tabItem(tabName = "clinicalenrichmentanalysis",
                           boxPlus(title = "Gene Cloud", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("enrichment")))),
                   
                   tabItem(tabName = "signallingpathways",
                           boxPlus(title = "Signalling Pathways", status = "primary", solidHeader = TRUE, collapsible = TRUE, closable = FALSE, width = "auto",
                                   withSpinner(plotOutput("oncopathway"))))
                 )))

server <- function(input, output) {
  ## MAF File summary
  # Read user file and clinical data
  usrvar <- reactive({
    validate(need(input$file1$datapath != "", "Please upload an MAF file."))
    
    if (is.null(input$file2$datapath)){
      read.maf(input$file1$datapath)}
    else {
      read.maf(maf = input$file1$datapath , clinicalData = input$file2$datapath)
    }})
  
  
  aff <- reactive({
    laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
    laml = read.maf(maf = laml.maf)
    faf <- list(getFields(laml))
    return(faf)
  })
  
  output$bugg <- renderUI({
    selectInput("bugg", label = "Choose:", choices = aff())
  })
  
  
  # Render summary tables
  output$samplesummary <- renderDataTable({
    getSampleSummary(usrvar())}, options = list(scrollX = TRUE))
  
  output$genesummary <- renderDataTable({
    getGeneSummary(usrvar())}, options = list(scrollX = TRUE))
  
  output$fields<- renderTable({
    getFields(usrvar())})
  
  output$clinicalsummary <- renderDataTable({
    getClinicalData(usrvar())}, options = list(scrollX = TRUE))
  
  # Plot visualizations
  output$mafsummary <- renderPlot({
    plotmafSummary(maf = usrvar(), rmOutlier = input$rmoutlier , addStat = input$addstat, dashboard = TRUE, titvRaw = TRUE)
  })
  
  output$oncoplots <- renderPlot({
    oncoplot(maf = usrvar() , top = input$oncoslider)
  })
  
  output$oncostrip <- renderPlot({
    oncostrip(maf = usrvar() , genes = unlist(strsplit(input$oncotext,",")))
  })
  
  output$transitionandtransversion <- renderPlot({
    usrvar.titv = titv(maf = usrvar(), plot = FALSE, useSyn = TRUE)
    plotTiTv(res = usrvar.titv)
  })
  
  output$rainfallplot <- renderPlot({
    rainfallPlot(maf = usrvar(), detectChangePoints = input$rainradio , pointSize = input$rainslider)
  })
  
  output$compare <- renderPlot({
    usrvar.mutload = tcgaCompare(maf = usrvar())
  })
  
  output$vafplot <- renderPlot({
    plotVaf(maf = usrvar(), vafCol = Variant_Type)
  })
  
  output$genecloud <- renderPlot({
    geneCloud(input = usrvar(), minMut = 3)
  })
  
  output$lollipopplot <- renderPlot({
    lollipopPlot(maf = usrvar(), gene = input$lollipopgene , showMutationRate = TRUE)
  })
  
  ## Copy number data
  
  # Read gistic files
  
  usrgis <- reactive({
    validate(need(input$file3$datapath != "", "Please upload all four GISTIC files."))
    validate(need(input$file4$datapath != "", "Please upload all four GISTIC files."))
    validate(need(input$file5$datapath != "", "Please upload all four GISTIC files."))
    validate(need(input$file6$datapath != "", "Please upload all four GISTIC files."))
    
    readGistic(gisticAllLesionsFile = input$file3$datapath, gisticAmpGenesFile = input$file4$datapath, gisticDelGenesFile = input$file5$datapath, gisticScoresFile = input$file6$datapath, isTCGA = TRUE)
  })
  
  # Summarize copy number data
  
  output$gissamplesum <- renderDataTable({
    getSampleSummary(usrgis())}, options = list(scrollX = TRUE))
  
  output$gisgenesum <- renderDataTable({
    getGeneSummary(usrgis())}, options = list(scrollX = TRUE))
  
  output$cytobandsum<- renderDataTable({
    getCytobandSummary(usrgis())})
  
  # Plot gistic plots
  
  output$genomeplot <- renderPlot({
    gisticChromPlot(gistic = usrgis(), markBands = "all")
    
  })
  
  output$bubbleplot <- renderPlot({
    gisticBubblePlot(gistic = usrgis())
    
  })
  
  output$gisticoncoplot <- renderPlot({
    gisticOncoPlot(gistic = usrgis(), clinicalData = getClinicalData(x = usrvar()), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)
    
  })
  
  output$cbssegments <- renderPlot({
    tcga.ab.009.seg <- system.file("extdata", "TCGA.AB.3009.hg19.seg.txt", package = "maftools")
    plotCBSsegments(cbsFile = tcga.ab.009.seg)
    
  })
  
  # Plot analysis
  
  output$somaticinteractions <- renderPlot({
    somaticInteractions(maf = usrvar(), top = 25, pvalue = c(0.05, 0.1))
    
  })
  
  output$oncodrive <- renderPlot({
    usr.sig = oncodrive(maf = usrvar(), AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
    head(usr.sig)
    
    plotOncodrive(res = usr.sig, fdrCutOff = 0.1, useFraction = TRUE)
    
  })
  
  output$pfamdomains <- renderPlot({
    laml.pfam = pfamDomains(maf = usrvar(), AACol = 'Protein_Change', top = 10)
    
  })
  
  output$pancancer <- renderPlot({
    usr.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
    pancanComparison(mutsigResults = usr.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1)
    
  })
  
  output$survivalanalysis <- renderPlot({
    mafSurvival(maf = usrvar(), genes = 'DNMT3A', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = TRUE)
    
  })
  
  output$enrichment <- renderPlot({
    fab.ce = clinicalEnrichment(maf = usrvar(), clinicalFeature = 'FAB_classification')
    plotEnrichmentResults(enrich_res = fab.ce, pVal = 0.05)
    
  })
  
  output$drug <- renderPlot({
    dgi = drugInteractions(maf = usrvar(), fontSize = 0.75)
    
  })
  
  output$oncopathway <- renderPlot({
    PlotOncogenicPathways(maf = usrvar(), pathways = "RTK-RAS")
    
  })
}

shinyApp(ui = ui, server = server) 
