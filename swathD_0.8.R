library(shinydashboard)
library(shiny)
library(DT)
library(readxl) # read_excel(),
library(readr) # read_tsv
library(dplyr)  # %>%, filter(), arrange(), select(), group_by(), left_join()
library(ggpubr)
library(biomaRt) # useMart()
library(impute) #impute.knn
library(pheatmap) # pheatmap
library(circlize)
library(grDevices) # pdf(), dev.off()
library(ggplot2)
library(reshape2)
library(ggrepel)
library(FactoMineR)
library(corrplot)
library(tidyr)
library(enrichR)
library(stringr)
library(purrr)
library(visNetwork)
#library(plotly)


options(shiny.maxRequestSize=100*1024^2)

sidebar <- dashboardSidebar(
  sidebarUserPanel("MSBoss",
                   subtitle = a(href = "#", icon("circle", class = "text-success"), "Online") #,
                   # Image file should be in www/ subdir
                   #image = "myimage.png"
  ),
  
  #sidebarSearchForm(label = "Enter a number", "searchText", "searchButton"),
  sidebarMenu(
    # Setting id makes input$tabs give the tabName of currently-selected tab
    id = "tabs",
    menuItem("1-Click SVVATH-D", tabName = "upload", icon = icon("import", lib = "glyphicon")), 
   # menuItem("Normalizations", tabName = "norm_compare", icon = icon("balance-scale", lib = "font-awesome")),
    menuItem("1-Click Results", tabName = "first_check", icon = icon("bar-chart-o")), 
    menuItem("Volcano2Enrichr", tabName = "re_analysis", icon = icon("list-alt", lib = "glyphicon")),
    menuItem("Heatmap2Enrichr", tabName = "re_analysis_2", icon = icon("list-alt", lib = "glyphicon")),
    menuItem("ML2Enrichr", tabName = "re_analysis_3", icon = icon("list-alt", lib = "glyphicon")),
    menuItem("Advance settings", tabName = "settings", icon = icon("cog", lib = "glyphicon")) #, badgeLabel = "new", badgeColor = "green")
    
    
    
    #menuItem("Widgets", icon = icon("th"), tabName = "widgets", badgeLabel = "new",
    #        badgeColor = "green"),
    #menuItem("Charts", icon = icon("bar-chart-o"),
    #        menuSubItem("Sub-item 1", tabName = "subitem1"),
    #       menuSubItem("Sub-item 2", tabName = "subitem2")
  )
)
#)

#pair_all <- NULL

#load("./lib/ensembl_dbs")
load("lib/mmusculus_ensembl_012219")
load("lib/rnorvegicus_ensembl_030319")
load("lib/hsapiens_ensembl_012219")
source("lib/swathD_idConv.R", local = TRUE)
source("lib/swathD_imputeFirst.R", local = TRUE)
source("lib/compare_norm.R", local = TRUE)
source("lib/qc_box.R", local = TRUE)
source("lib/qc_density.R", local = TRUE)
source("lib/swathD_naImpute.R", local = TRUE)
source("lib/swathD_adjBatch.R", local = TRUE) 
source("lib/swathD_volcano.R", local = TRUE)
source("lib/swathD_hm_var.R", local = TRUE)
source("lib/swathD_ml2enrichr.R", local = TRUE)
source("lib/swathD_ml2net.R", local = TRUE)
wd <- getwd()

body <- dashboardBody(
  tags$style(HTML("
                  .content-wrapper {
                  background-color: white !important;
                  }
                  .main-sidebar {
                  background-color: grey !important;
                  }
                  ")),
  tabItems(
    tabItem("upload", 
        fluidPage(
          column(3, 
            #sidebarLayout(
            wellPanel(
              fileInput("file1", "Choose metadata file (.csv)",
                        multiple = FALSE,
                        accept = c(".csv")
              ),
              fileInput("file2", "Choose Peakview SWATH full report file (.xlsx)",
                        multiple = FALSE,
                        accept = c(".xlsx")
              ),
              sliderInput("cutFC", label = p("Absolute fold change threshold"), min = 0, max = 4, value = 2, step = 0.25),
              #numericInput("cutFC", label = h4("Absolute fold change threshold"), min = 1, max = 10, value = 2, step = 0.25, width = "50%"),
              #sliderInput("cutP", label = h4("P-value threshold"), min = 0, max = 0.1, value = 0.05),
              #numericInput("cutP", label = h4("P-value threshold"), max = 0.99, value = 0.05),
              br(),
              textInput("cutP", label = p("P-value threshold"), value = 0.05),
              br(),
              selectInput("outplot", label = "Select 1-Click output", c("QC-sig", "QC-sig-enrichPlot"), selected = "QC-sig"),
              #checkboxInput("onlyQC", label = h4("Only QC plots"), value = FALSE, width = NULL),
              br(),
              actionButton("go", "RUN", icon = icon("play-circle"), class = "butt"),
                   tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                         .butt{box-shadow:
                                              0 1px 0 rgba(255, 255, 255, 0.25),
                                              0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                              0 0 0 rgba(0, 0, 0, 0.5) inset,
                                              0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                              0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                              0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
            )),
            mainPanel(
              tabsetPanel(type = "tabs",
                          tabPanel("Input", 
                                   column(width = 12, 
                                                   box(
                                                     title = "1. Upload metadata", width = 12, status = "primary",
                                                     div(style = 'overflow-x: scroll', DT::dataTableOutput("contents1"))
                                                     #div(style = 'overflow-x: scroll', tableOutput("contents1"))
                                                   )
                                    ),
                                    column(width = 12, 
                                               box(
                                                 title = "2. Upload SWATH full report", width = 12, status = "primary",
                                                 div(style = 'overflow-x: scroll', DT::dataTableOutput("contents2"))
                                               )
                                    ), 
                                    column(width = 12, 
                                           box(
                                             title = "3. Indicate thresholds and click RUN", width = 12, status = "primary")
                                            )
                                    ),
                          tabPanel("Summary", verbatimTextOutput("summary")),
                          #tabPanel("Table", tableOutput("table")),
                          tabPanel("Distribution", fluidPage(
                            fluidRow(column(12, h2("Distribution of normalized data"), 
                                            actionButton("normPlot", label="View Box and Density Plots"))),
                            br(), br(),
                            fluidRow(
                              column(6, 
                                     plotOutput("plot16", height = 350),
                                     plotOutput("plot17", height = 350),
                                     plotOutput("plot18", height = 350),
                                     plotOutput("plot19", height = 350),
                                     plotOutput("plot20", height = 350)),
                              column(6,
                                     plotOutput("plot21", height = 350),
                                     plotOutput("plot22", height = 350),
                                     plotOutput("plot23", height = 350),
                                     plotOutput("plot24", height = 350),
                                     plotOutput("plot25", height = 350))
                            )
                            
                          )),
                          tabPanel("Variance", fluidPage(
                            fluidRow(column(12, h2("Variance of normalization methods"), 
                                            actionButton("variancePlot", label="View Variance Plots"))),
                            br(), br(),
                            fluidRow(
                              column(8, plotOutput("plot26", height = 350)),
                              column(8, plotOutput("plotPCV", height = 350)),
                              column(8, plotOutput("plotPEV", height = 350))
                            )
                          ))
            ))
          )
    ),
#    tabItem("norm_compare", 
 #           tabsetPanel(type = "tabs" 
  #                      
   #         ) 
    #),
    tabItem("first_check", 
            tabsetPanel(type = "tabs", 
                        tabPanel("Quality control", fluidPage(
                          fluidRow(column(12, h2("Evaluation of Proteomic Data Quality"), 
                             actionButton("qcPlot", label="View Quality Plots"))),
                          br(), br(),
                            #actionButton("qcPlot", "", icon = icon("play-circle")), br(),
                            column(width = 8, plotOutput("plot3", height = 350)),
                            column(width = 4, plotOutput("plot1" ,height = 350)),
                            column(width = 4, plotOutput("plot2" ,height = 350)),
                            column(width = 4, plotOutput("plot5" ,height = 350)),
                            column(width = 4, plotOutput("plot4" ,height = 350))
                            )), # 
                       tabPanel("Significant proteins", fluidPage(
                         fluidRow(column(12, h2("Visualization of Significant Results: Volcano Plot and Self-Clusterred Heatmap"),
                                  actionButton("volcanoPlot", label="View Significant Plots"))),
                            #actionButton("volcanoPlot", "", icon = icon("play-circle")), 
                            br(), br(),
                            column(width = 8, plotOutput("plot6")), 
                           # br(), br(),
                            column(width = 4, plotOutput("plot7")) 
                            )),
                       tabPanel("Functional enrichment", fluidPage(
                         fluidRow(column(12, h2("Overview of Biological Interpretation"),
                                  actionButton("goPlot", label="View Enrich Plots")), #),
                                  #useShinyjs(),
                           #actionButton("goPlot", "", icon = icon("play-circle")), br(),
                           column(width = 6, plotOutput("plot8", height = 600)),
                           column(width = 6, plotOutput("plot9", height = 600)),
                           column(width = 6, plotOutput("plot10", height = 600)),
                #            )),
                 #      tabPanel("Pathways-enrichment", fluidPage(
                  #       fluidRow(column(12, h2("Overview of Biological Interpretation: Top10-matched Pathway"),
                   #               actionButton("pathwayPlot", label="View Pathway-Enrich Plots"))), 
                           #actionButton("pathwayPlot", "", icon = icon("play-circle")), br(),
                           column(width = 6, plotOutput("plot11", height = 600)),
                           column(width = 6, plotOutput("plot12", height = 600)),
                           column(width = 6, plotOutput("plot13", height = 600))
                           ))
                        ),
                       tabPanel("Enrich Table", fluidPage(
                           h2("Gene enrichment analysis: All matched results"),
                           actionButton("masterTable", "", icon = icon("play-circle")), br(),
                           column(width = 12, 
                                  box(width = 12, status = "primary", 
                                  div(style = 'overflow-x: scroll', DT::dataTableOutput("contents3"))
                                  ))
                       ))
            )
    ),
    
    tabItem("re_analysis", 
            tabsetPanel(type = "tabs",
                    tabPanel("Volcano Plot", fluidPage(
                        column(3, 
                            wellPanel(
                             #sidebarPanel(
                              selectInput("pair_select", label = "Pairwise subset", choices = NULL),
                              sliderInput("re_cutFC", label = p("Absolute fold change threshold"), min = 0, max = 4, value = 2, step = 0.05),
                              textInput("re_cutP", label = p("P-value threshold"), value = 0.05),
                              sliderInput("cutFC_label", label = p("Label: X-cutpoint (abs. log2)"), min = 0, max = 10, value = 1, step = 0.1), 
                              sliderInput("cutP_label", label = p("Label: Y-cutpoint (minus log10)"), min = 0, max = 10, value = 1.3, step = 0.1),
                              hr(),
                              sliderInput("xAxis", label = p("View: X-range"), min = -20, max = 20, value = c(-5, 5), step = 0.5),
                              sliderInput("yAxis", label = p("View: Y-range"), min = 0, max = 20, value = c(0, 5), step = 0.5),
                              actionButton("re_volcanoPlot", "Run", icon = icon("play-circle"), class = "butt"),  
                                  tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                        .butt{box-shadow:
                                                         0 1px 0 rgba(255, 255, 255, 0.25),
                                                         0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                         0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                         0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                         0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                         0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                            )),
                            mainPanel(
                              column(width = 12, plotOutput("plot14"),
                                     actionButton("select_protVolcano", "1. Generate and check new inputs", class = "butt2"), 
                                        tags$head(tags$style(".butt2{display: inline-block;
                                                             padding: .3em .75em;
                                                             background-color: #f6f6f6;
                                                             background-image: linear-gradient(rgba(0,0,0,0), rgba(0,0,0,.15));
                                                             text-decoration: none;
                                                             color: rgba(0,0,0,.8);
                                                             text-shadow: 0 1px rgba(255,255,255,.3);
                                                             font-weight: 700;
                                                             border-radius: 4px;
                                                             border: 1px solid #aaa;
                                                             box-shadow: inset 0 0 0 1px rgba(255,255,255,.1), inset 0 1px rgba(255,255,255,.3);}")),
                                     actionButton("submit_ReEnrich", "2. Submit to re_enrichment analysis", class = "butt2"),
                                     actionButton("export_reVolcano", "Optional: Download figure (pdf)", class = "butt2"),
                                     box(width = 12, status = "primary", 
                                         div(style = 'overflow-y: scroll', htmlOutput("contents4")))
                                     
                              ))
                    )),
                    tabPanel("Enrich Plot", fluidPage(
                      column(3, 
                             wellPanel(
                            #sidebarPanel(
                            #fluidRow(box(width = 12, wellPanel(
                              #column(2, 
                                     selectInput("pair_select2", label = "Pairwise", choices = NULL),
                              #),
                              #column(3, 
                                     selectInput("reg_select", label = "Input for re-enrichment", choices = NULL),
                              #),
                              #column(3, 
                                     selectInput("dbs_select", label = "Databases", choices = c("GO_Biological_Process_2015", "GO_Cellular_Component_2015", "GO_Molecular_Function_2015", 
                                                                                                "KEGG_2016", "Reactome_2016", "WikiPathways_2016", 
                                                                                                "ENCODE_TF_ChIP-seq_2015", "ChEA_2016", "KEA_2015",
                                                                                                "Transcription_Factor_PPIs", "CORUM", "PPI_Hub_Proteins",
                                                                                                "Jensen_DISEASES", "Human_Phenotype_Ontology", "dbGaP")),
                              #),
                              #column(2, 
                                     selectInput("top_n", label = "Top", choices = c(5, 10, 15, 20, 25), selected = 10),
                              #),
                              br(),
                              actionButton("re_enrichAnal", "Run", icon = icon("play-circle"), class = "butt"),  
                                    tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                         .butt{box-shadow:
                                                         0 1px 0 rgba(255, 255, 255, 0.25),
                                                         0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                         0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                         0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                         0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                         0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                            )),
                     # ), 
                            #   mainPanel(
                            #    column(width = 12, 
                            #            box(title = "ReEnrich", width = 12, status = "primary",
                            #             div(style = 'overflow-x: scroll', DT::dataTableOutput("contents5"))
                            #          )
                            #  )),
                            #br(),
                            br(),
                            br(),
                     mainPanel( 
                       fluidRow(
                         column(width = 2),
                         column(width = 10, plotOutput("plot15", width = "500px", height = "600px")) #, height = 500)
                            )) 
            ))
            )
  ),
  tabItem("re_analysis_2", 
          tabsetPanel(type = "tabs",
                      tabPanel("Heatmap", fluidPage(
                               column(3, 
                                      wellPanel(
                                        selectInput("hm_set", "Choose data type", 
                                                    choices = c("Significant features", "All features"), selected = "Significant features"),
                                        sliderInput("varTopN", label = p("Top variance feature (%)"), min = 1, max = 100, value = 10, step = 1),
                                        sliderInput("n_cluster", label = p("Numbers of cluster"), min = 1, max = 10, value = 2, step = 1),
                                        actionButton("re_hm_var", "Run", icon = icon("play-circle"), class = "butt"),  
                                                tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                                     .butt{box-shadow:
                                                                     0 1px 0 rgba(255, 255, 255, 0.25),
                                                                     0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                                     0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                                     0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                                     0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                                     0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                                        )
                                      ),
                              # ),
                               mainPanel( 
                                   fluidRow(
                                  column(width = 12, 
                                                plotOutput("plot30", height = 900))
                                  ))
                      )),
                      tabPanel("Enrich Plot", fluidPage(
                          #  fluidRow(box(width = 12, 
                            column(3,  
                                wellPanel(
                                    selectInput("cluster_select", label = "Selected cluster", choices = NULL),
                                #),
                                #column(3,
                                    selectInput("dbs_select2", label = "Databases", choices = c("GO_Biological_Process_2015", "GO_Cellular_Component_2015", "GO_Molecular_Function_2015", 
                                                                                                  "KEGG_2016", "Reactome_2016", "WikiPathways_2016", 
                                                                                                  "ENCODE_TF_ChIP-seq_2015", "ChEA_2016", "KEA_2015",
                                                                                                  "Transcription_Factor_PPIs", "CORUM", "PPI_Hub_Proteins",
                                                                                                  "Jensen_DISEASES", "Human_Phenotype_Ontology", "dbGaP")),
                              #  ),
                               # column(2, 
                                    selectInput("top_n2", label = "Top", choices = c(5, 10, 15, 20, 25), selected = 10),
                              #  ),
                                    br(),
                                    actionButton("cluster_enrichAnal", "Run", icon = icon("play-circle"), class = "butt"),  
                                    tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                           .butt{box-shadow:
                                                           0 1px 0 rgba(255, 255, 255, 0.25),
                                                           0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                           0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                           0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                           0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                           0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                            )),
                          #  br(),
                            br(),
                            br(),
                            mainPanel( 
                              fluidRow(
                                column(width = 2),
                                column(width = 10, plotOutput("plot31", width = "500px", height = "600px"))
                      ))
                      ))
  )),
  tabItem("re_analysis_3", 
          tabsetPanel(type = "tabs",
                      tabPanel("Machine Learning", fluidPage(
                        column(3, 
                               wellPanel(
                                 selectInput("ml", "Choose method", 
                                             choices = c("pca", "sipca", "nmf"), selected = "pca"),
                                 actionButton("ml_var", "Start machine learning", icon = icon("play-circle"), class = "butt"),  
                                       tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                            .butt{box-shadow:
                                                            0 1px 0 rgba(255, 255, 255, 0.25),
                                                            0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                            0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                            0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                            0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                            0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                                 ),
                               wellPanel(
                                 selectInput("chose_com_type", "Choose data type",
                                             choices = c("Significant features", "All features"), selected = "All features"),
                                 selectInput("choose_component_ml2hm", "Choose component",
                                             choices = c(1, 2, 3), selected = 1),
                                 actionButton("ml2hm", "Chosen dataset", icon = icon("play-circle"), class = "butt"),  
                                       tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                            .butt{box-shadow:
                                                            0 1px 0 rgba(255, 255, 255, 0.25),
                                                            0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                            0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                            0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                            0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                            0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                               )
                        ),
                        mainPanel( 
                          fluidRow(
                            column(width = 12, 
                                   textOutput("warning_ml2hm"),
                                   plotOutput("plot33", height = 900))
                          ))
                      )),
                      tabPanel("Enrich Plots", 
                            fluidPage(  
                                column(3, 
                                  wellPanel(
                                       selectInput("ml2", "Chosen method", choices = NULL),
                                       selectInput("choose_component", "Choose component", 
                                                     choices = c(1, 2, 3), selected = 1),
                                       selectInput("com_set", "Choose data type", 
                                                   choices = c("Significant features", "All features"), selected = "Significant features"),
                                       selectInput("dbs_select3", label = "Databases", choices = c("GO_Biological_Process_2015", "GO_Cellular_Component_2015", "GO_Molecular_Function_2015", 
                                                                                                     "KEGG_2016", "Reactome_2016", "WikiPathways_2016", 
                                                                                                     "ENCODE_TF_ChIP-seq_2015", "ChEA_2016", "KEA_2015",
                                                                                                     "Transcription_Factor_PPIs", "CORUM", "PPI_Hub_Proteins",
                                                                                                     "Jensen_DISEASES", "Human_Phenotype_Ontology", "dbGaP")),
                                        selectInput("top_n3", label = "Top", choices = c(5, 10, 15, 20, 25), selected = 10),
                                        br(),
                                        actionButton("ml_enrichAnal", "Run", icon = icon("play-circle"), class = "butt"),  
                                        tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                           .butt{box-shadow:
                                                           0 1px 0 rgba(255, 255, 255, 0.25),
                                                           0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                           0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                           0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                           0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                           0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                                  )),
                                  br(),
                                  br(),
                                 mainPanel( 
                                   fluidRow(
                                     column(width = 2),
                                     column(width = 10, plotOutput("plot32",  width = "500px", height = "600px"))
                                  ))
                    )),
                    tabPanel("EnrichNet Plots", fluidPage( 
                      column(3, 
                           wellPanel(
                          # selectInput("ml3", "Chosen method", choices = NULL),
                            selectInput("chose_ml2net_com", "Choose component", 
                                        choices = c(1, 2, 3), selected = 1),
                            tags$div(title = "Adjusted P-value of gene enrichment analysis",
                            textInput("apval", label = p("Apply adjusted P-value cutoff"), value = 0.05)
                            ),
                            tags$div(title = "%geneRatio * -log10(Adj.Pval))",
                            sliderInput("rscore", label = p("RankScore threshold (0 = not applied)"), min = 0, max = 10, value = 2, step = 0.5)
                            ),
                            tags$div(title = "Numbers of the input gene linked to particular functions/processes",
                            sliderInput("nlink", label = p("Connection degree threshold"), min = 1, max = 10, value = 4, step = 1)
                            ),
                            checkboxGroupInput("ml2net_dbs", label = "Choose databases", choices = list("GO_Biological_Process_2015" = "GO_Biological_Process_2015", 
                                                                                                "GO_Cellular_Component_2015" = "GO_Cellular_Component_2015", 
                                                                                                "GO_Molecular_Function_2015" = "GO_Molecular_Function_2015", 
                                                                                                "KEGG_2016" = "KEGG_2016", 
                                                                                                "Reactome_2016" = "Reactome_2016", 
                                                                                                "WikiPathways_2016" = "WikiPathways_2016", 
                                                                                                "ENCODE_TF_ChIP-seq_2015" = "ENCODE_TF_ChIP-seq_2015", 
                                                                                                "ChEA_2016" = "ChEA_2016", 
                                                                                                "KEA_2015" = "KEA_2015",
                                                                                                "Transcription_Factor_PPIs" = "Transcription_Factor_PPIs", 
                                                                                                "CORUM" = "CORUM", 
                                                                                                "PPI_Hub_Proteins" = "PPI_Hub_Proteins",
                                                                                                "Jensen_DISEASES" = "Jensen_DISEASES", 
                                                                                                "Human_Phenotype_Ontology" = "Human_Phenotype_Ontology", 
                                                                                                "dbGaP" = "dbGaP"),
                                              selected = "GO_Biological_Process_2015"),
                            actionButton("ml_enrichNet", "Run", icon = icon("play-circle"), class = "butt"),  
                                  tags$head(tags$style(".butt{background-color: ;} .butt{color: ;} .butt{border-color: black;} 
                                                       .butt{box-shadow:
                                                       0 1px 0 rgba(255, 255, 255, 0.25),
                                                       0 1px 0 rgba(255, 255, 255, 0.25) inset,
                                                       0 0 0 rgba(0, 0, 0, 0.5) inset,
                                                       0 1.25rem 0 rgba(255, 255, 255, 0.08) inset,
                                                       0 -1.25rem 1.25rem rgba(0, 0, 0, 0.3) inset,
                                                       0 1.25rem 1.25rem rgba(255, 255, 255, 0.1) inset;}"))
                            ),
                          wellPanel(
                            selectInput("layout", "Interactive layout", 
                                        choices = c("layout_nicely", "layout_in_circle", "layout_as_star", "layout_with_mds", "layout_randomly"), selected = "layout_nicely"),
                            checkboxInput("hierarchical", "Hierarchical layout", value = FALSE, width = '400px')
                            #actionButton("ml_reNet", "", icon = icon("play-circle"))
                          )
                   ),
                   mainPanel( 
                     fluidRow(
                       column(width = 12, 
                            textOutput("warning_ml2net"),  
                            visNetworkOutput("ml2net", width = "1000px", height = "900px")
                         )
                   ))
          ))
        )
  ),
  tabItem("settings",
    fluidPage( 
      fluidRow(
               column(3, 
                      wellPanel( h2("Data loading"),
                      h3("Input format"),              
                      selectInput("dataFormat", "Choose data format", 
                           choices = c("Peakview SWATH full report", "Processed dataset")), #"up_proteins (csv)", "up_peptides (csv)", "gn_proteins (csv)", "gn_peptides (csv)"), selected = "Peakview SWATH full report"), # add later:  "MRM (Sciex TripleTOF)", "iTRAQ-4plex", "iTRAQ-8plex", "SILAC",
                      checkboxInput("logdat", "Please check this box if the processed dataset is already log-transformed", value = FALSE, width = '400px')
                      #hr(),
                      #column(4, 
                      ),
                      wellPanel( h2("Data Preprocessing"),
                      h3("Normalization"),
                      selectInput("normMet", "Choose method", 
                                  choices = c("TAS", "MLR", "Quantile", "VSN", "none"), selected = "Quantile"), #), 
                      hr(),
                      #column(4, 
                      h3("Missing value imputation"),
                      selectInput("naImpute", "Choose method", 
                           choices = c("Zero replacement", "Minimum", "Median", "Average", "kNN"), selected = "Zero replacement"),
                      checkboxInput("imputeBeforeNorm", "Impute missing values before normalization", value = FALSE, width = '400px'),
              
         #      column(4, h2("Data transformation"),
          #            checkboxInput("transform", "Log2 transformation (only unchecks if the input data have already log2 transformed)", value = TRUE, width = '400px')),
              # column(4, 
                      hr(),
                      h3("Batch effect"),
                      checkboxInput("sampleBatch", "Remove batch effects (please indicate batch numbers in batch column of the metadata file)", value = FALSE, width = '400px'))
              ),
         #      column(4, h2("Organisms"),
          #            selectInput("organisms", "Choose databases",
           #                 choices = ensembl_dbs)) #c("H.sapiens or M.Musculus", "Others")))
               
                column(3, 
                      wellPanel( h2("UniProt-Gene mapping"),
                                  selectInput("species", "Choose species", choices = c("Human", "Mouse", "Rat"), selected = "Human")
                      ),
                      wellPanel( h2("Correlation"),
                      selectInput("chose_cor", "Choose method", 
                                  choices = c("spearman", "pearson"), selected = "spearman"),
                      selectInput("corrStyle", "Choose style", 
                                  choices = c('simple', 'complex'), selected = 'simple'),
                      h4("Choose options for complex style"),
                      selectInput("corrOrd", "Order", 
                                   choices = c('original', 'hclust', 'alphabet'), selected = 'original'),
                      selectInput("corrUp", "Method: Upper", 
                                  choices = c('circle', 'square', 'ellipse', 'number', 'pie', 'color'), selected = 'number'),
                      selectInput("corrLo", "Method: Lower", 
                                  choices = c('circle', 'square', 'ellipse', 'number', 'pie', 'color'), selected = 'color'),
      #                ),
       #        column(4, 
                      hr(),
                      h2("Heatmap"),
                      selectInput("chose_dist", "Clustering distance", 
                                  choices = c('euclidean', 'correlation', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'), selected = "correlation"), #),
              # column(4, h2("Heatmap"),
                      selectInput("chose_link", "Clustering linkage", 
                                  choices = c('average', 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'mcquitty', 'median', 'centroid'), selected = "average"))
                ),
                column(3, 
                       wellPanel( h2("Enrich plots"),
                                  h4("Select databases for 1-Click enrichment analysis"),
                       checkboxInput("ep_gobp", "GO-biological process", value = TRUE, width = '400px'),
                       checkboxInput("ep_gocc", "GO-cellular component", value = TRUE, width = '400px'),
                       checkboxInput("ep_gomf", "GO-molecular function", value = TRUE, width = '400px'),
                       checkboxInput("ep_kegg", "KEGG pathways", value = TRUE, width = '400px'),
                       checkboxInput("ep_reactome", "Reactome pathways", value = TRUE, width = '400px'),
                       checkboxInput("ep_wikipathways", "WikiPathways", value = TRUE, width = '400px')
                       ))
         
              )))
  
            )
  ) 


ui <- dashboardPage(
  dashboardHeader(title = "SVVATH-D", 
                  dropdownMenuOutput("messageMenu")),
  sidebar,
  body,
  dashboardSidebar(collapsed = TRUE),
  skin = 'purple'
  #setBackgroundColor("ghostwhite")
)


server <- function(input, output, session) {
  
  output$messageMenu <- renderMenu({
    dropdownMenu(type = "messages", headerText = "Need any help?", messageItem(
      from = "Contact me",
      message = "Email: chutipsi@ucmail.uc.edu"
    ), badgeStatus = NULL)
  })
  
  output$contents1 <- DT::renderDataTable({
    #output$contents1 <- renderTable({
    req(input$file1)
    input1 <- input$file1
    group_label_path <- input1$datapath
    #group_label_path <<- group_label_path
    assign("group_label_path", group_label_path, envir = .GlobalEnv) 
    
    group_label_filename <- input$file1$name
    #group_label_filename <<- group_label_filename
    assign("group_label_filename", group_label_filename, envir = .GlobalEnv) 
    
    read.csv(input1$datapath, row.names=1)
  },  options = list(pageLength = 3)
  )
  
  output$contents2 <- DT::renderDataTable({
    req(input$file2)
    input2 <- input$file2
    data_path <- input2$datapath
    data_path <<- data_path
    filename <-  input$file2$name  
    filename <<- filename
    read_excel(data_path, sheet = "Area - peptides")
  },  options = list(pageLength = 25)
  )
  
  observeEvent(input$cutFC, {
    cutFC <- input$cutFC
    cutFC <<- cutFC
  })
  
  observeEvent(input$cutP, {
    cutP <- as.numeric(input$cutP)
    cutP <<- cutP
  })
  
  
  
  ###### 1-Click ############
  observeEvent(input$go, {
    if(!exists('filename') | !exists('group_label_filename')){ 
      print("no file upload properly!") 
    } else {    # debug1
    
    withProgress(message = "Running", value = 0, {
      Sys.sleep(0.25)
      
      #dest <- paste0("~/Desktop/swathD2/", unlist(str_split(filename, "\\."))[1], "_", format(Sys.time(), "%Y%m%d%H%M"), "/")
      dest <- paste0("./", unlist(str_split(filename, "\\."))[1], "_", format(Sys.time(), "%Y%m%d%H%M"), "/")
      dir.create(dest)  
      setwd(dest)
      log_session <- paste("input =", filename, "; label =", group_label_filename, "; threshold FC=", cutFC, ", p-value =", cutP)
      log_session <<- log_session
      #write.csv(log_session, file = "log_session.csv")
      write.csv(log_session, file = "log_session.csv")
      
      group_label <- read.csv(group_label_path, row.names=1)
      assign("group_label", group_label, envir = .GlobalEnv)
      
      group <- group_label$group
      group <<- group
      sample_label <- as.character(group_label$sample_label)
      sample_label <<- sample_label

#-----UniProt accession ID convert to gene symbol -----------------------          
      
    #  source("lib/swathD_idConv.R", local = TRUE)
      dataformat <<- input$dataFormat  
      species <- as.character(input$species)
      assign("species", species, envir = .GlobalEnv)
      swathD_idConv(dataformat)
      
      idConvert_exprRaw <- cbind(id_all, expr_all)
      
      dir.create( paste0(dest, "datasets") )
    #  setwd( paste0(dest, "datasets") ) 
      write.csv(id_all, file = "idConvert.csv")
      #write.csv(id_all, file = paste0(dest, "datasets/idConvert.csv"))
      write.csv(idConvert_exprRaw, file = "idConvert_exprRaw.csv")
     # write.csv(idConvert_exprRaw, file = paste0(dest, "datasets/idConvert_exprRaw.csv"))
      print("Get raw data with gene sympbols")
      
      incProgress(0.1, message = "Data cleaning and preprocessing")
      
#------------------------------------------------------------------------

#----- Missing value imputation: if users want to impute missing values before normalization (e.g., for datasets with lots of NA) --------------------------------------
      
      if(input$imputeBeforeNorm){
      incProgress(0.1, message = "Impute missing values before normalization")  
      imputeBeforeNorm <<- input$imputeBeforeNorm
      naImpute <<- input$naImpute
     # source("lib/swathD_imputeFirst.R", local = TRUE)
      swathD_imputeFirst()
      }

#-----Expression data normalization --------------------------------      
      #source("lib/compare_norm.R", local = TRUE)
      compare_norm()
      choose_norm <<- input$normMet
      
#--------------------------------------------------------------------      
       
#----- Missing value imputation--------------------------------------
      #source("lib/swathD_naImpute.R", local = TRUE)
      naImpute <<- input$naImpute
      swathD_naImpute()
      
#----- Remove batch effect-------------------------------------------     
       
       if(input$sampleBatch){
      # source("lib/swathD_adjBatch.R", local = TRUE)  
         incProgress(0.1, message = "Remove batch effects")
         swathD_adjBatch(df3)
         print("Finish - Remove batch effects") 
       } else if(!(input$sampleBatch)){
         expr_processed <- df3
         assign("expr_processed", expr_processed, envir = .GlobalEnv)
       }
       
#--------------------------------------------------------------------       
      
       
      process_ds <- cbind(id_all, expr_processed) 
      process_ds <<- process_ds
      
      write.csv(process_ds, file = "process_ds.csv") 
     # write.csv(process_ds, file = paste0(dest, "datasets/process_ds.csv"))

      print("Get processed dataset!") 
      
      log_tf <- t(expr_processed)
      colnames(log_tf) <- id_all$gene.SYMBOL
      log_ds <- data.frame(group, log_tf)
      log_ds <<- log_ds 
      
      print("Finish - Data preprocessing")

######################################################################################      
      
 ###################################      
# Data analysis and visualization #
##################################
      
#-----Calculate fold changes and p-values ----------------------------------------------    
      
      incProgress(0.1, message = "Calculate fold changes and p-values")
      
      tmp <- data.frame(group = log_ds[ , 1], 2^log_ds[ , 2:length(log_ds)]) %>% 
        gather(gene.SYMBOL, expression, -group) %>%
        dplyr::group_by(group, gene.SYMBOL) %>% 
        dplyr::summarize(group_mean = mean(expression)) %>%
        spread(gene.SYMBOL, group_mean)

      
      gr_avr <- as.data.frame(tmp[ , 2:length(tmp)]) # for CV estimation later
      rownames(gr_avr) <- tmp$group
      
      ## generate group-pairs using combn function
      gr_pair <- combn(unique(tmp$group), 2)
      
      fc <- (gr_avr[gr_pair[1, ], ] / gr_avr[gr_pair[2, ], ]) %>% log2()
      
      rownames(fc) <- paste0('log2', '(', gr_pair[1, ], '/', gr_pair[2, ], ')')
      
      log2fc_ds <<- fc # for the record
      
      ## fit ANOVA
      proteins <- as.matrix(log_ds[, 2:length(log_ds)])
      fit.aov <- aov(proteins ~ group) # call group object
      output.aov <- summary.aov(fit.aov) 
      
      ## get ANOVA p-value for all proteins
      anova.pVal <- numeric(length = ncol(proteins))
      for (i in 1:length(output.aov)){
        anova.pVal[i] <- output.aov[[i]][1, 5]
      }
      
      ## get adjusted p-value from multiple pair-wise comparisons using Tukey-HSD post hoc
      ### calculate the first variable to get numbers of pairwise comparison
      tukey.dummy <- TukeyHSD(aov(proteins[ ,1] ~ group))
      
      ## generate a blank matrix; row = number of proteins, col = number of pairwise comparison obtained from tukey.dummy
      adj.pVal <- matrix(nrow = ncol(proteins), ncol = nrow(tukey.dummy[[1]]))
      
      ## label the blank column
      colnames(adj.pVal) <- c(rownames(tukey.dummy[[1]]))
      rownames(adj.pVal) <- colnames(proteins)
      
      ## get adjusted p-value after pairwise comparison from the column 4 of 'group' object
      for (i in 1:ncol(proteins)){ 
        adj.pVal[i, ] <- (TukeyHSD((aov(proteins[, i] ~ group))))[[1]][ ,4] 
      }
      
      ## merge ANOVA p-value with Tukey-HSD for each protein
      anova_ds <- cbind(anova.pVal, adj.pVal)
      anova_ds <<- anova_ds # for the record
      write.csv(anova_ds, file = "anova_ds.csv")
      #write.csv(anova_ds, file = paste0(dest, "datasets/anova_ds.csv"))
      
      print("Finish - ANOVA with Tukey's posthoc")
      
#------------------------------------------------------------------------------
      
      
      
      ################
      # Quality check #
      ################
      # calculate %CV using a formula: %CV = 100*sd/mean
      ## already has gr_avr for mean value, therefore just need to have gr_sd to fulfill the equation
      #gr_sd <- aggregate(log_ds[ , 2:length(log_ds)], list(group = log_ds$group), sd)
      
      incProgress(0.2, message = "Check data quality")
      
      #observe
    #  tmp <- data.frame(log_ds) %>%
      tmp <-  data.frame(group = log_ds[ , 1], 2^log_ds[ , 2:length(log_ds)]) %>%
        gather(gene.SYMBOL, expression, -group) %>%
        dplyr::group_by(group, gene.SYMBOL) %>% 
        dplyr::summarize(group_sd = sd(expression)) %>%
        spread(gene.SYMBOL, group_sd)
      gr_sd <- as.data.frame(tmp[ , 2:length(tmp)])
      rownames(gr_sd) <- tmp$group
      qc <- 100 *gr_sd/gr_avr # convert to non-log transformed data before CV calculation
      qc <- data.frame(group = tmp$group, qc)
      
      #qc <- 100 *(gr_sd[, 2:length(gr_sd)]/gr_avr[, 2:length(gr_avr)])
      #qc <- data.frame(group = gr_sd$group, qc)
      
      qc_ds <<- qc # for the record
      
      ## export anova dataset as .csv file
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
      #write.csv(qc_ds, file = "qc_ds.csv")
      
      ## prepare data for violin plot
      #qc_vector <- as.numeric(t(qc_ds[-1]))
      #gr_vector <- rep(qc_ds$group, each = length(qc_ds[-1]))
      
      #QC <- data.frame(group = gr_vector, CV = qc_vector)
      #QC <<- QC
      #m <- aggregate(CV~group, QC, function(i) round(median(i), 2))
      
      QC <- qc_ds %>% gather(gene, CV, -group)
      # calculate median CV of each group
      m <- QC %>% dplyr::group_by(group) %>% dplyr::summarise(CV = round(median(CV), 2))
      
      ## Violin plot 
      plot.qc <- ggplot(QC, aes(x=group, y=CV)) + 
        #scale_fill_brewer(palette="Blues")+
        #scale_fill_manual(values=c("#afceff", "#afceff"))+
        #scale_fill_manual(values=c("#7cb25b", "#db3f3f"))+ 
        geom_violin(aes(fill = group), trim=FALSE, width = 0.8, na.rm = TRUE, position = "dodge")+
        geom_boxplot(width=0.1, fill = 'white', outlier.size = 0.5, na.rm = TRUE, position = "dodge")+
        geom_text(data = m, aes(label = CV), position = position_dodge(width = 1), hjust = -0.5, vjust = -0.5, size = 5) +
        xlab("") + ylab("% Coefficient of Variation") +
        scale_y_continuous(breaks=c(0, 10, 20, 50, ceiling(max(QC$CV, na.rm=TRUE)))) +
        theme_light(base_size = 12) 
      #theme(axis.text.x = element_text(face="bold", color = "black", size=12),
      #     axis.text.y = element_text(face="bold", color = "black", size=12), 
      #    axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))
      #plot.qc <<- plot.qc
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
      plot.qc <<- plot.qc
    #  dir.create( paste0(dest, "quality") )
     # setwd( paste0(dest, "quality") ) 
      pdf("QC_violinPlot.pdf", width = 5, height = 4)
      print(plot.qc)
      dev.off()
      #dev.print(pdf, 'QC_violin.pdf')
      print("Get QC_violinPlot: Showing a distribution of the coefficient of variation of all variables in each sample group")
      
      
      ########################
      # correlation heatmap #
      ######################
      
      chose_cor <- as.character(input$chose_cor)
      #chose_cor <- 'spearman'
      #corr <- tas_ds[ , 5:length(tas_ds)]
      #corr <- 2^process_ds[, 5:length(process_ds)]
      corr <- 2^expr_processed
      # corr <- round(cor(corr, method = 'pearson'),3)
      corr <- round(cor(corr, method = chose_cor),3)
      
      if(input$corrStyle == 'simple'){
      corr[lower.tri(corr)] <- NA
      melted_corr <- melt(corr, na.rm = TRUE)
      corrHM <- ggplot(data = melted_corr, aes(x = Var2, y = Var1, fill = value))+  
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "white", high = "#ed4c07", mid = "yellow",   
                             midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                             name= paste(chose_cor, "\ncorrelation") ) +
        labs(x = "", y = "") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                         size = 12, hjust = 1)) +
        coord_fixed()
      
      plot_corrHM <- corrHM + 
        geom_text(aes(label = value), color = "black", size = 1) + #Var2, Var1, 
        theme(
        #  axis.title.x = element_blank(),
         # axis.title.y = element_blank(),
          #axis.text.x = element_text(face="bold", color = "black", size=12),
          axis.text.y = element_text(color = "black", size=11),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
        #  axis.ticks = element_blank(), 
       #   legend.justification = c(1, 0),
        #  legend.position = c(0.5, 0.7),
         # legend.direction = "horizontal")+
       # guides(fill = guide_colorbar(barwidth = 8, barheight = 1,
        #                             title.position = "top", title.hjust = 0.5))
      plot_corrHM <<- plot_corrHM
      } else if(input$corrStyle == 'complex'){
      #col <- colorRampPalette(c("#053061", "#2166AC", "#2166AC", "#92C5DE", "#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))
      #col2 <- colorRampPalette(c('#003d66', '#006bb3', '#0099ff', '#4db8ff', 'white', "#ff9999", "#ff4d4d", "#ff0000", '#cc0000'))
      corrplot(corr, type= 'upper', order = as.character(input$corrOrd), method = as.character(input$corrUp), #col = col2(200), 
                              tl.col="magenta", tl.srt=45, tl.pos = "d")
      corrplot(corr, add = TRUE, type = "lower", method = as.character(input$corrLo), order = as.character(input$corrOrd), #col = col2(200),
                        diag = FALSE, tl.pos = "n", cl.pos = 'n')
      plot_corrHM <- recordPlot()
      plot_corrHM <<- plot_corrHM
      }
      
      plot_corrHM <<- plot_corrHM
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
      pdf("corrHeatmap.pdf", width = 9, height = 9)
      print(plot_corrHM)
      dev.off()
      #dev.print(pdf, 'corrHM.pdf')
      print("Get corrHeatmap: A plot showing a correlation of expression profiles among all samples")
      
      ###################################
      # Number of peptides per protein #
      #################################
      
      n_pept_prot <- areaPept %>% 
        dplyr::group_by(Protein) %>% 
        dplyr::summarize(n_pept = n()) %>% 
        arrange(desc(n_pept))
      
      nPP <- data.frame(n_pept = c("1", "2-5", "6-10"), 
                        n_prot = rbind(n_pept_prot %>% filter(n_pept ==1) %>% nrow(),
                                       n_pept_prot %>% filter(n_pept >=2 & n_pept <= 5) %>% nrow(), 
                                       n_pept_prot %>% filter(n_pept >=6) %>% nrow()))
      
      nPP_plot <- ggplot(nPP, aes(x = n_pept, y= n_prot)) + 
        geom_bar(stat = "identity", fill = "steelblue") + 
        ylim(0, max(nPP$n_prot)+50) +
        geom_text(aes(label= n_prot), vjust=-0.3, color="black", size=4.5) +
        geom_text(aes(label= paste0(round(100*n_prot/sum(n_prot), 1), "%")), vjust=1.6, color="white", size=4.5) +
        xlab("Numbers of peptide") + ylab("Numbers of protein") +
        theme_light(base_size = 12)
      nPP_plot <<- nPP_plot
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
      pdf("number_pept_prot.pdf", width = 4, height = 3)
      print(nPP_plot)
      dev.off()
      print("Get number_pept_prot: A plot showing the numbers of peptides per protein")
      
      #######
      # PCA #
      ######
      
      #df_pca <- zScore_ds
      #fit_pca = PCA(df_pca[ , 2:length(df_pca)], graph = FALSE, scale.unit = FALSE)
      
      fit_pca <- PCA(log_ds[ , 2:length(log_ds)], graph = FALSE, scale.unit = TRUE, ncp = 5)
      percentage <- fit_pca$eig[ , 2]
      
      PCs <- data.frame(fit_pca$ind$coord)
      #PCs <- data.frame(fit_pca$ind$cos2)
      PCs$group <- group
      
      plotPCA <- ggplot(data = PCs, aes(x = Dim.1, y = Dim.2)) +
        geom_point(aes(colour = group), size = 3) +
        xlab(paste0('PC1', ' ', '(', round(percentage[1], 2), '%)')) + 
        ylab(paste0('PC2', ' ', '(', round(percentage[2], 2), '%)')) +
        scale_fill_hue(l=40) + 
        coord_fixed(ratio=1, xlim=range(PCs$Dim.1), ylim=range(PCs$Dim.2)) +
        geom_text_repel(label = rownames(PCs)) +
        theme_light()
      
      assign("plotPCA", plotPCA, envir = .GlobalEnv) 
      #plotPCA <<- plotPCA
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
      pdf("PCA.pdf", width = 5, height = 3.5)
      print(plotPCA)
      dev.off()
      print("Get PCA: Showing the unsupervised classification of sample groups")
      
      ##############################
      # protein abundance heatmap #
      ############################
      
      # use tas_ds as the input
     # tmp <- tas_ds[ , 4:length(tas_ds)] 
      
     # tmp <- data.frame(process_ds$gene.SYMBOL , 2^process_ds[, 5:length(process_ds)])
      tmp <- data.frame(process_ds$gene.SYMBOL , 2^expr_processed)
      
      qc_hm <- as.matrix(tmp[-1])
      rownames(qc_hm) <- tmp$gene.SYMBOL
      
      # keep missing value = 0 after log10 transform
      for(i in seq_along(qc_hm)){
        if(qc_hm[i] != 0){
          qc_hm[i] <- log10(qc_hm[i])
        } else {
          qc_hm[i] <- 0
        }}
      
      n_missing <- sum(qc_hm == 0)
      n_total <- dim(qc_hm)[1] * dim(qc_hm)[2]
      
      qc_hm_plot <- pheatmap(t(qc_hm),  breaks = seq(0, max(qc_hm), length.out=101), 
                             legend_breaks=seq(0, round(max(qc_hm), 0), length.out=8), 
                             legend_labels = c("1e+00", "1e+01", "1e+02", "1e+03", "1e+04", "1e+05", "1e+06", "1e+07"),
                             color = colorRampPalette(c("black", "#8ea1ff", "#14ff57", "yellow", "orange", "#ea4444"))(100), #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                             border_color = "gray",
                             #annotation_col = data.frame(group, row.names = sample_label),
                             #clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "maximum", #'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                             #clustering_method_rows = "average", 
                             clustering_method_columns = "complete",
                             cluster_rows = FALSE, #cluster_cols = TRUE,
                             #cellheight = 9, #cellwidth = 1.2, 
                             fontsize_row = 8, fontsize_col = 1, 
                             scale = "none",
                             main = paste0("Protein abundance heatmap (total ", n_total, " data points; ", n_missing, " missing values (", round(100*n_missing/n_total, 2), "%) showed in black)") )
      qc_hm_plot <<- qc_hm_plot
      pdf("qc_heatmap.pdf", width = 12, height = 4)
      print(qc_hm_plot)
      dev.off()
      print("Get qc_heatmap: Showing protein abundace across all samples")
      
      incProgress(0.2, message = "Prepare Volcano plot and heatmap")
      
      ################
      # Volcano plot #
      ###############
      # prepare data for volcano plot
      tmp <- data.frame(gene = rownames(anova_ds), anova_ds) 
      colnames(tmp) <- c("gene", "anova.pVal", paste0(gr_pair[1, ], "/", gr_pair[2, ]) )
      long_ano <- gather(tmp, compare, adj_pVal, -gene, -anova.pVal)
      
      #fc.wide <- log2fc_ds
      ##rownames(fc.wide) <- paste(gr_pair[1, ], "vs", gr_pair[2, ])
      #rownames(fc.wide) <- paste0(gr_pair[2, ], ".", gr_pair[1, ])
      #fc.vp <- data.frame(gene = colnames(fc.wide), t(fc.wide))
      
      fc.vp <- t(log2fc_ds)
      fc.vp <- data.frame(gene = colnames(log2fc_ds), fc.vp)
      colnames(fc.vp) <- c("gene", paste0(gr_pair[1, ], "/", gr_pair[2, ]) )
      long_fc <- gather(fc.vp, compare, log2FC, -gene)
      
      # left_join long_fc to long_ano using gene and compare as a key
      long_ano.fc <- long_ano %>% 
        left_join(long_fc, by = c("gene", "compare"))
      long_ano.fc$gene <- as.character(long_ano.fc$gene)
      long_ano.fc <<- long_ano.fc
      
      #### update the select input for replot: Volcano ######      
      #pair_all <- unique(long_ano.fc$compare)
      #pair_all <<- pair_all
      updateSelectInput(session, "pair_select", choices = unique(long_ano.fc$compare))
      
      ##########
      
      volcano_all <- ggplot(data = long_ano.fc, aes(x= log2FC, 
                                                    y=-log10(adj_pVal))) +
        #                     # geom_point(aes(color = as.factor(abs(log2FC) > log2(2) & anova.pVal < 0.05 & adj_pVal < 0.05)), size = 1) +
        geom_point(aes(color = as.factor(abs(log2FC) >= log2(cutFC) & anova.pVal < cutP & adj_pVal < cutP)), size = 1, alpha = 0.5, show.legend = FALSE) +
        scale_color_manual(values = c("grey", "red")) +
      #  theme(legend.position = "none") +
        # xlim(c(-6, 6)) + #+ ylim(c(0, 5)) +
        #scale_x_continuous(breaks = c(-6, -3, -1, 0, 1, 3, 6)) +
        xlab("log2 (fold change)") + ylab("-log10 (adjusted p-value)") +
        ggtitle(label = paste0("Volcano plot at ", cutFC, "x fold change and adj.P-value < ", cutP)) + 
        theme_grey(base_size = 15) +
       # theme(plot.title = element_text(lineheight= 0.8, face="bold")) +
        # add gene name labelling to points that their fc > absolute median + 1SD of all fc values and their p-value < 80% of the least p-value
   #     geom_text_repel(data = (subset(long_ano.fc, 
    #                                   abs(log2FC) > (median(abs(log2FC)) + sd(abs(log2FC))*0.5) & 
     #                                    log10(adj_pVal) < min(log10(adj_pVal))*0.8)),
      #                  aes(label = gene, size = 0.25),
       #                 box.padding = unit(0.35, "lines"),
        #                point.padding = unit(0.3, "lines")) +
        facet_wrap(~ compare)
      
      volcano_all <<- volcano_all
      
      #vp <- data.frame(t(log2fc_ds), anova_ds)
      #gene = as.character(rownames(vp))
      #vp <- data.frame(vp, gene = gene)
      
      ## export volcano plot
      #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
    #  dir.create( paste0(dest, "results") )
     # setwd( paste0(dest, "results") ) 
      pdf("volcano_all.pdf", width = 12, height = 4)
      print(volcano_all)
      dev.off()
      print("Get volcano_all: Volcano plot(s) of all pairwise comparisons - adjusted p-values from ANOVA with Tukey's posthoc")
      
      ########################
      # significant heatmap #
      ######################
      
      # use log_ds as the input, so heatmap represents data in log2 
      dfmat <- as.matrix(log_ds[ , 2:length(log_ds)])
      
      # calculate median for each row (protein) 
        # med <- apply(t(dfmat), 1, median)
      med <- apply(t(dfmat), 1, mean)
      # med <- rowMeans(t(dfmat))
      # scale to protein median
      dfmat_medScale <- (t(dfmat) - med)
      dfmat_medScale_all <<- dfmat_medScale # use later in re-Heatmap
      
      #dfvar <- dfmat_medScale
      
      print("Finish - Median scaling of each variable across all samples")
      
      tmp <- anova_ds[, 1]
      dfmat_medScale <- data.frame(dfmat_medScale, 
                                   anova_pVal = tmp, 
                                   gene = rownames(dfmat_medScale))
      
      ######## sig.peptides dataset ##########
      
      dfmat_medScale_sig <- dfmat_medScale %>% filter(anova_pVal < cutP)
      rownames(dfmat_medScale_sig) <- dfmat_medScale_sig$gene
      dfmat_medScale_sig <- dfmat_medScale_sig[, 1: (length(dfmat_medScale_sig) - 2)]
      
      dfmat_medScale_sig <<- dfmat_medScale_sig # use later in re-Heatmap
      
      nprot_sig <- nrow(dfmat_medScale_sig)
      
      chose_dist <- as.character(input$chose_dist)
      chose_dist <<- chose_dist
      chose_link <- as.character(input$chose_link)
      chose_link <<- chose_link
      
      hm_sig <- pheatmap(dfmat_medScale_sig, breaks = seq(-(max(round(dfmat_medScale_sig, 0))), max(round(dfmat_medScale_sig, 0)), length.out=101), #breaks = seq(-3, 3, length.out=101), 
                         legend_breaks=seq(-(max(round(dfmat_medScale_sig, 0))), max(round(dfmat_medScale_sig, 0)), length.out=5),
                         #legend_labels = c("-1", "1e+02", "0", "1e+04", "1"),
                         color = colorRampPalette(c("darkblue", "blue", "white", "orangered", "red"))(100),  #color = colorRampPalette(c("darkblue", "cornflowerblue", "white", "orangered", "darkred"))(100),
                         border_color = NA,
                         # annotation_col = data.frame(group, row.names = sample_label),
                         annotation_col = data.frame(group = factor(group), row.names = sample_label),
                         clustering_distance_rows = chose_dist,
                         clustering_distance_cols = chose_dist, #'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
                         clustering_method = chose_link, # 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
                         #clustering_method_columns = "ward",
                         fontsize_row = 2, fontsize_col = 10, 
                         scale = "none",
                         main = paste0(nprot_sig, " proteins with anova_pVal < ", cutP, ' (color in Log2-scale)'))
      hm_sig <<- hm_sig
      
      pdf("heatmap_sig.pdf", width = 9, height = 9)
      print(hm_sig)
      dev.off()
      print("Get heatmap_sig: A heatmap of significant variables based on ANOVA with Tukey's posthoc; Clustering by euclidean distance and average linkage")
      
      incProgress(1, message = "Finish QC-significant analysis")
      Sys.sleep(0.25)
    })
    
    ### output choice: run ORA when select QC-sig-enrichPlot ####      
    observeEvent(input$go, {
      if(input$outplot == "QC-sig-enrichPlot") {
        withProgress(message = "Start enrichr", value = 0, {
          Sys.sleep(0.1)
          
          #######
          # ORA #
          ######
          print("Overrepresentation gene enrichment analysis using Enrichr")
          print("Computational time depends on the complexity of input dataset and the responsiveness of Enrichr server")
          
          # subsetting of long.ano.fc to get genes with 2x foldchange and adj.p < 0.05
          tmp <- long_ano.fc %>% 
            filter(abs(log2FC) >= log2(2) & 
                     adj_pVal < 0.05) %>% 
            arrange(compare, adj_pVal)  
          
          # get a list of sig.genes for each pairwise comparison
          tmp <- split(tmp, factor(tmp$compare)) %>% 
            lapply(dplyr::select, "gene")  
          
          pair <<- names(tmp) # keep names of each pairwise compaison for later use
          
          # get all the first column (gene) of every df as the character vectors
          sigProtList <- sapply(tmp, `[`, 1) 
          print("Preparing inputs (a list of significant proteins) from all pairwise comparisons")
          print("Start Enrichr processing")
          
          # objects for adding a row to each dataframe in the lists; preventing errors due to ennrichr 0 matched 
          enrich_row_dummy <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)
          enrich_df_colnames <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Z.score", "Combined.Score", "Genes")
          
          incProgress(0.1, message = "Enrichment analysis: submit data")
          
          # get GOMF, the last row is the dummy NA 
          enrich_mf <- lapply(sigProtList, enrichr, "GO_Molecular_Function_2015") 
          enrich_mf <- reduce(enrich_mf, union) 
          names(enrich_mf) <- pair
          enrich_mf <- lapply(enrich_mf, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_mf <- mapply(cbind, enrich_mf, "class" = "mf", "compare" = pair,  SIMPLIFY=F)
          
          incProgress(0.15, message = "Enrichment analysis: GO-terms")
          # get GOCC, the last row is the dummy NA 
          enrich_cc <- lapply(sigProtList, enrichr, "GO_Cellular_Component_2015") 
          enrich_cc <- reduce(enrich_cc, union) 
          names(enrich_cc) <- pair
          enrich_cc <- lapply(enrich_cc, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_cc <- mapply(cbind, enrich_cc, "class" = "cc", "compare" = pair, SIMPLIFY=F)
          
          incProgress(0.15, message = "Enrichment analysis: GO-terms")
          
          # get GOBP, the last row is the dummy NA 
          enrich_bp <- lapply(sigProtList, enrichr, "GO_Biological_Process_2015") 
          enrich_bp <- reduce(enrich_bp, union) 
          names(enrich_bp) <- pair
          enrich_bp <- lapply(enrich_bp, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_bp <- mapply(cbind, enrich_bp, "class" = "bp", "compare" = pair,  SIMPLIFY=F)
          
          incProgress(0.1, message = "Enrichment analysis: pathways")
          
          # get KEGG, the last row is the dummy NA 
          enrich_kegg <- lapply(sigProtList, enrichr, "KEGG_2016")
          enrich_kegg <- reduce(enrich_kegg, union) 
          names(enrich_kegg) <- pair
          enrich_kegg <- lapply(enrich_kegg, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_kegg <- mapply(cbind, enrich_kegg, "class" = "kegg", "compare" = pair,  SIMPLIFY=F)
          
          incProgress(0.1, message = "Enrichment analysis: pathways")
          
          # get Reactome, the last row is the dummy NA 
          enrich_react <- lapply(sigProtList, enrichr, "Reactome_2016")
          enrich_react <- reduce(enrich_react, union) 
          names(enrich_react) <- pair
          enrich_react <- lapply(enrich_react, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_react <- mapply(cbind, enrich_react, "class" = "reactome", "compare" = pair,  SIMPLIFY=F)
          
          incProgress(0.1, message = "Enrichment analysis: pathways")
          
          # get WikiPathway, the last row is the dummy NA 
          enrich_wiki <- lapply(sigProtList, enrichr, "WikiPathways_2016")
          enrich_wiki <- reduce(enrich_wiki, union) 
          names(enrich_wiki) <- pair
          enrich_wiki <- lapply(enrich_wiki, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          enrich_wiki <- mapply(cbind, enrich_wiki, "class" = "wikiPathway", "compare" = pair,  SIMPLIFY=F)
          
          ###### for additional protein enrichment plot: enrich #####    
          # get ENCODE_TF, the last row is the dummy NA 
          #    enrich_encode <- lapply(sigProtList, enrichr, "ENCODE_TF_ChIP-seq_2015")
          #    enrich_encode <- reduce(enrich_encode, union) 
          #    names(enrich_encode) <- pair
          #    enrich_encode <- lapply(enrich_encode, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_encode <- mapply(cbind, enrich_encode, "class" = "ENCODE_TF", "compare" = pair,  SIMPLIFY=F)
          
          # get ChEA, the last row is the dummy NA 
          #    enrich_chea <- lapply(sigProtList, enrichr, "ChEA_2016")
          #    enrich_chea <- reduce(enrich_chea, union) 
          #    names(enrich_chea) <- pair
          #    enrich_chea <- lapply(enrich_chea, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_chea <- mapply(cbind, enrich_chea, "class" = "ChEA", "compare" = pair,  SIMPLIFY=F)
          
          # get KEA, the last row is the dummy NA 
          #    enrich_kea <- lapply(sigProtList, enrichr, "KEA_2015")
          #    enrich_kea <- reduce(enrich_kea, union) 
          #    names(enrich_kea) <- pair
          #    enrich_kea <- lapply(enrich_kea, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_kea <- mapply(cbind, enrich_kea, "class" = "KEA", "compare" = pair,  SIMPLIFY=F)
          
          # get TF_PPI (Transcription_Factor_Protein_Protein Interactions), the last row is the dummy NA 
          #    enrich_tfpi <- lapply(sigProtList, enrichr, "Transcription_Factor_PPIs")
          #    enrich_tfpi <- reduce(enrich_tfpi, union) 
          #    names(enrich_tfpi) <- pair
          #    enrich_tfpi <- lapply(enrich_tfpi, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_tfpi <- mapply(cbind, enrich_tfpi, "class" = "TFPI", "compare" = pair,  SIMPLIFY=F)
          
          # get CORUM (The comprehensive resource of mammalian protein complexes), the last row is the dummy NA 
          #    enrich_corum <- lapply(sigProtList, enrichr, "CORUM")
          #    enrich_corum <- reduce(enrich_corum, union) 
          #    names(enrich_corum) <- pair
          #    enrich_corum <- lapply(enrich_corum, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_corum <- mapply(cbind, enrich_corum, "class" = "CORUM", "compare" = pair,  SIMPLIFY=F)
          
          # get PPI (PPI_Hub_Proteins), the last row is the dummy NA 
          #    enrich_ppi <- lapply(sigProtList, enrichr, "PPI_Hub_Proteins")
          #    enrich_ppi <- reduce(enrich_ppi, union) 
          #    names(enrich_ppi) <- pair
          #    enrich_ppi <- lapply(enrich_ppi, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_ppi <- mapply(cbind, enrich_ppi, "class" = "PPI", "compare" = pair,  SIMPLIFY=F)
          
          # get DGA (Disease-gene associations mined from literature), the last row is the dummy NA 
          #    enrich_dga <- lapply(sigProtList, enrichr, "Jensen_DISEASES")
          #    enrich_dga <- reduce(enrich_dga, union) 
          #    names(enrich_dga) <- pair
          #    enrich_dga <- lapply(enrich_dga, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_dga <- mapply(cbind, enrich_dga, "class" = "DGA", "compare" = pair,  SIMPLIFY=F)
          
          # get HPO (Human_Phenotype_Ontology), the last row is the dummy NA 
          #    enrich_hpo <- lapply(sigProtList, enrichr, "Human_Phenotype_Ontology")
          #    enrich_hpo <- reduce(enrich_hpo, union) 
          #    names(enrich_hpo) <- pair
          #    enrich_hpo <- lapply(enrich_hpo, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_hpo <- mapply(cbind, enrich_hpo, "class" = "HPO", "compare" = pair,  SIMPLIFY=F)
          
          # get dbGaP (database for Genotype and Phenotype), the last row is the dummy NA 
          #    enrich_gap <- lapply(sigProtList, enrichr, "dbGaP")
          #    enrich_gap <- reduce(enrich_gap, union) 
          #    names(enrich_gap) <- pair
          #    enrich_gap <- lapply(enrich_gap, rbind, enrich_row_dummy) %>% lapply(setNames, enrich_df_colnames)
          #    enrich_gap <- mapply(cbind, enrich_gap, "class" = "dbGaP", "compare" = pair,  SIMPLIFY=F)
          
          ###########       
          
          ########## tmp 6 class ###########
          incProgress(0.1, message = "Prepare enrich table")
          
          tmp <- list(enrich_bp, enrich_cc, enrich_mf, 
                      enrich_kegg, enrich_react, enrich_wiki)  
          #     enrich_encode, enrich_chea, enrich_kea,
          #    enrich_tfpi, enrich_corum, enrich_ppi,
          #   enrich_dga, enrich_hpo, enrich_gap)
          names(tmp) <- c("GOBP", "GOCC", "GOMF", 
                          "KEGG", "Reactome", "WikiPathway") 
          #        "ENCODE_TF", "ChEA", "KEA",
          #       "TFPI", "CORUM", "PPI",
          #      "DGA", "HPO", "dbGaP")
          tmp <- unlist(tmp, recursive = FALSE, use.names = TRUE)
          enrich_master <- suppressWarnings(bind_rows(tmp)) %>% filter(!is.na(Term))
          enrich_master <<- enrich_master
          enrich_top10 <- lapply(tmp, arrange, Adjusted.P.value) %>% lapply(head, 10)
          enrich_top10 <- suppressWarnings(bind_rows(enrich_top10)) %>% filter(!is.na(Term))
          
          enrich_top10 <- separate(enrich_top10, Overlap, c("geneInput", "geneTotal"), sep = "/") 
          enrich_top10$geneInput <- as.numeric(enrich_top10$geneInput)
          enrich_top10$geneTotal <- as.numeric(enrich_top10$geneTotal)
          enrich_top10 <- enrich_top10 %>% mutate("geneRatio" = (geneInput / geneTotal)) %>% arrange(compare, class, Adjusted.P.value)
          enrich_top10 <<- enrich_top10
          
          write.csv(enrich_top10, file = "enrich_top10.csv")
          print("Get top10 gene enrichment analysis table")
          
          incProgress(0.1, message = "Prepare enrich plots")
          
          # get BP 
          BP <- enrich_top10 %>% filter(class == "bp")
          BP_plot <- ggplot(BP, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=BP$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "GOBP") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          BP_plot <<- BP_plot
          
          # get CC
          CC <- enrich_top10 %>% filter(class == "cc")
          CC_plot <- ggplot(CC, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=CC$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "GOCC") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          CC_plot <<- CC_plot
          
          # get MF
          MF <- enrich_top10 %>% filter(class == "mf")
          MF_plot <- ggplot(MF, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=CC$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "GOMF") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          MF_plot <<- MF_plot
          
          # get KEGG
          KEGG <- enrich_top10 %>% filter(class == "kegg")
          KEGG_plot <- ggplot(KEGG, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "#5656f7", high = "#63eaff", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=CC$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "KEGG") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          KEGG_plot <<- KEGG_plot
          
          # get Reactome
          Reactome <- enrich_top10 %>% filter(class == "reactome")
          Reactome_plot <- ggplot(Reactome, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=CC$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "Reactome") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          Reactome_plot <<- Reactome_plot
          
          # get WikiPathway
          WikiPathway <- enrich_top10 %>% filter(class == "wikiPathway")
          WikiPathway_plot <- ggplot(WikiPathway, aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 30), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
            #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
            scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
            #scale_x_discrete(limits=CC$Term) +
            labs(y = "geneRatio") +
            labs(x = "") +
            ggtitle(label = "WikiPathway") + 
            coord_flip() +
            facet_wrap(vars(compare), nrow = 1)
          WikiPathway_plot <<- WikiPathway_plot
          
          ##### for additional protein enrichment plot: Plot  ####    
          # get ENCODE_TF
          #    ENCODE_TF <- enrich_top10 %>% filter(class == "ENCODE_TF")
          #    ENCODE_TF_plot <- ggplot(ENCODE_TF, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "ENCODE_TF") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    ENCODE_TF_plot <<- ENCODE_TF_plot
          
          # get ChEA
          #    ChEA <- enrich_top10 %>% filter(class == "ChEA")
          #    ChEA_plot <- ggplot(ChEA, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "ChEA") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    ChEA_plot <<- ChEA_plot
          
          # get KEA
          #    KEA <- enrich_top10 %>% filter(class == "KEA")
          #    KEA_plot <- ggplot(KEA, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "KEA") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    KEA_plot <<- KEA_plot
          
          # get TFPI
          #    TFPI <- enrich_top10 %>% filter(class == "TFPI")
          #    TFPI_plot <- ggplot(TFPI, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "TFPI") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    TFPI_plot <<- TFPI_plot
          
          # get CORUM
          #    CORUM <- enrich_top10 %>% filter(class == "CORUM")
          #    CORUM_plot <- ggplot(CORUM, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "CORUM") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    CORUM_plot <<- CORUM_plot
          
          # get PPI
          #    PPI <- enrich_top10 %>% filter(class == "PPI")
          #    PPI_plot <- ggplot(PPI, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "PPI") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    PPI_plot <<- PPI_plot
          
          # get DGA
          #    DGA <- enrich_top10 %>% filter(class == "DGA")
          #    DGA_plot <- ggplot(DGA, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "DGA") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    DGA_plot <<- DGA_plot
          
          # get HPO
          #    HPO <- enrich_top10 %>% filter(class == "HPO")
          #    HPO_plot <- ggplot(HPO, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "HPO") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    HPO_plot <<- HPO_plot
          
          # get dbGAP
          #    dbGaP <- enrich_top10 %>% filter(class == "dbGaP")
          #    dbGaP_plot <- ggplot(dbGaP, aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
          #      geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #      #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          #      scale_colour_gradient(low = "#5656f7", high = "#63eaff", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #      #scale_x_discrete(limits=CC$Term) +
          #      labs(y = "geneRatio") +
          #      ggtitle(label = "dbGaP") + 
          #      coord_flip() +
          #      facet_wrap(vars(compare), nrow = 1)
          #    dbGaP_plot <<- dbGaP_plot
          
          #########
          
          #setwd("~/Desktop/swath_output") # replace Desktop with your destination path
          pdf("GOBP.pdf", width = 20, height = 6)
          print(BP_plot)
          dev.off()
          print("Get GOBP: GO-term Biological Processes")
          
          pdf("GOCC.pdf", width = 20, height = 6)
          print(CC_plot)
          dev.off()
          print("Get GOCC: GO-term Cellular Components")
          
          pdf("GOMF.pdf", width = 20, height = 6)
          print(MF_plot)
          dev.off()
          print("Get GOBP: GO-term Molecular Functions")
          
          pdf("KEGG.pdf", width = 20, height = 6)
          print(KEGG_plot)
          dev.off()
          print("Get KEGG: Pathways based on Kyoto Encyclopedia of Genes and Genomes database (https://www.genome.jp/kegg)")
          
          pdf("Reactome.pdf", width = 20, height = 6)
          print(Reactome_plot)
          dev.off()
          print("Get Reactome: a curated and peer-reviewed pathway database (https://reactome.org)")
          
          pdf("WikiPathway.pdf", width = 20, height = 6)
          print(WikiPathway_plot)
          dev.off()
          print("Get WikiPathway: a database of biological pathways maintained by and for the scientific community (https://www.wikipathways.org/index.php/WikiPathways)")
          
          ##### for additional protein enrichment plot: print #######        
          #    pdf("ENCODE_TF.pdf", width = 20, height = 6)
          #    print(ENCODE_TF_plot)
          #    dev.off()
          #    print("Get ENCODE_TF")
          
          #    pdf("ChEA.pdf", width = 20, height = 6)
          #    print(ChEA_plot)
          #    dev.off()
          #    print("Get ChEA")
          
          #    pdf("KEA.pdf", width = 20, height = 6)
          #    print(KEA_plot)
          #    dev.off()
          #    print("Get KEA")
          
          #    pdf("TFPI.pdf", width = 20, height = 6)
          #    print(TFPI_plot)
          #    dev.off()
          #    print("Get TFPI")
          
          #    pdf("CORUM.pdf", width = 20, height = 6)
          #    print(CORUM_plot)
          #    dev.off()
          #    print("Get CORUM")
          
          #    pdf("PPI.pdf", width = 20, height = 6)
          #    print(PPI_plot)
          #    dev.off()
          #    print("Get PPI")
          
          #    pdf("dbGaP.pdf", width = 20, height = 6)
          #    print(dbGaP_plot)
          #    dev.off()
          #    print("Get dbGaP")
          
          #    pdf("DGA.pdf", width = 20, height = 6)
          #    print(DGA_plot)
          #    dev.off()
          #    print("Get DGA")
          
          #    pdf("HPO.pdf", width = 20, height = 6)
          #    print(HPO_plot)
          #    dev.off()
          #    print("Get HPO")
          
          ######
          
          #setwd("~/Desktop/")
        #  setwd(wd)
          print('Finish smath function - please check your results in the swath_output folder')
          
          incProgress(1, message = "Complete enrichPlot")
          Sys.sleep(0.25)
        })
      } 
    })
    # })      
    } # debug1
  })  
  ##############
  
  observeEvent(input$normPlot, {
    output$plot16 <- renderPlot(box_all[['none']])  
    output$plot17 <- renderPlot(box_all[['TAS']])
    output$plot18 <- renderPlot(box_all[['Quantile']]) 
    output$plot19 <- renderPlot(box_all[['MLR']]) 
    output$plot20 <- renderPlot(box_all[['VSN']])
    output$plot21 <- renderPlot(density_all[['none']])  
    output$plot22 <- renderPlot(density_all[['TAS']])
    output$plot23 <- renderPlot(density_all[['Quantile']]) 
    output$plot24 <- renderPlot(density_all[['MLR']]) 
    output$plot25 <- renderPlot(density_all[['VSN']]) 
  })
  
  observeEvent(input$variancePlot, {
    output$plot26 <- renderPlot(
                      graphics::boxplot(pmad, main = 'PMAD', outline = FALSE, las = 2,
                                        border = 'darkred', boxwex = 0.6, horizontal = TRUE,
                                        col = colorRampPalette(c('white', 'red'))(5),
                                        notch = TRUE, frame.plot = FALSE, cex.axis = 0.9))  
    output$plotPCV <- renderPlot(
                      graphics::boxplot(pcv, main = 'PCV', outline = FALSE, las = 2,
                                        border = 'darkblue', boxwex = 0.6, horizontal = TRUE,
                                        col = colorRampPalette(c('white', 'blue'))(5),
                                        notch = TRUE, frame.plot = FALSE, cex.axis = 0.9))
    output$plotPEV <- renderPlot(
                      graphics::boxplot(pev, main = 'PEV', outline = FALSE, las = 2,
                                        border = 'darkgreen', boxwex = 0.6, horizontal = TRUE,
                                        col = colorRampPalette(c('white', 'lightgreen'))(5),
                                        notch = TRUE, frame.plot = FALSE, cex.axis = 0.9))
                    
  })
  
  observeEvent(input$qcPlot, {
    output$plot1 <- renderPlot(nPP_plot)  
    output$plot2 <- renderPlot(plot.qc)
    output$plot3 <- renderPlot(qc_hm_plot) 
    output$plot4 <- renderPlot(plotPCA) 
    output$plot5 <- renderPlot(plot_corrHM) 
  })
  
  observeEvent(input$volcanoPlot, {
    output$plot7 <- renderPlot(hm_sig)  
    output$plot6 <- renderPlot(volcano_all)  
  })
  
  observeEvent(input$goPlot, {
    
    if(input$ep_gobp){
      show("plot8")
      output$plot8 <- renderPlot(BP_plot)
    } else {
      hide("plot8")
    }
    
    
    
    if(input$ep_gocc){
    output$plot9 <- renderPlot(CC_plot)
    } else {}
    
    if(input$ep_gomf){
    output$plot10 <- renderPlot(MF_plot)
    } else {}
    
#  })
  
#  observeEvent(input$pathwayPlot, {
    
    if(input$ep_kegg){
    output$plot11 <- renderPlot(KEGG_plot)
    } else {}
    
    if(input$ep_reactome){
    output$plot12 <- renderPlot(Reactome_plot)
    } else {}
    
    if(input$ep_wikipathways){
    output$plot13 <- renderPlot(WikiPathway_plot)
    } else {}
  })
  
  observeEvent(input$masterTable, {
    output$contents3 <- DT::renderDataTable(enrich_master)
  })
  
  
  ##### Re_plot: Volcano ########
  
  # input$pair_select has already been updated after input$go session
  
  observeEvent(input$re_volcanoPlot, {
  #  source("lib/swathD_volcano.R", local = TRUE)
    
    observeEvent(input$pair_select, {
      pair_select <- input$pair_select
      pair_select <<- pair_select
      long_p.fc <- long_ano.fc %>% filter(compare == input$pair_select)
      #long_p.fc <<- long_p.fc
      assign("long_p.fc", long_p.fc, envir = .GlobalEnv)
    })
    
    observeEvent(input$re_cutFC, {
      re_cutFC <- input$re_cutFC
      re_cutFC <<- re_cutFC
    })
    
    observeEvent(input$re_cutP, {
      re_cutP <- as.numeric(input$re_cutP)
      re_cutP <<- re_cutP
    })
    
    observeEvent(input$cutFC_label, {
      cutFC_label <- input$cutFC_label
      cutFC_label <<- cutFC_label
    })
    
    observeEvent(input$cutP_label, {
      cutP_label <- as.numeric(input$cutP_label)
      cutP_label <<- cutP_label
    })
    
    
    #   output$plot14 <- renderPlot(swathD_volcano(data = long_p.fc, 
    #                                             re_cutFC = input$re_cutFC, re_cutP = input$re_cutP, 
    #                                            cutFC_label = input$cutFC_label, cutP_label = input$cutP_label, 
    #                                           range_x = input$xAxis,
    #                                          range_y = input$yAxis
    #                                         ))
    
    output$plot14 <- renderPlot(swathD_volcano(data = long_p.fc, 
                                               re_cutFC = re_cutFC, re_cutP = re_cutP, 
                                               cutFC_label = cutFC_label, cutP_label = cutP_label, 
                                               range_x = input$xAxis,
                                               range_y = input$yAxis
    ))
    
    
  })
  
  
  # export re_volcano in PDF file
  observeEvent(input$export_reVolcano, {
    re_volcano <- swathD_volcano(data = long_p.fc, 
                                 re_cutFC = re_cutFC, re_cutP = re_cutP, 
                                 cutFC_label = cutFC_label, cutP_label = cutP_label, 
                                 range_x = input$xAxis,
                                 range_y = input$yAxis)
    assign("re_volcano", re_volcano, envir = .GlobalEnv) 
    re_volcano_name <- paste0(c("re_volcano_", pair_select, "_fc", re_cutFC, "x_p", re_cutP, ".pdf"), collapse = "")
    #pdf(paste0(c("re_volcano_", pair_select, "_fc", re_cutFC, "x_p", re_cutP, ".pdf")), width = 9, height = 9)
    pdf("re_volcano.pdf", width = 12, height = 8)
    print(re_volcano)
    dev.off()
  })
  
  # generate lists of significant proteins for re_enrichment analysis
  observeEvent(input$select_protVolcano, {
    
    select_sigProt_up <- long_p.fc %>% filter(anova.pVal < re_cutP & adj_pVal < re_cutP & log2FC > log2(re_cutFC)) %>% dplyr::select(gene) 
    select_sigProt_up <- as.character(select_sigProt_up$gene)
    # select_sigProt_up <<- select_sigProt_up
    assign("select_sigProt_up", select_sigProt_up, envir = .GlobalEnv)
    
    select_sigProt_down <- long_p.fc %>% filter(anova.pVal < re_cutP & adj_pVal < re_cutP & log2FC < log2(1/re_cutFC)) %>% dplyr::select(gene) 
    select_sigProt_down <- as.character(select_sigProt_down$gene)
    # select_sigProt_down <<- select_sigProt_down
    assign("select_sigProt_down", select_sigProt_down, envir = .GlobalEnv)
    
    
    output$contents4 <- renderUI({
      #n_up <- paste(c(length(select_sigProt_up), "upregulated proteins:"))
      list_up <- paste(select_sigProt_up, collapse = ", ")  
      list_down <- paste(select_sigProt_down, collapse = ", ")  
      HTML(paste(c("<p> <b> <font color =", '"red"', ">", length(select_sigProt_up), " </font> upregulated proteins: </b> <br>", list_up, "</p> <br> ",
                   "<p> <b> <font color =", '"red"', ">", length(select_sigProt_down), " </font> downregulated proteins: </b> <br>", list_down, " </p> <br>")
      ))
    })
    
  })  
  
  observeEvent(input$submit_ReEnrich, {  
    
    updateSelectInput(session, "pair_select2", choices = pair_select)
    select_sigProt_all <- c(select_sigProt_up, select_sigProt_down)
    assign("select_sigProt_all", select_sigProt_all, envir = .GlobalEnv)
    updateSelectInput(session, "reg_select", choices = c("all", "upregulation", "downregulation"))
  })
  
  observeEvent(input$re_enrichAnal, {
    withProgress(message = "Submit data for gene enrichment analysis", value = 0.1, {
      Sys.sleep(0.25)
      
      if (input$reg_select == "all"){
        geneInput <- select_sigProt_all
      } else if (input$reg_select == "upregulation"){
        geneInput <- select_sigProt_up
      } else if (input$reg_select == "downregulation"){
        geneInput <- select_sigProt_down
      }
      
      dbs_chosen <- input$dbs_select
      
      #top_chosen <- as.numeric(input$top_n)
      
      incProgress(0.2, message = "Start Enrichr")
      
      re_enrichr <- enrichr(genes = geneInput, databases = dbs_chosen)
      assign("re_enrichr", re_enrichr, envir = .GlobalEnv)
      
      incProgress(0.2, message = "Prepare EnrichPlot")
      
      re_enrich_table <- separate(as.data.frame(re_enrichr[[1]]), Overlap, c("geneInput", "geneTotal"), sep = "/")
      re_enrich_table$geneInput <- as.numeric(re_enrich_table$geneInput)
      re_enrich_table$geneTotal <- as.numeric(re_enrich_table$geneTotal)
      re_enrich_table <- re_enrich_table %>% mutate("geneRatio" = (geneInput / geneTotal)) %>% dplyr::arrange(Adjusted.P.value)
      #re_enrich_table <- re_enrich_table %>% mutate("rankScore" = (100*geneRatio) * (-log10(Adjusted.P.value)) )
      
      assign("re_enrich_table", re_enrich_table, envir = .GlobalEnv)
      
      output$plot15 <- renderPlot(
        ggplot(head(re_enrich_table, n = as.numeric(input$top_n)), aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 15), 45), geneRatio), geneRatio), na.rm = TRUE) +
            geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
        #  labs(y = "[rankScore = %geneRatio * -log10(adj.P.val)]") +
            labs(y = "geneRatio") +
          ggtitle(label = dbs_chosen) + 
          xlab("")+ 
          theme(legend.title = element_text(size=15),
                legend.text = element_text(size=12)) +
          coord_flip() 
      )
      
      #   assign("re_enrichPlot", re_enrichPlot, envir = .GlobalEnv)
      
      #   output$plot15 <- renderPlotly(
      #                        ggplotly(width = 1000, height = 500, 
      #                                ggplot(head(re_enrich_table, n = as.numeric(input$top_n)), aes(reorder(Term, geneRatio), geneRatio), na.rm = TRUE) +
      #                                 geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
      #                                #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
      #                               scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
      #                              labs(y = "geneRatio") +
      #                             ggtitle(label = dbs_chosen) + 
      #                            coord_flip()
      #                ))
      
      
      #    output$contents5 <- DT::renderDataTable({
      #     re_enrich_table
      #  }, options = list(pageLength = 5))
      
    })
  })
  
  observeEvent(input$re_hm_var, {
  #  source("lib/swathD_hm_var.R", local = TRUE)
    
    hm_set <<- input$hm_set  
    varTopN <<- input$varTopN
    n_cluster <<- input$n_cluster
    
    if(hm_set == 'Significant features'){
      dfvar <- dfmat_medScale_sig
      dfvar <<- dfvar
    } else if (hm_set == 'All features'){
      dfvar <- dfmat_medScale_all
      dfvar <<- dfvar
    }
    
    swathD_hm_var(dfvar)
    
    output$plot30 <- renderPlot(hm_var)
    
    updateSelectInput(session, "cluster_select", choices = (1:n_cluster) )
    
  })
  
  observeEvent(input$cluster_enrichAnal, {
    withProgress(message = "Submit gene cluster for gene enrichment analysis", value = 0.1, {
      Sys.sleep(0.25)
    
      dbs_chosen2 <- input$dbs_select2
      dbs_chosen2 <<- dbs_chosen2
      
      incProgress(0.2, message = "Start Enrichr")
      
      #selectedClust <<- as.numeric(input$cluster_select)
      
      geneClust <- names(gene_cluster[gene_cluster == as.numeric(input$cluster_select)])
      geneClust <<- geneClust
      
      re_enrichr2 <- enrichr(genes = geneClust, databases = dbs_chosen2)
      assign("re_enrichr2", re_enrichr2, envir = .GlobalEnv)
      
      incProgress(0.2, message = "Prepare EnrichPlot")
      
      re_enrich2_table <- separate(as.data.frame(re_enrichr2[[1]]), Overlap, c("geneInput", "geneTotal"), sep = "/")
      re_enrich2_table$geneInput <- as.numeric(re_enrich2_table$geneInput)
      re_enrich2_table$geneTotal <- as.numeric(re_enrich2_table$geneTotal)
      re_enrich2_table <- re_enrich2_table %>% mutate("geneRatio" = (geneInput / geneTotal)) %>% dplyr::arrange(Adjusted.P.value)
      
      assign("re_enrich2_table", re_enrich2_table, envir = .GlobalEnv)
      
      output$plot31 <- renderPlot(
        ggplot(head(re_enrich2_table, n = as.numeric(input$top_n2)), aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 15), 45), geneRatio), geneRatio), na.rm = TRUE) +
          geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
          #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
          scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
          #labs(y = "Enrichr.Combined.Score") +
          labs(y = "geneRatio") +
          ggtitle(label = dbs_chosen2) + 
          xlab("")+
          theme(legend.title = element_text(size=15),
                legend.text = element_text(size=12)) +
          coord_flip() 
      )
      
      
    })
  
  })
  
  observeEvent(input$ml_var, {
  #  source("lib/swathD_ml2enrichr.R", local = TRUE)
    withProgress(message = "Start machine learning", value = 0.1, {
      Sys.sleep(0.25)
    
    chose_ml <<- input$ml
    
    if(chose_ml == 'nmf'){
      incProgress(0.3, message = paste("Run NMF with 100 iterations ... slow"))
    } else if(chose_ml == 'pca'){
      incProgress(0.3, message = paste("Run PCA"))
    } else if(chose_ml == 'sipca'){
      incProgress(0.3, message = paste("Run sIPCA"))
    }
    
    swathD_ml2enrichr()  
    incProgress(0.3, message = "Prepare data for EnrichPlot")
    
    updateSelectInput(session, "ml2", choices = chose_ml)
    incProgress(0.3, message = "Finish")
    
    })
  })
  
  observeEvent(input$ml2hm, {
    chose_com_ml2hm <<- as.numeric(input$choose_component_ml2hm)
    
    chose_com_type <<- input$chose_com_type
    
    if(input$chose_com_type == "Significant features"){
      ml_sig_com <- sig_coms[[chose_com_ml2hm]]
      chose_ml_com <<- ml_sig_com 
    } else if(input$chose_com_type == "All features"){
      ml_mat_com <- mat_coms[[chose_com_ml2hm]]
      chose_ml_com <<- ml_mat_com 
    }
    if(nrow(chose_ml_com) == 0){
      output$warning_ml2hm <- renderText("no information for the selected component and data type")
      output$plot33 <- renderPlot(NULL)
    } else {
    hm_chose_ml_com <- pheatmap(chose_ml_com, breaks = seq(-(max(round(chose_ml_com, 0))), max(round(chose_ml_com, 0)), length.out=101),
                                legend_breaks=seq(-(max(round(chose_ml_com, 0))), max(round(chose_ml_com, 0)), length.out=5),
                                color = colorRampPalette(c("darkblue", "blue", "white", "orangered", "red"))(100), 
                                border_color = NA,
                                annotation_col = data.frame(group = factor(group), row.names = sample_label),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_cols = "euclidean", 
                                clustering_method = "average",
                                fontsize_row = 4, fontsize_col = 10,
                                scale = "none",
                                main = paste0("ML2HM plot (com#",chose_com_ml2hm, "; ", nrow(chose_ml_com), " ", tolower(chose_com_type), ")"))
    
    output$warning_ml2hm <- renderText(NULL)
    output$plot33 <- renderPlot(hm_chose_ml_com)
    }
    
  })
  
  
  observeEvent(input$ml_enrichAnal, {
    withProgress(message = "Submit genes in the chosen component for gene enrichment analysis", value = 0.2, {
    
    com_set <- input$com_set      
    if(com_set == 'Significant features'){
      com_input <- com_sig_genes
      com_input <<- com_input
    } else if (com_set == 'All features'){
      com_input <- com_genes
      com_input <<- com_input
    }
    
    chose_com <<- as.numeric(input$choose_component)
    dbs_chosen3 <<- input$dbs_select3
    
    incProgress(0.4, message = "Start Enrichr")
    
    re_enrichr3 <- enrichr(genes = com_input[[chose_com]], databases = dbs_chosen3)
    
    assign("re_enrichr3", re_enrichr3, envir = .GlobalEnv)
    
    incProgress(0.2, message = "Prepare EnrichPlot")
    
    re_enrich3_table <- separate(as.data.frame(re_enrichr3[[1]]), Overlap, c("geneInput", "geneTotal"), sep = "/")
    re_enrich3_table$geneInput <- as.numeric(re_enrich3_table$geneInput)
    re_enrich3_table$geneTotal <- as.numeric(re_enrich3_table$geneTotal)
    re_enrich3_table <- re_enrich3_table %>% mutate("geneRatio" = (geneInput / geneTotal)) %>% dplyr::arrange(Adjusted.P.value)
    
    assign("re_enrich3_table", re_enrich3_table, envir = .GlobalEnv)
    
    
    output$plot32 <- renderPlot(
      ggplot(head(re_enrich3_table, n = as.numeric(input$top_n3)), aes(reorder(stringr::str_trunc(stringr::str_wrap(Term, 15), 45), geneRatio), geneRatio), na.rm = TRUE) +
        geom_point(aes(size = geneInput, colour = Adjusted.P.value)) + #(stat="identity") + 
        #scale_colour_gradient(low = "red", high = "black", name="Adjusted p-value")+ 
        scale_colour_gradient(low = "red", high = "yellow", breaks = c(0.05, 5e-8), limits = c(0, 0.05)) +
        #labs(y = "Enrichr.Combined.Score") +
        labs(y = "geneRatio") +
        ggtitle(label = dbs_chosen3) + 
        xlab("")+
        theme(legend.title = element_text(size=15),
              legend.text = element_text(size=12)) +
        coord_flip() 
      )
    incProgress(0.2, message = "Finish")
    
    })
  })
  
  observeEvent(input$ml_enrichNet, {
    withProgress(message = "Submit genes in the chosen component for gene enrichment analysis", value = 0.2, {
      
  #  source("lib/swathD_ml2net.R", local = TRUE)
    chose_ml2net_com <<- as.numeric(input$chose_ml2net_com)
    
    if( !exists("com_input") ){
      print("bypass enrich_plot; passing chosen features from heatmap to network")
      if(input$chose_com_type == "Significant features"){
        com_input <- com_sig_genes
        com_input <<- com_input
        ml2net <- com_input[[chose_ml2net_com]]
      } else if(input$chose_com_type == "All features"){
        com_input <- com_genes
        com_input <<- com_input 
        ml2net <- com_input[[chose_ml2net_com]]
      }
    } else {
    ml2net <- com_input[[chose_ml2net_com]]
    }
    
    ml2net <<- ml2net
    ml2net_dbs <- input$ml2net_dbs
    ml2net_dbs <<- ml2net_dbs
    print(paste("select", ml2net_dbs))
    apval <- as.numeric(input$apval)
    rscore <- as.numeric(input$rscore)
    rscore <<- rscore
    nlink <- as.numeric(input$nlink)
    nlink <<- nlink
    
    
    swathD_ml2net(dat = ml2net, dbs = ml2net_dbs, cutPval = apval, cutScore = rscore, cutLinkage = nlink)
    
    if(exists("warning_rankScore") | exists("warning_n_linkage")){
    output$warning_ml2net <- renderText("No data. Please try lower the rankScore or the numbers of linkage")
    #ourput$warning_n_linkage <- renderText(warning_n_linkage)
    output$ml2net <- renderVisNetwork(NULL)
    } else {
    output$warning_rankScore <- renderText(NULL)
    #output$warning_n_linkage <- renderText(NULL)
    output$ml2net <- renderVisNetwork(ml2net_vis %>% 
                                        visGroups(groupname = "Gene involved", shape = "dot", #size = 1,
                                                               font =  list(color = "magenta", size = 20),
                                                               color = list(shadow = TRUE, background = "ddd", border = "black",
                                                                            highlight = list(background = 'yellow', #"feff8f", 
                                                                                             border = "red")) ) %>%
                                        visGroups(groupname = "Function", shape = "dot", #size = 25,
                                                  font =  list(color = "blue", size = 20),
                                                  color = list(background = "yellow", border = "orange",
                                                               highlight = list(background = "feff8f", 
                                                                                border = "red")) ) %>%
                                        visHierarchicalLayout(enabled = input$hierarchical, levelSeparation = 150, treeSpacing = 100, nodeSpacing = 100, blockShifting = FALSE, edgeMinimization = FALSE) %>%
                                        visIgraphLayout(type = "full", randomSeed = 9, layout = as.character(input$layout), physics = FALSE, smooth = TRUE) %>%
                                        #visIgraphLayout(type = "full", randomSeed = 9, layout = "layout_nicely", physics = FALSE, smooth = TRUE) %>%
                                        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE, selectedBy = "group") %>%
                                        visInteraction(navigationButtons = TRUE)
                             
            )
          }
    
    })
  })
}


shinyApp(ui = ui, server = server)


