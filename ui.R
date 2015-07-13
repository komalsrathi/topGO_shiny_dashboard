library(shinydashboard)
library(shinyIncubator)
library(shiny)

dashboardPage(
  dashboardHeader(title = 'TopGO Web Interface',titleWidth = 500),
  dashboardSidebar(width = 500,
    div(style="overflow-y: scroll"),
#     h3('-------Input Data-------'),
#     fluidRow(
#       column(6, fileInput(inputId = 'FileInput', label = 'Upload LIMMA Output:', accept = c('csv','tsv','txt'))),
#       column(2, checkboxInput(inputId = 'header', label = 'Header', value = FALSE)),
#       column(2, div(style = "height:12px"), radioButtons(inputId = 'sep', label = 'Separator', choices = c(comma=',',tab="\t",space=' '), selected = ","),offset = 1)
#     ),
#     fluidRow(column(6, fileInput(inputId = 'FileInput1', label = 'Upload ExpressionSet:'))),
    #br(),
    fluidRow(
      column(6, div(style = "height:10px"), fileInput(inputId = 'inputdataset', label = 'Input Datasets', accept = c('csv','tsv','txt')))
    ),
    fluidRow(
      column(6, div(style = "height:10px"), selectInput(inputId='projects', label = 'Projects',choices = c('none'), selected = 'none'))
    ),
    h3('-------Select Foreground-------'),
    fluidRow(
      column(4, radioButtons(inputId = 'cutoff', label = 'Selection', choices = c('Upregulated'='pos','Downregulated'='neg','Both'='both'))),
      br(),
      column(4, textInput(inputId = 'fc', label = "Fold Change", value = '0')),
      column(4, textInput(inputId = 'pvalue', label = "Adj. Pvalue",value = '0.05'))
    ),
    fluidRow(column(3, actionButton(inputId = 'select', label = "Select Data"))),
    tags$style(type='text/css', "#select { width:100%; margin-left: 10px;}"),
    br(),
    h3('-------TopGO Function Input-------'),
    fluidRow(
      column(4,checkboxGroupInput(inputId = 'cbxgrp1', label = 'Select Tests', 
                                                              choices = c("Classic Fisher" = 1,
                                                                          "Classic KS" = 2,
                                                                          "Elim Fisher" = 3,
                                                                          "Elim KS" = 4,
                                                                          "Weight01 Fisher" = 5,
                                                                          "Weight01 KS" = 6), selected = NULL, inline = FALSE)),
      column(3, selectInput('ontology','Ontology',c('BP','MF','CC'))),
      column(5, selectInput(inputId = 'genome', label='Chip',choices=c('Mouse Gene 1.0'='mogene10sttranscriptcluster.db',
                                                                                                   'Mouse Gene 2.0'='mogene20sttranscriptcluster.db',
                                                                                                   'Human Gene 1.0'='hugene10sttranscriptcluster.db',
                                                                                                   'Human Gene 2.0'='hugene20sttranscriptcluster.db',
                                                                                                   'Human Genome U133'='hgu133a.db',
                                                                                                   'Human Genome U133 Plus 2.0'='hgu133plus2.db',
                                                                                                   'Human Genome U133A 2.0'='hgu133a2.db'), selected = 'mogene20sttranscriptcluster.db'))
    ),
    fluidRow(column(3, actionButton(inputId = 'topgo', label = "Run topGO"))),
    tags$style(type='text/css', "#topgo { width:100%; margin-left: 10px;}"),
    br(),
    h3('---Heatmap of Enriched Genes---'),
    fluidRow(
      column(4, actionButton(inputId = 'makeheat', label = 'Create Heatmap')),
      column(5, downloadButton('downloadheatmap', 'Download Heatmap'))
    ),
    tags$style(type='text/css', "#downloadheatmap { color:#000;}"),
    tags$style(type='text/css', "#makeheat { width:100%; margin-left: 10px;}"),
    br(),
    h3('------TopGO Enrichment Plot------'),
    fluidRow(
      column(5, selectInput('selectresult','Algorithms:',c('None'),selected = 'none')),
      column(3, div(style = "height:40px"), actionButton(inputId = 'makeplot', label = 'Create Plot')),
      column(2, div(style = "height:40px"), downloadButton('downloadplot','Download Plot'))
    ),
    tags$style(type='text/css', "#downloadplot { color:#000;}"),
    h3('-----Pathway Enrichment-----'),
    fluidRow(
      column(3, actionButton(inputId = 'runspia', label = 'Run SPIA'))
      ),
    tags$style(type='text/css', "#runspia { width:100%; margin-left: 10px; margin-bottom:40px;}")
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    tabsetPanel(type="tabs", id = "tabvalue",
                tabPanel(title = "Input Table", value = 'tab1', DT::dataTableOutput('table')),
                tabPanel(title = "Expression Plot", value = 'tab7', plotOutput('dotplot',width = 800,height = 800)),
                tabPanel(title = "PCA Plot", value = 'tab8', 
                         fluidRow(
                           column(5,div(style = "height:10px"),
                                  sliderInput("obs", "Number of observations:", min = 20, max = 2000, value = 500, width = 800, step = 10)),
                           column(3,offset = 1,div(style = "height:20px"),selectInput(inputId = 'selectlabel',label='Label',choices=c('Group','Sample'))),
                           column(3,div(style = "height:20px"),selectInput(inputId = 'selectcolor',label='Color',choices=c('Group','Sample')))
                         ),
                         plotOutput("distPlot", width = 800, height = 800)),
                tabPanel(title = "Filtered Table (Foreground)", value = 'tab2', DT::dataTableOutput('table2')),
                tabPanel(title = "TopGO Enrichment Results", value = 'tab3', DT::dataTableOutput('results')),
                tabPanel(title = "GO:Genes Annotation", value = 'tab5', DT::dataTableOutput('results2')),
                tabPanel(title = 'Heatmap of Enriched Genes', value = 'tab6', plotOutput('heatmap',width = 800,height = 1500)),
                tabPanel(title = "TopGO Enrichment Plot", value = 'tab4', plotOutput('plot')),
                tabPanel(title = "SPIA Pathway Enrichment", value = 'tab9', DT::dataTableOutput('spiaout')),
                tabPanel(title = "Help Page", value='tab10', 
                         h4("1. Upload the output of",strong("LIMMA (performs differentially expressed genes)."),br(),br(),
                           "2. Upload the",strong("Normalized Expression Set."),br(),br(),
                           "3. Click on any gene to view its expression in the tab called",strong("Expression Plot."),br(),br(),
                           "4. Click on the tab called",strong("PCA Plot"),"to view a PCA plot. You can change the number of top variant genes to view, color and labels.",br(),br(),
                           "5. When you are ready to select the foreground, go to",strong("Select Foreground"),"section, choose the Fold change & Pvalue cutoff and hit Select Data.",br(),br(),
                           "6. You can view your filtered table in the",strong("Filtered Table (Foreground)")," tab.",br(),br(),
                           "7. When you have decided which tests to use, go to",strong("TopGO Function Input"),"section, and select the proper tests, the Ontology & your chip.",br(),br(),
                           "8. Hit",strong("Run topGO."),br(),br(),
                           "9. topGO results will appear under the tab called",strong("TopGO Enrichment Results."),br(),br(),
                           "10. You can select the GO Terms that you are interested in. The genes belonging to those GO Terms will appear under the tab",strong("GO:Genes Annotation"),br(),br(),
                           "11. When you have selected all the GO Terms that you are interested in, you can go to the section",strong("Heatmap of Enriched Genes")," and hit",strong("Create Heatmap"),".", 
                           "This will create a heatmap under",strong("Heatmap of Enriched Genes"),".",br(),br(),
                           "12. You can create a plot of the top 5 GO Terms by selecting the appropriate test & clicking",strong("Create Plot"), "under", strong("TopGO Enrichment Plot"),".",br(),br(),
                           "13. Finally, you can do a pathway analysis of your foreground genes by clicking on",strong("Run SPIA")," under the",strong("Pathway Enrichment"), "tab.",
                           style = "font-si18pt; color:#000"
                           ))
  )
))
