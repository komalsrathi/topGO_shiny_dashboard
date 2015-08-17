library(shiny)
library(reshape2)
library(ggplot2)
library(topGO)
library(mogene20sttranscriptcluster.db)
library(mogene10sttranscriptcluster.db)
library(hugene10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)
library(hugene20sttranscriptcluster.db)
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133a2.db)
library(Rgraphviz)
library(DT)
library(plyr)
library(GMD)
library(RColorBrewer)
library(pheatmap)
library(gplots)
library(NMF)
library(SPIA)
library(shinydashboard)
library(shinyIncubator)

options(shiny.maxRequestSize = 30*1024^2)

shinyServer(function(input, output, session){
  
  ###################################### global variables #####################################
  # plot height & width
  varH <- 20
  varW <- 12
  
  # pData
  pData <- data.frame()
  
  # theme object 
  tt <- theme_bw() + theme(axis.title.x = element_text(face = "bold", size = 18),
                           axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 16),
                           axis.title.y = element_text(face = "bold", size = 18),
                           axis.text.y  = element_text(angle = 0, vjust = 0.5, size = 16),
                           title = element_text(face = 'bold',size = 18))
  ###################################### global variables #####################################
  
  # read in initial data input.csv i.e. a file containing list of projects
  datainput <- reactive({
    infile <- input$inputdataset
    if(is.null(infile))
      return(NULL)
    read.csv(infile$datapath,check.names=F)
  })
  
  # output list of projects in projects tab
  # you can select any one of the projects
  output$projects <- DT::renderDataTable({
    DT::datatable(datainput(), selection = 'single')
  })
  
  # go to tab1 when a project is selected 
  observe({
    s <- input$projects_rows_selected
    if(length(s)){
      updateTabsetPanel(session = session, inputId = "tabvalue", selected = 'tab1')
    }
  })
  
  # when a project is selected, the corresponding limma file is opened
  datasetInput <- reactive({
    d <- datainput()
    s <- input$projects_rows_selected
    dat <- d[s, , drop = FALSE]
    limma <- as.character(dat$LIMMA)
    file <- paste('data/',limma,sep='')
    dd <- read.csv(file=file)
    dd
  })
  
  # output data for selected project in tab1
  output$table <- DT::renderDataTable({
    DT::datatable(datasetInput(),
                  extensions = c('ColVis','Scroller'), 
                  options = list(dom = 'RMDCT<"clear">lfrtip', 
                                 searchHighlight = TRUE,
                                 initComplete = JS("function(settings, json) {",
                                                   "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                                                   "}"),
                                 pageLength = 5,
                                 lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                                 scrollX = TRUE),
                  selection = 'single')
  })
  
  # subset data based on fold change & up or downregulated genes - foreground data
  datasetInput2 <- reactive({
    validate(
      need(input$fc >= 0, "Input a Fold Change Value >= 0")
    )
    d <- datasetInput()
    if(input$cutoff=='pos')
    {
      d <- d[which(d$fc > as.numeric(input$fc) & d$adj.P.Val < as.numeric(input$pvalue)),]
    }
    else if(input$cutoff=='neg')
    {
      fc = as.numeric(input$fc)*(-1)
      d <- d[which(d$fc < fc & d$adj.P.Val < as.numeric(input$pvalue)),]
    }
    else if(input$cutoff=='both')
    {
      d <- d[which(abs(d$fc) > as.numeric(input$fc) & d$adj.P.Val < as.numeric(input$pvalue)),]
    }
  })
  
  # print out filtered data in tab2
  output$table2 <- DT::renderDataTable({
    input$select
    isolate({
      DT::datatable(datasetInput2(),
                    extensions = c('TableTools','ColVis','Scroller'),
                    options = list(
                      dom = 'RMDCT<"clear">lfrtip',
                      searchHighlight = TRUE,
                      tableTools = list(sSwfPath = '//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls_pdf.swf'),
                      pageLength = 5,
                      lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                      initComplete = JS("function(settings, json) {",
                                        "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                                        "}"),
                      scrollX = TRUE
                    ))
    })
  })
  
  # switch to tab2 on clicking the select button
  observe({
    if(input$select>0){
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab2') 
    }
  })
  
  # get GOdata object 
  datasetInput4 <- reactive({
    limma <- datasetInput()
    fg <- datasetInput2()
    
    validate(
      need(fg$fc, "Check Input...!!")
    )
    
    if (input$cutoff=='pos') {
      limma <- limma[which(limma$fc>0),]
    } else if (input$cutoff=='neg') {
      limma <- limma[which(limma$fc<0),]
    } else if(input$cutoff=='both') {
      limma <- limma      
    } 
    
    all_genes <- limma$adj.P.Val
    names(all_genes) <- limma$ID
    interesting_genes <- as.character(fg$ID)
    geneList <- as.logical(names(all_genes) %in% interesting_genes)
    names(geneList) <- names(all_genes)
    
    topDiffGenes <- function(allScore) {
      return(geneList)
    }
    
    genome.chip = sub('.db','',input$genome)
    genome.chip = paste(genome.chip,'GO2PROBE',sep='')
    genome.chip = get(genome.chip)
    
    GOdata <- new("topGOdata", 
                  description = "topGO session",
                  ontology = input$ontology, 
                  allGenes = all_genes, 
                  geneSel = topDiffGenes,
                  annot = annFUN.GO2genes,
                  GO2genes = as.list(genome.chip),
                  nodeSize = 10)
  })
  
  # update max node size & run selected tests
  datasetInput3 <- reactive({
    validate(
      need(input$cbxgrp1, "Select one or more tests")
    )
    
    GOdata <- datasetInput4()
    
    test = as.numeric(input$cbxgrp1)
    
    algos = list(c('classic','fisher'),c('classic','ks'),c('elim','fisher'),c('elim','ks'),c('weight01','fisher'),c('weight01','ks'))
    
    tmp <- data.frame(rownames=c(1:length(usedGO(GOdata))))
    for(i in test)
    {
      res <- runTest(object = GOdata,algorithm = algos[[i]][[1]],statistic = algos[[i]][[2]])
      allRes <- GenTable(GOdata,res,topNodes = length(usedGO(object = GOdata)))
      names(allRes)[6] = paste(algos[[i]][[1]],algos[[i]][[2]],sep='.')
      allRes = allRes[order(allRes$GO.ID),]
      allRes = cbind(allRes, rank = as.numeric(rownames(allRes)))
      names(allRes)[7] = paste(algos[[i]][[1]],algos[[i]][[2]],'rank',sep='.')
      tmp <- cbind(tmp,allRes)
      tmp <- tmp[,unique(colnames(tmp))]
    }
    tmp[,-1]
  })
  
  # output topGO results in tab3
  output$results <- DT::renderDataTable({
    input$topgo
    withProgress(session = session, message = 'Calculating...',detail = 'This may take a while...',{
      isolate({
        DT::datatable(data = datasetInput3(),
                      extensions = c('TableTools','ColVis','Scroller'),
                      options = list(
                        dom = 'RMDCT<"clear">lfrtip',
                        searchHighlight = TRUE,
                        tableTools = list(sSwfPath = '//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls_pdf.swf'),
                        pageLength = 5,
                        lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                        initComplete = JS("function(settings, json) {",
                                          "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                                          "}"),
                        scrollX = TRUE
                      ))
      })
    })
  })
  
  # switch to tab3 for topGO results
  observe({
    if(input$topgo > 0){
      updateTabsetPanel(session = session, inputId = "tabvalue", selected = "tab3")
    }
  })
  
  # update dropdown with selected tests for graph generation
  observe({
    if(input$topgo > 0)
    {
      sel <- as.character(input$cbxgrp1)
      seltests = c('1'='classic fisher','2'='classic ks','3'='elim fisher','4'='elim ks','5'='weight01 fisher','6'='weight01 ks')
      seltests = seltests[names(seltests) %in% sel]
      updateSelectInput(session = session, inputId = "selectresult", choices = as.character(seltests))
    }
  })
  
  # make another data input - to get results of the selected test
  datasetInput5 <- reactive({
    test <- as.character(input$selectresult)
    test <- unlist(strsplit(test,' '))
    GOdata <- datasetInput4()
    res <- runTest(object = GOdata, algorithm = as.character(test[1]), statistic = as.character(test[2]))
  })
  
  # plot function to plot the topGO graph using the selected test
  plot_out <- function(){
    GOdata <- datasetInput4()
    res <- datasetInput5()
    showSigOfNodes(GOdata, score(res), firstSigNodes = 3, useInfo = "all")
  }
  
  # create topGO plot in tab4
  output$plot <- renderPlot({
    input$makeplot
    withProgress(session = session, message = 'Calculating...',detail = 'This may take a while...',{
      isolate({
        plot_out()
      })
    })
  })
  
  # update tab4 with topGO plot
  observe({
    if(input$makeplot > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab4')
    }
  })
  
  # download topGO plot
  output$downloadplot <- downloadHandler(filename = function(){ 
    paste('output', '.pdf', sep='')},content = function(file){
      pdf(file = file, width = 10,height = 10)
      plot_out()
      dev.off()
    })
  
  # get GO associated genes 
  datasetInput6 <- reactive({
    limma <- datasetInput() 
    d <- datasetInput3()
    s <- input$results_rows_selected
    goid <- d[s, , drop = FALSE]
    goid <- as.character(goid$GO.ID)
    GOdata <- datasetInput4()
    gogenes <- printGenes(GOdata, whichTerms = goid, geneCutOff = 500, chip = input$genome)
    if(length(goid)>1)
    {
      dd = ldply(.data = gogenes, data.frame)
    } else 
    {
      dd = as.data.frame(gogenes)
      vec <- c('Chip.ID','LL.id','Symbol.id','Gene.name','raw.p.value')
      validate(
        need(length(vec)==ncol(dd), "Check Input...!!")
      )
      colnames(dd) = vec
      rownames(dd) = NULL
      dd = cbind(.id=goid,dd)
    }
    dd <- merge(dd,d[,1:2],by.x='.id',by.y='GO.ID')
    
    dd <- dd[,c(1,2,ncol(dd))] 
    dd <- merge(limma, dd, by.x='ID', by.y='Chip.ID') 
    dd
  })
  
  # output genes in selected GO Terms
  output$results2 <- DT::renderDataTable({
    DT::datatable(datasetInput6(),
                  extensions = c('TableTools','ColVis','Scroller'),
                  options = list(
                    dom = 'RMDCT<"clear">lfrtip',
                    searchHighlight = TRUE,
                    tableTools = list(sSwfPath = '//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls_pdf.swf'),
                    pageLength = 5,
                    lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                    initComplete = JS(
                      "function(settings, json) {",
                      "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                      "}"),
                    scrollX = TRUE
                  ))
  })
  
  # get eset and convert to dataframe
  datasetInputHeatmap <- reactive({
    
    validate(
      need(input$inputdataset, "Input list of datasets")
    )
    
    d <- datainput()
    s <- input$projects_rows_selected
    dat <- d[s, , drop = FALSE]
    eset <- as.character(dat$ESET)
    file <- paste('data/',eset,sep='')
    dd = load(file)
    
    dd <- get(dd)
    pData <<- pData(dd)
    pData[] <<- lapply(pData, as.character)
    updateSelectInput(session = session, inputId = 'selectcolor', choices = names(pData))
    updateSelectInput(session = session, inputId = 'selectlabel', choices = names(pData))
    dd <- exprs(dd)
    dd <- as.data.frame(dd)
  })
  
  # make dotplot of gene
  dotplot_out <- reactive({
    dd <- datasetInputHeatmap()
    validate(
      need(dd,'Check Input..!!!')
    )
    s <- input$table_rows_selected
    dt <- datasetInput()
    dt <- dt[s, , drop=FALSE]
    probeid <- dt$ID
    probename <- dt$SYMBOL
    dd$id <- rownames(dd)
    dd.m <- melt(dd,id.vars = 'id')
    dd.m <- dd.m[grep(probeid,dd.m$id),]
    dd.m = merge(dd.m, pData, by.x='variable', by.y='Sample')
    
    ggplot(dd.m, aes(x=Group,y=value,color=Group)) + geom_point(size=5,position=position_jitter(w = 0.1)) + 
      geom_errorbar(stat = "hline", yintercept = "mean", width=0.8,aes(ymax=..y..,ymin=..y..)) + 
      ggtitle(paste('\nExpression Plot:',probename,'\n')) + tt + theme(legend.position='none') + ylab('Normalized Expression\n') + xlab('\nConditions')
    
  })
  
  # output dotplot
  output$dotplot <- renderPlot({
    dotplot_out()
  },width = 'auto',height = 'auto')
  
  # update tab7 with dotplot
  observe({
    s <- input$table_rows_selected
    if(length(s)){
      updateTabsetPanel(session = session, inputId = "tabvalue", selected = 'tab7')
    }
  })
  
  # make a heatmap of selected genes
  # plot function
  heatmap_out <- function(){
    dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
    col_ylgnbu <- colorRampPalette(rev(brewer.pal(n = 9, "YlGnBu")))
    dd <- datasetInputHeatmap() 
    dt <- datasetInput6() 
    dd <- merge(dd,dt[,c(1,2,13)],by.x='row.names',by.y='ID')
    dd <- dd[order(dd$Term),]
    dd <- unique(dd)
    nr <- nrow(dd)
    ids <- as.character(dd$Term)
    annot = data.frame(GO.Term = ids)
    syms = as.character(dd$SYMBOL)
    dd = dd[,-c(1,ncol(dd),ncol(dd)-1)]
    nc = ncol(dd)
    # varH <<- ifelse(nr>=100,varH,round(nr/5))
    # varH <<- round(nr/5)
    aheatmap(x = as.matrix(dd), scale = 'row', distfun = dist2, fontsize = 16,annRow = annot, labRow = syms, labCol = pData$Group,
             cexRow = min(0.2 + 1/log10(nr), 1.2), cexCol = min(0.2 + 1/log10(nc), 1.2), Rowv = NA, Colv=NA,
             color = "YlGnBu",
             main = 'Heatmap of Enriched Genes')
  }
  
  # make heatmap for genes associated with selected GO terms
  output$heatmap <- renderPlot({
    input$makeheat
    isolate({
      heatmap_out()
    })
  },width = 'auto',height = 'auto')
  
  # update tab6 with heatmap
  observe({
    if(input$makeheat > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab6')
    }
  })
  
  # download button for heatmap
  # download plot
  output$downloadheatmap <- downloadHandler(filename = function(){ 
    paste('heatmap', '.pdf', sep='')},content = function(file){
      pdf(file = file, onefile = FALSE)
      heatmap_out()
      dev.off()
    })
  
  # PCA plot
  # make pca plot
  pca_plot <- reactive({
    expr <- datasetInputHeatmap()
    vst.var = sort(apply(expr,1,var),decreasing = T)
    s <- as.numeric(input$obs)
    vst.var = expr[rownames(expr) %in% names(vst.var)[1:s],]
    pca = prcomp(t(vst.var))
    scores = data.frame(pData,pca$x)
    tmp = summary(pca)
    pc1 = round(tmp$importance[2,1]*100,2)
    pc2 = round(tmp$importance[2,2]*100,2)
    title = paste('\nPrincipal Component Analysis\nTop',nrow(vst.var),'Variant Genes\n')
    
    ggplot(data=scores,aes(x=PC1,y=PC2)) + 
      geom_text(aes_string(label=as.character(input$selectlabel),color=as.character(input$selectcolor)),cex=8) +
      ggtitle(paste(title,'\n')) + 
      xlab(paste("\nPC1",paste(pc1,"%",sep=""),sep=": ")) + ylab(paste("PC2",paste(pc2,"%\n",sep=""),sep=": ")) +
      tt + theme(legend.text=element_text(size=16),
                 legend.title=element_text(size=18))
  })
  
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$inputdataset, {
    v$doPlot <- input$inputdataset
  })
  
  output$distPlot <- renderPlot({
    if (v$doPlot[[1]] == FALSE) return()
    pca_plot()
  })
  
  # SPIA
  spia_out <- reactive({
    limma <- datasetInput()
    fg <- datasetInput2()
    
    validate(
      need(fg, "Check Input...!!")
    )
    
    genome.chip <- as.character(input$genome)
    if(length(grep('m',genome.chip))){
      org <- 'mmu'
    } else if(length(grep('h',genome.chip))){
      org <- 'hsa'
    }
    
    all_genes = as.numeric(limma$ENTREZID)
    sig_genes = fg$logFC
    names(sig_genes) = fg$ENTREZID 
    sig_genes = sig_genes[complete.cases(names(sig_genes))]
    sig_genes = sig_genes[unique(names(sig_genes))] 
    spia.res <- spia(de=sig_genes, all=all_genes, organism=org)
    spia.res$KEGGLINK <- paste0("<a href='",spia.res$KEGGLINK,"' target='_blank'>",spia.res$KEGGLINK,"</a>")
    spia.res
  })
  
  output$spiaout <- DT::renderDataTable({
    input$runspia
    withProgress(session = session, message = 'Calculating...',detail = 'This may take a while...',{
      isolate({
        DT::datatable(spia_out(),escape = FALSE,selection = 'none',
                      extensions = c('TableTools','ColVis','Scroller'),
                      options = list(
                        dom = 'RMDCT<"clear">lfrtip',
                        searchHighlight = TRUE,
                        tableTools = list(sSwfPath = '//cdnjs.cloudflare.com/ajax/libs/datatables-tabletools/2.1.5/swf/copy_csv_xls_pdf.swf'),
                        pageLength = 5,
                        lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
                        initComplete = JS("function(settings, json) {",
                                          "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                                          "}"),
                        scrollX = TRUE
                      ))
      })
    })
  })
  
  # update tab9 with spia results
  observe({
    if(input$runspia > 0)
    {
      updateTabsetPanel(session = session, inputId = 'tabvalue', selected = 'tab9')
    }
  })
  
  # input file reactive function
  datasetInput7 <- reactive({
    infile <- input$inputdataset
    if(is.null(infile))
      return(NULL)
    read.csv(infile$datapath,check.names=F)
  })
  
  # update the dropdown with the csv file
  observe({
    dd <- datasetInput7()
    updateSelectInput(session = session, inputId = "projects", choices = as.character(dd$Project), selected = "none")
  })
  
})
