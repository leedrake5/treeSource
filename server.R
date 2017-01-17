library(shiny)
library(ggplot2)
library(reshape2)
library(pbapply)
library(data.table)
library(DT)



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    


sourceList <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    
    files <- inFile$datapath
    
    temp = inFile$name
    
    source.list <- lapply(files, function(x) read.csv(x, sep=","))
    source.list <- lapply(source.list, data.frame)

    names <- gsub(".csv", "", temp)
    

    
    names(source.list) <- names
    
    source.list
    

    
    })


chooseModernDefinitions <- reactive({
    inFile <- input$moderndefs
    if (is.null(inFile)) return(NULL)
    
    files <- inFile$datapath
    
    raw <- pblapply(files, readRWLArima)
    raw
})

modernDefNames <- reactive({
    inFile <- input$moderndefs
    if (is.null(inFile)) return(NULL)
    
    names <- inFile$name
    names


})

treeDataframe <- reactive({
    
    inFile <- input$file2
    δενδρος <- read.csv(inFile$datapath, sep=",")
    δενδρος
    
})





reactiveJackKnife <- reactive({
    
    χρονοσαρχειν <- as.numeric(as.vector(input$timemin))
    χρονοστελειν <- as.numeric(as.vector(input$timemax))
    πηγηΛογοι <- sourceList()
    δενδροςΚαταστασις <- treeDataframe()
    γενος <- "t-value"
    
    συλλογισμοσμοι <- treeJackKnifeMultipleSourceSig(timemin=χρονοσαρχειν, timemax=χρονοστελειν, tree.dataframe=δενδροςΚαταστασις, tree.source.list = πηγηΛογοι, γενος)
    
    συλλογισμοσμοι.λεγειν <- colnames(συλλογισμοσμοι)
    δενδροςΚαταστασις.λεγειν <- colnames(δενδροςΚαταστασις)
    
    
    συλλογισμοσμοι <- data.frame(δενδροςΚαταστασις.λεγειν[2:length(δενδροςΚαταστασις.λεγειν)],συλλογισμοσμοι)
    colnames(συλλογισμοσμοι) <- c("Tree", συλλογισμοσμοι.λεγειν)
    
    συλλογισμοσμοι



})

reactivePopSig <- reactive({
    
    πηγηΛογοι <- sourceList()
    γακΜαχαιρα <- reactiveJackKnife()
    
    σθλλογισμος <- populationSigDefinitionMultiple(x=γακΜαχαιρα, source.hypotheses <- names(πηγηΛογοι))
    
    σθλλογισμος
    
    
    
})


sourceGridReactive <- reactive({
    λεγεινμοι <- gsub(".txt", "", modernDefNames() )
    πηγηΛογοι <- sourceList()
    ΡΩΛμοι <- chooseModernDefinitions()
    χρονοσαρχειν <- as.numeric(as.vector(input$timemin))
    χρονοστελειν <- as.numeric(as.vector(input$timemax))
    δενδροςΚαταστασις <- treeDataframe()
    
    γακΜαχαιραμοι <- lapply(ΡΩΛμοι,  function(x) treeJackKnifeMultipleSourceSigHist(timemin=χρονοσαρχειν, timemax=χρονοστελειν, x, tree.source.list=πηγηΛογοι, return="t-value"))
    μεγεθοσΛογος <- lapply(γακΜαχαιραμοι, function(x) populationSigDefinitionMultiple(x, names(πηγηΛογοι)))
    συλλογισμοσμοι <- ldply(μεγεθοσΛογος, data.frame)
    συλλογισμοσμοι.λεγειν <- colnames(συλλογισμοσμοι)
    συλλογισμοσμοι <- data.frame(λεγεινμοι,συλλογισμοσμοι)
    colnames(συλλογισμοσμοι) <- c("Location", συλλογισμοσμοι.λεγειν)
    
    συλλογισμοσμοι

})



sourceNames <- reactive({
    λεγειν <- names(sourceList())
    λεγειν
    
})

output$sourcenamesui <- renderUI({
    selectInput(inputId = "sourcenames", label = h4("Source"), choices =  sourceNames())
})

output$sourcenamesmap <- renderUI({
    selectInput(inputId = "sourcenames", label = h4("Source"), choices =  sourceNames())
})


logText <- reactive({
    values[["log"]] <- capture.output(data <- queryMagic())
    
    
})


output$console <- renderPrint({
    logText()
    return(print(values[["log"]]))
    # You could also use grep("Warning", values[["log"]]) to get warning messages and use shinyBS package
    # to create alert message
})


output$sourcegridTable <- DT::renderDataTable({
    input$sourcegrid
    df <- sourceGridReactive()
     DT::datatable(df)
})

output$jackknifetable <- DT::renderDataTable({
    input$myjacknife
    df <- reactiveJackKnife()
    DT::datatable(df)
})


output$sourcemap <- renderPlot({
    
    NM <- get_map(location=c(input$mapregion), zoom=input$zoom, maptype=input$plottype)

    results.chu <- data.frame(results$Long, results$Lat, results[input$sources])
    colnames(results.chu) <- c("Long", "Lat", "Source")
    
    results.int.chu <- with(results.chu, interp(x=Long, y=Lat, z=Source, duplicate="mean", nx=input$resolution, ny=input$resolution))
    results.int.chu.melt <- reshape2::melt(results.int.chu$z, na.rm=TRUE)
    colnames(results.int.chu.melt) <- c("x", "y", "z")
    
    results.int.chu.melt$x <- results.int.chu$x[results.int.chu.melt$x]
    results.int.chu.melt$y <- results.int.chu$y[results.int.chu.melt$y]
    
    results.int.chu.melt[is.na(results.int.chu.melt)] <- 0
    
    ChuskaSigPlot <- ggmap(NM, alpha=0.2) +
    #stat_contour(data=results.int.chu.melt, aes(x=x, y=y, z=z, fill=z, alpha=z))+
    #stat_contour(data=results.int.cib.melt, aes(x=x, y=y, z=z, fill=z, alpha=z))+
    geom_tile(data=results.int.chu.melt, aes(x=x, y=y, z=z, fill=z, alpha=z))+
    geom_point(aes(Long, Lat, size=Source), colour="black", alpha=0.2, data=results.chu) +
    geom_text(aes(x=-108.08, y=36.25, label="Chaco Canyon"))+
    geom_point(aes(x=-108.05, y=36.1), shape=17, size=5) +
    coord_equal() +
    scale_fill_gradient(name = "Source Correlation", low="white", high="red", na.value="white", labels = scales::percent, limits = c(0, 1)) +
    scale_alpha(guide="none") +
    scale_y_continuous("Latitude", limits = c(35, 38)) +
    scale_x_continuous("Longitude", limits = c(-110, -106)) +
    coord_map() +
    theme_light() +
    ggtitle(input$sources) +
    guides(size=FALSE)
    
    ChuskaSigPlot

})


})


