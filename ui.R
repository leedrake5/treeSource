library(shiny)


# Define UI for application that draws a histogram
shinyUI(navbarPage("Data Selecter", id="nav",
tabPanel("Load Data",
titlePanel("Uploading Files"),
sidebarLayout(
sidebarPanel(
fileInput(inputId='file1', label='Choose source sequences',
accept = c('.csv'), multiple=TRUE
),

fileInput(inputId='file2', label='Choose archaeological tree ring data',
accept = c('.csv')
),



tags$hr(),


actionButton(inputId='sourcegrid', label="Run Full Model"),


tags$hr(),

downloadButton(outputId="downloadtable", label="Download")

),
mainPanel(
tabsetPanel(
tabPanel('Single Source Table',
DT::dataTableOutput('jackknifetable'))
)
)
)
),




tabPanel("Modern Source Definitions",
titlePanel("Modern Data"),
sidebarLayout(
sidebarPanel(

fileInput(inputId='moderndefs', label='Choose modern sequences',
accept = c('.txt', '.rwl'), multiple=TRUE
),

numericInput(inputId='timemin', label ="Beginning", value=1600),
numericInput(inputId='timemax', label ="End", value=1900),

actionButton(inputId='myjacknife', label="Run One Source"),

uiOutput('sourcenamesui')

),

mainPanel(
tabsetPanel(
tabPanel('Console',
tableOutput('verbatimTextOutput("console")')),
tabPanel('Full Table',
DT::dataTableOutput('sourcegridTable'))

)
)
)
),


tabPanel("Plot Map",
titlePanel("Create Plots"),
sidebarLayout(
sidebarPanel(

#uiOutput('sourcenamesmap'),

selectInput(inputId="sources", label="Choose Source", choices=c(
"Chuska Mountains" = "CHU",
"Zuni Mountains" = "CIB",
"Mesa Verde" = "MVER",
"North Jamez Mountains" = "MEA",
"South Jemez Mountains" = "JEM",
"San Juan Mountains" = "DUR",
"Mount Taylor" = "CEB",
"Gobernador" = "GOB"), selected="Chuska Mountains"),

textInput(inputId="mapregion", label="Map Region", value="New Mexico"),
sliderInput(inputId="zoom", label="Zoom", value=6, min=3, max=20),
selectInput(inputId="plottype", label="Plot Type", choices=c(
"Satellite" = "satellite",
"Terrain" = "terrain",
"Terrain Background" = "terrain-background",
"Watercolor" = "watercolor",
"Roads" = "roadmap",
"Hybrid" = "hybrid",
"Terrain Labels" = "terrain-labels",
"Terrain Lines" = "terrain-lines",
"Toner 2010" = "toner-2010",
"Toner 2011" = "toner-2011",
"Toner Lines" = "toner-lines",
"Toner Lite" = "toner-lite",
"Toner Background" = "toner-background",
"Toner Hybrid" = "toner-hybrid",
"Toner Labels" = "toner-labels"), selected="Terrain Background"),



sliderInput(inputId='resolution', label="Resolution", value=50, min=25,max=300),

downloadButton(outputId="downloadmap", label="Download Map")




),
mainPanel(
tabsetPanel(
tabPanel("Source Map",
plotOutput("sourcemap"),width="900px"),
tabPanel("Source Info",
plotOutput("simplesigtable"))

))




)
)))
