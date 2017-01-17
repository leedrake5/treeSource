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

fileInput(inputId='moderndefs', label='Choose modern sequences',
accept = c('.txt', '.rwl'), multiple=TRUE
),

numericInput(inputId='timemin', label ="Beginning", value=1600),
numericInput(inputId='timemax', label ="End", value=1900),

tags$hr(),


actionButton(inputId='sourcegrid', label="Run Full Model"),
actionButton(inputId='myjacknife', label="Run One Source"),

uiOutput('sourcenamesui'),

tags$hr(),

downloadButton(outputId="downloadtable", label="Download")

),
mainPanel(
tabsetPanel(
tabPanel('Console',
tableOutput('verbatimTextOutput("console")')),
tabPanel('Single Source Table',
DT::dataTableOutput('jackknifetable')),
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

uiOutput('sourcenamesmap'),

numericInput(inputId='resolution', label="Resolution", value=50),

downloadButton(outputId="downloadmap", label="Download Map")



),
mainPanel(
tabsetPanel(
tabPanel("Source Map",
plotOutput("sourcemap")),
tabPanel("Source Info",
plotOutput("simplesigtable"))

))




)
)))
