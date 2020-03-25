library(shiny)
#library(shinyIncubator)

shinyUI(fluidPage(
    # progressInit(),
	titlePanel(HTML("<h1><strong>SSD/TTC Analysis Tool</strong></h1><h4>Version 0.93 BETA</h4>"), 
			 windowTitle = "SSD/TTC Analysis Tool"),
     sidebarLayout(
          sidebarPanel(
               radioButtons("ssdttc", "Analysis Type", choices = c("SSD", "TTC"), inline=T),
               p("Input data files should be in Excel (.xls or .xlsx) or tab-delimited (.txt).  All character formatted columns
                               will be available to choose as species/chemical/group variables.  All numeric columns will be available to choose
                               as the analysis variable."),
               fileInput('file1', '1.) Choose file to upload (.xlsx, .xls or tab-delimited .txt):',
                         accept = c('.xlsx','.xls','.txt')
               ),
               uiOutput("sheets"),
               uiOutput("charvars"),
               uiOutput("numvars"),
               checkboxInput("calcgeommeans", "Calculate Geometric means?"),
               uiOutput("geommeansUI"),
               uiOutput("slider"),hr(),
               strong("Plot Options:"),
               uiOutput("uclbox"),
               uiOutput("listinmarginsbox"),
               selectInput("units", "Units", choices = c("mg/L", "ug/L")),
               numericInput("textSize", "Relative Text Size (Default is 1,  +/- for larger/smaller labels)", value = 1, step=.1, min=.5, max=2),
               hr(),
               strong("Analysis Options:"),
               uiOutput("leaveonebox"),
               uiOutput("addonebox"),br(),
               actionButton("run", "Run Analysis"),
               hr(),
                p("This web tool was created by QS-Informatics. Contact Jesse Krailler (krailler.j) with any comments or questions."),
               img(src="../../images/QSstacked.png", height=100),
               #,
               width=4
               
               
          ),
          mainPanel(
               tabsetPanel(id="inTabset", selected="Table",

                    tabPanel("Table", 
                             br(),
                             uiOutput("nrowmessage"),
                             tableOutput("FullTable")
                             )
                    ,  
                    tabPanel("Output", 
                             #textOutput("text"),
                             uiOutput("PlotImage"),
                             #plotOutput("analysisinfo", width="500px", height="500px"),
                             #tableOutput("analysisinfo2"),
                             br(),
                             #downloadButton('downloadExcel', 'Download Excel File'),
                             #downloadButton('downloadPDF', 'Download PDF File')
                             
                             uiOutput("Excelbutton"),br(),
                             uiOutput("PDFbutton")
                             )
               )
               
          )
     )
))