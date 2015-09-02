#---------------------------------------------------------------------------------------------------#
# AUTHOR : Grace Chang, Benjamin Hsieh
# Program name : lb_sdtm
# R Version : 3.2.0
# Purpose : LAB data cleaning, summary, and visualization
# Dataset: C:/Users/bhsieh/Rwd/ib_sdtm/LB_SDTM_CHEM_HEM.txt
# Output : Shiny App -  Data Tools to look for potential points of interest
# Version : 2.0 rCharts
# Date: 7/17/15

# Reason for modification: Tab1: Group by: lab test 
# 
# Modificaiton date: 8/18
# Modified By: BH, GC

#--------------------------------------------------------------------------------------------------#

#==================================================================================================#
# global.R
#==================================================================================================#

# Set working directory:
setwd("C:/Users/bhsieh/Rwd/ib_sdtm")

#--------------------------------------------------------------------------------------------------#
## capture all the output to a file, please choose your tracker file
# zz <- file("C:/Users/bhsieh/Rwd/ib_sdtm/trackerlog.txt", open = "a+")
# sink(zz)
# sink(zz, type = "message")
# Sys.time()

#--------------------------------------------------------------------------------------------------#
# Clear All objects
rm(list=ls())
#--------------------------------------------------------------------------------------------------#
# Load Libraries
library(ggplot2)
library(shiny)
# library(outliers)
# library(mvoutlier)
library(reshape2)
library(DT)
library(rCharts)
library(sas7bdat)
# library(timeline)
# library(googleVis)
library(dplyr)
library(GGally)
library(Hmisc)
# 

#--------------------------------------------------------------------------------------------------#
# Read in SAS
#--------------------------------------------------------------------------------------------------#
# df.lb.sdtm.v0 <<- read.sas7bdat( file = "lb.sas7bdat" )

# Change Date Format
## df.lb.sdtm.v0$LBTM <- as.Date( df.lb.sdtm.v0$LBTM, origin = "1960-01-01", tz = "UTC" )
## df.lb.sdtm.v0$LBDT <- as.Date( df.lb.sdtm.v0$LBDT, origin = "1960-01-01", tz = "UTC" )

#--------------------------------------------------------------------------------------------------#
# Clean [This will become reactive as file upload fucntion is added]
#--------------------------------------------------------------------------------------------------#
# Read in Text File
df.lb.sdtm.v0 <- data.frame(read.table("lb_sdtm_CHEM_HEM.txt", header=T, sep=","))

# Get change from baseline
df.lb.sdtm.v1 <- arrange(df.lb.sdtm.v0, Unique.Subject.Identifier,Lab.Test.or.Examination.Name, desc(Baseline.Flag))
df.lb.sdtm.v1 <- group_by(df.lb.sdtm.v1,Lab.Test.or.Examination.Name, Unique.Subject.Identifier)
df.lb.sdtm.v1$Change.From.Baseline.Numeric.Result <- do.call("c", 
                                                             with(df.lb.sdtm.v1, 
                                                                  lapply(split(Numeric.Result.Finding.in.Standard.Units,
                                                                               list(Lab.Test.or.Examination.Name, Unique.Subject.Identifier)), 
                                                                         function(x){
                                                                           x - x[1]
                                                                           }
                                                                         )
                                                             )
)
# Unique Subjects
unique.sub.id <- as.vector(sort(unique(df.lb.sdtm.v0$Unique.Subject.Identifier)))
names(unique.sub.id) <- unique.sub.id

# Unique variables
unique.var <- names(df.lb.sdtm.v0)

# Count number of subjects
count.id <- length(unique.sub.id)

# Count Number of Variable
count.var <- ncol(df.lb.sdtm.v0)

# Summary table
sum.table <- summary(df.lb.sdtm.v0)

# Lab Test Category
unique.cat.lab.test <- as.vector(sort(unique(df.lb.sdtm.v0$Category.for.Lab.Test)))
names(unique.cat.lab.test) <- unique.cat.lab.test

# Lab Test Subcatgory
sub.list <- lapply(unique.cat.lab.test, function(x){
  as.vector(unique(filter(df.lb.sdtm.v0, Category.for.Lab.Test == x)$Subcategory.for.Lab.Test))
})
names(sub.list) <- unique.cat.lab.test

# Unique Lab testname
unique.labname <- as.vector(sort(unique(df.lb.sdtm.v0$Lab.Test.or.Examination.Name)))
names(unique.labname) <- unique.labname

# Render Missing Variables
missing.var <- names(df.lb.sdtm.v0[apply(df.lb.sdtm.v0,2, 
                                         function(x){
                                           all(is.na(x))
                                           }
                                         )
                                   ]
                     )


# Non missing variables
nonmiss.unique.var <- unique.var[-which(unique.var %in% missing.var)]

# Render Visit Table
vtable <- df.lb.sdtm.v0 %>% 
  select(Visit.Number, Visit.Name, Study.Day.of.Specimen.Collection) %>%
  group_by(Visit.Number, Visit.Name) %>% 
  summarise(Range = paste(Min = min(Study.Day.of.Specimen.Collection, na.rm = T), Max=max(Study.Day.of.Specimen.Collection, na.rm = T), sep = " ~ "))
vtable <- vtable[, c(2,1,3)]

# Reference Range DF
ref.df <- df.lb.sdtm.v0 %>%
  select(Lab.Test.or.Examination.Name, Reference.Range.Lower.Limit.Std.Units, Reference.Range.Upper.Limit.Std.Units) %>%
  group_by(Lab.Test.or.Examination.Name) %>%
  summarise(Ref.Lower = min(Reference.Range.Lower.Limit.Std.Units, na.rm = T), Ref.Upper = max(Reference.Range.Upper.Limit.Std.Units, na.rm = T))

#--------------------------------------------------------------------------------------------------#
# Functions
#--------------------------------------------------------------------------------------------------#
## is.null get TF values
getTFIfNotNull <- function(df, variable, value){
  # Computes out.v, a logical vector of length=nrows(df), that checks whether value is contained in the corresponding variable/factor
  #
  # Arg:
  #   df: dataframe from which variable and value are taken
  #   variable: variable/factor name
  #   value: vector subset corresponding to the variable/factor
  #
  # Returns:
  #   1. If value is not NULL, returns logical vector from a match between value and entire variable column  
  #   2. If value is NULL (no subset taken)
  #     a) FALSE if variable name is Unique.Subject.Identifier (DEFAULT)
  #     b) TRUE otherwise: takes the entire subset
  if (!is.null(value)) {
    out.v <- df[,variable] %in% value
    } else {
      if (variable == "Unique.Subject.Identifier") {
        out.v <- FALSE
        } else {
          out.v <- TRUE
        }
      }
  return(out.v)
}

## Get current dataset 
getDataset <- function(df, tab.num, id = NULL, cat = NULL, subcat = NULL, test = NULL){
  # Computes a row subset of df based on selected vectors of Subject ID, Category of Lab Test, Subcategory of Lab Test, and Lab Test or Examination Name
  # 
  # Functions called:
  # getTFIfNotNull()
  #
  # Arg:
  #   df: dataframe to be subseted
  #   tab.num: version of function based on tab number
  #   id: vector of Unique.Subject.Identifiers
  #   cat: vector of Category.for.Lab.Test
  #   subcat: vector of Subcategory.for.Lab.Test
  #   test: vector of Lab.Test.or.Examination.Name
  # 
  # Returns:
  #   Dataframe with row subsets
  if (tab.num == 2) {
    dfout <- df[getTFIfNotNull(df, "Unique.Subject.Identifier",id) & getTFIfNotNull(df, "Category.for.Lab.Test", cat) & 
                    getTFIfNotNull(df, "Subcategory.for.Lab.Test", subcat) & getTFIfNotNull(df, "Lab.Test.or.Examination.Name", test),]
  } else {
    if (tab.num == 3) {
      dfout <- df[df$Lab.Test.or.Examination.Name==test,]
    }
  }
  return(dfout)
 
}

# Close Tracker
closetracker <- function(){
  # function that closes the tracker
  sink()
  closeAllConnections()
}

#==================================================================================================#
# UI.R
#==================================================================================================#


ui <- navbarPage("Lab SDTM Tools",
  #------------------------------------------------------------------------------------------------#
  # Tab 1 Initial Summary
  #------------------------------------------------------------------------------------------------#
  tabPanel(title="Summary",
           sidebarLayout(
             #---------------------------------------------------------------------------------------#
             sidebarPanel(
               h2("Welcome"),
               p("Hi, this is the Shiny application used for data cleaning and visualization."),
               # Upload File [STILL NEEDS TO BE IMPLEMENTED IN SERVER]
               div(class = "option-group",
                   fileInput(inputId = 'file1', 
                             label  = "Please upload your sas/txt file below:", 
                             accept = c('text/csv','text/comma-separated-values,text/plain', '.csv')
                             ),
                   tags$hr(),
                   checkboxInput(inputId = 'header', 
                                 label   = 'Header', TRUE
                                 ),
                   radioButtons(inputId  = 'sep', 
                                label    = 'Separator', 
                                choices  = c(Comma = ',', 
                                            Semicolon = ';', 
                                            Tab = '\t'
                                            ), 
                                selected = ',', 
                                inline   = T
                                ),
                   radioButtons(inputId = 'quote', 
                                label   = 'Quote', 
                                choices = c(None = '', 
                                            'Double Quote' = '"', 
                                            'Single Quote' = "'"
                                            ), 
                                selected = '"', 
                                inline = T
                                )
                   ),
               # [Will become reactive as we add the file upload function IN SERVER]
               div(class = "option-group",
                   h5("There are", span(count.id, style = "color:blue"), "unique subjects and", span(count.var, style = "color:blue"), "parameters."),
                   h5("Variables with all observations missing:"),
                   strong(paste(missing.var, collapse = " | "))
                   ),
               div(class = "option-group",
                   selectizeInput(inputId = "t1.choose.test", 
                                  label = "Group by a lab test", 
                                  choices = c("All tests",unique.labname) 
                                  ),
                   selectizeInput(inputId = "choosevar", 
                                  label = "Please choose a variable for summary statistics", 
                                  choices = nonmiss.unique.var, 
                                  multiple = T, 
                                  options = list(placeholder = 'Please select subject(s) below')
                                  )
                   )
               ),
             #---------------------------------------------------------------------------------------#
             mainPanel(
               # Summary given choosen variable, from package hmisc: describe
               fluidRow(
                 column(12,
                        textOutput("testtext"),
                        verbatimTextOutput("describe.text")
                        )
               ), 
               # Data table summary from DT
               fluidRow(
                 column(12,
                        # Show only if length of choosen variables is 1
                        conditionalPanel(
                          condition = "length(input.choosevar) == 1", 
                          dataTableOutput("summary.table")
                          )
                        )
                 )
               )
             #---------------------------------------------------------------------------------------#
             )
           ),
  
  #------------------------------------------------------------------------------------------------#
  # Tab 2 Subsetting Plot
  #------------------------------------------------------------------------------------------------#  
  tabPanel(title = "Tool 1",
           # Create Custom CSS for grouping and headers
           tags$head(
             tags$style(HTML("
                             .option-group {
                             border: 1px solid #ccc;
                             border-radius: 4px;
                             padding: 3px 8px;
                             margin: 3px -8px;
                             background-color: #f5f5f5;
                             }
                             
                             .option-header {
                             color: #79d;
                             text-transform: uppercase;
                             margin-bottom: 6px;
                             }
                             .side-side{
                             display:inline-block;
                             max-height: 50px;
                             vertical-align: top;
                             }
                            
                             ")
                        ),
             tags$style(type = 'text/css',
                            "
                            .btn-small {
                            border: 1px solid #ccc;
                            border-radius: 4px;
                            padding: 6px 12px;
                            vertical-align: middle;
                            max-width: 170px;
                            max-height: 34px;
                            line-height: 20px;
                            font-size: 14px;
                            font-weight: normal;
                            background-color: #fff;
                            }
                            .shiny-input-checkboxgroup label ~ .shiny-options-group, .shiny-input-radiogroup label ~ .shiny-options-group {
                            margin: auto;
                            }

                            .checkbox{
                            margin: 0px;  
                            margin-bottom: 10px;
                            }

                            #testtext{
                                 font-size: 20px;
                                 font-style: bold;
                                 }
                            ")
             ),
           
           #---------------------------------------------------------------------------------------#
           sidebarPanel(width = 2,
                        
                        # Select All/Manual select
                        radioButtons(inputId = "selectType",
                                     label   = h4("Subject Selection Type"),
                                     choices = c(
                                       "Manual Select" = "mselect", 
                                       "Select All" = "allselect"
                                       ),
                                     selected = "mselect"
                        ),
                        # Dynamic Choose Subject depending on selectType
                        uiOutput("choose.subject"),
                        # Dynamic Subsetting based on selection of Unique.Subject.Identifier
                        # Category
                        uiOutput("choose.cat"),
                        # Lab test
                        uiOutput("choose.lab")
                        ),
           #---------------------------------------------------------------------------------------#
           mainPanel(width = 10,
                     fluidRow(
                              # div of class option-group allows grouping of widgets
                              div(class = "option-group",
                                  # inline-block allows widgets to be displayed side by side
                                  div(style = "display:inline-block; margin: 3px; vertical-align: middle;",
                                      div(class = "option-header",
                                          "Brush to Zoom, Click to Reset, Hover over right plot for tooltip"
                                          )
                                      ),
                                  div(style = "display:inline-block; margin: 3px; vertical-align: middle;",
                                      div(class = "btn-small",
                                      # Checkbox to turn on/off the display of plotting options
                                      checkboxInput(inputId = "plotoptions", 
                                                    label = "Show Plot Options?", FALSE)
                                      )
                                  ),
                                  div(style = "display:inline-block; margin: 3px; vertical-align: middle;",
                                      div(class = "btn-small",
                                          checkboxInput(inputId = "t2showline",
                                                        label = "Show Plot Lines?",
                                                        value = FALSE)
                                          )
                                      ),
                                  div(style = "display:inline-block; margin: 3px; vertical-align: middle;",
                                      div(class = "btn-small",
                                          checkboxInput(inputId = "t2showrefline",
                                                        label = "Show Ref Lines?",
                                                        value = FALSE)
                                      )
                                  ), 
                                  
                                  div(style = "display:inline-block; margin: 3px; vertical-align: middle;",
                                      # download main scatter plot button
                                      downloadButton("downloadplot", label = "Download Plot as PDF")
                                  )
                              ),
                              # Dynamic main plot output
                              column(width = 7, 
                                     uiOutput("plotui"),
                                     # WellPanel wraps table in a nice box
                                     # Data Table from sidepanel subsetting
                                     wellPanel(width = 6,
                                               h4("Data Table from Subset"), 
                                               dataTableOutput("table.subset")
                                               ),
                                     # Dynamic plot options that displays when "plotoptions" is TRUE
                                     column(6,
                                            uiOutput("plotoption1"),
                                            uiOutput("plotoption2")
                                            ),
                                     column(6,
                                            uiOutput("plotoption3")
                                            )
                                     ),
                              column(5,
                                     # Zoom plot output
                              uiOutput("zoomplotui"),
                              verbatimTextOutput("hoverpoints.text"),
                              wellPanel(width = 6,
                                        h4("Visit/Study Table"), 
                                        dataTableOutput("visit.table")
                              ),
                              # Data Table from brushing
                              wellPanel(width = 6, 
                                        h4("Points from Brushing"),
                                        dataTableOutput("table.brushed")
                                        ),
                              wellPanel(width = 6,
                                  radioButtons("select.download.type", 
                                               label    = "Download brushed table as:", 
                                               choices  = list("CSV" = "csv","Text" = "txt"), 
                                               selected = "csv", 
                                               inline   = T
                                               ),
                                  downloadButton(outputId = "downloadtable", 
                                                 label="Download Brushed Table"
                                                 )
                              )
                       )
                     )
           )# end main panel
           #---------------------------------------------------------------------------------------#
  ), # End Tab
  
  #------------------------------------------------------------------------------------------------#
  # Tab 3 Change from baseline plots
  #------------------------------------------------------------------------------------------------#
  tabPanel("Tool 2",
           sidebarLayout(           
             #-------------------------------------------------------------------------------------#
             sidebarPanel(width = 2,
               h3("Change From Baseline Plots"),
               # Choose lab test
               selectizeInput(inputId = "t3.choose.lab.testname", 
                              label   = "Please choose category of lab test", 
                              choices = unique.labname,
                              options = list(placeholder = 'Choose a lab test')
                              ),
               # Choose type of plot
               checkboxGroupInput(inputId  = "cfbPlotType", 
                                  label    = "Type of Plot", 
                                  choices  = list("Scatter"  = "scat", 
                                                "Line"       = "line", 
                                                "Box"        = "box", 
                                                "Histogram"  = "hist"
                                                ), 
                                  selected = "scat")
               ),
             
             #-------------------------------------------------------------------------------------#
             mainPanel(width = 10,
                       div(class = "option-header",
                           "Brush and Double-Click to Zoom, Double-Click to Reset"
                           ),
                       helpText("Please wait while rChart loads on the top right, then click a few times on the graph to enable tooltips."),
                       fluidRow(
                         column(6,
                                div(class="option-group",
                                    # Scatter Plot Output
                                    plotOutput("scatgout",
                                               dblclick = "t3scat.dblclick",
                                               brush = brushOpts(
                                                 id = "t3scat.brush",
                                                 resetOnNew = T)
                                               )
                                    )
                                ),
                         column(6,
                                div(class = "option-group",
                                    # Scatter Plot Output
                                    showOutput("scatcout","nvd3")
                                    )
                                )
                         ),
                       fluidRow(
                         column(12,
                                div(class="option-group",
                                    # Scatter plot table of brushed points
                                    dataTableOutput("t3.scatsubset")
                                    )
                                )
                         ),
                       fluidRow(
                         column(6,
                                div(class = "option-group",
                                    # Boxplot output
                                    plotOutput("box.out",
                                               dblclick = "t3box.dblclick",
                                               brush    = brushOpts(
                                                 id = "t3box.brush",
                                                 resetOnNew = T)
                                               )
                                    )
                                ),
                         column(6,
                                div(class = "option-group",
                                    # Histogram Output
                                    plotOutput("hist.out",
                                               dblclick = "t3hist.dblclick",
                                               brush    = brushOpts(
                                                 id = "t3hist.brush",
                                                 resetOnNew = T)
                                               )
                                    )
                                )
                         )
                       )
             #-------------------------------------------------------------------------------------#
             ) # End sidebarLayout
           ) # End Tab 3
) # End ui


#==================================================================================================#
# server.R
#==================================================================================================#

server <- function(input, output, session) {
  
  #------------------------------------------------------------------------------------------------#
  # Tab 1
  #------------------------------------------------------------------------------------------------#
  # Upload File as SAS for txt (sdtm format)
  # [CODE HERE]
  
  
  # Create Describe Object
  t.describe <- reactive({
    if (input$t1.choose.test == "All tests") {
      return(describe(df.lb.sdtm.v0))
    } else {
    return(describe(filter(df.lb.sdtm.v0, Lab.Test.or.Examination.Name == input$t1.choose.test)))
    }
  })
  
  # Render Chosen Test
  output$testtext <- renderText({
    paste("Summary below, grouped by: ", input$t1.choose.test)
  })
  
  # Render Describe Text  
  output$describe.text <- renderPrint({  
    # To get rid of describe formatting [the dividing lines are wider than the mainpanel], have to 
    # unlist the describe object then assign elements of each object as a seperate object, coerce to new list and return
    if (is.null(input$choosevar)) {
      return()
    }
    this.describe <- t.describe()
    # get index of choosen
    choosen.var.index <- which(names(this.describe) %in% input$choosevar)
    # create the describe subset
    temp.desc <- this.describe[input$choosevar]
    list.out  <- lapply(temp.desc,function(x){
      this.describe[[
        choosen.var.index[which(input$choosevar %in% x)]
                  ]]
    })
    list.out
  })
  # Render Table of Summary, return if length is not equal to 1
  output$summary.table <- renderDataTable({
     if (length(input$choosevar) != 1 ) {
       return()
     }
     if (input$t1.choose.test == "All tests") {
       return(data.frame(table(df.lb.sdtm.v0[input$choosevar])))
     } else {
       data.frame(table(filter(df.lb.sdtm.v0, Lab.Test.or.Examination.Name == input$t1.choose.test)[input$choosevar]))
     }
  })
  
  #------------------------------------------------------------------------------------------------#
  # Tab 2
  #------------------------------------------------------------------------------------------------#
  
  # Render Visit Table
  output$visit.table <- renderDataTable({
    datatable(data = vtable, 
              options = list(scrollX = TRUE,
                             pageLength = 5
                             )
    )
  })
  
  
  #------------------------------------------------------------------------------------------------#
  # Subsetting
  #------------------------------------------------------------------------------------------------#
  
  # Create Reactive function to get the current subject id, depending on selection type
  cursubjectid <- reactive({
    if (input$selectType == "mselect"){
      return(input$sub.id)
    }
    if (input$selectType == "allselect")
      return(c(names(unique.sub.id)))
  })
  
  # Reset subjects
  observe({
    # Observes the actionbuttion "reset.all", a counter that resets sub.id as it is pressed
    input$reset.all
    updateSelectizeInput(session  = session, 
                         inputId  = "sub.id",
                         choices  = unique.sub.id,
                         selected = c()
    )
  })
  
  # Get the current selected dataset/subset
  t2.curdata <- reactive({
    df.out <- getDataset(df = df.lb.sdtm.v0, tab.num = 2, id = cursubjectid(), cat = input$t2.cat.lab, subcat = input$t2.sub.lab, test = input$t2.lab.testname)
    df.out
  })
  #------------------------------------------------------------------------------------------------#
  
  # Render the subject selection based on selectType
  output$choose.subject <- renderUI({
    if (input$selectType == "mselect"){
      tagList(
        actionButton(inputId = "reset.all",
                     label = "Reset Inputs" 
                     ),
        # Select Unique Subject Identifier
        selectizeInput(inputId  = 'sub.id', 
                       label    = h4("Unique Subject Identifier"),
                       multiple = TRUE, 
                       choices  = unique.sub.id,
                       options  = list(placeholder = 'Please select subject(s)')
                       )
        )
    }
  })
  
  
  # Category dynamic to subject id
  output$choose.cat <- renderUI({
    subject<-cursubjectid()
    thisdata <- filter(df.lb.sdtm.v0, Unique.Subject.Identifier %in% subject)
    category <- levels(droplevels(thisdata$Category.for.Lab.Test))
    selectizeInput(inputId  = "t2.cat.lab", 
                   label    = h4("Category of Lab Test"), 
                   choices  = category, 
                   multiple = T, 
                   options  = list(placeholder = 'Please select Category')
                   )
    
    })
  
  # Subcategory, dynamic to category of lab test [defunct]
#   output$t2.choose.sub.lab <- renderUI({
#     # If missing category, return to avoid error
#     if (is.null(input$t2.cat.lab)) {
#       return()
#     }
#     # If missing subject, return to avoid error
#     if (is.null(cursubjectid())) {
#       return()
#     }
#     # If subject does not have a subcategory, return to avoid error
#     if (all(subset(df.lb.sdtm.v0, (Unique.Subject.Identifier == cursubjectid()) & (Category.for.Lab.Test == input$t2.cat.lab), select=Subcategory.for.Lab.Test)=="")) {
#       return()
#     }
#     # Get selected category
#     selectizeInput('t2.sub.lab', 'Subcategory', choices = sub.list[input$t2.cat.lab], multiple = TRUE)
#   })
#   
  # Lab name, dynamic to Category and Subject
  output$choose.lab <- renderUI({
    category <- input$t2.cat.lab
    subject <- cursubjectid()
    thisdata <- filter(df.lb.sdtm.v0, Unique.Subject.Identifier %in% subject, Category.for.Lab.Test %in% category)
    labtest <- levels(droplevels(thisdata$Lab.Test.or.Examination.Name))
    selectizeInput(inputId  = "t2.lab.testname", 
                     label    = h4("Lab Test or Name"), 
                     choices  = labtest, 
                     multiple = T, 
                     options  = list(placeholder = 'Please select test(s)')
    )
    
  })
  
  #-------------------------------------------------------------------------------------------------#
  # PLOTTING
  #-------------------------------------------------------------------------------------------------#
  
  # Add plot line
  plotlines <- reactive({
    if (input$t2showline==TRUE) {
      return(geom_line(aes(colour = Lab.Test.or.Examination.Name), size = 1))
    } else {
      return(NULL)
      }
  })
  
  # Add Reference Lines
  reflines <- reactive({
    if (input$t2showrefline == TRUE) {
      thisdata <- t2.curdata()
      ggout <- list(
        geom_hline(data = thisdata, aes(yintercept = Reference.Range.Upper.Limit.Std.Units, colour = Lab.Test.or.Examination.Name), linetype = 2),
        geom_hline(data = thisdata, aes(yintercept = Reference.Range.Lower.Limit.Std.Units, colour = Lab.Test.or.Examination.Name), linetype = 2)
      )
      return(ggout)
    } else {
      return(NULL)
      }
  })
  
  # Create ggplot object
  plot.it <- reactive({
    # We plot it with a call to current data function t2.curdata()
    # get current dataset
    data.in <- t2.curdata()
    # if subset data has no numeric results (i.e. categorical variable results), then return
    if (all(is.na(data.in$Numeric.Result.Finding.in.Standard.Units))) {
      return()
    } else {
      plot.out<-ggplot(data = data.in, aes(x = Study.Day.of.Specimen.Collection,
                                           y = Numeric.Result.Finding.in.Standard.Units,
                                           group = Lab.Test.or.Examination.Name)) +
        geom_point(colour = "grey50", size = 4, na.rm = T) +
        geom_point(aes(colour = factor(Lab.Test.or.Examination.Name)), na.rm = T) + 
        plotlines() + 
        reflines() +
        theme_set(theme_gray(16)) +
        theme(axis.title.y    = element_text(size = rel(0.7), angle = 90),
              axis.title.x    = element_text(size = rel(0.7)),
              legend.position = "bottom", legend.text = element_text(size = 11)) + 
        guides(colour = guide_legend(title = NULL, ncol = 3)) +
        ylab("Standard Units")+
        xlab("Study Day")
      return(plot.out)
    }
  })
  
  # Render Plot
  output$plot <- renderPlot({
    plot.it()  
  })
  
  
  #-------------------------------------------------------------------------------------------------#
  # Interactive Plot
  
  # Dynamically changing the plot
  output$plotui <- renderUI({
    plotOutput("plot",
               height   = 800,
               click    = "plot.click",
               dblclick = dblclickOpts(
                 id     = "plot.dblclick",
                 delay  = input$dblclick.delay
               ),
               brush    = brushOpts(
                 id         = "plot.brush",
                 delay      = input$brush.delay,
                 delayType  = input$brush.policy,
                 direction  = input$brush.dir,
                 resetOnNew = input$brush.reset
                 )
               )
    })
  # Reactive Function for brushing
  t2.getbrushset <- reactive({
    dat <- t2.curdata()
    brush.data <- brushedPoints(dat, input$plot.brush)
    # Replace NAs with empty string to ensure each entry has a value and prevents corruption
    brush.data[is.na(brush.data)] <- ""
    brush.data <- brush.data[complete.cases(brush.data),]
    return(brush.data)
  })

  # Reactive Function for Hovering
  output$hoverpoints.text <- renderPrint({
    cat("Info from Hover:\n")
    if (is.null(input$plot.hover)) {
      return()
    }
    hoverinfo <- input$plot.hover
    dat <- isolate({t2.curdata()})
    xval <- hoverinfo$x
    yval <- hoverinfo$y
    dist <- sqrt((xval - dat$Study.Day.of.Specimen.Collection) ^ 2 + (yval - dat$Numeric.Result.Finding.in.Standard.Units) ^ 2)
    dfinfo <- select(dat, Unique.Subject.Identifier, Visit.Number, Numeric.Result.Finding.in.Standard.Units)
    dfinfo <- rename(dfinfo, Numeric.Result = Numeric.Result.Finding.in.Standard.Units)
    hoverdata <- dfinfo[which(dist < 2),]
    if (all(is.na(hoverdata))) {
      return(NULL)
    } else {
      return(hoverdata)
    }
  })
  
  #-------------------------------------------------------------------------------------------------#
  # Zoom Plot
  
  # Reactive function that takes the range of zoomed plot
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  # Dynamic Plot (hovering)
  output$zoomplotui <- renderUI({
    plotOutput("zoom.plot", 
               hover    = hoverOpts(
                 id          = "plot.hover",
                 delay       = input$hover.delay,
                 delayType   = input$hover.policy,
                 nullOutside = input$hover.null.outside
                 ),
               height = 600)
  })
  
  # Render Plot
  output$zoom.plot <- renderPlot({
    if (is.null(plot.it())) {
      return()
      }
    plot.it() + 
      coord_cartesian(xlim = ranges$x, ylim = ranges$y) + 
      theme(legend.position = "none")
    })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observe({
    brush <- input$plot.brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
        }
    })
  
  #--------------------------------------------------------------------------------------------------#
  # Hide/Show Plot Options
  
  # Double Click Options
  output$plotoption1 <- renderUI({
    if (input$plotoptions == FALSE) {
      return()
    } else {
      div(class = "option-group",
          div(class = "option-header", "Double-click"),
          sliderInput(inputId = "dblclick.delay", 
                      label   = "Delay", 
                      min     = 100, 
                      max     = 1000, 
                      value   = 400, 
                      step    = 100)
      )
      }
    })
  
  # Hover Options
  output$plotoption2 <- renderUI({
    if (input$plotoptions == FALSE) {
      return()
      } else {
        div(class = "option-group",
            div(class = "option-header", "Hover"),
            radioButtons("hover.policy", "Input rate policy",
                         c("debounce", "throttle"), inline = TRUE),
            sliderInput(inputId = "hover.delay", 
                        label   = "Delay", 
                        min     = 100, 
                        max     = 1000, 
                        value   = 200, 
                        step    = 100),
            checkboxInput(inputId = "hover.null.outside", 
                          label   = "NULL when outside", 
                          value   = TRUE)
        )
        }
    })
  
  # Brush Options
  output$plotoption3 <- renderUI({
    if (input$plotoptions == FALSE) {
      return()
      } else {
        div(class = "option-group",
            div(class = "option-header", "Brush"),
            radioButtons(inputId = "brush.dir", 
                         label   = "Direction(s)",
                         choices = c("xy", "x", "y"), 
                         inline  = TRUE),
            radioButtons(inputId = "brush.policy", 
                         label   = "Input rate policy",
                         choices = c("debounce", "throttle"), 
                         inline  = TRUE),
            sliderInput(inputId = "brush.delay", 
                        label   = "Delay", 
                        min     = 100, 
                        max     = 1000, 
                        value   = 200, 
                        step    = 100),
            checkboxInput(inputId = "brush.reset", 
                          label   = "Reset on new image")
        )
        }
    })

  #-------------------------------------------------------------------------------------------------#
  # Tables
  #-------------------------------------------------------------------------------------------------#
  
  # Render Table from Subsetting
  output$table.subset <- renderDataTable({
    datatable(data = t2.curdata(), 
              options = list(scrollX = TRUE)
              )
    })
  # Render Table from Brushing
  output$table.brushed <- renderDataTable({
    datatable(data = t2.getbrushset(), 
              options = list(scrollX = TRUE)
              )
    })
  
  
  #-------------------------------------------------------------------------------------------------#
  # Download
  #-------------------------------------------------------------------------------------------------#
  
  # Download Plot
  output$downloadplot <- downloadHandler(
    # save ggplot as pdf
    filename <- paste("myplot", Sys.Date(), ".pdf", sep = ''),
    content  <- function(file) {
      device <- function(..., width, height) {
        pdf(..., width = width, height = height)
        }
      ggsave(file, plot = plot.it(), device = device)
      }
    )
  
  # Download table
  output$downloadtable <- downloadHandler(
    # Write as cvs/txt file
    # Only works in browser
    filename <- function(){
      paste("brushtable", input$select.download.type, sep = ".")
      },
    content <- function(file) {
      sep <- switch(input$select.download.type, 
                    "csv" = ",", 
                    "txt" = "\t")
      write.table(t2.getbrushset(), file = file, sep = sep, row.names = FALSE)
      }
    )
  
  
  #-------------------------------------------------------------------------------------------------#
  # Tab 3
  #-------------------------------------------------------------------------------------------------#
  
  #-------------------------------------------------------------------------------------------------#
  # PLOTTING
  #-------------------------------------------------------------------------------------------------#
  
  # Get current data
  t3.curdata <- reactive({
    getDataset(df = df.lb.sdtm.v1, 
               tab.num = 3,
               test = input$t3.choose.lab.testname
               )
    })
  
  #-------------------------------------------------------------------------------------------------#
  # Reactive Functions for brushing, see tab 2 for details
  
  # [Possible IMPLMENTATON: To get rid of delay when checkboxes are changed, have  "isolate" each checkbox, so use 4 checkboxinput instead of a checkbox group]
  # Scatter/line plot
  t3.getbrushset.scat <- reactive({
    dat <- t3.curdata()
    brush.data <- brushedPoints(dat, input$t3scat.brush)
    brush.data[is.na(brush.data)] <- ""
    brush.data <- brush.data[complete.cases(brush.data),]
    return(brush.data)
    })
  
  # Box plot
  t3.getbrushset.hist <- reactive({
    dat <- t3.curdata()
    brush.data <- brushedPoints(dat, input$t3hist.brush)
    brush.data[is.na(brush.data)] <- ""
    brush.data <- brush.data[complete.cases(brush.data),]
    return(brush.data)
    })
  
  # Histogram
  t3.getbrushset.box <- reactive({
    dat <- t3.curdata()
    brush.data <- brushedPoints(dat, input$t3box.brush)
    brush.data[is.na(brush.data)] <- ""
    brush.data <- brush.data[complete.cases(brush.data),]
    return(brush.data)
    })
  
  #-------------------------------------------------------------------------------------------------#
  # Get gg object for plots
  getChangeFromBaseLine.gg <- reactive({
    # If no lab test chosen, return to avoid error
    # If no numeric values under selected factor, return string in console
    # Else return ggplot object 
    if (is.null(input$t3.choose.lab.testname)) {
      return()
      }
    gg.df <- t3.curdata()
    if (all(is.na(gg.df$Change.From.Baseline.Numeric.Result))) {
      return("No Numeric Values")
      } else {
        t3.ggout <- ggplot(data = gg.df)
        return(t3.ggout)
        }
    })
  #-------------------------------------------------------------------------------------------------#
  # Plot it
  
  # Scatter ggplot
  output$scatgout <- renderPlot({
    # If scat is checked in sidePanel (inputId=cfb.plot.type) & call to getChangeFromBaseLine.gg()
    # does not return a character object, then returns a scatter plot with:
    #   1. x = Sequence.Number 
    #   2. y = Change.From.Baseline.Numeric.Result
    #   3. colour = by Unique.Subject.Identifier
    #   4. no legend
    #   5. zoom coordinates
    if (!(any(input$cfbPlotType %in% "scat"))) {
      return()
    }
    obj.out<-getChangeFromBaseLine.gg()
    if (is.character(obj.out)) {
      print(obj.out)
    } else {
      obj.out + 
        geom_point(aes(x = Sequence.Number, 
                       y = Change.From.Baseline.Numeric.Result, 
                       colour = factor(Unique.Subject.Identifier))) + 
        t3.line() +
        theme(legend.position = "none") + 
        coord_cartesian(xlim = t3scat.range$x, ylim = t3scat.range$y)+
        ylab("Change From Baseline") +
        xlab("Sequence Number")
    }
  })
  # Scatter Rchart Plo, uses the NVD3.js library
  output$scatcout <- renderChart({
    # If scat is checked in sidePanel (inputId=cfbPlotType) & call to getChangeFromBaseLine.gg()
    # does not return a character object, then returns a scatter plot with:
    #   1. x = Sequence.Number 
    #   2. y = Change.From.Baseline.Numeric.Result
    #   3. colour = by Unique.Subject.Identifier
    #   4. no legend
    #   5. zoom coordinates
    if (!(any(input$cfbPlotType %in% "scat"))) {
      return()
    } else {
      dftemp <- t3.curdata()
      names(dftemp) <- gsub("\\.", "", names(dftemp))
      chartout<-nPlot(ChangeFromBaselineNumericResult~SequenceNumber, data = dftemp, group = "UniqueSubjectIdentifier", type = "scatterChart")
      chartout$chart(tooltipContent = "#!function(key, x, y, e, graph){ return key + ' | ' + e.point.VisitNumber} !#",
                     showLegend = FALSE
                     )
      chartout$xAxis(axisLabel = "Sequence Number")
      chartout$yAxis(axisLabel = "Change from Baseline")
      chartout$set(legendPosition = "none")
      chartout$addParams(dom = "scatcout")
      return(chartout)
    }
    })
  
  # Line Plot: a reactive function that feeds into scatter ggplot when checkbox is true
  t3.line <- reactive({
    # If line is checked in sidePanel (inputId=cfbPlotType) & call to getChangeFromBaseLine.gg()
    # does not return a character object, then returns a line plot with:
    #   1. x = Sequence.Number 
    #   2. y = Change.From.Baseline.Numeric.Result
    #   3. colour = by Unique.Subject.Identifier
    #   4. no legend
    #   5. zoom coordinates
    if (!(any(input$cfbPlotType %in% "line"))) {
      return()
      }
    geom_line(aes(x = Sequence.Number, 
                  y = Change.From.Baseline.Numeric.Result, 
                  colour = factor(Unique.Subject.Identifier)))
  })
  
  # Boxplot
  output$box.out <- renderPlot({
    # If box is checked in sidePanel (inputId=cfbPlotType) & call to getChangeFromBaseLine.gg()
    # does not return a character object, then returns a box plot with:
    #   1. x = Lab.Test.or.Examination.Name 
    #   2. y = Change.From.Baseline.Numeric.Result
    #   3. colour/each boxplot corresponds to individual patients (Unique.Subject.Identifier)
    #   4. no legend
    #   5. zoom coordinates
    
    if (!(any(input$cfbPlotType %in% "box"))) {
      return()
      }
    obj.out <- getChangeFromBaseLine.gg()
    if (is.character(obj.out)) {
      print(obj.out)
      } else {
        obj.out + 
          geom_boxplot(aes(x    = Lab.Test.or.Examination.Name, 
                           y    = Change.From.Baseline.Numeric.Result, 
                           fill = (Unique.Subject.Identifier))) + 
          theme(legend.position = "none") + 
          coord_cartesian(xlim = t3box.range$x, ylim = t3box.range$y)+
          ylab("Change from Baseline")+
          xlab("Lab Test or Examination Name by Subject")
        }
    })
  
  # Histogram
  output$hist.out<- renderPlot({
    # If hist is checked in sidePanel (inputId=cfbPlotType) & call to getChangeFromBaseLine.gg()
    # does not return a character object, then returns a histogram with:
    #   1. x = Change.From.Baseline.Numeric.Result 
    #   2. y = Counts (grouped by numeric result)
    #   3. colour = corresponds to individual patients (Unique.Subject.Identifier)
    #   4. no legend
    #   5. zoom coordinates
    if (!(any(input$cfbPlotType %in% "hist"))) {
      return()
      }
    obj.out <- getChangeFromBaseLine.gg()
    if (is.character(obj.out)) {
      print(obj.out)
      } else {
        obj.out + 
          geom_histogram(aes(x    = Change.From.Baseline.Numeric.Result, 
                             fill = Unique.Subject.Identifier)) + 
          theme(legend.position = "none") + 
          coord_cartesian(xlim = t3hist.range$x, ylim = t3hist.range$y)+
          ylab("Count")+
          xlab("Change From Baseline by Subject")
        }
    })
  
  
  #------------------------------------------Zoom Plot----------------------------------------#
  
  # Initialize ranges as reactive values
  t3scat.range <- reactiveValues(x = NULL, y = NULL)
  t3line.range <- reactiveValues(x = NULL, y = NULL)
  t3box.range  <- reactiveValues(x = NULL, y = NULL)
  t3hist.range <- reactiveValues(x = NULL, y = NULL)
  
  # Set ranges on double click
  
  # Scatter
  observeEvent(input$t3scat.dblclick, {
    brush <- input$t3scat.brush
    if (!is.null(brush)) {
      t3scat.range$x <- c(brush$xmin, brush$xmax) 
      t3scat.range$y <- c(brush$ymin, brush$ymax) 
      } else { 
        t3scat.range$x <- NULL 
        t3scat.range$y <- NULL 
        }
    }) 
  
  # Boxplot
  observeEvent(input$t3box.dblclick, { 
    brush <- input$t3box.brush
    if (!is.null(brush)) {
      t3box.range$x <- c(brush$xmin, brush$xmax)
      t3box.range$y <- c(brush$ymin, brush$ymax)
      } else {
        t3box.range$x <- NULL
        t3box.range$y <- NULL
        }
    })
  
  # Histogram
  observeEvent(input$t3hist.dblclick, {
    brush <- input$t3hist.brush
    if (!is.null(brush)) {
      t3hist.range$x <- c(brush$xmin, brush$xmax)
      t3hist.range$y <- c(brush$ymin, brush$ymax)
      } else {
        t3hist.range$x <- NULL
        t3hist.range$y <- NULL
        }
    })
    
  
  #-------------------------------------------------------------------------------------------------#
  # Show Tables
  #-------------------------------------------------------------------------------------------------#
  
  # Show brushed points in a table for scatterplot and lineplot as 
  # functionality does not work with histogram and boxplot. 
  # [In future, have to implement rCharts for a better user-interactive plots that includes histogram and boxplot]
  
  # Scatter
  output$t3.scatsubset <- renderDataTable({
    datatable(data = t3.getbrushset.scat(), 
              options = list(scrollX = TRUE)
              )
    })
  
  #------------------------------------------Outlier Detection--------------------------------------#
  
  # Outlier detection [To be implemented]
  output$outlier.test.type <- renderUI({
    # Return if no lab.test choose OR scatterplot is not selected
    if (!any(input$cfbPlotType %in% "scat") || is.null(input$t3.choose.lab.testname)) {
      return()
    }
    # [OTHER CODE HERE]
  })
}


#==================================================================================================# 
# Run App
#==================================================================================================#

shinyApp(ui = ui, server = server)

#---------------------------------------------------------------------------------------------------#
# Close Tracker
# sink()
# closeAllConnections()

#---------------------------------------------------------------------------------------------------#
# END SCRIPT
