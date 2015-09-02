#--------------------------95% CI Interval Reference Sheet, Binomial Response----------------------#
# install.packages(c("shiny", "htmlwidgets"))
library(binom)
library(shiny)
delta<-0.1
level<-0.95

#--------------------------------------------------------------------------------------------------#
# UI.R
ui <- fluidPage(
  h3("95% CI Interval Reference Sheet"),
  fluidRow(
    column(5,
           numericInput("N", label=h4("N Input"), value=20),
           sliderInput("int.null", label = h4("Null Hypothesis"),
              min = 0, max = 1, value = 0.3)
    ),
    column(7,
           selectInput("int.type", label = h4("Interval Type"), 
              choices = list("exact"="exact", "ac"="ac", "asymptotic"="asymptotic", "wilson"="wilson",
                             "prop.test"="prop.test","bayes"="bayes", "logit"="logit", "cloglog"="cloglog",
                             "probit"="probit"), selected = "exact"),
           br(),
           textInput("xlabel", "X Label", value ="Observed count of responders"),
           br(),
           textInput("ylabel", "Y Label", value ="Estimated proportion of responders"))
    ),
  plotOutput("plot"),
  downloadButton('downloadplot', "Download as PDF")
)

#--------------------------------------------------------------------------------------------------#
# server.R
server <- function(input, output) {
  binom.plot <- function(){
      N<-input$N
      plotdat<-binom.confint(0:N,N,conf.level=level,method=input$int.type)
      # Plot means
      plot(plotdat$mean,plotdat$x,
         ylab=input$xlabel,
         xlab=input$ylabel,
         pch=16)
      # Find Minimal Response Needed to Reject Null
      min.numres<- min(plotdat$x[plotdat$lower>input$int.null])
      if(min.numres>N){
        min.numres<- "N/A"
      }
      # Title
      title(paste(input$int.type,level*100,"% CIs for study size N =",N,sep=" "))
      # Plot Grid
      abline(v=c(0:N)/N,lty=3,lwd=1,col=8)
      abline(h=c(0:N),lty=3,lwd=1,col=8)
      # Plot CI
      for(i in 1:length(plotdat$x)){
        segments(plotdat$lower[i],plotdat$x[i],plotdat$upper[i],plotdat$x[i],lwd=1.3)
        segments(plotdat$lower[i],plotdat$x[i]-delta,plotdat$lower[i],plotdat$x[i]+delta,lwd=1.3)
        segments(plotdat$upper[i],plotdat$x[i]-delta,plotdat$upper[i],plotdat$x[i]+delta,lwd=1.3)
      }
      # Vertical Null Line
      abline(v=input$int.null, col="red")
      # Labels of Null and Minimal Response to Reject Null
      text(max(plotdat$x)/10, label = paste("Null = ", input$int.null, "\nMinimal Number of Response to Reject Null = ", 
                                            min.numres,sep=""), pos=2)
  }
  # Render Plot
  output$plot <- renderPlot({binom.plot()})
  # Download Plot
  output$downloadplot <- downloadHandler(
    filename <- paste("myplot", Sys.Date(), ".pdf", sep=''),
    content <- function(file) {
      pdf("plot.pdf", width=12, height=8)
      binom.plot()
      dev.off()
      file.copy("plot.pdf", file)
    }
  )
}

#--------------------------------------------------------------------------------------------------#
# Run App
shinyApp(ui, server)