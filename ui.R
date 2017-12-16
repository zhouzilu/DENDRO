#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel(title=h4("clonal detection power analysis by scRNA-seq",align='center')),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("this is side bar panel"),
      sliderInput("subclone_freq",
                  "Subclone cell percentage:",
                  min = 0,
                  max = 100,
                  value = 30,step=10),
      sliderInput("subclone_SNV_freq",
                  "Subclone specific mutation percentage:",
                  min = 0,
                  max = 100,
                  value = 25,step=10),
      sliderInput("fp",
                  "Percentage of mutations from ancestral branch (false positives)",
                  min = 0,
                  max = 100,
                  value = 0,step=20),
      sliderInput("n",
                  "Total number of mutations",
                  min = 0,
                  max = 10000,
                  value = 1000,step=1000)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot")
    )
  )
))
