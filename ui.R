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
  titlePanel(title=h2("Clonal Detection Power Analysis by scRNA-seq",align='center')),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("Please slide to adjust parameters"),
      p(""),
      sliderInput("subclone_freq",
                  "Subclone cell percentage (SCP):",
                  min = 0,
                  max = 100,
                  value = 30,step=10),
      sliderInput("subclone_SNV_freq",
                  "Subclone specific mutation percentage (SMP):",
                  min = 0,
                  max = 100,
                  value = 25,step=10),
      sliderInput("fp",
                  "Percentage of mutations from ancestral branch (false positives,fp)",
                  min = 0,
                  max = 80,
                  value = 0,step=20),
      sliderInput("n",
                  "Total number of mutations",
                  min = 0,
                  max = 10000,
                  value = 200,step=50),
      sliderInput("m",
                  "Total number of Cells",
                  min = 50,
                  max = 150,
                  value = 100,step=50),
      sliderInput("r",
                  "Read depth (Rd)",
                  min = 10,
                  max = 100,
                  value = 20,step=10)
    ),
    # Show a plot of the generated distribution
    mainPanel(
       h4(p("This app is mean to guide researcher on experimental design using Smart-seq2
         or Fluidigm protocal. We assume we have a three-clades tree. The statistic is 
         calculated by ", strong("DENDRO. "),  "All the stored simulation result is based on ", 
         span("Deng et al. Science 343.6167 (2014): 193-196. ",style="color:blue"),
         'In tumor cases, the estimation may be a little biased. We always recommend ',
         strong('more cell and higher depth if available.'))),
       plotOutput("distPlot",width = "100%"),
       br(),
       br(),
       br(),
       h4('Statistic definition:'),
       h5(p(strong("Detection probability:")," Probability that you detect this subclone given 100 permutation."),
          p(strong("Sensitivity:")," Sensitivity of the variable green cluster, i.e. out of all the cells in the 
            green cluster, how many of them is detected"),
          p(strong("Specificiy:")," Specificiy of the variable green cluster, i.e. out of the \'green cluster\' 
            you detected, how many are actually from the true green cluster"),
          p(strong("Adjusted Random Index:")," Adjusted rand index is a measure of the similarity between two 
            data clusterings after adjusted for the chance grouping of elements. For detali,", 
            span("https://www.wikiwand.com/en/Rand_index#/Adjusted_Rand_index",style='color:blue')))
    )
  )
))
