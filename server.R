#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(dendextend)
dend=as.dendrogram(hclust(dist(rbind(c(1,0,0),c(0,1,0),c(0,1,1)),method='binary')))
labels(dend)=c('a','b','c')
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  load('Lratiodist_eval_fp.rda')
  kprob3=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  lprob3=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  fprate_v=c(0,0.2,0.4,0.6,0.8)
  n_v=c(100,1000,2000,5000,10000)
  output$distPlot <- renderPlot({
    # generate kprob and lprob based on input from ui.R
    f = which.min(abs(fprate_v-input$fp/100))
    fp=round(fprate_v[f],2)
    k = which.min(abs(kprob3-input$subclone_freq/100))
    # kprob = round((1-fp)*kprob3[k],2)
    kprob = round(kprob3[k],2)
    
    l = which.min(abs(lprob3-input$subclone_SNV_freq/100))
    # lprob = round((1-fp)*lprob3[l],2)
    lprob = round(lprob3[l],2)
    
    # kprobi=round(((1-fp)-kprob)/2,2)
    kprobi=round((1-kprob)/2,2)
    # lprobj=round(((1-fp)-lprob)/3,2)
    lprobj=round((1-lprob)/3,2)
    
    n = which.min(abs(n_v-input$n))
    
    ari = round(eval_final[k,l,f,n,1],2)
    # draw the dendrogram with the specified kprob and lprob
    dend %>% set("branches_k_color", value = c("red", "blue",'green'), k = 3) %>%
      set("leaves_pch", c(17, 18, 19)) %>%  # node point type
      set("leaves_cex", 2) %>%  # node point size
      set("leaves_col", c("blue", "red", "green")) %>%
      plot(main = "Customized colors",ylim=c(-0.5,1.1))
    text(1.1,0.5,lprobj)
    text(2.1,0.25,lprobj)
    text(3.1,0.25,lprob)
    text(2.6,0.75,lprobj)
    text(1.7,1.1,fp)
    text(1,-0.1,kprobi)
    text(2,-0.1,kprobi)
    text(3,-0.1,kprob)
    text(2.5,-0.25, paste0('Total number of mutations are ',n))
    text(2.5,-0.5, paste0('Adjusted Random Index: ',ari))
    
  })
  
})

#  ------------------------------------------------------------------------

