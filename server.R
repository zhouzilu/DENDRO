#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
dend=as.dendrogram(hclust(dist(rbind(c(1,0,0),c(0,1,0),c(0,1,1)),method='binary')))
# labels(dend)=c('a','b','c')
dend=dendextend::`labels<-`(dend,c('a','b','c'))
`%>%`=tidyr::`%>%`
set=dendextend::set
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  # load('C:/Users/zhouzilu/Desktop/research/single cell/SNV simulation/hpc_res/appSNV_Lratiodist_eval_all.rda')
  load('appSNV_Lratiodist_eval_all.rda')
  kprob3=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98)
  lprob3=c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98)
  fprate_v=c(0,0.2,0.4,0.6,0.8)
  n_v=c(50,100,150,200,250,300,350,400,500,1000,1500,2000,5000,10000)
  rd_v=c(10,20,30,40,50,60,70,80,90,100)
  m_v=c(50,100,130)
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
    
    rd = which.min(abs(rd_v-input$r))
    
    m = which.min(abs(m_v-input$m))
    
    ari = round(eval_final[k,l,f,n,rd,m,1],2)
    sen = round(eval_final[k,l,f,n,rd,m,2],2)
    spe = round(eval_final[k,l,f,n,rd,m,3],2)
    det = round(eval_final[k,l,f,n,rd,m,4],2)
    # draw the dendrogram with the specified kprob and lprob
    par(mar=c(7.1,0.1,0.1,1.1),cex=1.5)
    dend %>% set("branches_k_color", value = c("red", "blue",'green'), k = 3) %>%
      set("leaves_pch", c(17, 18, 19)) %>%  # node point type
      set("leaves_cex", 2) %>%  # node point size
      set("leaves_col", c("blue", "red", "green")) %>%
      plot(main = NULL,yaxt='n',ylim=c(-1,1.1),cex=1.5)
    par(cex=1)
    text(1.1,0.5,lprobj,cex=1.5)
    text(2.1,0.25,lprobj,cex=1.5)
    text(3.1,0.25,paste0('SMP\n',lprob),cex=1.5)
    text(2.6,0.75,lprobj,cex=1.5)
    text(1.7,1.1,paste0('fp ',fp),cex=1.5)
    text(3.1,1.1,paste0('Rd ',rd_v[rd]),cex=1.5)
    text(1,-0.25,kprobi,cex=1.5)
    text(2,-0.25,kprobi,cex=1.5)
    text(3,-0.20,paste0('SCP\n',kprob),cex=1.5)
    text(2,-0.45, paste0('Total number of mutations are ',n_v[n]),cex=1.5)
    text(2,-0.60, paste0('Total number of cells are ',m_v[m]),cex=1.5)
    text(2,-0.75, paste0('Detection probability: ',det),cex=1.5)
    text(2,-0.90, paste0('Sensitivity: ',sen,'    Specificity: ',spe),cex=1.5)
    text(2,-1.05, paste0('Adjusted Rand Index: ',ari),cex=1.5)
  },height=600)
})

#  ------------------------------------------------------------------------

