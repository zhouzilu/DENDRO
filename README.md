# clonal-info-scRNAseq
This an R shiny app for researchers to estimate of power detection given clonal structure, mutation burden and subclone mutation estimate.

To RUN this app, we need the following package pre-installed:
```r
install.packages('shiny')
install.packages('dendextend')
install.packages('tidyr')
install.packages("viridisLite") # Issue with R 3.4.3
install.packages("viridis") # Issue with R 3.4.3
```
If all the library are pre-installed, open your R and use code:
```r 
library(shiny)
runGitHub("DENDRO","zhouzilu")
```
