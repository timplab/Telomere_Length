---
title: "QQ Plots"
output: html_notebook
author: Sam Sholes
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 

```{r}
plot(cars)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

# Load in preferences 
```{r setup, eval=TRUE, include=FALSE, cache=F, message=F, warning=F, results="hide"}
rm(list=ls());gc()
#knitr::opts_chunk$set(fig.path='figs/')
knitr::opts_chunk$set(cache = FALSE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
my_plot_hook <- function(x, options)
  paste("\n", knitr::hook_plot_tex(x, options), "\n")
knitr::knit_hooks$set(plot = my_plot_hook)
```

#Load in libraries 
```{r functions, include=F}
library(tidyverse)
library(ggsci)
library(zoo)
library(ggpmisc)
library(car)
```

#Load in telomere lengths file
```{r}
#start with telomere length csv output 

telomere_lengths_sample1 <- read.csv("telomere_length_data_sample1.csv")
southern_densitometry_sample1 <- read.csv("southern_densitometry_sample1.csv")

```

#Generate Q-Q plots:
```{r}
##Generate a normal QQ plot of Nanopore Sample 1 data
Sample1_Nanopore_qq <- qqnorm(telomere_lengths_sample1$length)
qqline(telomere_lengths_sample1$length, col = "red", lwd = 2)

##Generate a normal QQ plot of Southern Sample 1 data
Sample1_Southern_qq <- qqnorm(southern_densitometry_sample1$length)
qqline(southern_densitometry_sample1$length, col = "green", lwd = 2)

##Generate a QQ plot of Nanopore vs. Southern data 
NanoSouthqq <- qqplot(telomere_lengths_sample1$length, southern_densitometry_sample1 , xlab = "Sample1 Nanopore", ylab = "Sample1 Southern Density", main="Q-Q Plot of Nanopore vs Southern Density Distribution")
abline(a = 0, b = 1, col = "green", lwd = 2)

```

