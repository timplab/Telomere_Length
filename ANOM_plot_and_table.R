---
title: "ANOM for telomere lengths"
output: html_notebook
author: Sam Sholes and Thomas Kelly
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
library(ggExtra)
library(ANOM)
library(multcomp)
```

#Load in telomere lengths file
```{r}
#start with telomere length csv output 

telomere_lengths <- read.csv("telomere_length_data.csv")

```

#Output an ANOM plot and statistics summry table with a Bonferroni correction
```{r}
#Put the telomeres in desired order 
telomere_lengths$telomere <- factor(telomere_lengths$telomere, c("2L", "3L", "4L", "4R", "5L", "5R", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R", "16L"))

#Generate the ANOM plot and table
linmod <- lm(length ~ telomere, data = telomere_lengths)
glm <- glht(linmod, mcp(telomere = "GrandMean"))
print(summary(glm, test = adjusted("bonferroni")))
ANOM(glm, printp = FALSE, bg="white", bgrid = FALSE)
```

