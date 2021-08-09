---
title: "WALTER output box plots"
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

#Load in libraries 
```{r functions, include=F}
library(tidyverse)
library(ggsci)
library(zoo)
library(ggpmisc)
library(car)
library(BSDA)
```

#Load in output summary statistics from WALTER (see text references)
```{r}
#start with WALTER summary stat output and the Nanopore telomere length output

WALTER <- read.csv("WALTER_summary_stats.csv")
telomere_lengths <- read.csv("telomere_length_data_sample1.csv")

```

#Generate a boxplot 
```{r}

WALTERplot <- ggplot(WALTER, aes(x = as.factor(group2), color = group2)) +
  geom_boxplot(aes(
      lower = one, 
      upper = three, 
      middle = median, 
      ymin = min, 
      ymax = max),
    stat = "identity") +
  geom_boxplot(data=telomere_lengths, aes(x="Southern 1", y=length), draw_quantiles = c(0.25, 0.5, 0.75), color='darkorange3') +
  ggtitle("Telomere Length Comparison") + 
  labs(y="Telomere length (bp)", x="Telomere") +
  theme_classic() +
  scale_color_manual(values = c("darkorange3", "mediumpurple3")) +
  theme(legend.position = "none") +
   theme(text = element_text(size = 20))
  
plot(WALTERplot)
```
#Save 
```{r}
ggsave("WALTER_compare", plot = WALTERplot, device = "pdf", scale = 1, width = 10, height = 4, units = c("in", "cm", "mm"))
```
