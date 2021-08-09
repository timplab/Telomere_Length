---
title: "Correlation plots"
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
library(ggExtra)
```

#Load in telomere lengths file
```{r}
#start with telomere length csv output 

telomere_lengths_sample1 <- read.csv("telomere_length_data_sample1.csv")
telomere_lengths_sample2 <- read.csv("telomere_length_data_sample2.csv")

```

#Generate correlation plots between each sameple
```{r}


subset_sample1 <- subset(telomere_lengths_sample1, telomere_lengths_sample1$telomere %in% c("2L", "3L", "4L", "4R", "5L", "5R", "6L", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R", "16L"))
subset_sample1 <- subset_sample1 %>%
  rename(Sample1_length = length)
sample1_mean <-aggregate(subset_sample1[, 1], list(subset_sample1$telomere), mean)

subset_sample2 <- subset(telomere_lengths_sample2, telomere_lengths_sample2$telomere %in% c("2L", "3L", "4L", "4R", "5L", "5R", "6L", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R", "16L"))
subset_sample2 <- subset_sample2 %>%
  rename(Sample2_length = length)
sample2_mean <-aggregate(subset_sample2[, 1], list(subset_sample2$telomere), mean)

allmeans <- full_join(sample1_mean, sample2_mean)


###Plot the means against each other 

scatter_sample1and2 <- ggplot(data = NULL) +
  geom_point(data = allmeans, aes(x=Sample1_length, y=Sample2_length, color= Group.1), size =3) +
  geom_smooth(data = allmeans, aes(x=Sample1_length, y=Sample2_length), method = 'lm', color = "black", linetype ="dotted", se=FALSE) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y="Sample 1 mean telomere length", x="Sample 2 mean telomere length") +
   theme(legend.title = element_blank()) +
  xlim(250, 500) +
  ylim(250, 500)+
theme(legend.position = "none") 
scatter_sample1and2mar <- ggMarginal(scatter_sample1and2, type="density") 
plot(scatter_sample1and2mar)

###Calculate the correlation coefficient between the datasets 
cor(allmeans$Sample1_length, allmeans$Sample2_length)


```

#Save the correlation plot 
```{r}
ggsave("correlation plot", plot = sscatter_sample1and2mar, device = "pdf", scale = 1, width = 8, height = 6, units = c("in", "cm", "mm"))
