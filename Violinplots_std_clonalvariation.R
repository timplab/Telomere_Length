---
title: "Violin plots, std cuttoff, and clonal variation"
output: hmtl_notebook 
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
```

#Load in telomere lengths file
```{r}
#start with telomere length csv output 

telomere_lengths <- read.csv("telomere_length_data_all.csv")

```

#Generate a violin plot of telomere length
```{r}
#chose the telomeres of interest (that have above N threhold)
subset <- subset(telomere_lengths, telomere_lengths$telomere %in% c("1L", "2L", "3L", "4L", "4R", "5L", "5R", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R"))

subset$telomere <- factor(subset$telomere, c("1L", "2L", "3L", "4L", "4R", "5L", "5R", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R"))

#add the total mean violin of all telomeres 
total <- subset %>%
  mutate(telomere = "all") %>%
  mutate(genotype = "all")

totalsubset <- full_join(subset, total)
totalsubset$telomere <- factor(totalsubset$telomere, c("all", "1L", "2L", "3L", "4L", "4R", "5L", "5R", "6R", "7R", "8L", "8R", "9R", "10L", "10R", "11L", "11R", "13L", "14L", "15L", "15R"))

#generate a violin plot with ggplot2
violin <- ggplot(data=NULL) +
  geom_violin(data = totalsubset, aes(x=telomere, y=length, fill=genotype, width=2), position=position_dodge(0.75), draw_quantiles = c(0.5), color='black', alpha=0.6)+
  scale_fill_manual(values=c("mediumpurple4", "deepskyblue2", "aquamarine3")) +
  ggtitle("WT Single Telomere Length") + 
  geom_hline(yintercept = c("overall mean telomere length here"), linetype = "dotted", size=1) +
  labs(y="Telomere length (bp)", x="Chromosome End") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ylim ("lower y limit of plot", "upper y limit of plot") +
  scale_fill_brewer(palette="Oranges") 

plot(violin)
```

#Save the violin plot of telomere length 
```{r}
ggsave("violin_plot", plot = violin, device = "pdf", scale = 1, width = 25, height = 6, units = c("in", "cm", "mm"))
```

#Add a standard deviation threhold to the violin plot 
```{r}
telomere_lengthsmean <- aggregate(telomere_lengths[, 1], list(telomere_lengths$telomere), mean) %>%
  mutate(mean=length) %>%
  mutate(telomere=Group.1) %>%
  select(telomere, mean) %>%
  mutate(doubling="35")
telomere_lengthsstd <- aggregate(telomere_lengths[, 1], list(telomere_lengths$telomere), sd) %>%
  mutate(std=length) %>%
  mutate(telomere=Group.1) %>%
  select(telomere, std) %>%
  mutate(doubling="35")
telomere_lengths_std_mean <- full_join(telomere_lengthsmean, telomere_lengthsstd) %>%
  mutate(threestd = (telomere_lengthsmean$mean + (telomere_lengthsstd$std*3))) %>%
  select(telomere, threestd, doubling)


telomere_lengths_std_subset <- subset(telomere_lengths_std, telomere_lengths_std$telomere %in% c("1L", "1R","2L", "2R", "3L", "3R", "4L", "5L", "5R", "6R", "7L", "7R", "8L", "8R", "9L", "9R", "10R", "11L", "11R", "13L", "13R", "14L", "14R", "15L", "15R", "16L"))

telomere_lengths_std_subset$telomere <- factor(telomere_lengths_std_subset$telomere, c("1L", "1R","2L", "2R", "3L", "3R", "4L", "5L", "5R", "6R", "7L", "7R", "8L", "8R", "9L", "9R", "10R", "11L", "11R", "13L", "13R", "14L", "14R", "15L", "15R", "16L"))

subset <- subset(telomere_lengths, telomere_lengths$telomere %in% c("1L", "1R","2L", "2R", "3L", "3R", "4L", "5L", "5R", "6R", "7L", "7R", "8L", "8R", "9L", "9R", "10R", "11L", "11R", "13L", "13R", "14L", "14R", "15L", "15R", "16L"))

subset$telomere <- factor(subset$telomere, c("1L", "1R","2L", "2R", "3L", "3R", "4L", "5L", "5R", "6R", "7L", "7R", "8L", "8R", "9L", "9R", "10R", "11L", "11R", "13L", "13R", "14L", "14R", "15L", "15R", "16L"))

violin_with_std <- ggplot(data=NULL) +
  geom_violin(data = subset, aes(x=telomere, y=length, fill=doubling, width=3), position=position_dodge(0.8), draw_quantiles = c(0.5), color='black', alpha=0.6)+
  scale_fill_manual(values=c("mediumpurple4", "deepskyblue2", "aquamarine3")) +
  geom_point(data=telomere_lengths_std_subset, aes(x=telomere, y=onestd, color=doubling), position=position_dodge(0.8),shape=23, size=3) +
  scale_color_manual(values=c("darkviolet", "darkviolet", "darkviolet")) +
  #ggtitle("Violin with standard deviation") + 
  labs(y="Telomere length (bp)", x="") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
 ylim ("lower y limit plot", "upper y limit plot") +
  scale_fill_brewer(palette="RdYlBu")

plot(violin_with_std)
```

#Save the violin plot of telomere length with 3 std 
```{r}
ggsave("violin_plot_std", plot = violin_with_std, device = "pdf", scale = 1, width = 25, height = 6, units = c("in", "cm", "mm"))
```


#Plot examples of clonal variation with statistics 
```{r}
#generate a subset with the telomeres you are interested in examining for clonal variation
subset_clonal <- subset(telomere_lengths, telomere_lengths$telomere %in% c("3L", "8R", "10R"))

library(rstatix)
library(ggpubr)

#Run a t-test to compare the clones at each telomere 
stat.test <- subsetclonal %>%
  group_by(telomere) %>%
  t_test(length ~ genotype) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

stat.test <- stat.test %>% add_xy_position(x = "telomere")

#Generate violin plots with the p-values across clones for each telomere 
violinclonal <- ggplot(data=NULL) +
  geom_violin(data = subset_clonal, aes(x=telomere, y=length, fill=genotype, width=2), position=position_dodge(0.8), draw_quantiles = c(0.25, 0.5, 0.75), color='black', alpha=0.6)+
  geom_jitter(data= subsetclonal, aes(x=telomere, y=length, fill=genotype), position=position_jitterdodge(0.1), size=0.5, alpha=0.4) +
 stat_pvalue_manual(stat.test, label = "p.adj.signif") +
 scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y="Telomere length (bp)", x="Chromosome End") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Oranges")  

plot(violinclonal)
```

#Save the clonal variation violin plot 
```{r}
ggsave("violin_clonal_variation", plot = violinclonal, device = "pdf", scale = 1, width = 14, height = 7, units = c("in", "cm", "mm"))
```

