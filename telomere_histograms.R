---
title: "Figure 2 histograms"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 

```{r}
plot(cars)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

# Load preferences
```{r setup, eval=TRUE, include=FALSE, cache=F, message=F, warning=F, results="hide"}
rm(list=ls());gc()
#knitr::opts_chunk$set(fig.path='figs/')
knitr::opts_chunk$set(cache = FALSE, warning = FALSE,
                      message = FALSE, cache.lazy = FALSE)
my_plot_hook <- function(x, options)
  paste("\n", knitr::hook_plot_tex(x, options), "\n")
knitr::knit_hooks$set(plot = my_plot_hook)
```

# Load libraries 
```{r functions, include=F}
library(tidyverse)
library(ggsci)
library(zoo)
library(fGarch)
```

# To generate a histogram plot of the telomere length for two genotypes: 
# Load telomere length data for each genotype (from the telomere start position) with column names length, telomere, and genotype
```{r}
genotype1 <- read.csv("genotype1_telo_lengths.csv")
genotype2 <- read.csv("genotype2_telo_lengths.csv")

```

#Plot histogram comparing the two genotypes 
```{r}
# Calculate the mean, med, std for each genotype to generate a density curve using dsnorm
mean <- mean(genotype1$length)
med<- median(genotype1$length)
std <- sd(genotype1$length)
xcurve = c(genotype1$length)
ycurve <- dsnorm(xcurve, mean = mean, sd = std)

mean2 <- mean(genotype2$length)
med2<- median(genotype2$length)
std2 <- sd(genotype2$length)
xcurve2 = c(genotype2$length)
ycurve2 <- dsnorm(xcurve2, mean = mean2, sd = std2)

#Plot a histogram of the data with the density curve for each genotype using the ggplot2 histogram function
histogram <- ggplot(data=NULL) + 
  geom_histogram(data=genotype1, aes(x= length, fill=length),binwidth = 10, fill = "color of choice", alpha=0.4) + 
  geom_histogram(data=genotype2, aes(x= length, fill=length),binwidth = 10, fill = "mediumpurple4", alpha=0.4) + 
  geom_vline(xintercept = c('number_lower_bound', 'number_upper_bound'), linetype = "color of choice") +
  geom_line(aes(x= xcurve, y= ycurve, col="color of choice"), size=1.5) + 
  geom_line(aes(x= xcurve2, y= ycurve2, col="color of choice"), size=1.5) + 
  ggtitle("Title") + 
  labs(y="Axis title", x="Telomere length (bp)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
   xlim("expected lower bound", "expected upper bound") +
  theme(text = element_text(size = 20))

print(histogram)

```

#Save the histogram 
```{r}
ggsave("Name", plot = histogram, device = "pdf", scale = 1, width = 8, height = 4, units = c("in", "cm", "mm"))
```

# To generate a histogram plot of the telomere restriction fragment length for two genotypes: 
# Load telomere length data for each genotype (from restriction enzyme cut site for comparison to Southern blot)
# Load densitometry data 
```{r}
#load in telomere length data 
genotype1_Yprime <- read.csv("genotype1_telo_lengths_Yprime.csv")
genotype2_Yprime <- read.csv("genotype2_telo_lengths_Yprime.csv")

#load in densitometry data
genotype1_Southern_1 <- read.csv("southern_density_genotype1_1.csv")
genotype1_Southern_2 <- read.csv("southern_density_genotype1_2.csv")
genotype2_Southern_1 <- read.csv("southern_density_genotype2_1.csv")
genotype2_Southern_2 <- read.csv("southern_density_genotype2_2.csv")

```

```{r}

#annotate columns for telomere length data 
genotype1_Yprime <- genotype1_Yprime %>%
  mutate(qname = "1") %>%
  select(length_telo, qname)

genotype2_Yprime <- genotype2_Yprime %>%
  mutate(qname = "2") %>%
  mutate(length_telo = length_telo2) %>%
  select(length_telo2, qname)

# Calculate the mean, med, std for each genotype to generate a density curve using dsnorm
mean <- mean(genotype1_Yprime$length_telo)
med<- median(wgenotype1_Yprime$length_telo)
std <- sd(genotype1_Yprime$length_telo)
xcurve3 = c(genotype1_Yprime$length_telo)
ycurve3 <- dsnorm(xcurve3, mean = mean, sd = std)

mean2 <- mean(genotype2_Yprime$length_telo)
med2<- median(genotype2_Yprime$length_telo)
std2 <- sd(genotype2_Yprime$length_telo)
xcurve4 = c(genotype2_Yprime$length_telo)
ycurve4 <- dsnorm(xcurve4, mean = mean2, sd = std2)

#Plot the overlay of the two genotypes, the density curves for the nanopore data, and the density cruves from the Southern blot densitometry  
histogram_restriction <- ggplot(data=NULL) + 
  geom_histogram(data=wgenotype1_Yprime, aes(x= length_telo, fill=length_telo*0.3),binwidth = 10, fill = "darkorange2", alpha=0.4) + 
  geom_histogram(data=genotype2_Yprime, aes(x= length_telo2, fill=length),binwidth = 10, fill = "mediumpurple4", alpha=0.4) + 
  geom_line(aes(x= xcurve3, y= ycurve3, col="Nanopore WT"), size=1, alpha = 0.3) +
  geom_line(aes(x= xcurve4, y= ycurve4, col="Nanopore rif1"), size=1, alpha = 0.3) +
  geom_vline(xintercept = c("lower_intercept", "upper_intercept"), linetype = "dotted") +
  scale_color_manual(values = c("darkorange3", "mediumpurple4")) + 
  geom_line(data=genotype1_Southern_1, aes(x=length_telo, y=qname, color = "name1"), size =1, linetype="longdash")+
  geom_line(data=genotype1_Southern_2, aes(x=length_telo, y=qname, color = "name2"), size=1, linetype="longdash")+
  geom_line(data=genotype2_Southern_1, aes(x=length_telo, y=qname, color = "name3"), size=1, linetype="longdash")+
  geom_line(data=genotype2_Southern_2, aes(x=length_telo, y=qname, color = "name4"), size=1, linetype="longdash")+
  ggtitle("Y' Telomere Fragment Length") + 
  labs(x="XhoI telomere fragment length (bp)") +
  scale_color_manual(values = c("mediumpurple3", "darkorange2", "navy", "red3", "navy")) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
   theme(legend.title = element_blank()) +
  ylim(400, 2500) +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(labels = waiver(), breaks=waiver(), sec.axis = dup_axis(labels = waiver(), breaks=waiver()))

print(histogram_restriction)
```
#Save 
```{r}
ggsave("Yprime_histogram", plot = histogram_restriction, device = "pdf", scale = 1, width = 10, height = 6, units = c("in", "cm", "mm"))
```

