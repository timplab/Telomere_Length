---
title: "Telomerase null shortening rate line plots"
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
library(readxl)
```

#Load excel file with the telomere length means calculated in excel
```{r}
#start with telomere length csv output 

Passages0123 <- read_excel("passages_0123_means.xlsx")

```

#Aggregate the means and separate out into Passages 0 to 1 and 1 to 3 
```{r}

lengthmeans <- aggregate(Passages0123[, 3], list(Passages0123$Doubling), mean)

WTP1means <- subset(lengthmeans, lengthmeans$Group.1 %in% c("1", "35")) %>%
   mutate(slope = round(lm(length ~ Group.1)$coefficients[2], 2))
P123means <- subset(lengthmeans, lengthmeans$Group.1 %in% c("35", "55", "75")) %>%
  mutate(slope = round(lm(length ~ Group.1)$coefficients[2], 2))

```

#Plot the slopes for telomeres 1 - 8 
```{r}

fP123subset <- subset(P123means, P123means$TEL %in% c("1L","1R","2L", "3L", "3R", "4L", "5L", "5R", "6R", "7L", "7R", "8L", "8R")) %>%
  group_by(TEL)%>%
  mutate(slope = round(lm(length ~ Doubling)$coefficients[2], 2))

slopes18 <- fP123est2subset %>%
  select(TEL, slope)%>%
  unique()

plotf <- ggplot(data=NULL)+
  geom_point(data=fP123est2subset, aes(Doubling, length, color=TEL, group=TEL)) +
  geom_line(data=fP123est2subset, aes(Doubling, length, color=TEL, group=TEL)) +
  geom_point(data=WTP1means, aes(Group.1, length), size =4) +
  geom_smooth(data=WTP1means, aes(Group.1, length), method='lm', se=FALSE, linetype = "dashed", color="black")+
  geom_point(data=P123means, aes(Group.1, length), size =4) +
  geom_smooth(data=P123means, aes(Group.1, length), method='lm', se=FALSE, linetype = "dashed", color="black")+
  ggtitle("P0 to P3 Chr 1-8") + 
  labs(y="Telomere length (bp)", x="Number of Population Doublings") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ylim (130, 520) +
  scale_x_continuous(breaks=c(1, 35, 55, 75), labels=(c("0", "35", "55", "75")))

plot(plotf)

```

#Save the plot of slopes for telomeres 1 - 8 
```{r}
ggsave("slopes_1-8", plot = plotf, device = "pdf", scale = 1, width = 8, height = 4, units = c("in", "cm", "mm"))
```

#Plot the slopes for telomeres 9 - 16 
```{r}

sP123subset <- subset(P123est2, P123est2$TEL %in% c("9L", "9R", "10L", "10R", "11L", "13L", "13R", "14L", "15L", "15R", "16L")) %>%
  group_by(TEL)%>%
  mutate(slope = round(lm(length ~ Doubling)$coefficients[2], 2))

sP123est2subset$TEL <- factor(sP123est2subset$TEL, c("9L", "9R", "10L", "10R", "11L", "13L", "13R", "14L", "15L", "15R", "16L"))

slopes916 <- sP123est2subset %>%
  select(TEL, slope)%>%
  unique()

slopes <- full_join(slopes18, slopes916)

plots <- ggplot(data=NULL)+
  geom_point(data=sP123subset, aes(Doubling, length, color=TEL, group=TEL)) +
  geom_line(data=ssP123subset, aes(Doubling, length, color=TEL, group=TEL)) +
   geom_point(data=WTP1means, aes(Group.1, length), size =4) +
  geom_smooth(data=WTP1means, aes(Group.1, length), method='lm', se=FALSE, linetype = "dashed", color="black")+
  geom_point(data=P123means, aes(Group.1, length), size =4) +
  geom_smooth(data=P123means, aes(Group.1, length), method='lm', se=FALSE, linetype = "dashed", color="black")+
  ggtitle("P0 to P3 Chr 9-16") + 
  labs(y="Telomere length (bp)", x="Number of Population Doublings") +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ylim (130, 520) +
  scale_x_continuous(breaks=c(1, 35, 55, 75), labels=(c("0", "35", "55", "75")))

plot(plots)


```

#Save the plot of slopes for telomeres 9 - 16
```{r}
ggsave("slopes_9-16", plot = plots, device = "pdf", scale = 1, width = 8, height = 4, units = c("in", "cm", "mm"))
```

#Save slopes 
```{r}
write_csv(slopes, "P0123_slopes.csv")
```

