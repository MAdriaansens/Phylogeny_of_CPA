---
title: "Untitled"
author: "Mick Adriaansens"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
completeA <- read.delim('/nesi/project/uc04105/CPA_hit_ArchaeaGTDB_24032025.tsv', sep = '\t')
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
completeA %>% separate(GTDB_taxonomy, c('Domain_Phyla', 'Lower_class'), ";c__") -> completeA
completeA %>% separate(Domain_Phyla, c('Domain', 'Phyla'), ';p__') -> completeA
#red for reduced
#complete_red now contains 107092 species


#this breaks it up in taxons
completeA %>% separate(Lower_class, c('Class', 'Lower_order'), ";o__") -> completeA
completeA %>% separate(Lower_order, c('Order', 'Lower_family'), ";f__") -> completeA
completeA %>% separate(Lower_family, c('Family', 'Lower_genera'), ";g__") -> completeA
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
```{r}
Abundance <- completeA |> group_by(Phyla) |> summarize(observations=n(), groups = 'drops') |> arrange(desc(observations))



Phyla_CPA_count <- aggregate(completeA$X.CPA, by=list(Category=completeA$Phyla), FUN=sum)
Phyla_CPA_count <- Phyla_CPA_count %>% rename(
    Phyla = Category,
    total_CPA=x)



completeA$CPA_binary[completeA$X.CPA == 0] <- 0
completeA$CPA_binary[completeA$X.CPA > 0] <- 1


Uniq_CPA_Phyla <- aggregate(completeA$CPA_binary, by=list(Category=completeA$Phyla), FUN=sum)
Uniq_CPA_Phyla <- Uniq_CPA_Phyla %>% rename(
    Phyla = Category,
    Species_with_CPA=x)
Atotal <- merge(Abundance, Uniq_CPA_Phyla, by='Phyla')
Atotal <- merge(Atotal, Phyla_CPA_count, by='Phyla')
Atotal <- Atotal[-3]
Atotal$percentage = Atotal$Species_with_CPA/Atotal$observations*100

Atotal$cpa_ratio_if_cpa_present <- Atotal$total_CPA/Atotal$Species_with_CPA
mean(total$cpa_ratio_if_cpa_present, na.rm=TRUE)
Atotal$cpa_total_ratio <- Atotal$total_CPA/total$observations
mean(Atotal$cpa_total_ratio, na.rm=TRUE)

```
```{r make plot of most abudant phyla, with mean}

m <- ggplot(data=completeA, mapping = aes(y=Phyla,
                                           x=X.CPA,
                                           color=Phyla)) + 
  geom_jitter() + theme_minimal() + 
  stat_summary(fun=mean, 
               geom = 'point',
               size = 2,
               colour = 'gray20') + geom_vline(xintercept = 1)
```

```{r}
x <- ggplot(data = Atotal, mapping = aes(y=total_CPA,
                                   x=Species_with_CPA,
                                   color=Phyla)) + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point() + theme_minimal()
#note cyanobaceriota
```

```{r}
Halobacteriota <- filter(completeA, Phyla== 'Halobacteriota')
Phyla_HaloClass_CPA_count <- aggregate(Halobacteriota$X.CPA, by=list(Category=Halobacteriota$Class), FUN=sum)

Phyla_HaloOrder_CPA_count <- aggregate(Halobacteriota$X.CPA, by=list(Category=Halobacteriota$Order), FUN=sum)
Phyla_HaloFamily_CPA_count <- aggregate(Halobacteriota$X.CPA, by=list(Category=Halobacteriota$Family), FUN=sum)
```

```{r}
library(networkD3)
install.packages('networkD3')
library(networkD3)
library(dplyr)
links <- data.frame(source=c('Total_CPA', 'Total_CPA', 'Total_CPA', 'Total_CPA', 'Total_CPA','Total_CPA', 'Halobacteriota', 'Halobacteriota', 'Halobacteriota', 'Halobacteriota', 'Halobacteria', 'Halobacteria', 'Halobacteria', 'Halobacteria'),
                    target=c('Halobacteriota', 'Thermoproteota', 'Thermoplasmatota', 'Nanoarchaeota', 'Micrarchaeota','Other', 'Halobacteria', 'Methanosarcina', 'Methanomicrobia', 'Other_class', 'Natrialbaceae', 'Haloferacaceae', 'Haloarculaceae','Other_family'),
                    value=c(3860,2436,1730,1394,554, 2226, 2221, 639,712, 288, 846, 496,461,213))
nodes <- data.frame(
  name=c(as.character(links$source), 
  as.character(links$target)) %>% unique()
)
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
 
# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE)
p

#total CPA count
```

