---
title: "Bacteria_stat"
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

```{r}
library(ggplot2)
library(tidyr)
library(dplyr)
```


```{r load data and check phyla without CPA}
completeB <- read.delim('/nesi/nobackup/uc04105/CPA_hitBacteria28032025.tsv', sep = '\t')
missing <- read.delim('/nesi/nobackup/uc04105/missing_phyla.tsv', sep = '\t')

missing_phyla <- missing$Phyla

missing_phyla
#totaling ~143 species
```


```{r make taxonomic ordering from GTDB}
#
completeB %>% separate(GTDB_taxonomy, c('Domain_Phyla', 'Lower_class'), ";c__") -> completeB
completeB %>% separate(Domain_Phyla, c('Domain', 'Phyla'), ';p__') -> completeB
#red for reduced
complete_red <- data.frame(filter(completeB, !completeB$Phyla %in% missing_phyla))
#complete_red now contains 107092 species


#this breaks it up in taxons
complete_red %>% separate(Lower_class, c('Class', 'Lower_order'), ";o__") -> complete_red
complete_red %>% separate(Lower_order, c('Order', 'Lower_family'), ";f__") -> complete_red
complete_red %>% separate(Lower_family, c('Family', 'Lower_genera'), ";g__") -> complete_red

```



```{r}
ggplot(data = completeB, mapping = aes(x=CheckM2_completeness,
                                   y=CheckM2_contamination,
                                   color=CPA_binary))+ geom_point()  + geom_smooth(model=lm, col='red', se = FALSE) 
```

```{r create abundance column per phyla and top10 abundance}
#top species in GTDB bacteria
Top10_most_abundant <- complete_red |> group_by(Phyla) |> summarize(observations=n(), groups = 'drops') |> arrange(desc(observations)) |> head(10)
Most_abundant_phyla <- Top10_most_abundant $Phyla
Abundance <- complete_red |> group_by(Phyla) |> summarize(observations=n(), groups = 'drops') |> arrange(desc(observations)) 
Abundance <- Abundance[-3]

```

```{r create # total CPA  column per phyla and top10 }
#number of CPAs intotal
Phyla_CPA_count <- aggregate(complete_red$CPA_count, by=list(Category=complete_red$Phyla), FUN=sum)
T10_CPAPhyla<- Phyla_CPA_count %>% arrange(desc(Phyla_CPA_count$x)) |> head(10)
T10_CPAPhyla <- T10_CPAPhyla$Category

Phyla_CPA_count <- Phyla_CPA_count %>% rename(
    Phyla = Category,
    total_CPA=x)
#most CPA
```

```{r create #CPA specific species per phyla and top10, merge data}
Uniq_CPA_Phyla <- aggregate(complete_red$CPA_binary, by=list(Category=complete_red$Phyla), FUN=sum)
T10_uniq_CPA_phyla <- Uniq_CPA_Phyla  %>% arrange(desc(Uniq_CPA_Phyla $x)) |> head(10) 
T10_uniq_CPA_phyla  <- T10_uniq_CPA_phyla$Category

Uniq_CPA_Phyla <- Uniq_CPA_Phyla %>% rename(
    Phyla = Category,
    species_with_CPA=x)


total <- merge(Abundance, Uniq_CPA_Phyla, by='Phyla')
total <- merge(total, Phyla_CPA_count, by='Phyla')
total$percentage = total$species_with_CPA/total$observations*100
total$cpa_ratio_if_cpa_present <- total$total_CPA/total$species_with_CPA
mean(total$cpa_ratio_if_cpa_present)
total$cpa_total_ratio <- total$total_CPA/total$observations
mean(total$cpa_total_ratio)
```


```{r map CPA by phyla}
ggplot(data = complete_red, mapping = aes(y=Phyla,
                                   x=CPA_count,
                                   color=Phyla)) + geom_jitter() + geom_boxplot() + theme_minimal() + theme(legend.position = "none")
```

```{r plot percentage/observations per phyla}
ggplot(data = total, mapping = aes(y=percentage,
                                   x=observations,
                                   color=Phyla)) + geom_point() + theme_minimal() + theme(legend.position = "none")
```
```{r}
d<- ggplot(data = total, mapping = aes(y=total_CPA,
                                   x=observations,
                                   color=Phyla)) + geom_point() + theme_minimal()
```

```{r closer analysis into Pseudomonadota, as it contains ~1/3 CPAs}
Pseudomonadota <- filter(complete_red, Phyla== 'Pseudomonadota')
Phyla_PseuClass_CPA_count <- aggregate(Pseudomonadota$CPA_count, by=list(Category=Pseudomonadota$Class), FUN=sum)

Phyla_PseuOrder_CPA_count <- aggregate(Pseudomonadota$CPA_count, by=list(Category=Pseudomonadota$Order), FUN=sum)
```

```{r Sankey plot, with focus on Pseudomonadota}
library(networkD3)
install.packages('networkD3')
library(networkD3)

library(dplyr)
links <- data.frame(source=c("Total CPA", "Total CPA", "Total CPA", "Total CPA", "Total CPA", "Total CPA","Total CPA", "Total CPA", "Pseudomonadota", "Pseudomonadota", "Pseudomonadota", "Pseudomonadota", "Gamma proteobacteriota", "Gamma proteobacteriota", "Alpha proteobacteriota", "Gamma proteobacteriota", "Alpha proteobacteriota"),
                    target=c("Pseudomonadota", "Actinomycetota", "Bacteroidota", "Bacillota_A", "Bacillota", "Cyanobacteriota", "Planctomycetota", "Other", "Gamma proteobacteriota", "Alpha proteobacteriota", "Zeta proteobaceteriota", "Magnetococcia", "Burkholderiales","Pseudomonadales", "Rhizobiales", "Enterobacterales", "Sphingomonadales"),
                    value=c(28174, 12273, 11632, 4747, 3286, 2635, 1611, 13922, 16846,11204,96,28, 5792,4123,3894,3153,2083))
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

```{r take top10 of species with CPA count, total CPA and Overall species abundance}
#append top 10 of all to one list
important_phyla_list = c(T10_CPAPhyla, T10_uniq_CPA_phyla, Most_abundant_phyla)
ImpPhyla <- unique(important_phyla_list)
ImpPhyla #totals 13 phyla

#Cyanobacteria not in top10 most abudant taxonomic species/species with CPA but is 6th on number of CPAs
```


```{r make plot of most common species}
par(mfrow=c(2,3))
total_top <-filter(total, Phyla %in% ImpPhyla)
ggplot(data = total_top, mapping = aes(y=total_CPA,
                                   x=species_with_CPA,
                                   color=Phyla, show.legend = FALSE)) +  theme_minimal() + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point(size=5) +  scale_x_log10() +  scale_y_log10()

ggplot(data = total_top, mapping = aes(y=total_CPA,
                                   x=observations,
                                   color=Phyla, show.legend = FALSE)) +  theme_minimal() + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point(size=5) +  scale_x_log10() +  scale_y_log10()

ggplot(data = total_top, mapping = aes(y=species_with_CPA,
                                   x=observations,
                                   color=Phyla, show.legend = FALSE)) +  theme_minimal() + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point(size=5) +  scale_x_log10() +  scale_y_log10()
#note cyanobaceriota
```

```{r make plot of most abudant phyla, with mean}
most_abudant <- filter(complete_red,Phyla %in% ImpPhyla)

m <- ggplot(data=most_abudant, mapping = aes(y=Phyla,
                                           x=CPA_count,
                                           color=Phyla)) + 
  geom_jitter() + theme_minimal() + 
  stat_summary(fun=mean, 
               geom = 'point',
               size = 2,
               colour = 'gray20') + geom_vline(xintercept = 1)
```


```{r subset into other phyla if needed}
#subset into more interesting and smaller phyla
#Cyanobacteriota<- filter(complete_red, Phyla == 'Cyanobacteriota')


#Bacteroidota <- filter(complete_red, Phyla== 'Bacteroidota')
```


```{r get covariance completeness/cpa count}
covariance <- cov(completeB$CheckM2_completeness, completeB$CPA_count)
covariance
```

```{r correlation test completeness and CPA count}
cor.test(completeB$CheckM2_completeness, completeB$CPA_count)

```

```{r}
#linear regression

lmCcomPCbCcon = lm(CheckM2_completeness~Phyla + CPA_binary + CheckM2_contamination, data = completeB)
summary(lmComplete)
```

```{r}
#lmCcomPCoCcon = lm(CheckM2_completeness~Phyla + CPA_count + CheckM2_contamination, data = completeB)
#summary(lmComplete)
```

```{r}
#lmComPCcom = lm(CPA_count~Phyla + CheckM2_completeness, data = completeB)
#summary(lmComplete)
```

```{r}
#lmCoPCcomCcon = lm(CPA_count~Phyla + CheckM2_completeness + CheckM2_contamination, data = completeB)
#summary(lmComplete)
```

```{r}
#lmCcP = lm(CPA_count~Phyla + CheckM2_contamination, data = completeB)
#summary(lmComplete)
```

```{r}
#lmCoCcom_star_Cb = lm(CPA_count~CheckM2_contamination*CPA_binary, data = completeB)
#summary(lmComplete)
```

