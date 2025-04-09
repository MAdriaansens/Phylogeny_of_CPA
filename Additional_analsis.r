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
complete_NhaD_edit <- read.delim('/nesi/nobackup/uc04105/CPA_hitBacteria28032025_NhaD_edit.tsv', sep = '\t')
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
complete_NhaD_edit %>% separate(GTDB_taxonomy, c('Domain_Phyla', 'Lower_class'), ";c__") -> complete_NhaD_edit
complete_NhaD_edit %>% separate(Domain_Phyla, c('Domain', 'Phyla'), ';p__') -> complete_NhaD_edit
#red for reduced
#complete_red now contains 107092 species


#this breaks it up in taxons
complete_NhaD_edit %>% separate(Lower_class, c('Class', 'Lower_order'), ";o__") -> complete_NhaD_edit
complete_NhaD_edit %>% separate(Lower_order, c('Order', 'Lower_family'), ";f__") -> complete_NhaD_edit
complete_NhaD_edit %>% separate(Lower_family, c('Family', 'Lower_genera'), ";g__") -> complete_NhaD_edit

complete_red <- complete_NhaD_edit
Phyla_NhaD_count <- aggregate(complete_red$NhaD_count, by=list(Category=complete_red$Phyla), FUN=sum)
Uniq_NhaD_Phyla <- aggregate(complete_red$NhaD_binary, by=list(Category=complete_red$Phyla), FUN=sum)

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

Phyla_NhaD_count <- Phyla_NhaD_count %>% rename(
  Phyla=Category,
  NhaD_count = x
)
Uniq_NhaD_Phyla<- Uniq_NhaD_Phyla %>% rename(
    Phyla = Category,
    species_with_NhaD=x)

total <- merge(Abundance, Uniq_CPA_Phyla, by='Phyla')
total <- merge(total, Phyla_CPA_count, by='Phyla')
total <- merge(total, Phyla_NhaD_count, by='Phyla')
total <- merge(total, Uniq_NhaD_Phyla, by='Phyla')
total$species_without_CPA <- total$observations-total$species_with_CPA
lmComplete = lm(CPA_count~ NhaD_binary, data = complete_red)
summary(lmComplete)
```
```{r}
important_phyla_list = c(T10_CPAPhyla, T10_uniq_CPA_phyla, Most_abundant_phyla)
ImpPhyla <- unique(important_phyla_list)
ImpPhyla #totals 13 phyla
total_top <-filter(total, Phyla %in% ImpPhyla)

```


```{r}

ggplot(data = total_top, mapping = aes(y=species_with_CPA,
                                   x=species_with_NhaD,
                                   color=Phyla, show.legend = FALSE)) +  theme_minimal() + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point(size=5)
```


```{r}
NhaD_Family_count <- aggregate(complete_red$NhaD_count, by=list(Category=complete_red$Family), FUN=sum)
NhaD_Family_count <- NhaD_Family_count %>% rename(
    Family = Category,
    total_NhaD=x)
CPA_Family_count <- aggregate(complete_red$CPA_count, by=list(Category=complete_red$Family), FUN=sum)
CPA_Family_count <- CPA_Family_count %>% rename(
    Family = Category,
    total_CPA=x)
total_family <- merge(NhaD_Family_count, CPA_Family_count, by='Family')
total_family <- filter(total_family, total_NhaD > 0) 

```

```{r}
lmComplete = lm(total_CPA~total_NhaD, data = total_family)
summary(lmComplete)
```
```{r}
ggplot(data = total_family, mapping = aes(x=total_NhaD,
                                   y=total_CPA,
                                   color=Family, show.legend = FALSE)) +  theme_minimal() + geom_smooth(method=lm, se=FALSE,col='gray40', linetype='dashed') + geom_point(size=5) + theme(legend.position = "none")
```
```{r}
ggplot(data = complete_red, mapping = aes(x=CPA_count,
                                   y=NhaD_count, show.legend = FALSE)) +  theme_minimal()  + geom_point()
```

```{r}
par(mfrow=c(1,2))

hist(complete_red$CPA_count,
     main="",
     xlab = 'CPA_count',
     ylab = 'Observations',
     col='yellow')

hist(complete_red$NhaD_count,
     main="",
     xlab = 'NhaD_count',
     ylab = 'Observations',
     col='red')
```

