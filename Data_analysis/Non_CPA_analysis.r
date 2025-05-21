IT_archaea <- read.delim('/nesi/nobackup/uc04105/new_databases_May/GTDB_226/IT_hitArchaea2025.tsv', sep='\t', stringsAsFactors = TRUE)

IT_archaea$Total_IT <- IT_archaea$NhaB_count + IT_archaea$NhaC_count + IT_archaea$NhaD_count

IT_archaea$IT_binary[IT_archaea$Total_IT == 0] <- 0
IT_archaea$IT_binary[IT_archaea$Total_IT > 0] <- 1
IT_archaea %>% separate(GTDB_taxonomy, c('Domain_Phyla', 'Lower_class'), ";c__") -> IT_archaea
IT_archaea %>% separate(Domain_Phyla, c('Domain', 'Phyla'), ';p__') ->IT_archaea
IT_archaea %>% separate(Lower_class, c('Class', 'Lower_order'), ";o__") -> IT_archaea
IT_archaea %>% separate(Lower_order, c('Order', 'Lower_family'), ";f__") -> IT_archaea
IT_archaea %>% separate(Lower_family, c('Family', 'Lower_genera'), ";g__") -> IT_archaea


sum(IT_archaea$IT_binary)
sum(IT_archaea$NhaD_binary)


Abundance <- IT_archaea |> group_by(Phyla) |> summarize(observations=n(), groups = 'drops') |> arrange(desc(observations))

Phyla_NhaD_count <- aggregate(IT_archaea$NhaD_count, by=list(Category=IT_archaea$Phyla), FUN=sum)
Phyla_NhaD_count <- Phyla_NhaD_count %>% rename(
    Phyla = Category,
    total_NhaD=x)

Uniq_NhaD_Phyla <- aggregate(IT_archaea$NhaD_binary, by=list(Category=IT_archaea$Phyla), FUN=sum)
Uniq_NhaD_Phyla <- Uniq_NhaD_Phyla %>% rename(
    Phyla = Category,
    Species_with_NhaD=x)

Phyla_NhaC_count <- aggregate(IT_archaea$NhaC_count, by=list(Category=IT_archaea$Phyla), FUN=sum)
Phyla_NhaC_count <- Phyla_NhaC_count %>% rename(
    Phyla = Category,
    total_NhaC=x)

Uniq_NhaC_Phyla <- aggregate(IT_archaea$NhaC_binary, by=list(Category=IT_archaea$Phyla), FUN=sum)
Uniq_NhaC_Phyla <- Uniq_NhaC_Phyla %>% rename(
    Phyla = Category,
    Species_with_NhaC=x)

Phyla_NhaB_count <- aggregate(IT_archaea$NhaB_count, by=list(Category=IT_archaea$Phyla), FUN=sum)
Phyla_NhaB_count <- Phyla_NhaB_count %>% rename(
    Phyla = Category,
    total_NhaB=x)

Uniq_NhaB_Phyla <- aggregate(IT_archaea$NhaB_binary, by=list(Category=IT_archaea$Phyla), FUN=sum)
Uniq_NhaB_Phyla <- Uniq_NhaB_Phyla %>% rename(
    Phyla = Category,
    Species_with_NhaB=x)

Atotal <- merge(Abundance, Uniq_NhaD_Phyla, by='Phyla')
Atotal <- merge(Atotal, Phyla_NhaD_count, by='Phyla')
Atotal <- merge(Atotal, Uniq_NhaC_Phyla, by='Phyla')
Atotal <- merge(Atotal, Phyla_NhaC_count, by='Phyla')
Atotal <- merge(Atotal, Uniq_NhaB_Phyla, by='Phyla')
Atotal <- merge(Atotal, Phyla_NhaB_count, by='Phyla')
Atotal <- Atotal[-3]

lm(log(Atotal$Species_with_NhaB+1) ~log(Atotal$observations), data=Atotal)
NhaB<- ggplot(Atotal, aes(x=log(Atotal$observations), y=log(Atotal$Species_with_NhaB+1), col=Phyla)) + geom_smooth(method="lm", col='black', se=FALSE, linetype="dashed") + geom_point(size=3) + theme(legend.position="none") + ggtitle('NhaB') + xlim(0,8) + ylim(0,8) 

lm(log(Atotal$Species_with_NhaC+1) ~log(Atotal$observations), data=Atotal)
NhaC<- ggplot(Atotal, aes(x=log(Atotal$observations), y=log(Atotal$Species_with_NhaC+1), col=Phyla)) + geom_smooth(method="lm", col='black', se=FALSE, linetype="dashed") + geom_point(size=3) + theme(legend.position="none") + ggtitle('NhaC') + xlim(0,8) + ylim(0,8) 
lm(log(Atotal$Species_with_NhaD+1) ~log(Atotal$observations), data=Atotal)
NhaD<- ggplot(Atotal, aes(x=log(Atotal$observations), y=log(Atotal$Species_with_NhaD+1), col=Phyla)) + geom_smooth(method="lm", col='black', se=FALSE, linetype="dashed")+ geom_point(size=3) + ggtitle('NhaD') + xlim(0,8) + ylim(0,8) 
grid.arrange(NhaB, NhaC, NhaD, ncol=3)


