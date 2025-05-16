

count_genomes_subtree_Halo_Arc4079404 <- read.delim('/nesi/nobackup/uc04105/redvals/subset/Archaea_subtree_count_genomes_Halo_Arc4079404.tsv', sep='\t')
count_genomes_nodupes_subtree_Halo_Arc4079404<- ddply(count_genomes_subtree_Halo_Arc4079404,.(GTDB_ID,GTDB_tax),nrow)

ggplot(data=count_genomes_nodupes_subtree_Halo_Arc4079404, aes(x=V1)) + geom_histogram() + ggtitle('number of CPA Halo_Arc4079404')

count_genomes_nodupes_subtree_Halo_Arc4079404 %>% separate(GTDB_tax, c('Domain_Phyla', 'Lower_class'), ";c__") -> count_genomes_nodupes_subtree_Halo_Arc4079404
count_genomes_nodupes_subtree_Halo_Arc4079404 %>% separate(Domain_Phyla, c('Domain', 'Phyla'), ';p__') -> count_genomes_nodupes_subtree_Halo_Arc4079404

count_genomes_nodupes_subtree_Halo_Arc4079404 %>% separate(Lower_class, c('Class', 'Lower_order'), ";o__") -> count_genomes_nodupes_subtree_Halo_Arc4079404
count_genomes_nodupes_subtree_Halo_Arc4079404%>% separate(Lower_order, c('Order', 'Lower_family'), ";f__") -> count_genomes_nodupes_subtree_Halo_Arc4079404
count_genomes_nodupes_subtree_Halo_Arc4079404 %>% separate(Lower_family, c('Family', 'Lower_genera'), ";g__") -> count_genomes_nodupes_subtree_Halo_Arc4079404


Values_subtree_subtree_Halo_Arc4079404<- read.delim('/nesi/nobackup/uc04105/redvals/subset/red_of_all_pairs_Archaea_full__Halo_Arc4079404.tsv', sep ='\t')
head(Values_subtree_subtree_Halo_Arc4079404)
ggplot(data=Values_subtree_subtree_Halo_Arc4079404, mapping=aes(x=RED,
                                                           y=Phylo,
                                                           col=status)) + geom_point(alpha=0.2)


ggplot(data=Values_subtree_subtree_Halo_Arc4079404, mapping=aes(x=RED,
                                                           y=Phylo)) + geom_point() + geom_smooth(method='lm')




ggplot(data=Values_subtree_subtree_Halo_Arc4079404, aes(x=RED)) + geom_histogram(fill='blue') + geom_hline(yintercept = 100, col='red')

clusteroids_subtree_Halo_Arc4079404 = filter(Values_subtree_subtree_Halo_Arc4079404, Phylo==0)
ggplot(clusteroids_subtree_Halo_Arc4079404, aes(x=RED)) + geom_histogram(fill='blue', col='black')
mean(clusteroids_subtree_Halo_Arc4079404$RED)
nrow(Values_subtree_subtree_Halo_Arc4079404)
