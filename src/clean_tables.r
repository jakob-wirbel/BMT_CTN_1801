# ##############################################################################
#
## Clean all tables
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggthemes")

colours <- yaml::read_yaml(here('files/colours.yml'))

# raw (already somewhat fixed, actually) metadata
df.meta <- read_tsv(here('data/meta_sample.tsv'),
                    col_types = cols())
dim(df.meta)
# 1] 2630   5

#' originally, we had 2631 samples, but I had removed one already after
#' Emmes mentioned that there was a potential mislabeling event for
#' the tube P11E07
#' 
#' Another issue that is fixed already in this file:
#' I mislabeled the tubes for sequencing for P24E04 and P24F04
#' Which one is which?
#' 
#' Qiagen concentration
#' P24E04 126.20
#' P24F04 15.99
#' 
#' Novogene concentrations for this batch:
#' P24E04_1 16.74 --> this one is the real P24F04
#' P24E04_2 177.28 --> this one is the real P24E04
#' 
#' Keep the new names (P24E04_1 & P24E04_2), since the sequencing data keeps
#' these identifiers

# ##############################################################################
# add 16S prediction

df.16S <- read_tsv(here('data/copies_16S.tsv'), col_types = cols())

table(df.16S$copies)
# measured predicted 
# 858      1772 

# this comes out of the prediction model we built earlier!
# see the `absolute_abundance_prediction` script
df.meta <- df.meta %>% 
  full_join(df.16S, by='Sample_ID')

# ##############################################################################
# add the sequencing stats as well

df.seq_stats <- read_tsv(here('data/sequencing_stats.tsv'), col_types = cols())

df.meta <- df.meta %>% 
  full_join(df.seq_stats, by='Sample_ID')

# ##############################################################################
# some accounting!

# df.meta %>%
#   group_by(Timepoint) %>%
#   tally() %>% View

# 1 EXTRA (this one we need to filter out)
# 56 NA (those are ZYMO or negative controls)
# 2631 - 56 - 1 (mislabel) - 1 (EXTRA) -> 2573 stool samples sequenced

table(df.meta$Sample_type)
# next_day same_day 
# 2324      250 
# the 1 EXTRA timepoint sample is same_day, not next_day

# ##############################################################################
# perform sequencing checks on the raw motus/metaphlan data

#' load both motus and metaphlan
feat.all <- read.table(here('data/motus_all.tsv'), 
                       sep = '\t', comment.char = '', quote = '', header = TRUE,
                       row.names = 1, skip = 2)
feat.mpa <- read.table(here('data/metaphlan_all.tsv'), 
                       sep = '\t', comment.char = '', quote = '', header = TRUE,
                       row.names = 1)
rownames(feat.mpa)[1] <- 'x|s__UNCLASSIFIED'

#' remove motus with all zeros
feat.all <- feat.all[rowSums(feat.all==0)!=ncol(feat.all),]
#' restrict metaphlan output to species only (we don't really care about 
#' strains?)
feat.mpa <- feat.mpa[str_detect(rownames(feat.mpa), '\\|s__'),]
feat.mpa <- feat.mpa[str_detect(rownames(feat.mpa), '\\|t__', negate = TRUE),]

#' remove some samples with low counts
#' check that those are the ones that we would have expected from the failed
#' sequencing!
#' 
#' Keep samples with more than 1000 motus counts or 
#' samples with more than 10M reads?
#' 
#' The samples with lots of reads, but low motus counts, seem to be dominated
#' by funghi (according to Metaphlan)
#' Maybe we need to loop this in as well?
#' 
#' Most of the samples filtered out this way were from the failed sequencing 
#' anyway, so it would be expected to a certain extent
#' 
#' New approach: Remove all samples with fewer than 1M reads for failed 
#' sequencing. This kicks out 27/28 negative controls, and the left over 
#' negative controls seems to be contaminated by the adjacent well (see below)
#' 
#' The number of real samples kicked out is also pretty small
#' There are still a couple of samples with very low bacterial sequencing, 
#' those dominated by funghi, but that should be okay.
#' Maybe we can combine the motus counts and metaphlan relative abundances 
#' into eukaryotic counts? --> the problem is that metaphlan does not output
#' counts, so alpha diversity comparisons are very hard to do properly
#' 
#' 
#' How about similarity between adjacent wells? Maybe there was some spillover?
#' some of the negative controls have quite a bit of fecal-like bacteria, 
#' for example P08G12 or P01G12
#' seems like P08G11 and P08G12 are not too unsimilar (spearman's rho = 0.5637)

feat.neg <- feat.all[,df.meta %>% filter(Timepoint=='negative_control') %>% 
                       pull(Sample_ID)]

df.filt <- enframe(colSums(feat.all), name = 'Sample_ID', 
                   value = 'mOTUs_count') %>% 
  full_join(df.meta, by='Sample_ID') %>% 
  mutate(filter=case_when((mOTUs_count > 100 | 
                             hostremoved_reads>1e6)~'keep',
                          TRUE~'drop')) 
# visualize
df.filt %>% 
  ggplot(aes(x=log10(hostremoved_reads), y=log10(mOTUs_count), col=filter)) + 
  geom_point() + 
  facet_grid(~batch) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = log10(100)) + 
  geom_vline(xintercept = log10(1e06))

# number of human reads?
df.filt %>% 
  mutate(frac=1-(hostremoved_reads/trimmed_reads)) %>% 
  ggplot(aes(x=log10(mOTUs_count), y=frac)) + 
  facet_grid(~batch) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  geom_point() + 
  ylab('Fraction of human reads in the raw sequencing')

# ##############################################################################
# ZYMO

#' also look at the Zymo sequencing
#' There should be a specific number of species in there
#' Taken mostly from the AWI-Gen2 code

feat.zymo <- feat.all[,df.meta %>% filter(Timepoint=='ZYMO') %>% 
                        pull(Sample_ID)]
zymo.cols <- c('Escherichia coli'="#FF0000", 
               'Salmonella enterica'="#00A08A", 
               'Bacillus subtilis'="#F2AD00", 
               'Pseudomonas aeruginosa'="#F98400", 
               'Enterococcus faecalis'="#5BBCD6", 
               'Lactobacillus fermentum'='#ECCBAE',
               'Limosilactobacillus fermentum'='#ECCBAE',
               'Staphylococcus aureus'="#046C9A", 
               'Listeria monocytogenes'="#D69C4E",
               'Saccharomyces cerevisiae'='#ABDDDE', 
               'Cryptococcus neoformans'='#8C1515',
               'unassigned'='#707273',
               'unclassified'='#707273',
               'UNCLASSIFIED'='#707273',
               'other'='#D3D3D3')

ref.motus <- tibble(
  motus=c('Escherichia coli', 'Salmonella enterica', 'Bacillus subtilis',
          'Pseudomonas aeruginosa', 'Enterococcus faecalis', 
          'Lactobacillus fermentum', 'Staphylococcus aureus', 
          'Listeria monocytogenes', 'Saccharomyces cerevisiae', 
          'Cryptococcus neoformans'), name='reference', 
  rel.ab=c(rep(0.12, 8), 0.02, 0.02))
motus.zymo <- prop.table(as.matrix(feat.zymo),2)
motus.zymo <- motus.zymo[rowMeans(motus.zymo==0)!=1,]

g.zymo.motus <- motus.zymo %>% 
  as_tibble(rownames = 'motus') %>% 
  pivot_longer(-motus) %>% 
  mutate(motus=str_remove(motus, ' \\[.*\\]')) %>% 
  mutate(motus=case_when(motus %in% names(zymo.cols)~motus, TRUE~'other')) %>% 
  group_by(motus, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.motus) %>% 
  mutate(motus=factor(motus, levels=names(zymo.cols))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=motus)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  xlab('') + ylab('Relative abundance') + 
  scale_fill_manual(values = zymo.cols) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ggtitle('mOTUs3')
g.zymo.motus
ggsave(g.zymo.motus, filename = here('figures', 'general', 'zymo_motus.pdf'),
       width = 6, height = 4, useDingbats=FALSE)


#' seems okay overall, even though some ZYMO controls have lots and lots of
#' unassigned abundance... what is that? Maybe it looks different in the 
#' metaphlan output?

feat.zymo.mpa <- feat.mpa[,df.meta %>% filter(Timepoint=='ZYMO') %>% 
                            pull(Sample_ID)]
mpa.zymo <- prop.table(as.matrix(feat.zymo.mpa),2)
mpa.zymo <- mpa.zymo[rowMeans(mpa.zymo==0)!=1,]

g.zymo.mpa <- mpa.zymo %>% 
  as_tibble(rownames = 'metaphlan') %>% 
  pivot_longer(-metaphlan) %>% 
  mutate(metaphlan=str_remove(metaphlan, '^.*\\|s__')) %>% 
  mutate(metaphlan=str_replace_all(metaphlan, '_', ' ')) %>% 
  mutate(metaphlan=case_when(metaphlan %in% names(zymo.cols)~metaphlan, 
                             TRUE~'other')) %>%
  group_by(metaphlan, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.motus %>%
              rename(metaphlan=motus) %>% 
              mutate(metaphlan=case_when(
                metaphlan=='Lactobacillus fermentum'~
                  'Limosilactobacillus fermentum',
                TRUE~metaphlan))) %>% 
  mutate(metaphlan=factor(metaphlan, levels=names(zymo.cols))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=metaphlan)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  xlab('') + ylab('Relative abundance') + 
  scale_fill_manual(values = zymo.cols) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ggtitle('metaphlan')
g.zymo.mpa
ggsave(g.zymo.mpa, filename = here('figures', 'general', 'zymo_mpa.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

#' Metaphlan doesn't find the Bacillus subtilis at all, again
#' Super weird, but whatevs...

# we could show a PCoA here?
x <- mpa.zymo %>% 
  as_tibble(rownames = 'species') %>% 
  pivot_longer(-species) %>% 
  mutate(species=str_remove(species, '^.*\\|s__')) %>% 
  mutate(species=str_replace_all(species, '_', ' ')) %>% 
  mutate(species=case_when(species %in% names(zymo.cols)~species, 
                           TRUE~'other')) %>%
  mutate(species=case_when(
    species=='Limosilactobacillus fermentum'~'Lactobacillus fermentum', 
    species=='UNCLASSIFIED'~'unassigned',
    TRUE~species)) %>% 
  group_by(species, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>%
  mutate(type='metaphlan') %>% 
  bind_rows(motus.zymo %>% 
              as_tibble(rownames = 'species') %>% 
              pivot_longer(-species) %>% 
              mutate(species=str_remove(species, ' \\[.*\\]')) %>% 
              mutate(species=case_when(species %in% names(zymo.cols)~species, 
                                       TRUE~'other')) %>% 
              group_by(species, name) %>% 
              summarise(rel.ab=sum(value), .groups='drop') %>% 
              mutate(type='motus')) %>% 
  bind_rows(ref.motus %>%
              rename(species=motus) %>% 
              mutate(type='reference')) %>% 
  pivot_wider(names_from = c(type, name), values_from = rel.ab,
              values_fill = 0) %>% 
  as.data.frame()
rownames(x) <- x$species
x$species <- NULL

labdsv::pco(vegan::vegdist(t(x)))$points %>% 
  as.data.frame() %>% 
  as_tibble(rownames='sample') %>% 
  separate(sample, into=c('profiling', 'id'), sep='_') %>% 
  ggplot(aes(x=V1, y=V2, col=profiling)) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

# ##############################################################################
# Combine Metaphlan rel ab for Eukarya and motus counts for alpha-diversity?

eukaryotes <- feat.mpa[str_detect(rownames(feat.mpa), 'k__Eu'),]
# how much percent are the other kingdoms?
eukaryotes <- eukaryotes/100
rest <- 1 - colSums(eukaryotes)

# re-scale the eukaryotes to be similar to motus counts?
motus.counts <- colSums(feat.all[,names(rest)])

euk.counts <- (motus.counts/rest) - motus.counts
eukaryotes <- prop.table(as.matrix(eukaryotes), 2)
eukaryotes[is.na(eukaryotes)] <- 0
fixed.euk <- round(t(euk.counts * t(eukaryotes)))
rownames(fixed.euk) <- str_remove(rownames(fixed.euk), '.*\\|s__')

feat.adjusted <- rbind(feat.all, fixed.euk[,colnames(feat.all)])

g.before <- df.filt %>% 
  ggplot(aes(x=log10(hostremoved_reads), y=log10(mOTUs_count), col=filter)) + 
  geom_point() + 
  facet_grid(~batch) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = log10(100)) + 
  geom_vline(xintercept = log10(1e06))
g.before

df.filt <- enframe(colSums(feat.adjusted), name = 'Sample_ID', 
                   value = 'mOTUs_count') %>% 
  full_join(df.meta, by='Sample_ID') %>% 
  mutate(filter=case_when((mOTUs_count > 100 | 
                             hostremoved_reads>1e6)~'keep',
                          TRUE~'drop')) 

g.after <- df.filt %>% 
  ggplot(aes(x=log10(hostremoved_reads), y=log10(mOTUs_count), col=filter)) + 
  geom_point() + 
  facet_grid(~batch) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  geom_hline(yintercept = log10(100)) + 
  geom_vline(xintercept = log10(1e06))
g.after

feat.euk.added <- feat.adjusted[,df.meta$Sample_ID]
feat.euk.added <- feat.euk.added[rowSums(feat.euk.added==0)!=
                                   ncol(feat.euk.added),]

# ##############################################################################
# look at the negative controls 

# most of them fall out by filtering
# 2 stay

# what kind of taxa do we get in there?
motus.neg <- feat.euk.added[,df.meta %>% 
                              filter(Timepoint=='negative_control') %>% 
                              pull(Sample_ID)]
# mostly some weird ones, like skin ones
# most of the taxa have super low abundnace, single digit motus counts -> 
# spurious, or low-abundance

colSums(motus.neg)
# only 2 have more than 100 motus counts
# these ones might be well contamination

sp.cor <- cor(feat.euk.added[,c('P01G12', 'P08G12')], feat.euk.added, 
              method='spearman')
# > tail(sort(sp.cor[1,]))
# P01A02    P01E02    P01G11    P01B12    P01A08    P01G12 
# 0.5087290 0.5135943 0.5147058 0.5151517 0.5533719 1.0000000 
# > tail(sort(sp.cor[2,]))
# P08E10    P08E12    P08A12    P08G11    P08D11    P08G12 
# 0.4211645 0.4588259 0.4796312 0.5637959 0.6231559 1.0000000

# very likely some contamination! P08D11 -> P08G12, P01G11 -> P01G12

# in metaphlan, most of the top species is always unclassified, but it's a bit
# harder to assess overall


# ##############################################################################
# compare same_day vs next_day sample handling

# number of samples
g <- df.meta %>% 
  filter(!is.na(Participant_ID)) %>% 
  group_by(Sample_type) %>% tally() %>% 
  ggplot(aes(x=Sample_type, y=n, fill=Sample_type)) + 
  geom_bar(stat='identity') + theme_bw() + 
  xlab('') + ylab("Number of samples") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=c('#009B76', '#FFBF00'), guide='none')
ggsave(g, filename=here('figures/sample_type/number_of_samples.pdf'),
       width = 3, height = 4, useDingbats=FALSE)

# copies per extraction
df.cop <- df.meta %>% 
  filter(!is.na(Participant_ID)) %>% 
  select(Sample_ID, Participant_ID, Timepoint, Sample_type, DNA_concentration,
         DNA_concentration, copies_16S, copies) %>% 
  group_by(Timepoint, Participant_ID) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  filter(n==2) %>% 
  mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
  select(type, Sample_type, copies_16S, copies, Sample_ID)
df.plot.copies <- df.cop %>% 
  select(-c(copies, Sample_ID)) %>% 
  pivot_wider(names_from = Sample_type, values_from = c(copies_16S)) %>% 
  full_join(df.cop %>% 
              group_by(type) %>% 
              summarise(sample=paste(sort(unique(copies)), 
                                     collapse = ',')),
            by='type')
spearman.rho <- cor(df.plot.copies$next_day,
                    df.plot.copies$same_day, method='spearman')
g <- df.plot.copies %>% 
  ggplot(aes(x=same_day, y=next_day)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_smooth(method='lm', formula='y~x', fill='#ffffff00') +
  geom_point(aes(col=sample), alpha=0.6, pch=16) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab("16S copies from 'same_day' sample") + 
  ylab("16S copies from 'next_day' sample") + 
  annotate(x=-Inf, y=-Inf, geom='text', 
           hjust=-2, vjust=-20,
           label=paste0('rho=', sprintf(fmt='%.2f', spearman.rho))) +
  scale_colour_manual(values=c("#8C1515", "#E98300", "#7F7776"), name='')
ggsave(g, filename=here('figures/sample_type/copies_correlation.pdf'),
       width = 6, height = 4, useDingbats=FALSE)


# 234 random next_day samples
set.seed(2112)
df.rand <- df.meta %>% 
  filter(!Sample_ID %in% df.cop$Sample_ID) %>% 
  filter(!is.na(copies_16S))
df.plot.rand <- tibble(a=df.rand %>% sample_n(234) %>% pull(copies_16S),
                       b=df.rand %>% sample_n(234) %>% pull(copies_16S))
spearman.rho <- cor(df.plot.rand$a,
                    df.plot.rand$b, method='spearman')
g <- df.plot.rand %>% 
  ggplot(aes(x=a, y=b)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(col='#7F7776', alpha=0.6, pch=16) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab("16S copies for 234 random sample") + 
  ylab("16S copies for 234 random sample") + 
  annotate(x=-Inf, y=-Inf, geom='text', 
           hjust=-2, vjust=-20,
           label=paste0('rho=', sprintf(fmt='%.2f', spearman.rho)))
ggsave(g, filename=here('figures/sample_type/copies_correlation_random.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

# correlation
feat.rel <- prop.table(as.matrix(feat.euk.added), 2)
df.sample_type <- df.meta %>% 
  filter(Sample_ID %in% colnames(feat.rel)) %>% 
  group_by(Participant_ID, Timepoint) %>% 
  mutate(n=n()) %>% 
  ungroup() %>%
  filter(n==2)

# log-pearson correlation
# look out, this takes around 20-30 minutes
if (!file.exists('./figures/sample_type/same_sample_correlation.pdf')){
  df.cor <- cor(log10(feat.rel + 1e-05))
  df.cor[lower.tri(df.cor)] <- NA
  diag(df.cor) <- NA
  
  df.cor.plot <- df.cor %>% 
    as_tibble(rownames='Sample_ID') %>% 
    pivot_longer(-Sample_ID, names_to = 'Sample_ID_2', 
                 values_to = 'correlation') %>% 
    filter(!is.na(correlation)) %>% 
    left_join(df.meta %>% 
                mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
                select(Sample_ID, Participant_ID, type), by='Sample_ID') %>% 
    left_join(df.meta %>% 
                mutate(type2=paste0(Participant_ID, '-', Timepoint)) %>% 
                transmute(Sample_ID_2=Sample_ID, 
                          Participant_ID_2=Participant_ID, type2), by='Sample_ID_2') %>% 
    mutate(fill=case_when(type==type2~'same_sample',
                          Participant_ID==Participant_ID_2~'same_patient',
                          TRUE~'other')) %>% 
    filter(!is.na(Participant_ID)) %>% 
    filter(!is.na(Participant_ID_2))
  g <- df.cor.plot %>% 
    ggplot(aes(x=fill, y=correlation, fill=fill)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, fill='white', outlier.shape = NA) + 
    xlab('') + ylab('log-Pearson correlation') + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=c('#7F777695', '#009B7695', '#8C151595'),
                      guide='none')
  ggsave(g, filename='./figures/sample_type/same_sample_correlation.pdf',
         width = 4, height = 4, useDingbats=FALSE)
}

# alpha
alpha <- enframe(vegan::diversity(vegan::rrarefy(
  t(feat.euk.added), 3000), index='shannon'), 
  value ='alpha', name ='Sample_ID')
df.alpha <- df.sample_type %>% 
  left_join(alpha, by='Sample_ID') %>% 
  select(Participant_ID, Timepoint, Sample_type, alpha) %>% 
  mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
  pivot_wider(names_from = Sample_type, values_from = alpha) 
spearman.rho <- cor(df.alpha$same_day,
                    df.alpha$next_day,
                    method='spearman')
g <- df.alpha %>% 
  ggplot(aes(x=same_day, y=next_day)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point(col='#7F7776', alpha=0.8, pch=16) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab("Alpha diversity [same day sample]") + 
  ylab("Alpha diversity [next day sample]") + 
  annotate(x=-Inf, y=-Inf, geom='text', 
           hjust=-2, vjust=-20,
           label=paste0('rho=', sprintf(fmt='%.2f', spearman.rho)))
ggsave(g, filename=here('figures/sample_type/alpha_diversity.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

# beta
df.bc <- vegan::vegdist(vegan::rrarefy(t(feat.euk.added), 3000), method='bray')
df.bc <- as.matrix(df.bc) 
diag(df.bc) <- NA
df.bc[lower.tri(df.bc)] <- NA
df.bc.plot <- df.bc %>% 
  as_tibble(rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, values_to = 'BC_dist', names_to = 'Sample_ID_2') %>% 
  filter(!is.na(BC_dist)) %>% 
  left_join(df.meta %>% select(Sample_ID, Participant_ID, Timepoint) %>% 
              mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
              select(Sample_ID, Participant_ID, type), by='Sample_ID') %>% 
  left_join(df.meta %>% select(Sample_ID, Participant_ID, Timepoint) %>% 
              mutate(type2=paste0(Participant_ID, '-', Timepoint)) %>% 
              transmute(Sample_ID_2=Sample_ID,
                        Participant_ID_2=Participant_ID, type2), 
            by='Sample_ID_2')  %>% 
  mutate(fill=case_when(type==type2~'same_sample',
                        Participant_ID==Participant_ID_2~'same_patient',
                        TRUE~'other')) %>% 
  filter(!is.na(Participant_ID)) %>% 
  filter(!is.na(Participant_ID_2))
g <- df.bc.plot %>% 
  ggplot(aes(x=fill, y=BC_dist, fill=fill)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, fill='white', outlier.shape = NA) + 
  xlab('') + ylab('Bray-Curtis distance') + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=c('#7F777695', '#009B7695', '#8C151595'),
                    guide='none')
ggsave(g, filename=here('figures/sample_type/same_sample_bray.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

# pcoa
df.bc <- vegan::vegdist(vegan::rrarefy(
  t(feat.euk.added[,df.sample_type$Sample_ID]), 3000), 
  method='bray')
pco.res <- labdsv::pco(df.bc)
df.pco <- as_tibble(as.data.frame(pco.res$points), rownames='Sample_ID') %>% 
  left_join(df.sample_type, by='Sample_ID')

df.pco %>% 
  mutate(sample=paste0(Participant_ID, '-', Timepoint)) %>% 
  ggplot(aes(x=V1, y=V2, col=Sample_type)) +
  geom_line(aes(group=sample), col='lightgrey') +
  geom_point() + 
  theme_bw() + theme(panel.grid=element_blank())

df.pco %>% 
  mutate(sample=paste0(Participant_ID, '-', Timepoint)) %>% 
  group_by(sample) %>% 
  summarise(dist=diff(V1) + diff(V2)) %>% 
  arrange(desc(dist))

as_tibble(feat.euk.added[,df.sample_type %>% 
                           filter(Participant_ID=='000014', Timepoint=='42') %>% 
                           pull(Sample_ID)], rownames='species') 
# ggplot(aes(x=log10(P20F04+1e-04), y=log10(P20A09+1e-04))) + 
# ggplot(aes(x=log10(P10D04+1e-04), y=log10(P10B06+1e-04))) +
  # geom_point()

# feat.all[,c('P10D04', 'P10B06')] 

# test
feat.filt <- feat.rel[,df.sample_type$Sample_ID] 
feat.filt <- feat.filt[rowMeans(feat.filt!=0) > 0.05,]
res <- map(rownames(feat.filt), .f=function(x){
  df <- df.sample_type %>% 
    select(Sample_ID, Participant_ID, Sample_type, Timepoint) %>% 
    left_join(enframe(feat.filt[x,], name='Sample_ID', value='rel.ab'),
              by='Sample_ID') %>% 
    mutate(rel.ab=log10(rel.ab+1e-05)) %>% 
    mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
    arrange(type, Sample_type)
  fit <- lmerTest::lmer(rel.ab~Sample_type+(1|type), data=df)
  tmp <- coefficients(summary(fit))
  tibble(species=x, p.val=tmp[2,5], coef=tmp[2,1])
}) %>% bind_rows()
g <- res %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coef, y=-log10(q.val))) + 
  geom_point() + 
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) + 
  xlab('Linear model coefficient') + 
  ylab('-log10(q-value)') + 
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), lty=2) + 
  xlim(-0.2, 0.2)
ggsave(g, file=here('figures/sample_type/volcano.pdf'),
       width = 4, height = 4, useDingbats=FALSE)


# there is only a single one that might be slightly interesting, 
# but not even very strong association
x <- res %>% 
  filter(p.val < 0.001) %>% pull(species)
df <- df.sample_type %>% 
  select(Sample_ID, Participant_ID, Sample_type, Timepoint) %>% 
  left_join(enframe(feat.filt[x,], name='Sample_ID', value='rel.ab'),
            by='Sample_ID') %>% 
  mutate(rel.ab=log10(rel.ab+1e-05)) %>% 
  mutate(type=paste0(Participant_ID, '-', Timepoint)) %>% 
  arrange(type, Sample_type)
cor(df %>% filter(Sample_type=='next_day') %>% pull(rel.ab), 
    df %>% filter(Sample_type=='same_day') %>% pull(rel.ab), 
    method='spearman') # 0.9193312
g <- df %>% 
  select(Sample_type, rel.ab, type) %>% 
  pivot_wider(values_from = rel.ab, names_from = Sample_type) %>% 
  mutate(diff=next_day-same_day) %>% 
  mutate(diff=case_when(diff==0~'no_change', diff < 0 ~ 'decrease', 
                        diff > 0 ~ 'increase')) %>% 
  ggplot(aes(x=same_day, y=next_day)) + 
  geom_abline(slope = 1, intercept = 0, lty=2) +
  geom_point(pch=16) +
  theme_bw() + 
  theme(panel.grid.minor=element_blank()) +
  ggtitle(x)
ggsave(g, file='./figures/sample_type/species_difference.pdf',
       width = 4, height = 4, useDingbats=FALSE)

# ##############################################################################
# some accounting, again

# 2631 samples
# 1 mislabel
# 28 Zymo
# 28 negative control
# 1 EXTRA timepoint (this one is in the same_day samples!)

# 250 same_day
# -----------------
# 2324 samples

keep <- df.filt %>% 
  filter(filter=='keep') %>% 
  pull(Sample_ID)

df.meta.clean <- df.meta %>% 
  filter(!is.na(Participant_ID)) %>% 
  filter(Timepoint!='EXTRA') %>% 
  filter(Sample_type!='same_day') %>% # dim 2324
  filter(Sample_ID %in% keep) # %>% # dim 2309 --> 15 dropped out 

# this makes sense! 42 samples drop out with this filter, most of them are 
# negative controls. 15 in next_day, 1 in same_day, leaves 26 negative controls
# meaning 2 would survive

# filtered out for not enough sequencing
# > 1 million reads and less than 100 in adjusted motus counts

# get the Timepoint into a useable format
df.meta.clean <- df.meta.clean %>% 
  mutate(Timepoint=factor(Timepoint, levels = c(
    'PCON', 'PINF', sort(as.numeric(unique(Timepoint))))))


# ##############################################################################
# Patient accounting

df.outcome <- read_tsv(here('data/meta_participants.tsv'), col_types = cols())

# some encoding
df.outcome <- df.outcome %>% 
  mutate(Relapse_or_death=factor(
    Relapse_or_death, levels=c('censored', 'Relapse', 'Death'))) %>% 
  mutate(cGVHD_or_death=factor(
    cGVHD_or_death, levels=c('censored', 'cGVHD', 'Death'))) %>% 
  mutate(aGVHD24_or_death=factor(
    aGVHD24_or_death, levels=c('censored', 'aGVHD24', 'Death'))) %>% 
  mutate(aGVHD34_or_death=factor(
    aGVHD34_or_death, levels=c('censored', 'aGVHD34', 'Death'))) %>% 
  mutate(cGVHD_ms_or_death=factor(
    cGVHD_ms_or_death, levels=c('censored', 'cGVHD_MS', 'Death')))


# patient info
dim(df.outcome)
# [1] 304  38

# some accounting again
# 324 (158 PTCy vs 147 Tac/MTX) co-enrolled in BMT CTN 1801

# 304 patients with any stool sample
# 157 PTCy vs 147 Tac/MTX
# same as in Mike's table

# ##############################################################################
# Metadata table 

for (x in c('Sex', 'Age', 'Race', 'Ethnicity', 'KLP_score', 'Primary_disease', 
            'CMV_status', 'HLA_matching', 'Conditioning', 'Comorbidity_score', 
            'Disease_risk')){
  if (is.character(df.outcome[[x]])){
    tmp <- df.outcome[,c('Treatment_group', 'Participant_ID', x)] %>% 
      filter(Participant_ID %in% df.meta$Participant_ID) 
    tmp %>% 
      group_by(Treatment_group, !!!syms(x)) %>% 
      tally() %>% 
      group_by(Treatment_group) %>% 
      mutate(n.all=sum(n)) %>% 
      mutate(freq=n/n.all) %>% 
      select(-n.all) %>% print
    fisher.test(tmp[[x]], tmp$Treatment_group) %>% print
  } else if (is.numeric(df.outcome[[x]])){
    tmp <- df.outcome[,c('Treatment_group', 'Participant_ID', x)] %>% 
      filter(Participant_ID %in% df.meta$Participant_ID)
    tmp %>% 
      group_by(Treatment_group) %>% 
      reframe(m=mean(!!!syms(x)), s=sd(!!!syms(x)),
              med=median(!!!syms(x)),
              min=min(!!!syms(x)),
              max=max(!!!syms(x))) %>% print
    wilcox.test(tmp[[x]]~tmp$Treatment_group) %>% print
  }
}

# which sample is closest to the engraftment day?
t.engraftment <- df.outcome %>% 
  select(Participant_ID, Treatment_group, Engraftment, Engraftment_day) %>% 
  full_join(df.meta.clean %>% 
              select(Participant_ID, Timepoint, Sample_ID), 
            by='Participant_ID') %>% 
  filter(Engraftment=='1') %>% 
  filter(Engraftment_day < 40) %>%
  filter(!Timepoint %in% c('PCON', 'PINF', 'EXTRA')) %>% 
  filter(!is.na(Timepoint)) %>% 
  mutate(diff=Engraftment_day - as.numeric(as.character(Timepoint))) %>% 
  group_by(Participant_ID) %>% 
  group_map(.f=function(.x, .y){
    if (any(.x$diff >= -3 & .x$diff < 0)){
      return(.x %>% filter(diff >= -3 & diff < 0) %>% mutate(.y))
    }
    # if (any(.x$diff <= 2 & .x$diff > 0)){
    # return(.x %>% filter(diff <= 3 & diff > 0) %>% mutate(.y))
    # }
    if (any(.x$diff <= 0)){
      return(.x %>% filter(diff <= 0) %>% filter(diff==max(diff)) %>% 
               mutate(.y))
    } else {
      return(tibble(Treatment_group=unique(.x$Treatment_group), 
                    Engraftment_day=unique(.x$Engraftment_day),
                    Sample_ID=NA_character_,
                    Timepoint=NA_character_, diff=NA_real_) %>%
               mutate(.y))
    }}) %>% 
  bind_rows() %>% 
  filter(!is.na(diff)) %>% 
  filter(abs(diff) < 14)

t.engraftment %>% View

# how many participants never nadired in their neutrophil counts?
t.engraftment %>% 
  filter(Engraftment_day==1) %>% 
  dim
# 34
t.engraftment %>% 
  filter(diff <= -9) %>% 
  dim
# 27
# 27/245 -> 0.1102

# fix the outcome table by adding the corresponding Engraftment sample
df.outcome <- df.outcome %>% 
  full_join(t.engraftment %>% 
              select(Sample_ID, Participant_ID) %>% 
              rename(Engraftment_sample=Sample_ID), by='Participant_ID')

# ##############################################################################
# number of samples plots

g <- df.meta.clean %>%
  group_by(Participant_ID) %>% 
  tally() %>% 
  left_join(df.outcome %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  ggplot(aes(x=n, fill=Treatment_group)) + 
    geom_bar() +
    xlab('Number of samples per participant') + 
    ylab('Participant count') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=unlist(colours$group.colours))
ggsave(g, width = 5, height = 4, useDingbats=FALSE,
       filename=here('./figures/general/number_of_samples_per_participant.pdf'))

# ##############################################################################
# export cleaned tables

# motus
# metadata
# patient outcome data
df.response <- df.outcome %>% 
  filter(Participant_ID %in% df.meta.clean$Participant_ID)
feat.motus <- feat.all[,df.meta.clean$Sample_ID]
feat.motus.plus <- feat.euk.added[,df.meta.clean$Sample_ID]
feat.metaphlan <- feat.mpa[,df.meta.clean$Sample_ID]
save(feat.motus,
     feat.motus.plus, 
     feat.metaphlan,
     df.meta.clean, 
     df.response, 
     file=here('data/all_data.RData'))
