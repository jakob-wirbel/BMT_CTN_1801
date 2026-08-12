# ##############################################################################
#
## Domination events
##  Are they associated with treatment groups?
##  Are they associated with outcome in any way?
#
# ##############################################################################

library("tidyverse")
library("here")
library("vegan")
library("survival")
library("ggfortify")

colours <- yaml::read_yaml('./files/colours.yml')

load('./data/all_data.RData')

feat.rel <- prop.table(as.matrix(feat.motus.plus),2)

# ##############################################################################
# Which sample is dominated, and by what?

df.domination <- as_tibble(feat.rel, rownames='species') %>% 
  pivot_longer(-species) %>% 
  filter(value > 0.3) 

g <- df.meta.clean %>% 
  mutate(dominated=Sample_ID %in% df.domination$name) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Timepoint, Treatment_group, dominated) %>% 
  tally() %>% 
  group_by(Timepoint, Treatment_group) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(g=paste0(Treatment_group, dominated)) %>% 
  ggplot(aes(x=Treatment_group, y=n, fill=g)) + 
  geom_bar(stat='identity') +
  facet_grid(~Timepoint) +
  scale_fill_manual(values=c('#ff9f25', '#E98300',
                             '#00bfe0', '#007C92'), guide='none') + 
  xlab('') + ylab('Number of dominated samples') + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank())
ggsave(g, filename=here('figures/domination/number_of_samples.pdf'),
       width = 10, height = 4, useDingbats=FALSE)

g <- df.meta.clean %>% 
  mutate(dominated=Sample_ID %in% df.domination$name) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Timepoint, Treatment_group, dominated) %>% 
  tally() %>% 
  group_by(Timepoint, Treatment_group) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  mutate(g=paste0(Treatment_group, dominated)) %>% 
  filter(dominated) %>% 
  ggplot(aes(x=Timepoint, y=freq, col=Treatment_group)) + 
  geom_point() + 
  geom_line(aes(group=Treatment_group)) +
  scale_colour_manual(values=unlist(colours$group.colours)) +
  xlab('') + ylab('Percentage of dominated samples') + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     axis.ticks.x = element_blank())
ggsave(g, filename=here('figures/domination/percentage_of_samples.pdf'),
       width = 10, height = 4, useDingbats=FALSE)

# export source data
g$data %>% 
  select(Timepoint, Treatment_group, n, n.all, freq) %>% 
  write_tsv('./figures/source_data/Fig3a.tsv')

# test group influence
df.meta.clean %>% 
  mutate(dominated=Sample_ID %in% df.domination$name) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Timepoint, Treatment_group, dominated) %>% 
  tally() %>% 
  group_by(Timepoint) %>% 
  group_map(.f=function(.x, .y){
    if (nrow(.x) == 4){
      tbl <- matrix(.x$n, ncol=2)
      t <- fisher.test(tbl)
      tibble(pval=t$p.value, odds=t$estimate, .y, n=sum(.x$n))
    } else {
      tibble(pval=NA, odds=NA, .y, n=sum(.x$n))
    }
  }) %>% bind_rows() %>% 
  mutate(qval=p.adjust(pval, method='fdr'))
# # A tibble: 19 × 5
# pval  odds Timepoint     n  qval
# <dbl> <dbl> <fct>     <int> <dbl>
# 1 0.151  1.75  PCON        230 0.574
# 2 0.442  1.32  PINF        236 0.763
# 3 0.254  1.37  7           258 0.604
# 4 0.417  0.782 14          224 0.763
# 5 0.777  0.906 21          199 0.869
# 6 0.654  0.850 28          193 0.777
# 7 0.565  1.26  35          111 0.767
# 8 0.0837 2.25  42          104 0.574
# 9 0.402  1.55  49          100 0.763
# 10 0.241  1.94  56          100 0.604
# 11 0.0736 2.36  63           87 0.574
# 12 0.235  1.84  70           91 0.604
# 13 0.151  2.14  77           86 0.574
# 14 1      1.24  84           39 1    
# 15 0.505  1.77  98           80 0.767
# 16 0.536  1.72  180          63 0.767
# 17 0.0946 7.55  270          44 0.574
# 18 1      0.978 365          42 1    
# 19 0.609  2.62  730          22 0.771

# domination by participant
df.dom.pat <- df.meta.clean %>% 
  mutate(dominated=Sample_ID %in% df.domination$name) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  select(Participant_ID, Treatment_group, Timepoint, dominated) %>% 
  arrange(Participant_ID, Timepoint) %>% 
  mutate(dominated=case_when(dominated==FALSE~0, TRUE~1)) %>% 
  group_by(Participant_ID) %>% 
  group_map(.f=function(.x, .y){
    .x %>% mutate(dominated=cumsum(dominated)) %>% 
      mutate(dominated=dominated > 0) %>% 
      mutate(.y)
  }) %>% bind_rows()

df.dom.pat %>% 
  mutate(n.all=length(unique(Participant_ID))) %>%
  group_by(Participant_ID) %>% 
  filter(dominated) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  arrange(Timepoint) %>%
  select(n.all, Timepoint) %>% 
  group_by(Timepoint) %>% 
  reframe(n.all=unique(n.all), s=n()) %>% 
  mutate(cumulative_sum=cumsum(s)) %>% 
  mutate(freq=cumulative_sum/n.all) %>% 
  ggplot(aes(x=Timepoint, y=freq)) + 
    geom_vline(xintercept = 0, lty=2) +
    geom_line(group=1) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    xlab('Day relative to infusion') + 
    ylab('Percentage of patients with domination') + 
    scale_colour_manual(values=unlist(colours$group.colours)) +
    ylim(0,1)

df.dom.pat %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=length(unique(Participant_ID))) %>% 
  group_map(.f=function(.x, .y){
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      filter(dominated) %>% 
      slice_head(n=1) %>% 
      ungroup() %>% 
      arrange(Timepoint) %>%
      select(n.all, Timepoint) %>% 
      group_by(Timepoint) %>% 
      reframe(n.all=unique(n.all), s=n()) %>% 
      mutate(cumulative_sum=cumsum(s)) %>% 
      mutate(freq=cumulative_sum/n.all) 
    tmp %>% mutate(.y)
  }) %>% bind_rows() %>% 
  ggplot(aes(x=Timepoint, y=freq, col=Treatment_group)) + 
  geom_vline(xintercept = 0, lty=2) +
  geom_line(aes(group=Treatment_group)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('Day relative to infusion') + 
  ylab('Percentage of patients with domination') + 
  scale_colour_manual(values=unlist(colours$group.colours)) +
  ylim(0,1)



# which species?
df.dom.meta <- df.meta.clean %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  left_join(df.domination %>% rename(Sample_ID=name), by='Sample_ID') 
named.species <- df.dom.meta %>% 
  group_by(species, Treatment_group) %>% 
  tally() %>% 
  group_by(species) %>% 
  mutate(n.all=sum(n)) %>% 
  arrange(desc(n.all)) %>% 
  filter(!is.na(species)) %>% 
  head(n=20) %>% 
  mutate(species_nice=str_remove(species,' \\[.*\\]')) %>% 
  mutate(species_nice=str_remove(species_nice,' incertae sedis')) %>% 
  mutate(species_nice=case_when(
    species_nice=='Clostridiales species'~'Blautia_A caecimuris',
    TRUE~species_nice))

g <- named.species %>% 
  mutate(species_nice=factor(species_nice, 
                             levels=rev(unique(named.species$species_nice)))) %>% 
  ggplot(aes(x=species_nice, y=n, fill=Treatment_group)) +
  geom_bar(stat='identity', position = position_dodge()) +
  coord_flip() + xlab('') + ylab('Number of dominated samples') + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.y = element_blank()) + 
  scale_fill_manual(values=unlist(colours$group.colours))
ggsave(g, filename=here('figures/domination/dominating_species_filled.pdf'),
       width = 8, height = 4, useDingbats=FALSE)

# export source data
g$data %>% 
  select(species, species_nice, Treatment_group, n, n.all) %>% 
  write_tsv('./figures/source_data/Fig3b.tsv')

pdf(here('figures/domination/single_species.pdf'),
    width = 8, height = 4, useDingbats = FALSE)
df.rates <- list()
p.vals <- rep(NA, length.out=length(unique(named.species$species)))
names(p.vals) <- unique(named.species$species)
for (x in unique(named.species$species)){
  g <- enframe(feat.rel[x,] > 0.3, name='Sample_ID') %>% 
    full_join(df.meta.clean %>% 
                select(Sample_ID, Timepoint, Participant_ID),
              by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') %>% 
    group_by(Timepoint, Treatment_group, value) %>% 
    tally() %>% 
    ungroup() %>% 
    pivot_wider(values_from = n, names_from = value, values_fill = 0) %>% 
    mutate(all=`FALSE` + `TRUE`) %>% 
    mutate(freq=`TRUE`/all) %>% 
    ggplot(aes(x=Timepoint, y=freq, col=Treatment_group)) + 
    geom_point() +
    xlab('') + ylab('Percentage of sampled dominated') + 
    ggtitle(x) +
    geom_line(aes(group=Treatment_group)) + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.y = element_blank()) + 
    scale_colour_manual(values=unlist(colours$group.colours))
  print(g)
  df.rates[[x]] <- enframe(feat.rel[x,] > 0.3, name='Sample_ID') %>% 
    full_join(df.meta.clean %>% 
                select(Sample_ID, Timepoint, Participant_ID),
              by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') %>% 
    group_by(Timepoint, value) %>% 
    tally() %>% 
    group_by(Timepoint) %>% 
    mutate(n.all=sum(n)) %>% 
    mutate(freq=n/n.all) %>% 
    filter(value) %>% 
    mutate(species=x)
  # any of them are enriched per group?
  df.test <- enframe(feat.rel[x,], name='Sample_ID') %>% 
    mutate(dominated=value > 0.3) %>% 
    left_join(df.meta.clean %>% select(Sample_ID, Participant_ID), 
              by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') 
  t <- fisher.test(table(df.test$dominated, df.test$Treatment_group))
  p.vals[x] <- t$p.value
}
dev.off()
p.adjust(p.vals, method='BH') %>% sort
# Clostridiales species incertae sedis [meta_mOTU_v3_13012]                 
# 0.0009805086
# [Ruminococcus] gnavus [ref_mOTU_v3_01594] 
# 0.3112238667 
# Enterococcus faecium [ref_mOTU_v3_00321]   
# 0.5397189790
# Akkermansia species incertae sedis [meta_mOTU_v3_12805] 
# 0.6736207935 
# Streptococcus thermophilus [ref_mOTU_v3_01348]               
# 0.8270733218
# Akkermansia muciniphila [ref_mOTU_v3_03591] 
# 0.8774740614 
# Escherichia coli [ref_mOTU_v3_00095]
# 0.8774740614
# Enterococcus faecalis [ref_mOTU_v3_00318] 
# 0.8774740614 
# Klebsiella pneumoniae [ref_mOTU_v3_00085]               
# 1.0000000000
# Lactobacillus rhamnosus [ref_mOTU_v3_00710] 
# 1.0000000000 

# are they actually dominating or only 'the last ones standing'
species.key <- named.species$species_nice
names(species.key) <- named.species$species
g <- as_tibble(t(feat.rel[unique(named.species$species),]), 
               rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, names_to = 'species') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(species %in% named.species$species) %>% 
  mutate(species_nice=species.key[species]) %>% 
  filter(value!=0) %>%
  ggplot(aes(x=log10(value+1e-05), col=species_nice)) + 
  # geom_histogram(bins=20) + 
  geom_density()+
  geom_vline(xintercept = log10(0.3)) +
  theme_bw() +
  scale_colour_viridis_d() +
  # coord_cartesian(ylim=c(0,250)) + 
  NULL
ggsave(g, filename=here('figures/domination/abundance_distribution.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

# export source data
# this is a bit annoying, since we need to access the underlying ggplot2 stuff

data <- g@data %>% 
  select(species_nice, value) %>% 
  mutate(weight=1, value=log10(value+1e-05))
range <- c(-4.4, 0)
map(unique(data$species_nice), .f = function(x){
  tmp <- data %>% 
    filter(species_nice==x)
  density <- ggplot2:::compute_density(tmp$value, tmp$weight, from = range[1], 
                                       to = range[2], bw = 'nrd0', adjust = 1, 
                                       kernel = "gaussian", 
                                       n = 100, bounds = c(-Inf, Inf))
  tibble(value=density$x, density=density$density, species=x)
}) %>% bind_rows() %>% write_tsv('./figures/source_data/Fig3c.tsv')



as_tibble(t(feat.rel[unique(named.species$species),]), 
          rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, names_to = 'species') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(species %in% named.species$species) %>% 
  mutate(species_nice=species.key[species]) %>% 
  filter(value!=0, value < 0.3) %>%
  group_by(species) %>%
  reframe(m=mean(value))
# it's complicated :D 
# some are just randomly a bit higher than 0.3 (R. gnavus, Akkermansia, Blautia)
# how about alpha diversity as well?

g <- as_tibble(t(feat.rel[unique(named.species$species),]), 
               rownames='Sample_ID') %>% 
  pivot_longer(-Sample_ID, names_to = 'species') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(species %in% named.species$species) %>% 
  mutate(species_nice=species.key[species]) %>% 
  filter(value!=0) %>%
  mutate(dom=value > 0.3) %>% 
  ggplot(aes(x=value, y=copies_16S)) + 
    geom_point(aes(alpha=dom), pch=16) +
    facet_wrap(~species_nice) +
    geom_smooth(data=. %>% filter(value > 0.3), method='lm') +
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    xlab('Relative abundance') + ylab('Sample log10(16S copy number)') + 
    geom_vline(xintercept = 0.3, lty=2) + 
    scale_alpha_manual(values=c(0.4, 1))
ggsave(g, filename=here('figures/domination/rel_ab_vs_copies.pdf'),
       width = 7, height = 7, useDingbats=FALSE)

# export source data
g@data %>% 
  select(Sample_ID, species_nice, value, copies_16S, dom) %>% 
  write_tsv('./figures/source_data/EDFig6a.tsv')

# this tells me that domination does not mean domination: some species have low
# absolute abundance, when dominating, others have high abundance even when
# they are the only one left
#
# e.g. R. gnavus is barely ever over 30%, same as Blautia 
# Akkermansia is over the whole relative abundance space, but lower absolute 
# abundance when it's the only one left
# E.coli or E. faecium are high even when they are dominating
#
# can we re-define domination? Either based on alpha diversity or on absolute
# abundance?
as_tibble(feat.rel, rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  filter(value!=0) %>%
  left_join(df.meta.clean %>% 
              select(Sample_ID, copies_16S), by='Sample_ID') %>% 
  mutate(value_ab=value*copies_16S) %>% 
  ggplot(aes(x=value, y=value_ab)) + 
  geom_point() 
# not really a clear cutoff


# very very strong domination is basically only E. faecium and K. pneumoniae

# okay, what to do now?
# the main points should be: domination is not really domination, if it's
# a species that is anyway super high in relative abundance (Blautia)
# second is: is it domination or relative abundance effect?

alpha <- enframe(vegan::diversity(vegan::rrarefy(t(feat.motus.plus), 3000), 
                                  index='shannon'), 
                 value ='alpha', name ='Sample_ID')

tmp <- as_tibble(feat.rel, rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  filter(value > 0) %>% 
  group_by(Sample_ID) %>% 
  mutate(dominated=any(value > 0.3)) %>% 
  group_map(.f=function(.x, .y){
    if (!unique(.x$dominated)){
      tibble(.y, dominated=FALSE, species='not_dominated')
    } else {
      .xred <- .x %>% 
        filter(value > 0.3)
      if (any(.xred$species %in% named.species$species)){
        if (nrow(.xred)==1){
          tibble(.y, dominated=TRUE, species=.xred$species)
        } else {
          .xred <- .xred %>% 
            filter(species %in% named.species$species)
          if (nrow(.xred) > 1){
            tibble(.y, dominated=TRUE, 
                   species=paste(.xred$species, collapse=';'))
          } else {
            tibble(.y, dominated=TRUE, species=paste0(.xred$species, '+'))
          }
        }
      } else {
        tibble(.y, dominated=TRUE, species='other')
      }
    }
  }) %>% 
  bind_rows() %>% 
  mutate(other=str_detect(species, '\\+')) %>% 
  mutate(species=str_remove(species, '\\+')) %>% 
  mutate(species=str_split(species, ';')) %>% 
  rowwise() %>% 
  mutate(number=length(species)) %>% 
  mutate(other=case_when(number > 1~TRUE, TRUE~other)) %>% 
  unnest(species) %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  select(Sample_ID, dominated, species, other, 
         copies_16S, copies, DNA_concentration) %>% 
  left_join(alpha, by='Sample_ID')

# n
tmp %>%
  group_by(species) %>%
  tally()


g <- tmp %>% 
  group_by(species) %>% 
  mutate(m_copies=median(copies_16S)) %>% 
  arrange(desc(m_copies)) %>% 
  mutate(species=factor(species, levels=unique(.$species))) %>% 
  ggplot(aes(x=species, y=copies_16S)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('') + ylab('value') + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(g, filename=here('./figures/domination/copies_16S.pdf'),
       width = 8, height = 5, useDingbats=FALSE)

# export source data
g@data %>% 
  select(Sample_ID, species, copies_16S) %>% 
  write_tsv('./figures/source_data/Fig3d.tsv')

# actual values
tmp %>% 
  group_by(species) %>% 
  reframe(m_copies=median(copies_16S)) %>% 
  arrange(desc(m_copies))

# test difference to non-dominated samples
pvals <- rep(NA, length.out=length(unique(named.species$species)))
names(pvals) <- unique(named.species$species)
for (i in unique(named.species$species)){
  t <- wilcox.test(copies_16S~species, data=tmp %>%
                     filter(species %in% c('not_dominated', i)))
  pvals[i] <- t$p.value
}
p.adjust(pvals, method='BH')
# Enterococcus faecium [ref_mOTU_v3_00321]               
# 9.458013e-19
# Akkermansia muciniphila [ref_mOTU_v3_03591] 
# 5.688451e-04 
# Clostridiales species incertae sedis [meta_mOTU_v3_13012]   
# 9.424122e-01 
# Akkermansia species incertae sedis [meta_mOTU_v3_12805] 
# 1.549430e-04
# Escherichia coli [ref_mOTU_v3_00095]                 
# 8.452388e-09
# Klebsiella pneumoniae [ref_mOTU_v3_00085] 
# 2.705617e-08 
# Enterococcus faecalis [ref_mOTU_v3_00318]
# 4.891122e-06
# Streptococcus thermophilus [ref_mOTU_v3_01348] 
# 3.489334e-13 
# Lactobacillus rhamnosus [ref_mOTU_v3_00710]                  
# 1.664637e-06
# [Ruminococcus] gnavus [ref_mOTU_v3_01594]
# 4.993669e-02 

# co-occurrence?
pheatmap::pheatmap(cor(t(feat.rel[unique(named.species$species),]), method='spearman'))
pheatmap::pheatmap(cor(log10(t(feat.rel[unique(named.species$species),])+1e-05)))

# no

# how about we look only at pathobionts?
pathos <- unique(named.species$species)[c(1, 5, 6, 8)]

as_tibble(feat.rel[pathos,], rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  reframe(dominated=any(value > 0.3)) %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  group_by(Timepoint, dominated) %>% 
  tally() %>% 
  group_by(Timepoint) %>% 
  mutate(n.all=sum(n)) %>% 
  filter(dominated) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=Timepoint, y=freq, group=1)) + 
  geom_line()

g <- as_tibble(feat.rel[pathos,], rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  filter(value > 0.3) %>% 
  full_join(df.meta.clean, by='Sample_ID') %>% 
  select(species, Timepoint, Sample_ID) %>% 
  group_by(Timepoint) %>% 
  mutate(n.all=length(unique(Sample_ID))) %>% 
  group_by(species, Timepoint) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  mutate(freq=n/n.all) %>% 
  filter(!is.na(species)) %>% 
  select(species, Timepoint, freq) %>% 
  distinct() %>% 
  pivot_wider(names_from = species, values_from = freq, values_fill = 0) %>% 
  pivot_longer(-Timepoint, names_to = 'species', values_to = 'freq') %>% 
  ggplot(aes(x=Timepoint, y=freq, group=species, colour=species)) + 
  geom_line() +
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('') + ylab('Percentage of samples')
ggsave(g, filename=here('figures/domination/freq_domination_pathogens.pdf'),
       width = 7, height = 4, useDingbats=FALSE)

# export source data
g@data %>% 
  arrange(species, Timepoint) %>% 
  write_tsv('./figures/source_data/EDFig6b.tsv')

# association with outcome?
source(here('src/utils.r'))
# any dominated at day 14
df.pred.14 <- as_tibble(feat.rel[,df.meta.clean %>% 
                                   filter(Timepoint=='14') %>% 
                                   pull(Sample_ID)], 
                        rownames='species') %>% 
  pivot_longer(-species, names_to='Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  reframe(m=as.numeric(any(value > 0.3))) %>% 
  left_join(df.meta.clean %>% select(Sample_ID, Participant_ID),
            by='Sample_ID') %>% 
  left_join(df.response, by='Participant_ID') %>% 
  mutate(m=as.logical(m))
assoc.14 <- .f_test_outcome(df.pred.14, landmarked = 14)

# any dominated at day 21
df.pred.21 <- as_tibble(feat.rel[,df.meta.clean %>% 
                                   filter(Timepoint=='21') %>% 
                                   pull(Sample_ID)], 
                        rownames='species') %>% 
  pivot_longer(-species, names_to='Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  reframe(m=as.numeric(any(value > 0.3))) %>% 
  left_join(df.meta.clean %>% select(Sample_ID, Participant_ID),
            by='Sample_ID') %>% 
  left_join(df.response, by='Participant_ID') %>% 
  mutate(m=as.logical(m))
assoc.21 <- .f_test_outcome(df.pred.21, landmarked = 21)

# domination with enteric pathobionts at day 21
df.pred.21.p <- as_tibble(feat.rel[pathos,df.meta.clean %>% 
                                     filter(Timepoint=='21') %>% 
                                     pull(Sample_ID)], 
                          rownames='species') %>% 
  pivot_longer(-species, names_to='Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  reframe(m=as.numeric(any(value > 0.3))) %>% 
  left_join(df.meta.clean %>% select(Sample_ID, Participant_ID),
            by='Sample_ID') %>% 
  left_join(df.response, by='Participant_ID') %>% 
  mutate(m=as.logical(m))
assoc.21p <- .f_test_outcome(df.pred.21.p, landmarked = 21)

df.hr <- assoc.14$associations %>% 
  filter(set=='all') %>% 
  mutate(type='day14') %>% 
  bind_rows(assoc.21$associations %>% 
              filter(set=='all') %>% 
              mutate(type='day21')) %>% 
  bind_rows(assoc.21p$associations %>% 
              filter(set=='all') %>% 
              mutate(type='day21p'))
  

g <- df.hr %>% 
  ggplot(aes(x=outcome, y=estimate, fill=type)) + 
  geom_hline(yintercept = 1) +
  geom_linerange(aes(ymin=low, ymax=high, col=type), 
                 position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2), pch=23, size=1.5) +
  coord_flip() + 
  theme_bw() + 
  ylab('Hazard ratio') + xlab('') + 
  scale_fill_manual(values = c('#F98C0AFF', '#BB3754FF', '#56106EFF')) + 
  scale_colour_manual(values = c('#F98C0AFF', '#BB3754FF', '#56106EFF'))
ggsave(g, filename=here('figures/domination/domination_associations.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
# nothing here is significant! (except for one tiny thing, but who knows with
# multiple testing)

df.hr %>% 
  mutate(qval=p.adjust(pval, method='BH')) %>% 
  arrange(pval)

# this might be chance?

# export source data
df.hr %>% 
  write_tsv('./figures/source_data/Fig3e.tsv')
