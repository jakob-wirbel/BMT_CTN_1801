# ##############################################################################
#
## Explore differences in the microbiome composition
#
# ##############################################################################

library("tidyverse")
library("here")

load('./data/all_data.RData')
colours <- yaml::read_yaml('./files/colours.yml')

feat.rel <- prop.table(as.matrix(feat.motus), 2)

# function to test
.f_test <- function(.x, .y){
  if (mean(.x$value!=log10(5e-05)) > 0.1){
    fit <- lmerTest::lmer(value~Treatment_group + drug + (1|Participant_ID), 
                          data=.x)
    coefs <- coefficients(summary(fit))
    fit.abs <- lmerTest::lmer(log10(value.abs+1)~Treatment_group + 
                                (1|Participant_ID), data=.x)
    coefs.abs <- coefficients(summary(fit.abs))
    tibble(pval=coefs[2,5], coef=coefs[2,1],
           pval.abs=coefs.abs[2,5], coef.abs=coefs.abs[2,1], .y)
  } else {
    tibble(pval=NA, coef=NA, pval.abs=NA, coef.abs=NA, .y)
  }
}

.f_test_single <- function(.x, .y){
  # browser()
  if (mean(.x$value!=-4) > 0.1){
    fit <- lm(value~Treatment_group + drug, data=.x)
    coefs <- coefficients(summary(fit))
    fit.abs <- lm(log10(value.abs+1)~Treatment_group + drug, data=.x)
    coefs.abs <- coefficients(summary(fit.abs))
    tibble(pval=coefs[2,4], coef=coefs[2,1],
           pval.abs=coefs.abs[2,4], coef.abs=coefs.abs[2,1], .y)
  } else {
    tibble(pval=NA, coef=NA, pval.abs=NA, coef.abs=NA, .y)
  }
}

# get drug exposure data as well
df.abx.exp <- read_tsv('./files/abx_exposure.tsv')

# ##############################################################################
# Differences between arms at baseline

meta.baseline <- df.meta.clean %>% 
  filter(Timepoint=='PCON')
feat.baseline <- feat.rel[,meta.baseline$Sample_ID]
feat.baseline <- feat.baseline[rowMeans(feat.baseline!=0) > 0.1,]

# volcano
df.test <- as_tibble(feat.baseline,
                     rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(meta.baseline, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
  mutate(value.abs=value*10^copies_16S) %>% 
  mutate(value=log10(value + 5e-05)) %>% 
  group_by(species) %>% 
  group_map(.f=.f_test_single) %>% bind_rows() %>% 
  filter(!is.na(pval))

# baseline differences?
df.test %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  mutate(label=case_when(qval < 0.01~species, TRUE~'')) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) +
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  ggrepel::geom_text_repel(aes(label=label), size=3) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  geom_hline(yintercept = -log10(0.01), lty=2) + 
  ggtitle('Baseline differences') 
# none!
# same for absolute abundance
g <- df.test %>% 
  filter(species!='unassigned') %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  select(species, coef, coef.abs, qval, qval.abs) %>% 
  pivot_longer(-species) %>% 
  mutate(abs=case_when(str_detect(name, 'abs')~'Absolute abundance', 
                       TRUE~'Relative abundance')) %>% 
  mutate(name=str_remove(name, '.abs')) %>% 
  pivot_wider(names_from=name, values_from = value) %>% 
  mutate(label=case_when(qval < 0.05~species, TRUE~'')) %>% 
  mutate(label=str_remove(label, ' \\[.*\\]$')) %>% 
  mutate(highlight=qval < 0.05) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) + 
  geom_hline(yintercept = -log10(c(0.05, 0.01)), colour='darkgrey', lty=3) +
  geom_point(aes(col=highlight)) + 
  facet_grid(~abs, scales = 'free_x') + 
  ggrepel::geom_text_repel(aes(label=label), size=2) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     strip.background = element_blank()) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  scale_colour_manual(values=c('darkgrey', '#8C1515'), guide='none')
ggsave(g, filename=here('figures/composition/volcano_differences_baseline.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# export source data
df.test %>% 
  select(pval, coef, species) %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  write_tsv('./figures/source_data/EDFig7a.tsv')

# ##############################################################################
# Differences between arms overall

# remove all samples from single patients
feat.test <- feat.rel[,df.meta.clean$Sample_ID]
feat.test <- feat.test[rowMeans(feat.test!=0) > 0.1,]

# volcano
df.test <- as_tibble(feat.test,
                     rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  mutate(value.abs=value*10^copies_16S) %>% 
  mutate(value.abs=log10(value.abs+1000)) %>% 
  mutate(value=log10(value + 5e-05)) %>% 
  group_by(species) %>% 
  group_map(.f=.f_test) %>% bind_rows() %>% 
  filter(!is.na(pval))

# baseline differences?
df.test %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  mutate(label=case_when(qval < 0.01~species, TRUE~'')) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) +
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  ggrepel::geom_text_repel(aes(label=label), size=3) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  geom_hline(yintercept = -log10(0.01), lty=2) + 
  ggtitle('Baseline differences') 
# none!
# same for absolute abundance, only Streptococcus mutans
g <- df.test %>% 
  filter(species!='unassigned') %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  select(species, coef, coef.abs, qval, qval.abs) %>% 
  pivot_longer(-species) %>% 
  mutate(abs=case_when(str_detect(name, 'abs')~'Absolute abundance', 
                       TRUE~'Relative abundance')) %>% 
  mutate(name=str_remove(name, '.abs')) %>% 
  pivot_wider(names_from=name, values_from = value) %>% 
  mutate(label=case_when(qval < 0.05~species, TRUE~'')) %>% 
  mutate(label=str_remove(label, ' \\[.*\\]$')) %>% 
  mutate(highlight=qval < 0.05) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) + 
  geom_hline(yintercept = -log10(c(0.05, 0.01)), colour='darkgrey', lty=3) +
  geom_point(aes(col=highlight)) + 
  facet_grid(~abs, scales = 'free_x') + 
  ggrepel::geom_text_repel(aes(label=label), size=2) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     strip.background = element_blank()) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  scale_colour_manual(values=c('darkgrey', '#8C1515'), guide='none')
ggsave(g, filename=here('figures/composition/volcano_differences_overall.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# export source data
df.test %>% 
  select(pval, coef, species) %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  write_tsv('./figures/source_data/EDFig7b.tsv')

df.sm <- enframe(feat.test['Streptococcus mutans [ref_mOTU_v3_01605]',],
                 name = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') 
df.sm %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=log10(value+5e-05),
             col=Treatment_group)) + 
  geom_point() + 
  geom_smooth(method='loess')
df.sm %>% 
  mutate(value.ab = value*10^copies_16S) %>% 
  mutate(value.ab=log10(value.ab+1e3)) %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=value.ab,
             col=Treatment_group)) + 
  geom_point() + 
  geom_smooth(method='loess')
df.sm %>% 
  mutate(present=value > 0) %>% 
  group_by(Treatment_group, Timepoint, present) %>% 
  tally() %>% 
  group_by(Treatment_group, Timepoint) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  filter(present) %>% 
  ggplot(aes(x=Timepoint, y=freq, col=Treatment_group, group=Treatment_group)) + 
  geom_line()


# ##############################################################################
# Differences between arms at days 14-28 post-HCT

df.test <- as_tibble(feat.rel,
                     rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  mutate(value.abs=value*10^copies_16S) %>% 
  mutate(value.abs=log10(value.abs+1e3)) %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(species) %>% 
  mutate(prev=mean(value!=0)) %>% 
  filter(prev > 0.1) %>% 
  mutate(value=log10(value + 5e-05)) %>% 
  left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
  group_map(.f=.f_test) %>% bind_rows() %>% 
  filter(!is.na(pval))

df.test %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  mutate(label=case_when(qval < 0.01~species, TRUE~'')) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) +
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  ggrepel::geom_text_repel(aes(label=label), size=3) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  geom_hline(yintercept = -log10(0.01), lty=2)

g <- df.test %>% 
  filter(species!='unassigned') %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  mutate(qval.abs=p.adjust(pval.abs, method='fdr')) %>% 
  select(species, coef, coef.abs, qval, qval.abs) %>% 
  pivot_longer(-species) %>% 
  mutate(abs=case_when(str_detect(name, 'abs')~'Absolute abundance', 
                       TRUE~'Relative abundance')) %>% 
  mutate(name=str_remove(name, '.abs')) %>% 
  pivot_wider(names_from=name, values_from = value) %>% 
  mutate(label=case_when(qval < 0.05~species, TRUE~'')) %>% 
  mutate(label=str_remove(label, ' \\[.*\\]$')) %>% 
  mutate(highlight=qval < 0.05) %>% 
  ggplot(aes(x=coef, y=-log10(qval))) + 
  geom_hline(yintercept = -log10(0.05), colour='darkgrey', lty=3) +
  geom_point(aes(col=highlight)) + 
  facet_grid(~abs, scales = 'free_x') + 
  ggrepel::geom_text_repel(aes(label=label), size=2) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     strip.background = element_blank()) + 
  xlab('Linear model effect size') + ylab('-log10(Q-value)') + 
  scale_colour_manual(values=c('darkgrey', '#8C1515'), guide='none')
ggsave(g, filename=here('figures/composition/volcano_differences_7_28.pdf'),
       width = 7, height = 4, useDingbats=FALSE)

# export
df.test %>% 
  filter(species!='unassigned') %>% 
  mutate(qval=p.adjust(pval, method='fdr')) %>% 
  select(species, coef, pval, qval) %>% 
  arrange(coef) %>% 
  write_tsv('./files/enrichment_results.tsv')

# it's the same with Metaphlan, so that's good/encouraging

# ##############################################################################
# associations with outcome for C. scindens (top hit)
df.prediction <- as_tibble(feat.rel[626,,drop=FALSE],
          rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(Participant_ID) %>% 
  reframe(m=mean(value > 0) > 0.5) %>% 
  full_join(df.response, by='Participant_ID') %>% 
  filter(!is.na(m))
source('./src/utils.r')
df.res <- .f_test_outcome(df.prediction, landmarked = 28)
df.res$associations %>% filter(set=='all') %>% 
  mutate(eff=paste0(sprintf(fmt='%.3f', estimate), ' (', 
                    sprintf('%.3f', low), ', ', 
                    sprintf('%.3f', high), ')')) %>% 
  mutate(pval=sprintf('%.3f', pval))  %>% 
  select(outcome, pval, eff) %>% 
  mutate(outcome=factor(
    outcome, levels = c('NRM', 'Relapse', 'aGVHD24', 'aGVHD34', 
                        'cGVHD', 'cGVHD_MS', 'OS', 'GRFS'))) %>% 
  arrange(outcome)

# prevalence of C. scindens over time
g <- as_tibble(feat.rel[626,,drop=FALSE],
          rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  group_by(Timepoint, Treatment_group) %>% 
  reframe(prevalence=mean(value>0)) %>% 
  ggplot(aes(x=Timepoint, y=prevalence, colour=Treatment_group)) + 
    geom_point() + 
    geom_line(aes(group=Treatment_group)) + 
    scale_colour_manual(values=unlist(colours$group.colours)) +
    theme_bw() + theme(panel.grid.minor = element_blank())
ggsave(g, filename='./figures/composition/c_scindens_prevalence.pdf',
       width = 6, height = 3, useDingbats=FALSE)

# export source data
g@data %>% 
  write_tsv('./figures/source_data/Fig4b.tsv')

# ##############################################################################
# let's make this cute with the phylum?
df.motus.tax <- read_tsv('./files/mOTUs_3.0.0_GTDB_tax.tsv',
                         col_names = c('motus_ID', 'domain', 'phylum', 'class', 
                                       'order', 'family', 'genus', 'species'),
                         col_types = cols())

g <- df.test %>% 
  mutate(motus_ID=str_extract(species, '(ref|meta|ext)_mOTU_v3_[0-9]{5}')) %>% 
  left_join(df.motus.tax %>% select(motus_ID, class), by='motus_ID') %>%
  mutate(qval=p.adjust(pval, method='BH')) %>% 
  group_by(class) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  mutate(class=case_when(n < 15~'other', TRUE~class)) %>% 
  mutate(highlight=qval < 0.05) %>% 
  mutate(label=case_when(qval < 0.05~species, TRUE~'')) %>% 
  mutate(label=str_remove(label, ' \\[.*\\]$')) %>% 
  ggplot(aes(x=coef, y=-log10(qval), col=class)) + 
    geom_hline(yintercept = -log10(0.05), colour='darkgrey', lty=3) +
    geom_point(aes(alpha=highlight), pch=16) + 
    scale_alpha_manual(values=c(1, 1), guide='none') +
    ggrepel::geom_text_repel(aes(label=label), size=2) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       strip.background = element_blank()) +
    scale_colour_manual(values=c('#CC4678FF', '#E16462FF', '#900DA4FF', 
                                 '#FCCE25FF', '#42049EFF', 'grey'), name='Class')
  
    # scale_colour_manual(values=c('#8F993E', '#279989', '#175E54', '#E04F39', 
                                 # '#6FA287', 'darkgrey'), name='Class')
ggsave(g, filename=here('figures/composition/volcano_differences_7_28_class.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# export source data
g@data %>% 
  select(pval, coef, qval, species, class) %>% 
  write_tsv('./figures/source_data/Fig4a.tsv')

# check the actual genome names in the GTDB table?
df.genome.metadata <- read_tsv('./files/mOTUs3.1.0.genome_metadata.tsv', 
                               col_types = cols()) %>% 
  select(GENOME_SOURCE, ENVIRONMENT)

df.test %>% 
  mutate(qval=p.adjust(pval, method='BH')) %>% 
  filter(qval < 0.05) %>% 
  mutate(GENOME_SOURCE=str_extract(species, 
                                   '(ref|meta|ext)_mOTU_v3_[0-9]{5}')) %>% 
  mutate(GENOME_SOURCE=str_replace(GENOME_SOURCE, 'v3_', 'v31_')) %>% 
  left_join(df.genome.metadata, by='GENOME_SOURCE') %>% 
  select(coef, qval, species, GENOME_SOURCE, ENVIRONMENT) %>% 
  group_by(coef, qval, species, GENOME_SOURCE, ENVIRONMENT) %>% 
  tally() %>% 
  group_by(GENOME_SOURCE) %>% 
  arrange(GENOME_SOURCE, desc(n)) %>% View

# abundance distribution for all the ones with differential abundance
g <- as_tibble(feat.rel[df.test %>% mutate(qval=p.adjust(pval, method='BH')) %>% 
                     filter(qval < 0.05) %>% pull(species),], 
          rownames='species') %>% 
  pivot_longer(-species, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  ggplot(aes(x=Treatment_group, y=log10(value + 5e-05))) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.1) +
    facet_wrap(~species)
ggsave(g, filename=here('figures/composition/boxplot_7_28.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# export source data
g@data %>% 
  select(species, Sample_ID, value, Timepoint, 
         Participant_ID, Treatment_group) %>% 
  write_tsv('./figures/source_data/Fig4c.tsv')
