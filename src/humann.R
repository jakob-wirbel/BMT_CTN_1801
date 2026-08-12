# ##############################################################################
#
## Functional enrichment data
#
# ##############################################################################

library("tidyverse")
library("here")

source(here('src/utils.r'))

colours <- yaml::read_yaml('./files/colours.yml')

# metadata
load('./data/all_data.RData')
df.abx.exp <- read_tsv('./files/abx_exposure.tsv')

df.humann.raw <- read_tsv(here('data/humann.tsv'), col_types = cols())

# the detailed, species-resolved pathway abundances do not sum up to the
# community-level pathway abundances (and apparently, this is fine and expected
# https://github.com/biobakery/humann?tab=readme-ov-file#3-pathway-abundance-file)

# UNMAPPED/UNINTEGRATED make up most of the abundance
# maybe we should filter them out for convenience?
# I guess it should be relatively consistent/constant across samples
# maybe we should check this hunch, though

# how should we deal with the species-level contributions to each pathway?
# we could get relative contributions for pathways of interest? 
# we could check contributions for specific species of interest?

# df.humann.raw <- map(all.files, .f = function(x){
#   sample <- str_split(x, '_')[[1]][1]
#   tmp <- read_tsv(paste0(humann.loc, x), skip=1, 
#                   col_names = c('pathway', 'abundance'),
#                   col_types=cols()) %>% 
#     mutate(sample=sample)}) %>% 
#   bind_rows() %>% 
#   pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)
# 
# write_tsv(df.humann.raw, file = './data/humann_all.tsv')
# this takes too long...
# we can combine everything on the cluster and just load the combined
# abundance file

# check the rate of UNMAPPED/UNINTEGRATED
g <- df.humann.raw %>% 
  filter(!str_detect(pathway, '\\|')) %>% 
  pivot_longer(-pathway) %>% 
  mutate(pathway=case_when(pathway %in% c('UNMAPPED', 'UNINTEGRATED')~pathway,
                           TRUE~'INTEGRATED')) %>% 
  group_by(name, pathway) %>%
  reframe(value=sum(value)) %>% 
  group_by(name) %>% 
  mutate(rel_value=value/sum(value)) %>% 
   ggplot(aes(x=pathway, y=rel_value, fill=pathway)) + 
    geom_boxplot() + 
    theme_bw() + 
    xlab('') + ylab('Relative abundance')
ggsave(g, filename=here('figures/humann/humann_basic_stats.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

# most things are not integrated or unmapped... nothing to do about this, tbh

df.humann <- df.humann.raw %>% 
  filter(!str_detect(pathway, '\\|')) %>% 
  filter(!pathway%in%c('UNMAPPED', 'UNINTEGRATED')) %>% 
  pivot_longer(-pathway) %>% 
  group_by(name) %>% 
  mutate(rel_value=value/sum(value)) %>% 
  select(-value) %>% 
  pivot_wider(names_from = name, values_from = rel_value)

df.humann <- df.humann %>% 
  as.data.frame()
rownames(df.humann) <- df.humann$pathway
df.humann$pathway <- NULL
df.humann <- as.matrix(df.humann)


df.humann <- df.humann[,df.meta.clean$Sample_ID]

# general stuff
hist(log10(df.humann), 100) # log.n0 <- 1e-06
hist(rowMeans(df.humann!=0), 100) # 5%? 

df.humann <- df.humann[rowMeans(df.humann!=0) > 0.05,]

# PCoA
df.pco <- labdsv::pco(vegan::vegdist(t(log10(df.humann + 1e-07)), 
                                     method='euclidean'))
head(df.pco$eig/sum(df.pco$eig[df.pco$eig > 0]))

as_tibble(as.data.frame(df.pco$points), rownames='Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  ggplot(aes(x=V1, y=V2, col=Treatment_group)) + 
    geom_point()
# nothing, really

# ##############################################################################
# differential abundance?

# relative abundance difference between arms at baseline?
df.baseline <- df.humann[,df.meta.clean %>% 
                           filter(Timepoint=='PCON') %>% 
                           pull(Sample_ID)]
df.baseline <- df.baseline[rowMeans(df.baseline!=0) > 0.05,]
df.res <- map(rownames(df.baseline), .f=function(x){
  tmp <- enframe(df.baseline[x,], name = 'Sample_ID') %>% 
    left_join(df.meta.clean, by='Sample_ID') %>% 
    left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') %>% 
    mutate(value=log10(value+5e-06))
  fit <- lm(value~Treatment_group+drug, data=tmp)
  f <- coefficients(summary(fit))
  tibble(pathway=x, coef=f[2,1], p.val=f[2,4])}) %>% 
  bind_rows()

g <- df.res %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coef, y=-log10(q.val))) + 
  geom_point() + 
  theme_bw() + xlab('Linear model coefficient') + 
  ylab('-log10(q-value)') + 
  theme(panel.grid.minor=element_blank()) + 
  geom_hline(yintercept = -log10(c(0.05, 0.001)), lty=2)
ggsave(g, filename=here('figures/humann/volcano_baseline.pdf'),
       width = 4, height = 4, useDingbats=FALSE)
# nothing really

# relative abundance difference between arms overall?
df.res <- map(rownames(df.humann), .f=function(x){
  tmp <- enframe(df.humann[x,], name = 'Sample_ID') %>% 
    left_join(df.meta.clean, by='Sample_ID') %>% 
    left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') %>% 
    mutate(value=log10(value+5e-06))
  fit <- lmerTest::lmer(value~Treatment_group+drug+(1|Participant_ID), data=tmp)
  f <- coefficients(summary(fit))
  tibble(pathway=x, coef=f[2,1], p.val=f[2,5])}) %>% 
  bind_rows()

g <- df.res %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coef, y=-log10(q.val))) + 
  geom_point() + 
  theme_bw() + xlab('Linear model coefficient') + 
  ylab('-log10(q-value)') + 
  theme(panel.grid.minor=element_blank()) + 
  geom_hline(yintercept = -log10(c(0.05, 0.001)), lty=2)
ggsave(g, filename=here('figures/humann/volcano_overall.pdf'),
       width = 4, height = 4, useDingbats=FALSE)
# nothing really, either

# relative abundance difference between arms 14-28?
df.14_28 <- df.humann[,df.meta.clean %>% 
                           filter(Timepoint%in%c('14', '21', '28')) %>% 
                           pull(Sample_ID)]
df.14_28 <- df.14_28[rowMeans(df.14_28!=0) > 0.05,]
df.res <- map(rownames(df.14_28), .f=function(x){
  tmp <- enframe(df.14_28[x,], name = 'Sample_ID') %>% 
    left_join(df.meta.clean, by='Sample_ID') %>% 
    left_join(df.abx.exp %>% select(Sample_ID, drug), by='Sample_ID') %>% 
    left_join(df.response %>% select(Participant_ID, Treatment_group),
              by='Participant_ID') %>% 
    mutate(value=log10(value+5e-06))
  fit <- lmerTest::lmer(value~Treatment_group+drug+(1|Participant_ID), data=tmp)
  f <- coefficients(summary(fit))
  tibble(pathway=x, coef=f[2,1], p.val=f[2,5])}) %>% 
  bind_rows()

g <- df.res %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coef, y=-log10(q.val))) + 
  geom_point() + 
  theme_bw() + xlab('Linear model coefficient') + 
  ylab('-log10(q-value)') + 
  theme(panel.grid.minor=element_blank()) + 
  geom_hline(yintercept = -log10(c(0.05, 0.001)), lty=2)
ggsave(g, filename='./figures/humann/volcano_14_to_28.pdf',
       width = 4, height = 4, useDingbats=FALSE)

# wow, bile acids is the stronges differentially abudnant pathway between
# treatment arms at this time
# there must be something there!

# second is starch biosynthesis, last one is bile acid epimerization

# bile acids # PWY-7754/PWY-8134 is only in C. scindens
# starch stuff # PWY-622 is in Bifido, Klebsiella, Fuso, but mostly unclassified
# bile acid epimerization # PWY-6518, only unclassified

df.humann.raw %>% 
  filter(str_detect(pathway, 'PWY-7754')) %>% 
  filter(str_detect(pathway, '\\|')) %>% 
  pivot_longer(-pathway, names_to='Sample_ID') %>% 
  inner_join(df.meta.clean %>% 
               filter(Timepoint%in%c('14', '21', '28')), 
             by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Sample_ID) %>% 
  ggplot(aes(x=Treatment_group, y=value)) + 
    geom_boxplot() +
    facet_grid(~pathway)

# ##############################################################################
# visualize bile acids

as_tibble(df.humann[str_detect(rownames(df.humann), 'PWY-7754|PWY-8134'),],
          rownames='pathway') %>% 
  pivot_longer(-pathway, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=log10(value+5e-05), 
             col=Treatment_group))+ 
    geom_jitter() +
    geom_smooth(method='loess') +
    facet_grid(~pathway) 
    
df.bile <- as_tibble(df.humann[str_detect(rownames(df.humann), 
                                          'PWY-7754|PWY-8134'),],
          rownames='pathway') %>% 
  pivot_longer(-pathway, names_to = 'Sample_ID') %>% 
  group_by(Sample_ID) %>% 
  reframe(m=mean(value)) %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) 
g <- df.bile %>% 
  ggplot(aes(x=Treatment_group, fill=Treatment_group, y=log10(m + 1e-05))) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.2) +
    theme_bw() +
    theme(panel.grid.minor=element_blank()) + 
    xlab('') + ylab('Relative abundance') 
ggsave(g, filename=here('figures/humann/bile_acid_14_28.pdf'),
       width = 4, height = 4, useDingbats=FALSE)
summary(lmerTest::lmer(m~Treatment_group+(1|Participant_ID), 
                       data=df.bile %>% mutate(m=log10(m+1e-05))))
# Fixed effects:
# Estimate Std. Error        df  t value Pr(>|t|)    
# (Intercept)             -4.85492    0.03713 262.35782 -130.765  < 2e-16 ***
# Treatment_groupTac/MTX   0.27850    0.05320 267.56177    5.235 3.33e-07 ***

# maybe account for sequencing depth?
summary(lmerTest::lmer(m~Treatment_group+hostremoved_reads+(1|Participant_ID), 
                       data=df.bile %>% mutate(m=log10(m+1e-05)) %>% 
                         mutate(hostremoved_reads=log10(hostremoved_reads)) %>% 
                         mutate(depth=hostremoved_reads < 1e7)))
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)             -6.37926    0.71958 608.13260  -8.865  < 2e-16 ***
# Treatment_groupTac/MTX   0.27798    0.05281 266.83277   5.264  2.9e-07 ***
# hostremoved_reads        0.20524    0.09676 607.91339   2.121   0.0343 *  
# still significant!

# maybe account for antibiotic treatment?
summary(lmerTest::lmer(m~Treatment_group+hostremoved_reads+(1|Participant_ID), 
                       data=df.bile %>% mutate(m=log10(m+1e-05)) %>% 
                         mutate(hostremoved_reads=log10(hostremoved_reads)) %>% 
                         mutate(depth=hostremoved_reads < 1e7)))
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)             -6.37926    0.71958 608.13260  -8.865  < 2e-16 ***
# Treatment_groupTac/MTX   0.27798    0.05281 266.83277   5.264  2.9e-07 ***
# hostremoved_reads        0.20524    0.09676 607.91339   2.121   0.0343 *  


# ##############################################################################
# visualize bile acids
df.bile.acid <- as_tibble(df.humann[
  str_detect(rownames(df.humann), 
             'PWY-7754|PWY-8134'),
  df.meta.clean %>% 
    pull(Sample_ID)], rownames='pathway') %>%
  pivot_longer(-pathway, names_to = 'Sample_ID') 

df.bile.acid %>% 
  mutate(pathway=str_remove(pathway, '\\:.*')) %>% 
  pivot_wider(names_from = pathway, values_from = value) %>% 
  ggplot(aes(x=log10(`PWY-7754`+1e-06), y=log10(`PWY-8134`+1e-06))) + 
    geom_point()
# basically, very much co-occurent 

g <- df.bile.acid %>% 
  left_join(as_tibble(prop.table(as.matrix(feat.motus), 2), 
                      rownames='species') %>% 
              filter(str_detect(species, 'scindens')) %>% 
              pivot_longer(-species, names_to = 'Sample_ID', 
                           values_to = 'C_scindens'),
            by='Sample_ID') %>% 
  ggplot(aes(x=log10(C_scindens+5e-05), y=log10(value+5e-06))) + 
    geom_point() + 
    facet_grid(species~pathway) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    xlab('log10(C. scindens relative abundance)') + 
    ylab('log10(Metabolic pathway relative abundance)')
# not all C. scindens seem to be corresponding with bile acid pathways here
# interesting!
ggsave(g, filename=here('figures/humann/bile_acids_vs_scindens.pdf'),
       width = 7, height = 4.5, useDingbats=FALSE)

# ##############################################################################
# association between bile acids and outcome?

df.prediction <- df.bile.acid %>% 
  filter(str_detect(pathway, '8134')) %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(Participant_ID) %>% 
  reframe(m=mean(value > 0) > 0.5) %>% 
  full_join(df.response, by='Participant_ID') %>% 
  filter(!is.na(m))

df.res <- .f_test_outcome(df.prediction, landmarked = 28)
# export plots 
for (x in c('cGVHD', 'aGVHD24', 'GRFS')){
  g <- df.res$plots[[x]] + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    scale_colour_manual(values=c('#E98300', '#007C92')) + 
    scale_y_continuous(limits=c(0,1)) 
  g %>% 
    ggsave(filename=here('figures/humann/', 
                         paste0('bile_acids_', x, '_association.pdf')),
           width = 5, height = 4, useDingbats=FALSE)
}

df.res$associations %>% 
  filter(set=='all') %>% View
# # A tibble: 8 × 8
# outcome  estimate  high   low   pval set   n_present n_absent
# <chr>       <dbl> <dbl> <dbl>  <dbl> <chr>     <dbl>    <dbl>
# 1 NRM         1.07  2.60  0.443 0.877  all          58      215
# 2 Relapse     0.876 1.65  0.464 0.682  all          58      215
# 3 cGVHD       0.540 0.982 0.296 0.0436 all          59      215
# 4 cGVHD_MS    0.755 1.60  0.356 0.464  all          59      215
# 5 aGVHD24     0.588 1.18  0.292 0.137  all          42      130
# 6 aGVHD34     1.03  2.78  0.385 0.948  all          58      204
# 7 GRFS        0.671 1.03  0.436 0.0698 all          58      211
# 8 OS          0.794 1.53  0.412 0.491  all          59      215

# ##############################################################################
# any associations for other pathways?

tmp <- as_tibble(df.14_28, rownames='pathway') %>% 
  pivot_longer(-pathway, names_to = 'Sample_ID') %>% 
  mutate(present=value > 0) %>% 
  left_join(df.meta.clean %>% select(Sample_ID, Participant_ID), 
            by='Sample_ID') %>% 
  group_by(Participant_ID, pathway) %>% 
  reframe(presence=mean(present) > 0.4)
tmp %>% 
  group_by(pathway) %>% 
  reframe(freq=mean(presence)) %>% 
  pull(freq) %>% hist(100)
df.associations.all <- tmp %>% 
  group_by(pathway) %>% 
  mutate(freq=mean(presence)) %>% 
  filter(freq > 0.2 & freq < 0.8) %>% 
  group_map(.f=function(.x, .y){
    df.pred <- .x %>% 
      left_join(df.response, by='Participant_ID') %>% 
      mutate(m=presence)
    df.x <- .f_test_outcome(df.pred, landmarked = 28)
    df.x$association %>% mutate(.y)
  }) %>% bind_rows()

# some of the pathways are not found very often at all

# how about the cGVHD associations?
View(df.associations.all %>% filter(set=='all', outcome=='cGVHD'))

# bile acids the second-lowest p-value, and the lowest for protective effects


# examples
# PWY-702: L-methionine biosynthesis II
x <- 'PWY-702: L-methionine biosynthesis II'
x <- 'PWY-5130: 2-oxobutanoate degradation I'
x <- 'PWY-5677: succinate fermentation to butanoate'
x <- 'PWY-7456: &beta;-(1,4)-mannan degradation'
x <- 'PWY-622: starch biosynthesis'
x <- 'PWY-7312: dTDP-&beta;-D-fucofuranose biosynthesis'
x <- 'TRPSYN-PWY: L-tryptophan biosynthesis'
df.pred <- enframe(df.14_28[x,], name='Sample_ID') %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(Participant_ID) %>% 
  reframe(m=mean(value > 0) > 0.4) %>% 
  left_join(df.response, by='Participant_ID')
# df.pred %>% ggplot(aes(x=m, y=log10(m2+5e-05)))+geom_boxplot()
# table(df.pred$Treatment_group, df.pred$m)
df.pred %>% group_by(Treatment_group, m) %>% tally() %>% 
  ggplot(aes(x=Treatment_group, y=n, fill=m)) + 
  geom_bar(stat='identity')
a <- .f_test_outcome(df.pred, landmarked = 28)


enframe(df.humann[x,], name='Sample_ID') %>% 
  full_join(df.meta.clean, 
             by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=log10(value+5e-05), 
             col=Treatment_group))+ 
  geom_jitter() +
  geom_smooth(method='loess')


# ##############################################################################
# specifically how does it look for butyrate production?

df.associations.all %>% 
  filter(set=='all') %>% 
  filter(str_detect(pathway, 'butanoate')) %>% View

scfa.pathways <- c(
  'CENTFERM-PWY: pyruvate fermentation to butanoate',
  'PWY-5676: acetyl-CoA fermentation to butanoate II',
  'P163-PWY: L-lysine fermentation to acetate and butanoate',
  'PWY-5677: succinate fermentation to butanoate',
  
  'P108-PWY: pyruvate fermentation to propanoate I',
  'PWY-5494: pyruvate fermentation to propanoate II (acrylate pathway)',
  'PWY-5088: L-glutamate degradation VIII (to propanoate)'
)
df.butyrate <- as_tibble(df.humann[scfa.pathways,], rownames='pathway') %>% 
  pivot_longer(-pathway, names_to = 'Sample_ID') %>% 
  left_join(df.meta.clean %>% select(Sample_ID, Participant_ID, Timepoint), 
            by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID')

tmp %>% 
  group_by(pathway) %>% 
  reframe(freq=mean(presence)) %>% 
  filter(pathway %in% scfa.pathways)


df.associations.all %>% 
  filter(set=='all') %>% 
  filter(pathway %in% scfa.pathways) %>% 
  mutate(eff=paste0(sprintf(fmt='%.3f', estimate), ' (', 
                    sprintf('%.3f', low), ', ', 
                    sprintf('%.3f', high), ')')) %>% 
  mutate(pval=sprintf('%.3f', pval)) %>% 
  select(pathway, outcome, pval, eff) %>% 
  mutate(outcome=factor(
    outcome, levels = c('NRM', 'Relapse', 'aGVHD24', 'aGVHD34', 
                        'cGVHD', 'cGVHD_MS', 'OS', 'GRFS'))) %>% 
  arrange(pathway, outcome) %>% 
  View()




# where do the pathways come from?
# df.humann.raw %>% 
#   filter(str_detect(pathway, 'PWY-5088')) %>% 
#   filter(str_detect(pathway, '\\|')) %>% 
#   pivot_longer(-pathway, names_to='Sample_ID') %>% 
#   inner_join(df.meta.clean %>% 
#                filter(Timepoint%in%c('14', '21', '28')), 
#              by='Sample_ID') %>% 
#   left_join(df.response %>% select(Participant_ID, Treatment_group),
#             by='Participant_ID') %>% 
#   group_by(Sample_ID) %>% 
#   ggplot(aes(x=pathway, y=value, fill=Treatment_group)) + 
#   geom_boxplot() +
#   coord_flip()

g <- df.butyrate %>% 
  ggplot(aes(x=as.numeric(Timepoint), 
             y=log10(value + 5e-05), col=Treatment_group)) + 
  geom_jitter(alpha=0.3, pch=16) +
  xlab('Sampling timepoint') + ylab('Relative pathway abundance') + 
  scale_colour_manual(values=unlist(colours$group.colours)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  geom_smooth(method='loess') + 
  scale_x_continuous(breaks=seq_len(length(levels(df.butyrate$Timepoint))),
                     labels=levels(df.butyrate$Timepoint)) +
  facet_wrap(~pathway)
ggsave(g, filename=here('figures/humann/butyrate_pathways_time.pdf'),
       width = 10, height = 7, useDingbats=FALSE)

g <- df.butyrate %>% 
  group_by(Sample_ID) %>% 
  reframe(value=sum(value), Timepoint, Treatment_group) %>% 
  distinct() %>% 
  ggplot(aes(x=as.numeric(Timepoint), 
             y=log10(value + 5e-05), col=Treatment_group)) + 
  geom_jitter(alpha=0.3, pch=16) +
  xlab('Sampling timepoint') + ylab('Relative pathway abundance') + 
  scale_colour_manual(values=unlist(colours$group.colours)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  geom_smooth(method='loess') + 
  scale_x_continuous(breaks=seq_len(length(levels(df.butyrate$Timepoint))),
                     labels=levels(df.butyrate$Timepoint)) 
ggsave(g, filename='./figures/humann/butyrate_time.pdf',
       width = 8, height = 6, useDingbats=FALSE)

# any difference?
lmerTest::lmer(value~Treatment_group + (1|Participant_ID), 
               data=df.butyrate %>% 
                 group_by(Sample_ID) %>% 
                 reframe(value=sum(value), Timepoint, 
                         Treatment_group, Participant_ID) %>% 
                 distinct() %>% 
                 mutate(value=log10(value+5e-05)) %>% 
                 filter(Timepoint %in% c('14', '21', '28'))) %>% 
  summary
  
# associations?
df.pred <- df.butyrate %>% 
  group_by(Sample_ID) %>% 
  reframe(value=sum(value)) %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(Participant_ID) %>% 
  reframe(m2=mean(value)) %>% 
  full_join(df.response, by='Participant_ID') %>% 
  filter(!is.na(m2)) %>% 
  mutate(m=m2>median(m2)) 
a <- .f_test_outcome(df.pred, landmarked = 28)
a$associations %>% filter(set=='all') %>% 
  mutate(eff=paste0(sprintf(fmt='%.3f', estimate), ' (', 
                    sprintf('%.3f', low), ', ', 
                    sprintf('%.3f', high), ')')) %>% 
  mutate(pval=sprintf('%.3f', pval))  %>% 
  select(outcome, pval, eff) %>% 
  mutate(outcome=factor(
    outcome, levels = c('NRM', 'Relapse', 'aGVHD24', 'aGVHD34', 
                        'cGVHD', 'cGVHD_MS', 'OS', 'GRFS'))) %>% 
  arrange(outcome) %>% 
  View()
# no associations

