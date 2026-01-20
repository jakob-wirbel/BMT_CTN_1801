# ##############################################################################
#
## Microbial diversity and absolute abundance
#
# ##############################################################################

library("tidyverse")
library("here")
library("vegan")
library("survival")
library("ggfortify")

set.seed(1234)
colours <- yaml::read_yaml('./files/colours.yml')

load('./data/all_data.RData')

# ##############################################################################
# alpha diversity
feat.rar <- rrarefy(t(feat.motus), 3000)
alpha <- enframe(diversity(feat.rar, index='shannon'), 
                 value='shannon', name='Sample_ID') %>% 
  full_join(enframe(diversity(feat.rar, index='invsimpson'), 
                    value='invsimpson', name='Sample_ID'), by='Sample_ID') %>% 
  full_join(enframe(diversity(feat.rar, index='simpson'), 
                    value='simpson', name='Sample_ID'), by='Sample_ID')

df.alpha <- alpha %>% 
  left_join(df.meta.clean, by='Sample_ID') %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group), 
            by='Participant_ID')

# ##############################################################################
# alpha over time
g <- df.alpha %>% 
  select(Timepoint, shannon, invsimpson, simpson, 
         Sample_ID, Treatment_group) %>% 
  pivot_longer(c(shannon, invsimpson, simpson)) %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=value, col=Treatment_group)) + 
    geom_jitter(alpha=0.3, pch=16) +
    xlab('Sampling timepoint') + ylab('Prokaryotic alpha diversity') + 
    scale_colour_manual(values=unlist(colours$group.colours)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    geom_smooth(method='loess') + 
    scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                       labels=levels(df.alpha$Timepoint)) + 
    facet_grid(name~., scale='free')
ggsave(g, filename='./figures/alpha/alpha_time.pdf',
       width = 8, height = 8, useDingbats=FALSE)

 
# test each time point
p.vals.alpha <- df.alpha %>% 
  select(Timepoint, shannon, invsimpson, simpson, 
         Sample_ID, Treatment_group) %>% 
  pivot_longer(c(shannon, invsimpson, simpson)) %>%  
  group_by(Timepoint, name) %>% 
  group_map(.f=function(.x, .y){
    fit <- lm(value~Treatment_group, data=.x)
    tibble(p=coefficients(summary(fit))[2,4],
           ef=coefficients(summary(fit))[2,1]) %>% 
      mutate(.y) %>% mutate(n=nrow(.x))
  }) %>% bind_rows()

sel.points <- c('PCON', 'PINF', '7', '14', '21', '28', '35', '42')
df.alpha %>% 
  filter(Timepoint %in% sel.points) %>%
  select(Timepoint, shannon, invsimpson, simpson, 
         Sample_ID, Treatment_group) %>% 
  pivot_longer(c(shannon, invsimpson, simpson)) %>% 
  ggplot(aes(x=Timepoint, y=value, fill=Treatment_group)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=unlist(colours$group.colours)) + 
    geom_text(aes(label=p_nice, y=-1), data=p.vals.alpha %>% 
                filter(Timepoint %in% sel.points) %>% 
                mutate(p_nice=paste0('p = ', sprintf('%.3f', p)),
                       Treatment_group=NA)) + 
    xlab('Sampling timepoint') + 
    ylab('Prokaryotic alpha diversity') + 
    facet_grid(name~., scales='free_y')

df.lmer <- df.alpha %>% 
  filter(Timepoint %in% c(14, 21, 28)) %>% 
  select(Timepoint, shannon, invsimpson, simpson, 
         Sample_ID, Treatment_group, Participant_ID) %>% 
  pivot_longer(c(shannon, invsimpson, simpson)) %>% 
  group_by(name) %>% 
  group_map(.f=function(.x, .y){
    fit <- lmerTest::lmer(value~Treatment_group + (1|Participant_ID),
                          data=.x)
    p <- coefficients(summary(fit))
    as_tibble(as.data.frame(p), rownames='id') %>% mutate(.y)
    }) %>% 
  bind_rows()
df.lmer
# # A tibble: 6 × 7
# id                     Estimate `Std. Error`    df `t value` `Pr(>|t|)` name      
# <chr>                     <dbl>        <dbl> <dbl>     <dbl>      <dbl> <chr>     
# 1 (Intercept)              8.47         0.577   267.     14.7   3.71e- 36 invsimpson
# 2 Treatment_groupTac/MTX   2.11         0.827   273.      2.55  1.13e-  2 invsimpson
# 3 (Intercept)              2.37         0.0714  256.     33.2   1.48e- 94 shannon   
# 4 Treatment_groupTac/MTX   0.213        0.102   261.      2.08  3.88e-  2 shannon   
# 5 (Intercept)              0.755        0.0153  243.     49.2   2.67e-128 simpson   
# 6 Treatment_groupTac/MTX   0.0339       0.0220  249.      1.54  1.25e-  1 simpson 


# In the Peled NEJM paper, they show invSimpson with a log scale
g.alpha.peled <- df.alpha %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=invsimpson, col=Treatment_group)) + 
  geom_jitter(alpha=0.3, pch=16) +
  xlab('Sampling timepoint') + ylab('Prokaryotic alpha diversity') + 
  scale_colour_manual(values=unlist(colours$group.colours)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  geom_smooth(method='loess') + 
  scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                     labels=levels(df.alpha$Timepoint)) + 
  scale_y_log10(guide = "axis_logticks")
# looks very similar, overall
ggsave(g.alpha.peled, filename='./figures/alpha/alpha_time_peled.pdf',
       width = 8, height = 5, useDingbats=FALSE)


# ##############################################################################
# extraction copies 

g <- df.alpha %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=copies_16S, col=Treatment_group)) + 
    geom_jitter(alpha=0.3, pch=16) +
    xlab('Sampling timepoint') + ylab('log10(16S copies per extraction)') + 
    scale_colour_manual(values=unlist(colours$group.colours)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    geom_smooth(method='loess') + 
    scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                       labels=levels(df.alpha$Timepoint)) + 
    NULL
ggsave(g, filename='./figures/ddPCR/copies_time.pdf',
       width = 8, height = 4, useDingbats=FALSE)

# only measured ones
g <- df.alpha %>% 
  filter(copies=='measured') %>% 
  ggplot(aes(x=as.numeric(Timepoint), y=copies_16S, col=Treatment_group)) + 
  geom_jitter(alpha=0.3, pch=16) +
  xlab('Sampling timepoint') + ylab('log10(16S copies per extraction)') + 
  scale_colour_manual(values=unlist(colours$group.colours)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  geom_smooth(method='loess') + 
  scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                     labels=levels(df.alpha$Timepoint)) + 
  NULL
ggsave(g, filename='./figures/ddPCR/copies_time_measured_only.pdf',
       width = 8, height = 4, useDingbats=FALSE)

# DNA concentration
g <- df.alpha %>% 
  mutate(all=paste0(Treatment_group, '-', copies)) %>% 
  mutate(DNA_concentration=case_when(
    DNA_concentration < 0.1~0.1, TRUE~DNA_concentration)) %>% 
  ggplot(aes(x=as.numeric(Timepoint), 
             y=log10(DNA_concentration), col=Treatment_group)) + 
  geom_jitter(alpha=0.3, pch=16) +
  xlab('Sampling timepoint') + ylab('log10(16S copies per extraction)') + 
  scale_colour_manual(values=unlist(colours$group.colours)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  geom_smooth(method='loess') + 
  scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                     labels=levels(df.alpha$Timepoint)) + 
  # scale_shape_manual(values=c(4, 16))  +
  NULL
ggsave(g, filename='./figures/ddPCR/concentrations_time.pdf',
       width = 8, height = 4, useDingbats=FALSE)

# test each time point
p.vals.abs <- df.alpha %>% 
  # filter(Timepoint %in% c('14', '21', '28')) %>% 
  group_by(Timepoint) %>% 
  group_map(.f=function(.x, .y){
    fit <- lm(copies_16S~Treatment_group, data=.x)
    mean.value <- .x %>% group_by(Treatment_group) %>% 
      reframe(m=mean(copies_16S))
    tibble(p=coefficients(summary(fit))[2,4],
           ef=coefficients(summary(fit))[2,1],
           m.ptyc=mean.value %>% filter(Treatment_group=='PTCy/Tac/MMF') %>% 
             pull(m),
           m.tac=mean.value %>% filter(Treatment_group=='Tac/MTX') %>% 
             pull(m)) %>% 
      mutate(.y) %>% 
      mutate(n=nrow(.x))
  }) %>% bind_rows()

g <- df.alpha %>% 
  filter(Timepoint %in% sel.points) %>%
  ggplot(aes(x=Timepoint, y=copies_16S, fill=Treatment_group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) + 
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=unlist(colours$group.colours)) + 
  geom_text(aes(label=p_nice, y=12), data=p.vals.abs %>% 
              filter(Timepoint %in% sel.points) %>% 
              mutate(p_nice=paste0('p = ', sprintf('%.3f', p)),
                     Treatment_group=NA)) + 
  xlab('Sampling timepoint') + 
  ylab('log10(16S copies per extraction)')
ggsave(g, filename=here('figures/ddPCR/copies_time_test.pdf'),
       width = 10, height = 4, useDingbats=FALSE)

p <- lmerTest::lmer(copies_16S ~ Treatment_group+(1|Participant_ID),
                    data=df.alpha %>% 
                      filter(Timepoint %in% c('14', '21', '28')))
p <- coefficients(summary(p))[2,5] # 2.984152e-05
g <- df.alpha %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>%
  ggplot(aes(x=Treatment_group, y=copies_16S, fill=Treatment_group)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) + 
    geom_boxplot(outlier.shape = NA, alpha=0.8) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=unlist(colours$group.colours), guide='none') + 
    annotate(geom='text', label=paste0('p = ', sprintf(fmt='%.2e', p)), 
             y=10.5, x=1.5) +
    xlab('') + ylab('log10(16S copies) across days 14 to 28')
ggsave(g, filename=here('figures/ddPCR/copies_days.pdf'),
       width = 3, height = 4, useDingbats=FALSE)

# copies vs alpha
spearman.rho <- cor(df.alpha$shannon, df.alpha$copies_16S, 
                    method='spearman')
g <- df.alpha %>% 
  ggplot(aes(y=copies_16S, x=shannon)) + 
    geom_point() + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    ylab('log10(16S copies per extraction)') + 
    xlab("Prokaryotic alpha diversity") + 
    geom_smooth(method='lm', formula = 'y~x') + 
    annotate(x=-Inf, y=-Inf, geom='text', 
             hjust=-6, vjust=-20,
             label=paste0('rho=', sprintf(fmt='%.2f', spearman.rho)))
ggsave(g, filename=here('figures/alpha/alpha_vs_abs_ab.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

# ##############################################################################
# How about host DNA in the stool (as a marker for inflammation?)

g <- df.alpha %>% 
  ggplot(aes(x=as.numeric(Timepoint), 
             y=(1-hostremoved_frac)*100, 
             col=Treatment_group)) + 
    geom_jitter(alpha=0.3, pch=16) +
    xlab('Sampling timepoint') + ylab('Percentage of human reads') + 
    scale_colour_manual(values=unlist(colours$group.colours)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    geom_smooth(method='loess') + 
    scale_x_continuous(breaks=seq_len(length(levels(df.alpha$Timepoint))),
                       labels=levels(df.alpha$Timepoint)) + 
    NULL
ggsave(g, filename=here('figures/general/host_fraction_time.pdf'),
       width = 6, height = 4, useDingbats=FALSE)


df.alpha %>% 
  filter(Timepoint %in% c('14', '21', '28')) %>% 
  ggplot(aes(x=Treatment_group, y=1-hostremoved_frac)) + 
    geom_boxplot()

fit <- lmerTest::lmer(hostremoved_frac~Treatment_group+(1|Participant_ID), 
                      data=df.alpha %>% 
                        filter(Timepoint %in% c('14', '21', '28')))
summary(fit)
# p-value = 0.815
