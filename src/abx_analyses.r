# ##############################################################################
#
## ABX analyses
#
# ##############################################################################

library("tidyverse")
library("here")
library("vegan")
library("survival")

colours <- yaml::read_yaml('./files/colours.yml')

load('./data/all_data.RData')

# raw antibiotics data
abx.data <- read_tsv(here('./data/antibiotics_data.tsv'), 
                     col_types = cols()) 

# ##############################################################################
# Exposure to antimicrobial

# overall exposure to antimicrobials
g <- abx.data %>% 
  filter(!Class %in% c('antifungal', 'antiprotist', 'antiviral')) %>% 
  full_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=length(unique(Participant_ID))) %>% 
  group_map(.f=function(.x, .y){
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      filter(Start_rel==min(Start_rel)) %>%
      slice_head(n=1) %>% 
      ungroup() %>% 
      arrange(Start_rel) %>% 
      select(n.all, Start_rel) %>% 
      group_by(Start_rel) %>% 
      reframe(n.all=unique(n.all), s=n()) %>% 
      mutate(cumulative_sum=cumsum(s)) %>% 
      mutate(freq=cumulative_sum/n.all) 
    tmp %>% 
      add_row(Start_rel=356, freq=tmp$freq[nrow(tmp)]) %>% 
      mutate(.y)
  }) %>% 
  bind_rows() %>% 
  ggplot(aes(x=Start_rel, y=freq, col=Treatment_group)) + 
    geom_vline(xintercept = 0, lty=2) +
    geom_line(aes(group=Treatment_group)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    xlab('Day relative to infusion') + 
    ylab('Percentage of patients exposed to antibiotics') + 
    scale_colour_manual(values=unlist(colours$group.colours)) +
    ylim(0,1)
ggsave(g, filename='./figures/abx/cumulative_exposure.pdf',
       width = 5, height = 4, useDingbats=FALSE)    

ggsave(g + coord_cartesian(xlim=c(-20, 35)), 
       filename='./figures/abx/cumulative_exposure_35.pdf',
       width = 5, height = 4, useDingbats=FALSE)    


# exposure by class
g <- abx.data %>% 
  full_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=length(unique(Participant_ID))) %>% 
  group_by(Treatment_group, Class) %>% 
  group_map(.f=function(.x, .y){
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      filter(Start_rel==min(Start_rel)) %>%
      slice_head(n=1) %>% 
      ungroup() %>% 
      arrange(Start_rel) %>% 
      select(n.all, Start_rel) %>% 
      group_by(Start_rel) %>% 
      reframe(n.all=unique(n.all), s=n()) %>% 
      mutate(cumulative_sum=cumsum(s)) %>% 
      mutate(freq=cumulative_sum/n.all) 
    tmp %>% 
      add_row(Start_rel=365, freq=tmp$freq[nrow(tmp)]) %>% 
      mutate(.y)
  }) %>% 
  bind_rows() %>% 
  ggplot(aes(x=Start_rel, y=freq, col=Treatment_group)) + 
  geom_vline(xintercept = 0, lty=2) +
  geom_line(aes(group=Treatment_group)) + 
  facet_wrap(~Class) + theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  xlab('Day relative to infusion') + 
  ylab('Percentage of patients exposed') + 
  scale_colour_manual(values=unlist(colours$group.colours)) +
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave(g, filename=here('figures/abx/cumulative_exposure_detail.pdf'),
       width = 8, height = 5, useDingbats=FALSE)    

ggsave(g + coord_cartesian(xlim=c(-20, 35)), 
       filename=here('figures/abx/cumulative_exposure_detail_35.pdf'),
       width = 8, height = 5, useDingbats=FALSE)   

#
df.test.logrank <- abx.data %>% 
  filter(!Class %in% c('antifungal', 'antiprotist', 'antiviral')) %>% 
  full_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=length(unique(Participant_ID))) %>% 
  group_by(Participant_ID) %>% 
  reframe(m=min(Start_rel), Treatment_group=unique(Treatment_group)) %>% 
  mutate(status=case_when(is.na(m)~0,
                          m > 35 ~ 0,
                          TRUE~1)) %>% 
  mutate(m=case_when(is.na(m)~35, m > 35~35, TRUE~m))

surv_object <- Surv(df.test.logrank$m, df.test.logrank$status)

logrank_test <- survdiff(surv_object ~ df.test.logrank$Treatment_group)
# Call:
# survdiff(formula = surv_object ~ df.test.logrank$Treatment_group)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# df.test.logrank$Treatment_group=PTCy/Tac/MMF 157      135      149      1.36      3.46
# df.test.logrank$Treatment_group=Tac/MTX      147      140      126      1.62      3.46
# 
# Chisq= 3.5  on 1 degrees of freedom, p= 0.06
logrank_test$pvalue
# 0.06275168


abx.data %>% 
  filter(!Class %in% c('antifungal', 'antiprotist', 'antiviral')) %>% 
  full_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=length(unique(Participant_ID))) %>% 
  group_by(Class) %>% 
  group_map(.f=function(.x, .y){
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      reframe(m=min(Start_rel)) %>% 
      full_join(df.response %>% select(Participant_ID, Treatment_group),
                by='Participant_ID') %>% 
      mutate(status=case_when(is.na(m)~0,
                              m > 35 ~ 0,
                              TRUE~1)) %>% 
      mutate(m=case_when(is.na(m)~35, m > 35~35, TRUE~m))
    surv_object <- Surv(tmp$m, tmp$status)
    logrank_test <- survdiff(surv_object ~ tmp$Treatment_group)
    tibble(pval=logrank_test$pvalue, .y)
  }) %>% 
  bind_rows() 
# # A tibble: 15 × 2
# pval Class         
# <dbl> <chr>         
# 1   0.454  aminoglycoside
# 2 NaN      ansamycins    
# 3   0.241  carbapenem    
# 4   0.512  cephalosporin 
# 5   0.591  clindamycin   
# 6   0.232  daptomycin    
# 7   0.663  glycopeptide  
# 8   0.408  linezolid     
# 9   0.348  macrolide     
# 10   0.182  monobactam    
# 11 NaN      nitrofuran    
# 12   0.460  other         
# 13   0.104  penicillin    
# 14   0.301  polypeptide   
# 15   0.691  quinolones    
# 16   0.0397 sulfonamide   
# 17   0.699  tetracycline  

# ##############################################################################
# Get the number of doses per class/participant/day


abx.others <- abx.data %>% filter(is.na(Abx)) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') 
abx.condensed <- abx.data %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  filter(!is.na(Abx)) %>% 
  group_by(Participant_ID, Treatment_group, Abx, Route, Class, Class_2) %>% 
  group_map(.f=function(.x, .y){
    tmp <- IRanges::IRanges(start=as.numeric(.x$Start_rel), 
                            end=as.numeric(.x$End_rel))
    as_tibble(as.data.frame(IRanges::reduce(tmp))) %>% 
      mutate(.y)
  }) %>% 
  bind_rows()

# count the number of doses (full day prescribed an abx) in a way that they
# can be combined by class/over-all class

# up to day 35
df.doses <- abx.condensed %>% 
  group_by(Participant_ID, Treatment_group, Abx, Route, Class, Class_2) %>% 
  group_map(.f=function(.x, .y){
    q <- IRanges::IRanges(start=.x$start, end=.x$end)
    ref <- IRanges::IRanges(start=-7, end=35)
    ov <- as.data.frame(IRanges::pintersect(IRanges::findOverlapPairs(q, ref)))
    if (nrow(ov) > 0){
      doses <- sum(ov$width)
    } else {
      doses <- 0
    }
    tibble(doses=doses) %>% mutate(.y)
  }) %>% 
  bind_rows() %>% 
  group_by(Participant_ID, Class) %>% 
  reframe(m=sum(doses), Class_2, Treatment_group) %>% 
  distinct() %>% 
  select(-Class_2) %>% 
  bind_rows(tibble(Participant_ID=abx.others$Participant_ID,
                   Treatment_group=abx.others$Treatment_group)) %>% 
  pivot_wider(names_from = Class, values_from = m, values_fill = 0) %>% 
  pivot_longer(-c(Participant_ID, Treatment_group), 
               names_to = 'class', values_to = 'm') %>% 
  filter(class!='NA')

g <- df.doses %>% 
  ggplot(aes(x=Treatment_group, y=m, fill=Treatment_group)) + 
  geom_boxplot() + 
  facet_wrap(~class, scales='free') + 
  xlab('') + ylab('Number of antibiotic doses') + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=unlist(colours$group.colours))
ggsave(g, filename=here('figures/abx/number_of_doses_detail.pdf'),
       width = 12, height = 8, useDingbats=FALSE)

df.doses %>% 
  group_by(class) %>% 
  group_map(.f=function(.x, .y){
    t <- wilcox.test(.x$m~.x$Treatment_group)
    t.t <- t.test(.x$m~.x$Treatment_group)
    tibble(p=t$p.value, mean.ptcy=t.t$estimate['mean in group PTCy/Tac/MMF'],
           mean.tacmtx=t.t$estimate['mean in group Tac/MTX']) %>% 
      mutate(diff=mean.tacmtx-mean.ptcy) %>% 
      mutate(.y)
  }) %>% bind_rows()
# # A tibble: 20 × 5
# p mean.ptcy mean.tacmtx    diff class         
# <dbl>     <dbl>       <dbl>   <dbl> <chr>         
# 1   0.456      0.134      0.0340 -0.0997 aminoglycoside
# 2 NaN          0          0       0      ansamycins    
# 3   0.660     32.8       34.3     1.44   antifungal    
# 4   0.456      4.60       5.01    0.408  antiprotist   
# 5   0.770     46.6       45.5    -1.08   antiviral     
# 6   0.203      1.66       0.639  -1.02   carbapenem    
# 7   0.174      7.72       6.29   -1.43   cephalosporin 
# 8   0.584      0.140      0.0680 -0.0721 clindamycin   
# 9   0.225      0.592      0.0816 -0.511  daptomycin    
# 10   0.536      5.09       4.48   -0.606  glycopeptide  
# 11   0.605      0.376      0.177  -0.199  linezolid     
# 12   0.348      0.401      0.0816 -0.320  macrolide     
# 13   0.174      0.268      0.0272 -0.240  monobactam    
# 14 NaN          0          0       0      nitrofuran    
# 15   0.464      1.36       1.69    0.324  other         
# 16   0.0844     1.55       2.54    0.996  penicillin    
# 17   0.304      0          0.0408  0.0408 polypeptide   
# 18   0.881     12.6       12.0    -0.613  quinolones    
# 19   0.0346     4.65       6.68    2.03   sulfonamide   
# 20   0.669      0.682      1.01    0.332  tetracycline  


g <- df.doses %>% 
  filter(!class %in% c('antifungal', 'antiprotist', 'antiviral')) %>% 
  group_by(Treatment_group, Participant_ID) %>% 
  reframe(s=sum(m)) %>% 
  ggplot(aes(x=Treatment_group, y=s, fill=Treatment_group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1) + 
  xlab('') + ylab('Number of antibiotic doses') + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=unlist(colours$group.colours))
ggsave(g, filename=here('figures/abx/number_of_doses.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
df.doses %>% 
  filter(!class %in% c('antifungal', 'antiprotist', 'antiviral')) %>% 
  group_by(Treatment_group, Participant_ID) %>% 
  reframe(s=sum(m)) %>% 
  wilcox.test(data=., s~Treatment_group)
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  s by Treatment_group
# W = 11908, p-value = 0.6302
# alternative hypothesis: true location shift is not equal to 0


# ##############################################################################
# how about prophylaxis vs non-prophylaxis abx

#' this is a bit harder to assess
#' talked with Nathan about the different Abx
#' some of them are clearly non-prophylaxis, some others are clearly prophylaxis
#' again for others it might be either, depending on the type of administration
#' Nathan's rule of thumb was that prophylaxis is probably not given IV, but
#' that might be different in some cases
#' 
#' Okay, how to approach this?
#' Look at the most common abx, just by sheer number of doses
#' classify those as prophylaxis or not, hopefully that will catch the most
#' important cases
#' 
#' Okay, problem with doses alone: some patients might get the abx all the time,
#' but not that many patients overall. So number of doses alone is not the
#' way forward, i guess.
#' 
#' abx.data %>%
#' left_join(df.response %>% select(Participant_ID, Treatment_group),
#'  by='Participant_ID') %>% 
#' filter(Abx=='Amoxicillin') %>%
#' ggplot(aes(x=Start_rel, xend=End_rel, y=Participant_ID, 
#'  col=Treatment_group)) +
#'  geom_segment() +
#'  facet_grid(Treatment_group~., scales='free_y') +
#'  coord_cartesian(xlim=c(-10, 35), clip = 'off')
#' 
#' Any difference in exposure for the most common abx across arm? that makes 
#' maybe more sense
#' 
#' The other thing we should look at is the effect of abx treatment directly 
#' before a sample. Meaning, was the patient exposed to this abx before sampling
#' occurred? That might be best

abx.select <- abx.condensed %>% 
  group_by(Abx, Route, Class, Class_2) %>% 
  group_map(.f=function(.x, .y){
    q <- IRanges::IRanges(start=.x$start, end=.x$end)
    ref <- IRanges::IRanges(start=-7, end=35)
    ov <- as.data.frame(IRanges::pintersect(IRanges::findOverlapPairs(q, ref)))
    if (nrow(ov) > 0){
      doses <- sum(ov$width)
    } else {
      doses <- 0
    }
    tibble(doses=doses) %>% mutate(.y)
  }) %>% 
  bind_rows() %>% 
  group_by(Abx, Route) %>% 
  reframe(m=sum(doses), Class, Class_2) %>% 
  distinct() %>% 
  arrange(desc(m)) %>% 
  filter(!Class %in% c('antiviral', 'antifungal', 'antiprotist'))

#View(abx.select) 100 seems like a good cutoff, maybe?
abx.select <- abx.select %>% 
  filter(m > 10) %>% pull(Abx)

# exposure across arms for single abx
n.all <- df.response %>% 
  select(Treatment_group, Participant_ID) %>% 
  distinct() %>% 
  group_by(Treatment_group) %>% 
  reframe(n.all=n())

df.exposure.most.common <- abx.data %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  filter(Abx %in% abx.select) %>% 
  filter(!is.na(Treatment_group)) %>% 
  mutate(Abx=paste0(Abx, '-', Route)) %>% 
  group_by(Treatment_group, Abx) %>% 
  left_join(n.all, by='Treatment_group') %>% 
  group_map(.f=function(.x, .y){
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      filter(Start_rel==min(Start_rel)) %>%
      slice_head(n=1) %>% 
      ungroup() %>% 
      arrange(Start_rel) %>% 
      select(n.all, Start_rel) %>% 
      group_by(Start_rel) %>% 
      reframe(n.all=unique(n.all), s=n()) %>% 
      mutate(cumulative_sum=cumsum(s)) %>% 
      mutate(freq=cumulative_sum/n.all) 
    tmp %>% 
      add_row(Start_rel=356, freq=tmp$freq[nrow(tmp)]) %>% 
      mutate(.y)
  }) %>% 
  bind_rows() 

g <- df.exposure.most.common %>% 
  ggplot(aes(x=Start_rel, y=freq, col=Treatment_group)) + 
  geom_vline(xintercept = 0, lty=2) +
  geom_line(aes(group=Treatment_group)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('Day relative to infusion') + 
  ylab('Percentage of patients exposed to antibiotics') + 
  scale_colour_manual(values=unlist(colours$group.colours)) +
  ylim(0,1) + coord_cartesian(xlim=c(-20, 35)) +
  facet_wrap(~Abx)
ggsave(g, filename=here('figures/abx/exposure_most_common.pdf'),
       width = 8, height = 7, useDingbats=FALSE)

# test this as well
abx.data %>% 
  filter(Abx %in% abx.select) %>% 
  mutate(Abx=paste0(Abx, '-', Route)) %>%
  group_by(Abx) %>% 
  group_map(.f=function(.x, .y){ 
    # browser()
    tmp <- .x %>% 
      group_by(Participant_ID) %>% 
      reframe(m=min(Start_rel)) %>% 
      full_join(df.response %>% select(Participant_ID, Treatment_group),
                by='Participant_ID') %>% 
      mutate(status=case_when(is.na(m)~0,
                              m > 35 ~ 0,
                              TRUE~1)) %>% 
      mutate(m=case_when(is.na(m)~35, m > 35~35, TRUE~m))
    surv_object <- Surv(tmp$m, tmp$status)
    logrank_test <- survdiff(surv_object ~ tmp$Treatment_group)
    tibble(pval=logrank_test$pvalue, .y)
  }) %>% bind_rows() %>% 
  filter(!is.na(pval)) %>% 
  mutate(q.val=p.adjust(pval, method='BH')) %>% 
  arrange(pval)
# # A tibble: 27 × 3
# pval Abx                                             q.val
# <dbl> <chr>                                           <dbl>
# 1 0.0433 Trimethoprim-Sulfamethoxazole/Bactrim/Septra-PO 0.600
# 2 0.0725 Doxycycline/Doryx/Doxy, adoxa-IV                0.600
# 3 0.0835 PenicillinVK/Apo-pen-VK/Novo-pen-VK-PO          0.600
# 4 0.0927 Cefepime/Maxipime-PO                            0.600
# 5 0.124  Dapsone-PO                                      0.600
# 6 0.148  Cefepime/Maxipime-IV                            0.600
# 7 0.160  Ceftriaxone/Rocephin-IV                         0.600
# 8 0.179  Levofloxacin/Levaquin-PO                        0.600
# 9 0.232  Daptomycin/Cubicin-IV                           0.600
# 10 0.248  Pipercillin-Tazobactam/Zosyn-IV                 0.600
# 11 0.271  Amoxicillin-PO                                  0.600
# 12 0.290  Metronidazole/Flagyl-PO                         0.600
# 13 0.301  Pipercillin-Tazobactam/Zosyn-PO                 0.600
# 14 0.333  Ceftriaxone/Rocephin-PO                         0.600
# 15 0.333  Metronidazole/Flagyl-NA                         0.600
# 16 0.372  Meropenem/Merrem-IV                             0.627
# 17 0.522  Levofloxacin/Levaquin-IV                        0.828
# 18 0.633  Vancomycin-IV                                   0.949
# 19 0.681  Doxycycline/Doryx/Doxy, adoxa-PO                0.968
# 20 0.776  Cefazolin-IV                                    0.975
# 21 0.823  Ciprofloxacin/Cipro-PO                          0.975
# 22 0.863  Trimethoprim-Sulfamethoxazole/Bactrim/Septra-IV 0.975
# 23 0.865  Cefdinir-PO                                     0.975
# 24 0.930  Amoxicillin Clavulanate/Augmentin-PO            0.975
# 25 0.941  Vancomycin-PO                                   0.975
# 26 0.965  Ciprofloxacin/Cipro-IV                          0.975
# 27 0.975  Metronidazole/Flagyl-IV                         0.975

# doses
df.doses.most.common <- abx.condensed %>% 
  group_by(Participant_ID, Treatment_group, Abx, Route) %>% 
  group_map(.f=function(.x, .y){
    q <- IRanges::IRanges(start=.x$start, end=.x$end)
    ref <- IRanges::IRanges(start=-7, end=35)
    ov <- as.data.frame(IRanges::pintersect(IRanges::findOverlapPairs(q, ref)))
    if (nrow(ov) > 0){
      doses <- sum(ov$width)
    } else {
      doses <- 0
    }
    tibble(doses=doses) %>% mutate(.y)
  }) %>% 
  bind_rows() %>% 
  group_by(Participant_ID, Abx, Route) %>% 
  reframe(m=sum(doses), Treatment_group) %>% 
  distinct() %>% 
  filter(Abx%in%abx.select) %>% 
  mutate(Abx=paste0(Abx, '-', Route)) %>% 
  select(-Route) %>% 
  pivot_wider(names_from = Abx, values_from = m, values_fill = 0) %>% 
  pivot_longer(-c(Participant_ID, Treatment_group), 
               names_to = 'Abx', values_to = 'm')
df.doses.most.common %>% 
  group_by(Abx, Treatment_group) %>% reframe(m=mean(m)) %>% 
  pivot_wider(names_from = Treatment_group, values_from = m) %>% 
  mutate(diff=`Tac/MTX`-`PTCy/Tac/MMF`) %>% 
  arrange(desc(abs(diff)))
g <- df.doses.most.common %>% 
  ggplot(aes(x=Treatment_group, y=m, fill=Treatment_group)) + 
  geom_boxplot() + 
  facet_wrap(~Abx, scales='free') + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('') + ylab('Number of doses') +
  scale_fill_manual(values=unlist(colours$group.colours))
ggsave(g, filename=here('figures/abx/doses_most_common.pdf'),
       width = 8, height = 8, useDingbats=FALSE)

df.doses.most.common %>% 
  group_by(Abx) %>% 
  group_map(.f=function(.x, .y){
    if (length(unique(.x$Treatment_group))==1){
      tibble(p.val=NA, .y)
    } else {
      t <- wilcox.test(m~Treatment_group, data=.x)
      tibble(p.val=t$p.value, .y)
    }
  }) %>% bind_rows() %>% 
  filter(!is.na(p.val)) %>% print(n=27)

# # A tibble: 27 × 2
# p.val Abx                                            
# <dbl> <chr>                                          
# 1 0.682  Amoxicillin Clavulanate/Augmentin-PO           
# 2 0.309  Amoxicillin-PO                                 
# 3 0.813  Cefazolin-IV                                   
# 4 0.944  Cefdinir-PO                                    
# 5 0.0240 Cefepime/Maxipime-IV                           
# 6 0.0862 Cefepime/Maxipime-PO                           
# 7 0.141  Ceftriaxone/Rocephin-IV                        
# 8 0.326  Ceftriaxone/Rocephin-PO                        
# 9 0.996  Ciprofloxacin/Cipro-IV                         
# 10 0.950  Ciprofloxacin/Cipro-PO                         
# 11 0.145  Dapsone-PO                                     
# 12 0.198  Daptomycin/Cubicin-IV                          
# 13 0.0797 Doxycycline/Doryx/Doxy, adoxa-IV               
# 14 0.719  Doxycycline/Doryx/Doxy, adoxa-PO               
# 15 0.436  Levofloxacin/Levaquin-IV                       
# 16 0.720  Levofloxacin/Levaquin-PO                       
# 17 0.277  Meropenem/Merrem-IV                            
# 18 0.891  Metronidazole/Flagyl-IV                        
# 19 0.326  Metronidazole/Flagyl-NA                        
# 20 0.261  Metronidazole/Flagyl-PO                        
# 21 0.0953 PenicillinVK/Apo-pen-VK/Novo-pen-VK-PO         
# 22 0.285  Pipercillin-Tazobactam/Zosyn-IV                
# 23 0.316  Pipercillin-Tazobactam/Zosyn-PO                
# 24 0.910  Trimethoprim-Sulfamethoxazole/Bactrim/Septra-IV
# 25 0.0596 Trimethoprim-Sulfamethoxazole/Bactrim/Septra-PO
# 26 0.440  Vancomycin-IV                                  
# 27 0.886  Vancomycin-PO  

# how to systematize this? Like how can we know that these minor differences
# in a couple of patients explain the difference in absolute abundance?
# let's be honest, Ceftriaxone is the most significant one, but this is in like
# 12 patients in PTCy


# okay, let's go with correlation beween previous abx exposure and absolute 
# abundance

df.cor <- df.meta.clean %>% 
  select(Sample_ID, Participant_ID, Timepoint, copies_16S) %>% 
  left_join(df.response %>% select(Participant_ID, Treatment_group),
            by='Participant_ID') %>% 
  mutate(Timepoint=case_when(Timepoint=='PINF'~0, 
                             Timepoint=='PCON'~-5,
                             TRUE~as.numeric(as.character(Timepoint)))) %>% 
  arrange(Participant_ID, Timepoint)

df.drugs <- abx.data %>% 
  filter(!Class_2 %in% c('antiviral', 'antifungal', 'antiprotist')) %>% 
  select(Abx, Route) %>% 
  distinct() %>% 
  filter(!is.na(Abx))

df.res <- tibble(drug=character(0), route=character(0),
                 n=integer(0), p.val=double(0), coefficent=double(0),
                 enrichment=double(0), odds=double(0), n_interest=double(0))
for (i in seq_len(nrow(df.drugs))){
  tmp <- abx.data %>% 
    filter(Abx==df.drugs$Abx[i], Route==df.drugs$Route[i]) %>% 
    select(Participant_ID, Start_rel, End_rel) %>% 
    arrange(Participant_ID) 
  
  df.test <- tmp %>% 
    full_join(df.cor, by='Participant_ID', relationship = 'many-to-many') %>% 
    mutate(drug=case_when(
      ((Start_rel < (Timepoint-7)) & (End_rel > (Timepoint-7)))~'yes',
      ((Start_rel < (Timepoint)) & (End_rel > (Timepoint)))~'yes',
      TRUE~'no')) %>% 
    group_by(Participant_ID, Timepoint, Treatment_group, copies_16S) %>% 
    reframe(drug=any(drug=='yes'))
  n.samples <- df.test %>% 
    pull(drug) %>% sum
  if (n.samples > 10){
    red <- df.test %>% 
      filter(Timepoint < 30) %>% 
      filter(Timepoint > 10)
    if (sum(red$drug) < 2){
      t <- list(p.value=NA, estimate=NA)
    } else {
      t <- fisher.test(table(red$Treatment_group, red$drug))
      red %>% 
        group_by(Treatment_group, drug) %>% 
        tally() %>% 
        group_by(Treatment_group) %>% 
        mutate(n.all=sum(n)) %>% 
        mutate(freq=n/n.all) %>% 
        ggplot(aes(x=Treatment_group, y=freq, fill=drug)) + 
        geom_bar(stat='identity')
      fit.red <- coefficients(summary(lmerTest::lmer(
        copies_16S~Treatment_group+drug+(1|Participant_ID), data=red)))
      
    }
    fit <- lm(copies_16S~drug, data=df.test)
    res <- coefficients(summary(fit))
    df.res <- df.res %>% 
      add_row(drug=df.drugs$Abx[i], route=df.drugs$Route[i],
              p.val=res[2,4], coefficent=res[2,1], n=n.samples,
              enrichment=t$p.value, odds=t$estimate, n_interest=sum(red$drug))
  }
}

df.res %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coefficent, y=-log10(q.val), size=n)) + 
  geom_point() + 
  theme_bw() + 
  xlab('Effect on absolute abundance') + ylab('-log10(q-value)')

# okay, interesting
# let's look at the most common non-prophylaxis abx in the interesting
# timeframe

drugs <- df.res %>% 
  mutate(comb=paste0(drug, '-', route)) %>% 
  filter(n_interest > 10) %>%
  pull(comb)
# remove prophylaxis abx
abx.pro <- c('Ciprofloxacin/Cipro-PO', 'PenicillinVK/Apo-pen-VK/Novo-pen-VK-PO',
             'Trimethoprim-Sulfamethoxazole/Bactrim/Septra-PO',
             'Levofloxacin/Levaquin-PO', 'Dapsone-PO', 'Cefdinir-PO')
abx.treat <- setdiff(drugs, abx.pro)

g <- df.res %>% 
  mutate(comb=paste0(drug, '-', route)) %>% 
  mutate(prophylaxis=comb%in%abx.pro) %>% 
  mutate(q.val=p.adjust(p.val, method='BH')) %>% 
  ggplot(aes(x=coefficent, y=-log10(q.val), size=n, col=prophylaxis)) + 
  geom_point() + 
  theme_bw() + 
  xlab('Effect on absolute abundance') + ylab('-log10(q-value)') + 
  ggthemes::scale_colour_tableau() +
  ggrepel::geom_text_repel(aes(label=comb), size=2)
ggsave(g, filename=here('figures/abx/volcano_drug_effect.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

tmp <- abx.data %>% 
  mutate(comb=paste0(Abx, '-', Route)) %>% 
  filter(comb %in% abx.treat) %>% 
  select(Participant_ID, Start_rel, End_rel) %>% 
  arrange(Participant_ID) 

tmp2 <- tmp %>% 
  full_join(df.cor, by='Participant_ID', relationship = 'many-to-many') %>% 
  mutate(period=Timepoint-7) %>% 
  filter(!is.na(Start_rel))

subject <- IRanges::IRanges(tmp2$Start_rel, tmp2$End_rel)
query <- IRanges::IRanges(tmp2$period, tmp2$Timepoint)
p <- IRanges::pintersect(subject, query, resolve.empty='max.start')

df.test <- tmp2 %>% 
  mutate(doses=p@width) %>% 
  bind_rows(tmp %>% 
              full_join(df.cor, by='Participant_ID', 
                        relationship = 'many-to-many') %>% 
              mutate(period=Timepoint-7) %>% 
              filter(is.na(Start_rel)) %>% 
              mutate(doses=0)) %>% 
  mutate(drug=case_when(
    ((Start_rel < (Timepoint-7)) & (End_rel > (Timepoint-7)))~'yes',
    ((Start_rel < (Timepoint)) & (End_rel > (Timepoint)))~'yes',
    TRUE~'no')) %>% 
  group_by(Participant_ID, Timepoint, Treatment_group, 
           copies_16S, Sample_ID) %>% 
  reframe(drug=any(drug=='yes'), doses=sum(doses))

g <- df.test %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  group_by(Treatment_group, drug) %>% 
  tally() %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=Treatment_group, y=freq, fill=drug)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  scale_fill_manual(values=c('#007B53', '#563D82')) +
  xlab('') + ylab('Proportion')
ggsave(g, filename=here('figures/abx/proportion_non_prophylaxis_14_28.pdf'),
       width = 4, height = 4, useDingbats=FALSE)


fisher.test(table(df.test %>% 
                    filter(Timepoint > 10, Timepoint < 30) %>% 
                    pull(Treatment_group),
                  df.test %>% filter(Timepoint > 10, Timepoint < 30) %>% 
                    pull(drug)))
# Fisher's Exact Test for Count Data
# 
# data:  table(df.test %>% filter(Timepoint > 10, Timepoint < 30) %>% pull(Treatment_group), df.test %>% filter(Timepoint > 10, Timepoint < 30) %>% pull(drug))
# p-value = 0.02828
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.4973898 0.9681318
# sample estimates:
# odds ratio 
#  0.6944621 

# number of doses?

g <- df.test %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  ggplot(aes(x=Treatment_group, y=doses)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1) +
    theme_bw() + 
    xlab('') 
ggsave(g, filename=here('figures/abx/doses_non_prophylaxis_14_28.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

wilcox.test(doses~Treatment_group, 
            data=df.test %>% filter(Timepoint > 10, Timepoint < 30))

# data:  doses by Treatment_group
# W = 52762, p-value = 0.008665
# alternative hypothesis: true location shift is not equal to 0

write_tsv(df.test, file='./files/abx_exposure.tsv')


# ##############################################################################
# Test everything, while taking the antibiotic exposure into account


### alpha diversity
feat.rar <- rrarefy(t(feat.motus), 3000)
alpha <- enframe(diversity(feat.rar, index='shannon'), 
                 value='shannon', name='Sample_ID') %>% 
  full_join(enframe(diversity(feat.rar, index='invsimpson'), 
                    value='invsimpson', name='Sample_ID'), by='Sample_ID') %>% 
  full_join(enframe(diversity(feat.rar, index='simpson'), 
                    value='simpson', name='Sample_ID'), by='Sample_ID')
df.test <- df.test %>% 
  full_join(alpha, by='Sample_ID')

for (var in c('shannon', 'invsimpson', 'simpson', 'copies_16S')){
  # with/without abx
  fit.without <- lmerTest::lmer(
    paste0(var, '~Treatment_group+(1|Participant_ID)'), 
    data=df.test %>% 
      filter(Timepoint> 10, Timepoint < 30))
  fit.with <- lmerTest::lmer(
    paste0(var, '~Treatment_group+drug+(1|Participant_ID)'), 
    data=df.test %>% 
      filter(Timepoint> 10, Timepoint < 30))
  fit.with.c <- lmerTest::lmer(
    paste0(var, '~Treatment_group+doses+(1|Participant_ID)'), 
    data=df.test %>% 
      filter(Timepoint> 10, Timepoint < 30))
  p.wo <- coefficients(summary(fit.without))['Treatment_groupTac/MTX','Pr(>|t|)']
  p.wi <- coefficients(summary(fit.with))['Treatment_groupTac/MTX','Pr(>|t|)']
  p.wi.c <- coefficients(summary(fit.with.c))['Treatment_groupTac/MTX','Pr(>|t|)']
  
  message(var)
  message("P-value without Abx: ", sprintf(fmt='%.5f',p.wo))
  message("P-value with Abx: ", sprintf(fmt='%.5f',p.wi))
  message("P-value with Abx (continuous): ", sprintf(fmt='%.5f',p.wi.c))
}

# shannon
# P-value without Abx: 0.03793
# P-value with Abx: 0.08043
# P-value with Abx (continuous): 0.08736
# invsimpson
# P-value without Abx: 0.01095
# P-value with Abx: 0.02372
# P-value with Abx (continuous): 0.02318
# simpson
# P-value without Abx: 0.12696
# P-value with Abx: 0.26405
# P-value with Abx (continuous): 0.28978
# copies_16S
# P-value without Abx: 0.00003
# P-value with Abx: 0.00007
# P-value with Abx (continuous): 0.00007


#### visualize the difference in microbial load

df.m <- df.test %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  group_by(drug, Treatment_group) %>% 
  reframe(m=median(copies_16S), s=sd(copies_16S), 
          q25=quantile(copies_16S, probs=0.25),
          q75=quantile(copies_16S, probs=0.75))

g <- df.test %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  ggplot(aes(x=Treatment_group, fill=Treatment_group, y=copies_16S)) + 
  geom_boxplot(alpha=0.85, outlier.colour = NA) +
  geom_jitter(aes(colour=drug), 
              position = position_jitterdodge(jitter.width = 0.1)) +
  scale_colour_manual(values=c('#007B53', '#563D82')) +
  scale_fill_manual(values=c('#E98300', '#007C92')) +
  theme_bw() + theme(panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank()) +
  xlab('') + ylab('log10(16S copies)') + 
  geom_linerange(aes(ymin=q25, ymax=q75, y=m, x=Treatment_group,
                     col=drug), data=df.m) +
  geom_point(aes(y=m, x=Treatment_group, col=drug), data=df.m)
ggsave(g, filename=here('figures/abx/load_non_prophylaxis_14_28.pdf'),
       width = 5, height = 4, useDingbats=FALSE)


# ##############################################################################
#' Let's try this differently
#'
#' Okay, talked with Nathan again and he thinks we should probably just use
#' PO vs IV as prophylaxis vs non-prophylaxis. This might be the easier way to
#' classify antibiotics and less fraught with differences between hospitals and
#' so on. Still not perfect, since maybe participants are switched to IV 
#' because they cannot take oral medication? So nothing is perfect here, i guess
#'
#' If both approaches end up showing a similar pictures, that would be 
#' pretty strong and then it would convince myself as well

#' all the IV abx are considered to be treatment!
#' the rest of the code is the same as before
tmp <- abx.data %>% 
  filter(!Class%in%c('antiviral', 'antifungal')) %>%
  filter(Route=='IV') %>% 
  select(Participant_ID, Start_rel, End_rel) %>% 
  arrange(Participant_ID)

tmp2 <- tmp %>% 
  full_join(df.cor, by='Participant_ID', relationship = 'many-to-many') %>% 
  mutate(period=Timepoint-7) %>% 
  filter(!is.na(Start_rel))

subject <- IRanges::IRanges(tmp2$Start_rel, tmp2$End_rel)
query <- IRanges::IRanges(tmp2$period, tmp2$Timepoint)
p <- IRanges::pintersect(subject, query, resolve.empty='max.start')

df.test.iv <- tmp2 %>% 
  mutate(doses=p@width) %>% 
  bind_rows(tmp %>% 
              full_join(df.cor, by='Participant_ID', 
                        relationship = 'many-to-many') %>% 
              mutate(period=Timepoint-7) %>% 
              filter(is.na(Start_rel)) %>% 
              mutate(doses=0)) %>% 
  mutate(drug=case_when(
    ((Start_rel < (Timepoint-7)) & (End_rel > (Timepoint-7)))~'yes',
    ((Start_rel < (Timepoint)) & (End_rel > (Timepoint)))~'yes',
    TRUE~'no')) %>% 
  group_by(Participant_ID, Timepoint, Treatment_group, 
           copies_16S, Sample_ID) %>% 
  reframe(drug=any(drug=='yes'), doses=sum(doses))

g <- df.test.iv %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  group_by(Treatment_group, drug) %>% 
  tally() %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=Treatment_group, y=freq, fill=drug)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  scale_fill_manual(values=c('#007B53', '#563D82')) +
  xlab('') + ylab('Proportion')
ggsave(g, filename=here('figures/abx/proportion_non_prophylaxis_14_28_IV.pdf'),
       width = 4, height = 4, useDingbats=FALSE)


fisher.test(table(df.test.iv %>% 
                    filter(Timepoint > 10, Timepoint < 30) %>% 
                    pull(Treatment_group),
                  df.test.iv %>% filter(Timepoint > 10, Timepoint < 30) %>% 
                    pull(drug)))
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(df.test %>% filter(Timepoint > 10, Timepoint < 30) %>% pull(Treatment_group), df.test %>% filter(Timepoint > 10, Timepoint < 30) %>% pull(drug))
# p-value = 0.001283
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.412189 0.815084
# sample estimates:
# odds ratio 
#   0.580387 

#' the enrichment is stronger in the PTCy group with this approach!

# number of doses?

g <- df.test.iv %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  ggplot(aes(x=Treatment_group, y=doses)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_bw() + 
  xlab('') 
ggsave(g, filename=here('figures/abx/doses_non_prophylaxis_14_28_IV.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

wilcox.test(doses~Treatment_group, 
            data=df.test.iv %>% filter(Timepoint > 10, Timepoint < 30))

# data:  doses by Treatment_group
# W = 54012, p-value = 0.0009636
# alternative hypothesis: true location shift is not equal to 0

#' same here

df.m <- df.test.iv %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  group_by(drug, Treatment_group) %>% 
  reframe(m=median(copies_16S), s=sd(copies_16S), 
          q25=quantile(copies_16S, probs=0.25),
          q75=quantile(copies_16S, probs=0.75))

g <- df.test.iv %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  ggplot(aes(x=Treatment_group, fill=Treatment_group, y=copies_16S)) + 
  geom_boxplot(alpha=0.85, outlier.colour = NA) +
  geom_jitter(aes(colour=drug), 
              position = position_jitterdodge(jitter.width = 0.1)) +
  scale_colour_manual(values=c('#007B53', '#563D82')) +
  scale_fill_manual(values=c('#E98300', '#007C92')) +
  theme_bw() + theme(panel.grid.major.x = element_blank(),
                     panel.grid.minor = element_blank()) +
  xlab('') + ylab('log10(16S copies)') + 
  geom_linerange(aes(ymin=q25, ymax=q75, y=m, x=Treatment_group,
                     col=drug), data=df.m) +
  geom_point(aes(y=m, x=Treatment_group, col=drug), data=df.m)
ggsave(g, filename=here('figures/abx/load_non_prophylaxis_14_28_IV.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

fit <- lmerTest::lmer(copies_16S~Treatment_group+drug+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))
#                         Estimate Std. Error       df    t value      Pr(>|t|)
# (Intercept)             8.6997659 0.06230675 328.1369 139.627985 2.621489e-294
# Treatment_groupTac/MTX  0.3039697 0.07845102 266.2966   3.874643  1.344993e-04
# drugTRUE               -0.3194158 0.06759541 586.8225  -4.725407  2.877572e-06

fit <- lmerTest::lmer(copies_16S~Treatment_group+doses+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))
#                         Estimate  Std. Error       df    t value      Pr(>|t|)
# (Intercept)             8.73302871 0.059695700 317.0186 146.292424 4.694565e-293
# Treatment_groupTac/MTX  0.29295224 0.077222117 266.4189   3.793631  1.837642e-04
# doses                  -0.04187209 0.006219033 605.7389  -6.732894  3.867229e-11

fit <- lmerTest::lmer(copies_16S~Treatment_group+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))
#                        Estimate Std. Error       df    t value      Pr(>|t|)
# (Intercept)            8.5551982 0.05601974 257.1058 152.717578 2.649049e-254
# Treatment_groupTac/MTX 0.3414158 0.08035317 263.0165   4.248941  2.984152e-05

#' the same conclusion still holds!!!
#' the difference in microbial load between arms is not due to IV-given 
#' antibiotics, even though there is a over-representation of IV antibiotics 
#' in the PTCy group
#' I would think that this is a pretty strong indication that Abx are not the
#' driving force here, at least not alone


#' Should probably also test alpha diversity again
df.test.iv <- df.test.iv %>% full_join(alpha, by='Sample_ID')
fit <- lmerTest::lmer(shannon~Treatment_group+drug+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # no
fit <- lmerTest::lmer(shannon~Treatment_group+doses+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # no 
fit <- lmerTest::lmer(invsimpson~Treatment_group+drug+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # yes
fit <- lmerTest::lmer(invsimpson~Treatment_group+doses+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # yes
fit <- lmerTest::lmer(simpson~Treatment_group+drug+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # no
fit <- lmerTest::lmer(simpson~Treatment_group+doses+(1|Participant_ID), 
                      data=df.test.iv %>% filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))['Treatment_groupTac/MTX','Pr(>|t|)'] # no

#' same picture as above: the difference in Inverse Simpson remains significant
#' also when including antibiotics, not so for the other diversity metrics
#' This shows again that maybe inverse simpson might have been the better
#' metric to use, but we had to use Shannon because of the prespecified 
#' analyses plan. However, Peled et al. also used Inverse Simpson, so for
#' comparison it makes sense to include it as well



# ##############################################################################
#' Edit:
#' Reviewer 4 mentioned that infections might be additional confounders
#' 
#' we could add this here?
#' 
#' not really, since the infection data is not yet fully adjudicated
#' however, if we add the preliminary infection data, do they have a significant
#' contribution?

# get infection data?
df.inf <- read_csv('./files/inf.csv') %>%  # This is data from Mike
  rename(Participant_ID=PATID, infection_day=binfday)

df.response %>% 
  select(Participant_ID, Treatment_group) %>% 
  left_join(df.inf) %>% 
  ggplot(aes(x=infection_day, fill=Treatment_group)) + 
  geom_histogram(bins=100) + 
  facet_grid(Treatment_group~.)

tmp.2 <- df.inf %>% 
  mutate(Start_rel=infection_day, End_rel=infection_day+7) %>% 
  arrange(Participant_ID)

df.test.2 <- tmp.2 %>% 
  full_join(df.cor, by='Participant_ID', relationship = 'many-to-many') %>% 
  mutate(infection=case_when(
    ((Start_rel < (Timepoint-7)) & (End_rel > (Timepoint-7)))~'yes',
    ((Start_rel < (Timepoint)) & (End_rel > (Timepoint)))~'yes',
    TRUE~'no')) %>% 
  group_by(Participant_ID, Timepoint, Treatment_group, 
           copies_16S, Sample_ID) %>% 
  reframe(infection=any(infection=='yes'))

g <- df.test.2 %>% 
  filter(Timepoint > 10, Timepoint < 30) %>% 
  group_by(Treatment_group, infection) %>% 
  tally() %>% 
  group_by(Treatment_group) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=Treatment_group, y=freq, fill=infection)) + 
  geom_bar(stat='identity') + 
  theme_bw() +
  scale_fill_manual(values=c('#202864', '#D23264')) +
  xlab('') + ylab('Proportion') + 
  coord_flip() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank())
ggsave(g, filename='./figures/revision/infections_bias.png',
       width = 6, height = 2)

# also more common in PTCy

fisher.test(table(df.test.2 %>% 
                     filter(Timepoint > 10, Timepoint < 30) %>% 
                     pull(Treatment_group),
                   df.test.2 %>% filter(Timepoint > 10, Timepoint < 30) %>% 
                     pull(infection)))
# p-value = 0.001685
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.2450301 0.7591764
# sample estimates:
#   odds ratio 
# 0.4373416
# yes, this is significant

fit <- lmerTest::lmer(copies_16S~Treatment_group+drug+infection+(1|Participant_ID), 
                      data=df.test %>% 
                        full_join(df.test.2 %>% 
                                    select(Sample_ID, infection), 
                                  by='Sample_ID') %>% 
                        filter(Timepoint> 10, Timepoint < 30))
coefficients(summary(fit))
# Estimate Std. Error       df    t value      Pr(>|t|)
# (Intercept)             8.7282552 0.06266124 329.6833 139.292736 5.126462e-295
# Treatment_groupTac/MTX  0.3138135 0.07770614 265.3585   4.038464  7.047726e-05
# drugTRUE               -0.3593372 0.06628659 579.7957  -5.420964  8.708413e-08

# are those two entangled
df.test %>% 
  full_join(df.test.2 %>% 
              select(Sample_ID, infection), 
            by='Sample_ID') %>% 
  filter(Timepoint> 10, Timepoint < 30) %>% 
  group_by(drug, infection) %>% 
  tally() %>% 
  ggplot(aes(x=drug, y=n, fill=infection)) +
  geom_bar(stat='identity')
fisher.test(table(df.test %>% 
                    full_join(df.test.2 %>% 
                                select(Sample_ID, infection), 
                              by='Sample_ID') %>% 
                    filter(Timepoint> 10, Timepoint < 30) %>% pull(drug),
                  df.test %>% 
                    full_join(df.test.2 %>% 
                                select(Sample_ID, infection), 
                              by='Sample_ID') %>% 
                    filter(Timepoint> 10, Timepoint < 30) %>% pull(infection)))
#' yes, of course
#' it would have been weird if not
#' so I guess it's fine to not include infections, since antibiotics and 
#' infections are anyway correlated, so much of the possible 'confounding' from
#' infections is already covered by abx treatment
