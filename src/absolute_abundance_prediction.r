# ##############################################################################
#
## Script to predict 16S copy number from measured samples
#
# ##############################################################################

library("tidyverse")
library("mlr3")
library("mlr3learners")
library("mlr3tuning")
library("here")

set.seed(2024)

# expected number of copies for the NIST controls
# from Boryana's paper
nist.exp <- tibble(Sample_ID=c("PosMix", "PosO", 'PosE'), 
                   expected_copies=c(77160000, 61440000, 82450000))

# ##############################################################################
# load all needed data

# sample metadata
df.meta <- read_tsv(here('data/meta_sample.tsv'), 
                    col_types = cols())

# sequencing information
df.sequencing <- read_tsv(here('data/sequencing_stats.tsv'),
                          col_types = cols())

# 16S raw data
df.ddPCR <- read_tsv(here('data/ddPCR_results_batch1.tsv'), col_types=cols()) %>% 
  mutate(batch='batch1') %>% 
  bind_rows(read_tsv(here('data/ddPCR_results_batch2.tsv'), col_types=cols()) %>% 
              mutate(batch='batch2'))
  
# 16S plate layout
df.ddPCR.layout <- read_tsv(here('data/ddPCR_layout_batch1.tsv'), 
                            col_types=cols()) %>% 
  bind_rows(read_tsv('./data/ddPCR_layout_batch2.tsv', 
                     col_types=cols())) 

# metaphlan data
feat.mpa <- read.table(here('data/metaphlan_all.tsv'), 
                       sep = '\t', comment.char = '', 
                       quote = '', header = TRUE,
                       row.names = 1)
rownames(feat.mpa)[1] <- 'x|s__UNCLASSIFIED'

feat.mpa.species <- feat.mpa[str_detect(rownames(feat.mpa), 's__'),]
feat.mpa.species <- feat.mpa.species[str_detect(rownames(feat.mpa.species), 
                                                't__', negate = TRUE),]
alpha <- vegan::diversity(t(feat.mpa.species), index = 'shannon')

# ##############################################################################
# a bit of data exploration

# join everything 
df.table <- full_join(df.ddPCR, df.ddPCR.layout, by=c('Well', 'Plate')) %>% 
  left_join(df.meta, by='Sample_ID') %>% 
  mutate(copies_reaction=log(No_droplets/No_negatives) * 1.0/0.795 * 
           1000 * 22) %>% 
  filter(Sample_ID != 'Pipetting_error')

# filter based on the no-template control
copies_lob <- df.table %>% 
  group_by(Plate) %>% 
  filter(Sample_ID=='Neg') %>% 
  reframe(neg_max=max(copies_reaction))
copies_lob_min <- df.table %>% 
  group_by(Plate) %>% 
  filter(Sample_ID=='Neg') %>% 
  reframe(neg_min=min(copies_reaction))
copies_rxn_NTC_span <- copies_lob %>% 
  full_join(copies_lob_min, by='Plate') %>% 
  mutate(span=case_when(neg_min==0~0, 
                        TRUE~neg_max/neg_min)) %>% 
  mutate(loq=neg_max*4)

high.neg <- copies_rxn_NTC_span %>% 
  filter(neg_max > 25)
if (nrow(high.neg) > 0){
  message("Something went wrong with plate(s): ", 
          paste(high.neg$Plate, sep = ', '))
}
high.noise <- copies_rxn_NTC_span %>% 
  filter(span > 10)
if (nrow(high.noise) > 0){
  message("Something went wrong with plate(s): ", 
          paste(high.noise$Plate, sep = ', '))
}

# compute the number of copies
df.table <- df.table %>% 
  filter(!is.na(No_droplets)) %>% 
  mutate(Not_enough_droplets=case_when(No_droplets < 10000 ~ "not_enough", 
                                       TRUE~"pass")) %>% 
  mutate(Too_concentrated=case_when(No_negatives < 10~'not_concentrated',
                                    TRUE~"pass")) %>% 
  mutate(Too_dilute=case_when(
    copies_reaction < max(copies_rxn_NTC_span$loq)~"too_dilute",
    TRUE~'pass')) %>% 
  left_join(nist.exp, by='Sample_ID') %>% 
  filter(!is.na(dilution)) %>% 
  mutate(copies_microL=copies_reaction * dilution * 1/6) %>% 
  mutate(measured_vs_expected=copies_microL/expected_copies) %>% 
  mutate(flag=case_when(measured_vs_expected < 1/5~'too_low',
                        measured_vs_expected > 5 ~ 'too high',
                        TRUE~'within_range'))

# negative/positive controls
df.table %>% 
  filter(Not_enough_droplets=='pass') %>% 
  filter(str_detect(Sample_ID, 'Pos')) %>% 
  ggplot(aes(x=Plate, y=measured_vs_expected, fill=Sample_ID)) + 
    geom_boxplot()
# everything within the expected range! Yay!

df.table %>% 
  filter(Not_enough_droplets=='pass') %>% 
  filter(str_detect(Sample_ID, 'P[0-9]', negate=TRUE)) %>% 
  ggplot(aes(x=Sample_ID, y=copies_reaction, fill=Sample_ID)) + 
  geom_boxplot() + 
  facet_grid(~Plate, scales = 'free')
# Same here: everything looks good! Yay!

# filter out everything that failed

df.all <- df.table %>% 
  filter(str_detect(Sample_ID, 'P[0-9]')) %>% 
  filter(Not_enough_droplets=='pass') %>% 
  mutate(copies_reaction=case_when(Too_dilute=='pass'~copies_reaction,
                                   TRUE~50)) %>% 
  select(-c(expected_copies, copies_microL, measured_vs_expected, flag)) %>% 
  mutate(copies_per_extraction=copies_reaction*dilution*1/6*100) 

# correlations with other measurements
# DNA concentration from Novogene?
spearman.rho <- cor(df.all$copies_per_extraction, 
                    df.all$DNA_concentration, 
                    method='spearman')
x <- cor.test(df.all$copies_per_extraction, 
              df.all$DNA_concentration, 
              method='spearman', exact=FALSE)

g <- df.all %>%
  ggplot(aes(x=DNA_concentration,
             y=log10(copies_per_extraction))) + 
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  annotate(x=-Inf, y=-Inf, geom='text', 
           hjust=-2, vjust=-20,
           label=paste0('rho=', sprintf(fmt='%.2f', spearman.rho))) + 
  xlab('DNA concentration [ng/µL]') +
  ylab('16S copies per extraction') + 
  geom_smooth(method='lm', formula='y~log(x)')
ggsave(g, filename=here('figures/ddPCR', 'correlation_concentration.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

ggsave(g + facet_grid(~batch), 
       filename=here('figures/ddPCR', 'correlation_concentration_batch.pdf'),
       width = 9, height = 5, useDingbats=FALSE)

# why did we run the second batch?
x <- df.meta %>% 
  left_join(df.table %>% select(Sample_ID, batch), by='Sample_ID') %>% 
  filter(Timepoint!='EXTRA') %>% 
  mutate(batch=case_when(is.na(batch)~'not_measured', TRUE~batch)) %>% 
  mutate(first_comp=case_when(batch=='batch2'~'not_measured', TRUE~batch)) %>%
  mutate(second_comp=case_when(batch=='batch2'~'batch1', TRUE~batch)) 
g <- x %>% 
  ggplot(aes(x=batch, y=log10(DNA_concentration))) + 
    geom_boxplot(outlier.shape = NA) +
    theme_bw() + xlab('') + ylab('log10(DNA concentration)') + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank()) + 
    geom_jitter(width = 0.1) +
    NULL
t.test(log10(DNA_concentration)~first_comp, data=x)  # 1.621e-08
t.test(log10(DNA_concentration)~second_comp, data=x)  # 9.864e-05 # still different
t.test(log10(DNA_concentration)~batch, data=x %>% filter(batch!='batch1'))  # 0.9204 # batch two is very representative
ggsave(g, filename=here('figures/ddPCR', 'batch_bias.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

# ##############################################################################
# train the model

.f_train_model <- function(data){
  task.pred <- as_task_regr(data, target='copies_per_extraction')
  lrn.regr <- lrn('regr.ranger', importance='impurity')
  rcv10 <- rsmp("repeated_cv", repeats = 10, folds = 10)
  rr <- resample(task.pred, lrn.regr, rcv10, store_models = TRUE)
  rsq <- rr$aggregate(msr('regr.rsq'))
  x <- rr$prediction()
  df.plot <- tibble(row_ids=x$row_ids,
                    truth=x$truth,
                    response=x$response) %>% 
    group_by(row_ids) %>% 
    summarise(t=mean(truth), r=mean(response))
  pearson.r <- cor(df.plot$t, df.plot$r)
  x <- epiR::epi.ccc(df.plot$t, df.plot$r)
  g <- df.plot %>% 
    ggplot(aes(x=t, y=r)) + 
    geom_point() + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    annotate(x=-Inf, y=-Inf, geom='text', 
             hjust=-0.7, vjust=-30,
             label=paste0('R^2=', sprintf(fmt='%.2f', rsq), "; Pearson's r=",
                          sprintf(fmt='%.2f', pearson.r))) +
    ylab('Predicted 16S copies per extraction') + 
    xlab('Real 16S copies per extraction')
  return(list('plot'=g, 'models'=rr, 'ccc'=x))
}

# prep the feature table

df.feat <- df.all %>% 
  left_join(as_tibble(feat.mpa[str_detect(rownames(feat.mpa), 
                                          pattern = 'p__', negate = TRUE),], 
                      rownames='kingdom') %>% 
              pivot_longer(-kingdom, names_to = 'Sample_ID') %>% 
              pivot_wider(values_from = value, names_from = 'kingdom'), 
            by='Sample_ID') %>% 
  left_join(df.sequencing %>% select(Sample_ID, hostremoved_frac),
            by='Sample_ID') %>% 
  select(Sample_ID, DNA_concentration,
         k__Bacteria, k__Archaea, k__Eukaryota, copies_per_extraction,
         `x|s__UNCLASSIFIED`, hostremoved_frac, Sample_type, Timepoint) %>% 
  mutate(Sample_type=as.factor(Sample_type)) %>% 
  mutate(Timepoint=factor(Timepoint, levels=c(
    'PCON', 'PINF', sort(unique(as.numeric(Timepoint))))))
colnames(df.feat) <- make.names(colnames(df.feat))


# train for DNA concentration only
rr.conc <- .f_train_model(
  df.feat %>% 
    mutate(copies_per_extraction=log10(copies_per_extraction)) %>% 
    select(copies_per_extraction, DNA_concentration))
ggsave(rr.conc$plot, filename=here('figures/ddPCR', 'prediction_dna_only.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

# train the full model
rr.full <- .f_train_model(
  df.feat %>% 
    mutate(copies_per_extraction=log10(copies_per_extraction)) %>% 
    select(-Sample_ID))
ggsave(rr.full$plot, filename=here('figures/ddPCR', 'prediction_full.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

rr.full$ccc$rho.c
#         est     lower     upper
# 1 0.9045751 0.8926906 0.9152023
# ##############################################################################
# predict new data

df.ext <- df.meta %>% 
  select(-Participant_ID) %>% 
  left_join(as_tibble(feat.mpa[str_detect(rownames(feat.mpa), 
                                          pattern = 'p__', negate = TRUE),], 
                      rownames='kingdom') %>% 
              pivot_longer(-kingdom, names_to = 'Sample_ID') %>% 
              pivot_wider(values_from = value, names_from = 'kingdom'), 
            by='Sample_ID') %>% 
  left_join(df.sequencing %>% select(Sample_ID, hostremoved_frac),
            by='Sample_ID') %>% 
  select(Sample_ID, DNA_concentration,
         k__Bacteria, k__Archaea, k__Eukaryota, 
         `x|s__UNCLASSIFIED`, hostremoved_frac, Sample_type, Timepoint) %>% 
  mutate(Sample_type=as.factor(Sample_type)) %>% 
  mutate(Timepoint=factor(Timepoint, levels=c(
    'PCON', 'PINF', sort(unique(as.numeric(Timepoint)))))) %>% 
  filter(!Sample_ID %in% df.feat$Sample_ID) %>% 
  mutate(copies_per_extraction=0)
colnames(df.ext) <- make.names(colnames(df.ext))

task.pred <- as_task_regr(df.ext, target='copies_per_extraction')
pred.matrix <- matrix(NA, nrow = nrow(df.ext), 
                      ncol=length(rr.full$models$learners))
rownames(pred.matrix) <- df.ext$Sample_ID
colnames(pred.matrix) <- paste0('model_', seq(ncol(pred.matrix)))
for (i in seq(ncol(pred.matrix))){
  pred <- rr.full$models$learners[[i]]$predict(task=task.pred)
  pred.matrix[,i] <- pred$response
}

copies <- enframe(rowMeans(pred.matrix), name='Sample_ID',
                  value='copies_16S') %>% 
  mutate(copies='predicted') %>%
  bind_rows(
    df.all %>% 
      transmute(Sample_ID, copies_16S=log10(copies_per_extraction), 
                copies='measured'))

# export
write_tsv(copies, file = './data/copies_16S.tsv')
