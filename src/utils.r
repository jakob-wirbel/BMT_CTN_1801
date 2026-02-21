# from the survival package
smat <- function(x) {
  # the rest of the routine is simpler if everything is a matrix
  dd <- dim(x)
  if (is.null(dd)) as.matrix(x)
  else if (length(dd) ==2) x
  else matrix(x, nrow=dd[1])
}

# plot the outcome of an Aalen-Johanson estimator
.f_plot_ajfit <- function(tac, pty){
  
  df.plot.tac <- tibble(time=tac$time, cumhaz=smat(tac$cumhaz)[,1], 
                        type=c(rep('absent', tac$strata[1]), 
                               rep('present', tac$strata[2]))) %>% 
    mutate(change=c(1, diff(cumhaz))) %>% 
    filter(change!=0) %>% 
    select(-change)
  
  x.a <- df.plot.tac %>% 
    filter(type=='absent')
  x.b <- df.plot.tac %>% 
    filter(type!='absent')
  
  df.plot.tac <- df.plot.tac %>% 
    bind_rows(tibble(time=x.a$time[-1],
                     cumhaz=x.a$cumhaz[-nrow(x.a)],
                     type='absent')) %>% 
    bind_rows(tibble(time=x.b$time[-1],
                     cumhaz=x.b$cumhaz[-nrow(x.b)],
                     type='present')) %>% 
    bind_rows(tibble(time=365, cumhaz=c(x.a$cumhaz[nrow(x.a)], 
                                        x.b$cumhaz[nrow(x.b)]),
                     type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=0, cumhaz=0, type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=c(x.a$time[1], x.b$time[1]), cumhaz=c(0, 0), 
                     type=c('absent', 'present'))) %>% 
    arrange(type, time, cumhaz) 
  
  df.plot.pty <- tibble(time=pty$time, cumhaz=smat(pty$cumhaz)[,1], 
                        type=c(rep('absent', pty$strata[1]), 
                               rep('present', pty$strata[2]))) %>% 
    mutate(change=c(1, diff(cumhaz))) %>% 
    filter(change!=0) %>% 
    select(-change)
  
  x.a <- df.plot.pty %>% 
    filter(type=='absent')
  x.b <- df.plot.pty %>% 
    filter(type!='absent')
  
  df.plot.pty <- df.plot.pty %>% 
    bind_rows(tibble(time=x.a$time[-1],
                     cumhaz=x.a$cumhaz[-nrow(x.a)],
                     type='absent')) %>% 
    bind_rows(tibble(time=x.b$time[-1],
                     cumhaz=x.b$cumhaz[-nrow(x.b)],
                     type='present')) %>% 
    bind_rows(tibble(time=365, cumhaz=c(x.a$cumhaz[nrow(x.a)], 
                                        x.b$cumhaz[nrow(x.b)]),
                     type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=0, cumhaz=0, type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=c(x.a$time[1], x.b$time[1]), cumhaz=c(0, 0), 
                     type=c('absent', 'present'))) %>% 
    arrange(type, time, cumhaz) 
  
  df.plot.all <- bind_rows(df.plot.tac %>% mutate(group='Tac'),
                           df.plot.pty %>% mutate(group='PTCy'))
  
  g <- df.plot.all %>% 
    ggplot(aes(x=time, y=cumhaz, colour=group)) + 
    geom_line(aes(lty=type))
  return(g)
}

# plot the outcome of cox ph model
.f_plot_coxph <- function(tac, pty){
  
  df.plot.tac <- tibble(time=tac$time, surv=smat(tac$surv)[,1], 
                        type=c(rep('absent', tac$strata[1]), 
                               rep('present', tac$strata[2]))) %>% 
    mutate(change=c(1, diff(surv))) %>% 
    filter(change!=0) %>% 
    select(-change)
  
  x.a <- df.plot.tac %>% 
    filter(type=='absent')
  x.b <- df.plot.tac %>% 
    filter(type!='absent')
  
  df.plot.tac <- df.plot.tac %>% 
    bind_rows(tibble(time=x.a$time[-1],
                     surv=x.a$surv[-nrow(x.a)],
                     type='absent')) %>% 
    bind_rows(tibble(time=x.b$time[-1],
                     surv=x.b$surv[-nrow(x.b)],
                     type='present')) %>% 
    bind_rows(tibble(time=365, surv=c(x.a$surv[nrow(x.a)], 
                                        x.b$surv[nrow(x.b)]),
                     type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=0, surv=1, type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=c(x.a$time[1], x.b$time[1]), surv=1, 
                     type=c('absent', 'present'))) %>% 
    arrange(type, time, desc(surv))
  
  df.plot.pty <- tibble(time=pty$time, surv=smat(pty$surv)[,1], 
                        type=c(rep('absent', pty$strata[1]), 
                               rep('present', pty$strata[2]))) %>% 
    mutate(change=c(1, diff(surv))) %>% 
    filter(change!=0) %>% 
    select(-change)
  
  x.a <- df.plot.pty %>% 
    filter(type=='absent')
  x.b <- df.plot.pty %>% 
    filter(type!='absent')
  
  df.plot.pty <- df.plot.pty %>% 
    bind_rows(tibble(time=x.a$time[-1],
                     surv=x.a$surv[-nrow(x.a)],
                     type='absent')) %>% 
    bind_rows(tibble(time=x.b$time[-1],
                     surv=x.b$surv[-nrow(x.b)],
                     type='present')) %>% 
    bind_rows(tibble(time=365, surv=c(x.a$surv[nrow(x.a)], 
                                        x.b$surv[nrow(x.b)]),
                     type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=0, surv=1, type=c('absent', 'present'))) %>% 
    bind_rows(tibble(time=c(x.a$time[1], x.b$time[1]), surv=1, 
                     type=c('absent', 'present'))) %>% 
    arrange(type, time, desc(surv)) 
  
  df.plot.all <- bind_rows(df.plot.tac %>% mutate(group='Tac'),
                           df.plot.pty %>% mutate(group='PTCy'))
  
  g <- df.plot.all %>% 
    ggplot(aes(x=time, y=surv, colour=group)) + 
    geom_line(aes(lty=type))
  return(g)
}

# test for associations
.f_test_outcome <- function(pred, landmarked=NA){
  require('survival')
  if (!is.na(landmarked)){
    # browser()
    if (!is.numeric(landmarked)){
      stop("Landmark needs to be numeric!")
    }
    pred$Relapse_days <- pred$Relapse_days - landmarked
    pred$cGVHD_days <- pred$cGVHD_days - landmarked
    pred$aGVHD24_days <- pred$aGVHD24_days - landmarked
    pred$aGVHD34_days <- pred$aGVHD34_days - landmarked
    pred$OS_days <- pred$OS_days - landmarked
    pred$GRFS_days <- pred$GRFS_days - landmarked
  }
  
  
  df.hr <- tibble(outcome=character(0),
                  estimate=double(0), high=double(0), low=double(0),
                  pval=double(0), set=character(0), 
                  n_present=double(0), n_absent=double(0))
  
  n.present <- pred %>% group_by(m) %>% tally()
  n.present.group <- pred %>% group_by(m, Treatment_group) %>% tally()
  
  # NRM
  fdata1 <- finegray(Surv(Relapse_days, Relapse_or_death) ~ ., 
                     data = pred %>% filter(Relapse_days > 0), 
                     etype = 'Death')
  fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ Treatment_group + m, 
                  data = fdata1, weight = fgwt)
  x <- summary(fgfit1)
  n.present <- pred %>% filter(Relapse_days > 0) %>% group_by(m) %>% tally()
  n.present.group <- pred %>% filter(Relapse_days > 0) %>% 
    group_by(m, Treatment_group) %>% tally()
  df.hr <- df.hr %>% 
    add_row(outcome='NRM', set='all',
            high=x$conf.int['mTRUE','upper .95'], 
            low=x$conf.int['mTRUE','lower .95'], 
            estimate=x$conf.int['mTRUE','exp(coef)'], 
            pval=x$coefficients['mTRUE','Pr(>|z|)'],
            n_present=n.present %>% filter(m) %>% pull(n),
            n_absent=n.present %>% filter(!m) %>% pull(n))
  fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ m, 
                  data = fdata1, weight = fgwt)
  x <- summary(fgfit1)
  df.hr <- df.hr %>% 
    add_row(outcome='NRM', set='no_group',
            high=x$conf.int['mTRUE','upper .95'], 
            low=x$conf.int['mTRUE','lower .95'], 
            estimate=x$conf.int['mTRUE','exp(coef)'], 
            pval=x$coefficients['mTRUE','Pr(>|z|)'],
            n_present=n.present %>% filter(m) %>% pull(n),
            n_absent=n.present %>% filter(!m) %>% pull(n))
  for (grouping in unique(pred$Treatment_group)){
    fdata1 <- finegray(Surv(Relapse_days, Relapse_or_death) ~ ., 
                       data = pred %>% filter(Treatment_group==grouping) %>% 
                         filter(Relapse_days > 0), 
                       etype = 'Death')
    fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ m, 
                    data = fdata1, weight = fgwt)
    x <- summary(fgfit1)
    df.hr <- df.hr %>% 
      add_row(outcome='NRM', set=grouping,
              high=x$conf.int['mTRUE','upper .95'], 
              low=x$conf.int['mTRUE','lower .95'], 
              estimate=x$conf.int['mTRUE','exp(coef)'], 
              pval=x$coefficients['mTRUE','Pr(>|z|)'],
              n_present=n.present.group %>% filter(m) %>% 
                filter(Treatment_group==grouping) %>%  pull(n),
              n_absent=n.present.group %>% filter(!m) %>% 
                filter(Treatment_group==grouping) %>% pull(n))
  }
  
  pred$cGVHD_MS_days <- pred$cGVHD_days
  pred$cGVHD_MS_or_death <- pred$cGVHD_ms_or_death
  all.plots <- list()
  # browser()
  for (var in c('Relapse', 'cGVHD', 'cGVHD_MS', 'aGVHD24', 'aGVHD34')){
    data <- pred %>% filter(!!sym(paste0(var, '_days')) > 0)
    n.present <- data %>% group_by(m) %>% tally()
    n.present.group <- data %>% group_by(m, Treatment_group) %>% tally()
    fdata1 <- finegray(Surv(data[[paste0(var, '_days')]], 
                            data[[paste0(var, '_or_death')]])~., 
                       data=data, etype = var)
    fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ Treatment_group + m, 
                    data = fdata1, weight = fgwt)
    x <- summary(fgfit1)
    
    # Aalen-Johansen curves
    tmp <- pred %>% 
      filter(Treatment_group=='PTCy/Tac/MMF')
    ajfit.pty <- survfit(Surv(tmp[[paste0(var, '_days')]], 
                              tmp[[paste0(var, '_or_death')]])~tmp$m)
    tmp <- pred %>% 
      filter(Treatment_group=='Tac/MTX')
    ajfit.tac <- survfit(Surv(tmp[[paste0(var, '_days')]], 
                              tmp[[paste0(var, '_or_death')]])~tmp$m)
    all.plots[[var]] <- .f_plot_ajfit(ajfit.tac, ajfit.pty) + ggtitle(var)
    df.hr <- df.hr %>% 
      add_row(outcome=var, set='all',
              high=x$conf.int['mTRUE','upper .95'], 
              low=x$conf.int['mTRUE','lower .95'], 
              estimate=x$conf.int['mTRUE','exp(coef)'], 
              pval=x$coefficients['mTRUE','Pr(>|z|)'],
              n_present=n.present %>% filter(m) %>% pull(n),
              n_absent=n.present %>% filter(!m) %>% pull(n))
    fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ m, 
                    data = fdata1, weight = fgwt)
    x <- summary(fgfit1)
    df.hr <- df.hr %>% 
      add_row(outcome=var, set='no_group',
              high=x$conf.int['mTRUE','upper .95'], 
              low=x$conf.int['mTRUE','lower .95'], 
              estimate=x$conf.int['mTRUE','exp(coef)'], 
              pval=x$coefficients['mTRUE','Pr(>|z|)'],
              n_present=n.present %>% filter(m) %>% pull(n),
              n_absent=n.present %>% filter(!m) %>% pull(n))
    for (grouping in unique(pred$Treatment_group)){
      pred.g <- pred %>% 
        filter(!!sym(paste0(var, '_days')) > 0) %>% 
        filter(Treatment_group==grouping)
      fdata1 <- finegray(Surv(pred.g[[paste0(var, '_days')]], 
                              pred.g[[paste0(var, '_or_death')]])~., 
                         data=pred.g, etype = var)
      fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~  m, 
                      data = fdata1, weight = fgwt)
      x <- summary(fgfit1)
      df.hr <- df.hr %>% 
        add_row(outcome=var, set=grouping,
                high=x$conf.int['mTRUE','upper .95'], 
                low=x$conf.int['mTRUE','lower .95'], 
                estimate=x$conf.int['mTRUE','exp(coef)'], 
                pval=x$coefficients['mTRUE','Pr(>|z|)'],
                n_present=n.present.group %>% filter(m) %>% 
                  filter(Treatment_group==grouping) %>%  pull(n),
                n_absent=n.present.group %>% filter(!m) %>% 
                  filter(Treatment_group==grouping) %>% pull(n))
    }
  }
  
  # GRFS
  # browser()
  for (var in c('GRFS', 'OS')){
    data <- pred %>% filter(!!sym(paste0(var, '_days')) > 0)
    n.present <- data %>% group_by(m) %>% tally()
    n.present.group <- data %>% group_by(m, Treatment_group) %>% tally()
    coxphfit <- coxph(Surv(data[[paste0(var, '_days')]], 
                           data[[paste0(var, '_event')]]) ~ 
                        data$Treatment_group + data$m)
    x <- summary(coxphfit)
    df.hr <- df.hr %>% 
      add_row(outcome=var, set='all',
              high=x$conf.int['data$mTRUE','upper .95'], 
              low=x$conf.int['data$mTRUE','lower .95'], 
              estimate=x$conf.int['data$mTRUE','exp(coef)'], 
              pval=x$coefficients['data$mTRUE','Pr(>|z|)'],
              n_present=n.present %>% filter(m) %>% pull(n),
              n_absent=n.present %>% filter(!m) %>% pull(n))
    coxphfit <- coxph(Surv(data[[paste0(var, '_days')]], 
                           data[[paste0(var, '_event')]]) ~ 
                        data$m)
    x <- summary(coxphfit)
    df.hr <- df.hr %>% 
      add_row(outcome=var, set='no_group',
              high=x$conf.int['data$mTRUE','upper .95'], 
              low=x$conf.int['data$mTRUE','lower .95'], 
              estimate=x$conf.int['data$mTRUE','exp(coef)'], 
              pval=x$coefficients['data$mTRUE','Pr(>|z|)'],
              n_present=n.present %>% filter(m) %>% pull(n),
              n_absent=n.present %>% filter(!m) %>% pull(n))
    
    tmp <- data %>% 
      filter(Treatment_group=='PTCy/Tac/MMF')
    coxfit.pty <- survfit(Surv(tmp[[paste0(var, '_days')]], 
                              tmp[[paste0(var, '_event')]])~tmp$m)
    tmp <- data %>% 
      filter(Treatment_group=='Tac/MTX')
    coxfit.tac <- survfit(Surv(tmp[[paste0(var, '_days')]], 
                              tmp[[paste0(var, '_event')]])~tmp$m)
    all.plots[[var]] <- .f_plot_coxph(coxfit.tac, coxfit.pty) + ggtitle(var)
    
    
    for (grouping in unique(pred$Treatment_group)){
      pred.g <- data %>% 
        filter(Treatment_group==grouping)
      coxphfit <- coxph(Surv(pred.g[[paste0(var, '_days')]], 
                             pred.g[[paste0(var, '_event')]]) ~ pred.g$m)
      x <- summary(coxphfit)
      df.hr <- df.hr %>% 
        add_row(outcome=var, set=grouping,
                high=x$conf.int['pred.g$mTRUE','upper .95'], 
                low=x$conf.int['pred.g$mTRUE','lower .95'], 
                estimate=x$conf.int['pred.g$mTRUE','exp(coef)'], 
                pval=x$coefficients['pred.g$mTRUE','Pr(>|z|)'],
                n_present=n.present.group %>% filter(m) %>% 
                  filter(Treatment_group==grouping) %>%  pull(n),
                n_absent=n.present.group %>% filter(!m) %>% 
                  filter(Treatment_group==grouping) %>% pull(n))
    }
  }
  return(list('associations'=df.hr, 'plots'=all.plots))
}

