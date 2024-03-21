
# Define the mean delta LEMS week 04 to recovery (week 52 with LOCF when week 52 is missing) below NLI
# compute the delta at the myotome level

# LOCF for individual myotomes
myotomes_EMSCI <- c('RMS_C5','RMS_C6','RMS_C7','RMS_C8', 'RMS_T1',
                    'LMS_C5','LMS_C6','LMS_C7','LMS_C8', 'LMS_T1',
                    'RMS_L2','RMS_L3', 'RMS_L4','RMS_L5', 'RMS_S1',
                    'LMS_L2','LMS_L3', 'LMS_L4','LMS_L5', 'LMS_S1',
                    'Patientennummer')

raw_EMSCI_all_myo_26 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute III') %>%
  dplyr::select(all_of(myotomes_EMSCI))
names(raw_EMSCI_all_myo_26)[1:20] <- paste0(names(raw_EMSCI_all_myo_26)[1:20], '_week26')
raw_EMSCI_all_myo_52 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'chronic') %>%
  dplyr::select(all_of(myotomes_EMSCI))
names(raw_EMSCI_all_myo_52)[1:20] <- paste0(names(raw_EMSCI_all_myo_52)[1:20], '_week52')
raw_EMSCI_all_myo_4 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute I') %>%
  dplyr::select(all_of(myotomes_EMSCI))
names(raw_EMSCI_all_myo_4)[1:20] <- paste0(names(raw_EMSCI_all_myo_4)[1:20], '_week4')

wide_EMSCI_myo_4_26 <- merge(raw_EMSCI_all_myo_4, raw_EMSCI_all_myo_26, by=c('Patientennummer'))
wide_EMSCI_myo_4_26_52 <- merge(wide_EMSCI_myo_4_26, raw_EMSCI_all_myo_52, by=c('Patientennummer'))
wide_EMSCI_myo_4_26_52[wide_EMSCI_myo_4_26_52 == ''] <- NA

for (col in myotomes_EMSCI[1:20]){
  wide_EMSCI_myo_4_26_52[paste0(col, '_locf')] <- NA
}

for (i in c(1:dim(wide_EMSCI_myo_4_26_52)[1])){
  for (col in myotomes_EMSCI[1:20]){
    if (!(is.na(wide_EMSCI_myo_4_26_52[i, paste0(col, '_week52')]))){
      wide_EMSCI_myo_4_26_52[i, paste0(col, '_locf')] <- wide_EMSCI_myo_4_26_52[i, paste0(col, '_week52')]
    } else {
      wide_EMSCI_myo_4_26_52[i, paste0(col, '_locf')] <- wide_EMSCI_myo_4_26_52[i, paste0(col, '_week26')]
    }
  }
}

for (var in myotomes_EMSCI[1:20]){
  var_diff <- paste0('delta_', var, '_04_recovery')
  var04 <- paste0(var, '_week4')
  var_recovery <- paste0(var, '_locf')
  #wide_EMSCI_myo_4_26_52[var_diff] <- NA
  wide_EMSCI_myo_4_26_52[var_diff] <- as.numeric(wide_EMSCI_myo_4_26_52[[var_recovery]]) - as.numeric(wide_EMSCI_myo_4_26_52[[var04]])
}

summary(wide_EMSCI_myo_4_26_52)

# compute mean delta below NLI
motor_vars_delta_EMSCI <- myotomes_EMSCI[1:20]
for (var in myotomes_EMSCI[1:20]){
  motor_vars_delta_EMSCI[motor_vars_delta_EMSCI == var] <- paste0('delta_', var, '_04_recovery')
}

wide_EMSCI_myo_4_26_52$mean_delta_blwNLI <- NA
wide_EMSCI_myo_4_26_52 <- wide_EMSCI_myo_4_26_52[rowSums(is.na(wide_EMSCI_myo_4_26_52)) != ncol(wide_EMSCI_myo_4_26_52), ]
for (i in c(1:dim(wide_EMSCI_myo_4_26_52)[1])){
  ID = wide_EMSCI_myo_4_26_52$Patientennummer[i]
  NLI = raw_data_EMSCI %>% filter(Patientennummer == ID, ExamStage == 'acute I') %>%
    select(NLI)
  if (is.na(NLI[[1]])) {
    print('NLI missing')
  } else if (match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI))) < 6){
    blwNLI <- motor_vars_delta_EMSCI
  } else if (match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI))) %in% c(6:9)) {
    blwNLI <- motor_vars_delta_EMSCI[-c(1:((match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI)))-4)*2))]
  } else {
    blwNLI <- motor_vars_delta_EMSCI[-c(1:10)]
  }
  if (!(is.na(NLI[[1]]))) {
    temp <- wide_EMSCI_myo_4_26_52 %>% dplyr::filter(Patientennummer == ID) %>%
      dplyr::select(blwNLI)
    wide_EMSCI_myo_4_26_52$mean_delta_blwNLI[wide_EMSCI_myo_4_26_52$Patientennummer == ID] <- rowMeans(temp, na.rm=T)[[1]]
  }
}

################################################################################
# Identify outliers...
################################################################################

# Define potential population based on inclusion criteria
demo_psm_EMSCI <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute I', NLI_level %in% c('cervical', 'thoracic')) %>%
  dplyr::select(all_of(c('Patientennummer', 'Sex', 'AgeAtDOI', 'LEMS', 'UEMS', 'AIS', 'NLI_level'))) %>%
  dplyr::filter(Patientennummer %in% unique(true_severe_EMSCI$Patientennummer)) # only define outliers in the true severe defined in the clinical definition
demo_psm_EMSCI[demo_psm_EMSCI == ''] <- NA

demo_delta_psm_EMSCI <- merge(demo_psm_EMSCI, wide_EMSCI_myo_4_26_52 %>% select('Patientennummer', 'mean_delta_blwNLI'))
  
EMSCI_LEMS_26 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute III') %>%
  dplyr::select(all_of(c('Patientennummer', 'LEMS')))
names(EMSCI_LEMS_26)[2] <- paste0(names(EMSCI_LEMS_26)[2], '_week26')
EMSCI_LEMS_52 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'chronic') %>%
  dplyr::select(all_of(c('Patientennummer', 'LEMS')))
names(EMSCI_LEMS_52)[2] <- paste0(names(EMSCI_LEMS_52)[2], '_week52')

EMSCI_LEMS_recovery <- merge(EMSCI_LEMS_26, EMSCI_LEMS_52)
EMSCI_LEMS_recovery[EMSCI_LEMS_recovery == ''] <- NA
EMSCI_LEMS_recovery$LEMS_locf <- NA
for (i in c(1:dim(EMSCI_LEMS_recovery)[1])){
    if (!(is.na(EMSCI_LEMS_recovery[i, 'LEMS_week52']))){
      EMSCI_LEMS_recovery[i, 'LEMS_locf'] <- EMSCI_LEMS_recovery[i, 'LEMS_week52']
    } else {
      EMSCI_LEMS_recovery[i, 'LEMS_locf'] <- EMSCI_LEMS_recovery[i, 'LEMS_week26']
    }
}

demo_delta_locf_psm_EMSCI <- merge(demo_delta_psm_EMSCI, EMSCI_LEMS_recovery)
demo_delta_locf_psm_EMSCI$UEMS <- as.numeric(demo_delta_locf_psm_EMSCI$UEMS)
demo_delta_locf_psm_EMSCI$LEMS <- as.numeric(demo_delta_locf_psm_EMSCI$LEMS)
demo_delta_locf_psm_EMSCI$LEMS_locf <- as.numeric(demo_delta_locf_psm_EMSCI$LEMS_locf)
demo_delta_locf_psm_EMSCI <- demo_delta_locf_psm_EMSCI %>% tidyr::drop_na(LEMS_locf, mean_delta_blwNLI, UEMS, AgeAtDOI)
  #demo_delta_locf_psm_EMSCI[complete.cases(demo_delta_locf_psm_EMSCI), ]


# ...95% percentile in terms of mean delta LEMS week 04 to recovery (week 52 with LOCF when week 52 is missing) below NLI
id_top5perc_diff_blwnli_EMSCI <- as.vector(names(table(demo_delta_locf_psm_EMSCI$Patientennummer[demo_delta_locf_psm_EMSCI$mean_delta_blwNLI > quantile(demo_delta_locf_psm_EMSCI$mean_delta_blwNLI, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
demo_delta_locf_psm_EMSCI$tail95_blwnli <- ifelse(demo_delta_locf_psm_EMSCI$Patientennummer %in% id_top5perc_diff_blwnli_EMSCI, 1, 0)


################################################################################
# Comparison of baseline characteristics of the top 5% in terms of mean delta MS below NLI
# compared to the PSM matched comparators

# Perform propensity score matching to obtain comparable groups in terms of baseline characteristics
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}
 
  print('Characteristics before PSM:')
  print(table1(formula(paste("~ factor(Sex) + AgeAtDOI + as.numeric(LEMS) + as.numeric(UEMS) + factor(AIS) + NLI_level | factor(", criteria, ")")), 
               data=demo_delta_locf_psm_EMSCI, overall=F, extra.col=list(`P-value`=pvalue)))
  
  mod_test1 <- glm(formula(paste('factor(', criteria, ") ~ factor(Sex) + AgeAtDOI + as.numeric(LEMS) + as.numeric(UEMS) + NLI_level")), #+ factor(AIS)
                   data=demo_delta_locf_psm_EMSCI, family='binomial')
  print(summary(mod_test1))
  
  psa_nn <- matchit(formula(paste0('factor(', criteria, ') ~ as.numeric(UEMS)  + factor(Sex) + AgeAtDOI + as.numeric(LEMS) + NLI_level')), #factor(AIS) + 
                    data=demo_delta_locf_psm_EMSCI, distance='glm', method = 'nearest', m.order = 'largest', replace = F, 
                    ratio=4)
  psa_nn
  print(summary(psa_nn))
  
  love.plot(bal.tab(psa_nn),
            stat=c('m', 'v'),
            grid=T,
            thresholds = c(m=.25,v=1.25))
  
  psa_matched <- match.data(psa_nn)
  
  print('Characteristics after PSM:')
  print(table1(formula(paste("~ factor(Sex) + AgeAtDOI + LEMS + UEMS + factor(AIS) + LEMS_locf + NLI_level | factor(", criteria, ")")),
               data=psa_matched, overall=F, extra.col=list(`P-value`=pvalue)))
  
  EMSCI_UEMS_26 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute III') %>%
    dplyr::select(all_of(c('Patientennummer', 'UEMS')))
  names(EMSCI_UEMS_26)[2] <- paste0(names(EMSCI_UEMS_26)[2], '_week26')
  EMSCI_UEMS_52 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'chronic') %>%
    dplyr::select(all_of(c('Patientennummer', 'UEMS')))
  names(EMSCI_UEMS_52)[2] <- paste0(names(EMSCI_UEMS_52)[2], '_week52')
  
  EMSCI_UEMS_recovery <- merge(EMSCI_UEMS_26, EMSCI_UEMS_52)
  EMSCI_UEMS_recovery[EMSCI_UEMS_recovery == ''] <- NA
  EMSCI_UEMS_recovery$UEMS_locf <- NA
  for (i in c(1:dim(EMSCI_UEMS_recovery)[1])){
    if (!(is.na(EMSCI_UEMS_recovery[i, 'UEMS_week52']))){
      EMSCI_UEMS_recovery[i, 'UEMS_locf'] <- EMSCI_UEMS_recovery[i, 'UEMS_week52']
    } else {
      EMSCI_UEMS_recovery[i, 'UEMS_locf'] <- EMSCI_UEMS_recovery[i, 'UEMS_week26']
    }
  }
  
  psa_matched <- merge(psa_matched, EMSCI_UEMS_recovery)
  
  psa_matched <- mutate(psa_matched, sum = as.numeric(UEMS_locf) - as.numeric(LEMS_locf))
  ## Definition 1 CCS: UEMS - LEMS < 0
  psa_matched$diff_ms_recovery_bin <- ifelse(psa_matched$sum<0, 1, 0)
  
  ## Definition 2 CCS: UEMS - LEMS < 4
  psa_matched$diff_ms_52_bin5 <- ifelse(psa_matched$sum<(-4), 1, 0)
  
  ## Definition 3 CCS: UEMS - LEMS < 9
  psa_matched$diff_ms_52_bin10 <- ifelse(psa_matched$sum<(-9), 1, 0)
  
  ## Definition 4 CCS: UEMS - LEMS < 18
  psa_matched$diff_ms_52_bin19 <- ifelse(psa_matched$sum<(-18), 1, 0)

  
  ## Definition 5 CCS: based on NLI --> (1- (mean_UEMS_blw_NLI/mean_LEMS))*100 > 0.1
  
  wide_EMSCI_myo_locf <- wide_EMSCI_myo_4_26_52 %>% select('Patientennummer', paste0(myotomes_EMSCI[1:20], '_locf'))
  wide_EMSCI_myo_locf[, 2:21] <- sapply(wide_EMSCI_myo_locf[, 2:21], as.numeric)
  
  wide_EMSCI_myo_locf$mean_LEMS_w52 <- rowMeans(wide_EMSCI_myo_locf[ , c(paste0(myotomes_EMSCI[11:20], '_locf'))], 
                                                na.rm=TRUE)
  
  motor_vars_UEMS_52 <- paste0(myotomes_EMSCI[1:10], '_locf')
  wide_EMSCI_myo_locf$mean_UEMS_blw_NLI_w52 <- NA
  for (i in c(1:dim(wide_EMSCI_myo_locf)[1])){
    ID = wide_EMSCI_myo_locf$Patientennummer[i]
    NLI = raw_data_EMSCI %>% filter(Patientennummer == ID, ExamStage == 'acute I') %>%
      select(NLI)
    if (is.na(NLI[[1]])) {
      print('NLI missing')
    } else if (match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI))) < 6){
      blwNLI <- motor_vars_UEMS_52
    } else if (match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI))) %in% c(6:9)) {
      blwNLI <- motor_vars_UEMS_52[-c(1:((match(NLI[[1]], levels(factor(raw_data_EMSCI$NLI)))-4)*2))]
    } else {
      blwNLI <- motor_vars_UEMS_52[-c(1:10)]
    }
    if (!(is.na(NLI[[1]]))) {
      temp <- wide_EMSCI_myo_locf %>% dplyr::filter(Patientennummer == ID) %>%
        dplyr::select(blwNLI)
      wide_EMSCI_myo_locf$mean_UEMS_blw_NLI_w52[wide_EMSCI_myo_locf$Patientennummer == ID] <- rowMeans(temp, na.rm=T)[[1]]
    }
  }
  
  wide_EMSCI_myo_locf$diff_ms_52_NLI <- (1 - (wide_EMSCI_myo_locf$mean_UEMS_blw_NLI_w52/
                                                wide_EMSCI_myo_locf$mean_LEMS_w52))*100
  wide_EMSCI_myo_locf$diff_ms_52_bin_NLI <- ifelse(wide_EMSCI_myo_locf$diff_ms_52_NLI>0.1, 1, 0)
  
  wide_EMSCI_myo_locf %>% filter(Patientennummer %in% psa_matched$Patientennummer)
  
  wide_EMSCI_myo_locf_sub <- wide_EMSCI_myo_locf %>% select(Patientennummer, diff_ms_52_bin_NLI)
  
  psa_matched <- merge(psa_matched, wide_EMSCI_myo_locf_sub)
  
  print(table1(formula(paste('~ factor(diff_ms_recovery_bin) + factor(diff_ms_52_bin5) + 
           factor(diff_ms_52_bin10) + factor(diff_ms_52_bin19) + factor(diff_ms_52_bin_NLI) + as.numeric(UEMS_locf) | factor(', 'tail95_blwnli', ')')), 
               data=psa_matched, overall=F, extra.col=list(`P-value`=pvalue)))
  
  chisq.test((table(psa_matched$diff_ms_52_bin_NLI, psa_matched$tail95_blwnli)))
  
  
  psa_matched_true_severe <- psa_matched
  
# Compare the IDs between clinical definition and stats def in AIS A

unique(psa_matched_aisA$Patientennummer)
unique(id_recoverers_EMSCI$Patientennummer)
intersect(unique(psa_matched_true_severe$Patientennummer[psa_matched_true_severe$tail95_blwnli == 1]), unique(id_recoverers_EMSCI$Patientennummer) )

raw_data_EMSCI %>% filter(Patientennummer %in% unique(psa_matched_true_severe$Patientennummer[psa_matched_true_severe$tail95_blwnli == 1]), 
                          ExamStage == 'acute I') %>% 
  select(Patientennummer, NLI)
