
pca_outlier <- read.csv('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/outliers_week4_full.csv')

sygen_raw_pca <- merge(raw_sygen, pca_outlier, by = "ptid")

print(table1(formula(paste("~ factor(sexcd) + age + lower04 + upper04 + factor(ais4) + NLI_calc_bin | factor(", 'svm_outlier', ")")),
             data=sygen_raw_pca, overall=F, extra.col=list(`P-value`=pvalue)))

print(table1(formula(paste('~ lower_locf + upper_locf + factor(ais52) + 
                       factor(anyana04) + factor(vaccd04) + 
                       factor(anyana52) + factor(vaccd52) + factor(tx1_r) | factor(', 'svm_outlier', ')')), 
             data=sygen_raw_pca, overall=F, extra.col=list(`P-value`= pvalue)))

print(table1(formula(paste('~ factor(diff_ms_52_bin) + factor(diff_ms_52_bin5) + 
           factor(diff_ms_52_bin10) + factor(diff_ms_52_bin19) +
           factor(diff_ms_52_bin_NLI) | factor(', 'svm_outlier', ')')), 
             data=sygen_raw_pca, overall=F, extra.col=list(`P-value`=pvalue)))



ID_tail <- raw_sygen %>% filter(tail95_blwnli == 1) %>% select(ptid)

c('129202-A', '129806-AB', '19206-ABCDEFGHI', '21302-ABCD', '219101-ABCD',
'23601-AB', '249902-AB', '261001-ABCDEFG', '269711-ABCDEFGH',
'29701-ABCDE', '4304-ABCDEFGH', '6105-ABCDE', '79201-ACD', '9203-ABCD',
'9602-ABCDEFGHI') %in% ID_tail['ptid'][[1]]

'9203-ABCD' %in% ID_tail['ptid'][[1]]

print(table1(formula(paste("~ factor(SEXCD) + AGE + Lower04 + Upper04 + factor(AIS.4) | factor(", 'bin_lems_improv_5', ")")),
             data=df_window.1.4[!is.na(df_window.1.4$bin_lems_improv_5),], overall=F, extra.col=list(`P-value`=pvalue)))

temp_pr_clinical <- raw_sygen %>% filter(ptid %in% df_window.1.4$new_PTID) %>% select(upper_locf, lower_locf, ptid)
temp_pr_clinical_window <- merge(df_window.1.4, temp_pr_clinical, by.x=c('new_PTID'), by.y=c('ptid'))

temp_pr_clinical_window_NLI <- merge(temp_pr_clinical_window, raw_sygen %>% select(NLI_calc_bin, ptid), by.x=c('new_PTID'), by.y=c('ptid'))

print(table1(formula(paste("~ factor(SEXCD) + AGE + Lower04 + Upper04 + factor(AIS.4) + upper_locf + lower_locf + NLI_calc_bin | factor(", 'bin_lems_improv_5', ")")),
             data=temp_pr_clinical_window_NLI[!is.na(temp_pr_clinical_window_NLI$bin_lems_improv_5),], overall=F, extra.col=list(`P-value`=pvalue)))


## Comparing CCS in PR versus comparator in Sygen clinical definition

ID_recover_sygen <- df_window.1.4 %>% filter(bin_lems_improv_5 ==1) %>% select(new_PTID)
ID_comparator_sygen <- df_window.1.4 %>% filter(bin_lems_improv_5 ==0) %>% select(new_PTID)

raw_sygen$PR_clinical <- NA
for (i in c(1:length(raw_sygen$ptid))){
  #print(i)
  id = raw_sygen$ptid[i]
  #print(id)
  if (id %in% ID_recover_sygen[[1]]){
    raw_sygen$PR_clinical[i] <- 'pr'
  } else if (id %in% ID_comparator_sygen[[1]]){
    raw_sygen$PR_clinical[i] <- 'comparator'
  }
}

print(table1(formula(paste('~ factor(diff_ms_52_bin) + factor(diff_ms_52_bin5) + 
           factor(diff_ms_52_bin10) + factor(diff_ms_52_bin19) +
           factor(diff_ms_52_bin_NLI) | factor(', 'PR_clinical', ')')), 
             data=raw_sygen[!is.na(raw_sygen$PR_clinical),], overall=F, extra.col=list(`P-value`=pvalue)))

#################################################################################
## Comparing CCS in PR versus comparator in EMSCI clinical definition

true_severe_EMSCI_complete <- mutate(true_severe_EMSCI_complete, sum = as.numeric(UEMS.recovery) - as.numeric(LEMS.recovery))
## Definition 1 CCS: UEMS - LEMS < 0
true_severe_EMSCI_complete$diff_ms_recovery_bin <- ifelse(true_severe_EMSCI_complete$sum<0, 1, 0)

## Definition 2 CCS: UEMS - LEMS < 4
true_severe_EMSCI_complete$diff_ms_52_bin5 <- ifelse(true_severe_EMSCI_complete$sum<(-4), 1, 0)

## Definition 3 CCS: UEMS - LEMS < 9
true_severe_EMSCI_complete$diff_ms_52_bin10 <- ifelse(true_severe_EMSCI_complete$sum<(-9), 1, 0)

## Definition 4 CCS: UEMS - LEMS < 18
true_severe_EMSCI_complete$diff_ms_52_bin19 <- ifelse(true_severe_EMSCI_complete$sum<(-18), 1, 0)

print(table1(formula(paste('~ factor(diff_ms_recovery_bin) + factor(diff_ms_52_bin5) + 
           factor(diff_ms_52_bin10) + factor(diff_ms_52_bin19) | factor(', 'recover', ')')), 
             data=true_severe_EMSCI_complete[!is.na(true_severe_EMSCI_complete$recover),], overall=F, extra.col=list(`P-value`=pvalue)))

## Definition 5 CCS: based on NLI --> (1- (mean_UEMS_blw_NLI/mean_LEMS))*100 > 0.1

raw_data_EMSCI_myo_rec <- raw_data_EMSCI %>% 
  filter(Patientennummer %in% true_severe_EMSCI_complete$Patientennummer,
         ExamStage %in% c('chronic') ) %>% 
  select(c(paste0(rep('RMS', length(levels_myotomes)), '_',
                   levels_myotomes),
           paste0(rep('LMS', length(levels_myotomes)), '_',
                  levels_myotomes),
           'Patientennummer'))

raw_data_EMSCI_myo_rec[raw_data_EMSCI_myo_rec == ''] <- NA
ID_NA_52 <- raw_data_EMSCI_myo_rec[is.na(raw_data_EMSCI_myo_rec$RMS_C5),]['Patientennummer'][[1]]
raw_data_EMSCI_myo_rec.complete <- raw_data_EMSCI_myo_rec[complete.cases(raw_data_EMSCI_myo_rec), ]

raw_data_EMSCI_myo_26 <- raw_data_EMSCI %>% 
  filter(Patientennummer %in% ID_NA_52,
         ExamStage %in% c('acute III') ) %>% 
  select(c(paste0(rep('RMS', length(levels_myotomes)), '_',
                  levels_myotomes),
           paste0(rep('LMS', length(levels_myotomes)), '_',
                  levels_myotomes),
           'Patientennummer'))

raw_data_EMSCI_myo_rec.complete.all <- rbind(raw_data_EMSCI_myo_rec.complete, raw_data_EMSCI_myo_26)
raw_data_EMSCI_myo_rec.complete.all[, 1:20] <- sapply(raw_data_EMSCI_myo_rec.complete.all[, 1:20], as.numeric)

raw_data_EMSCI_myo_rec.complete.all$mean_LEMS_w52 <- rowMeans(raw_data_EMSCI_myo_rec.complete.all[ , c('RMS_L2','RMS_L3',
                                                                'RMS_L4','RMS_L5',
                                                                'LMS_S1','LMS_L2','LMS_L3',
                                                                'LMS_L4','LMS_L5',
                                                                'LMS_S1')], 
                                    na.rm=TRUE)

motor_vars_UEMS_52 <- c('RMS_C5','RMS_C6','RMS_C7','RMS_C8', 'RMS_T1',
                        'LMS_C5','LMS_C6','LMS_C7','LMS_C8', 'LMS_T1')
raw_data_EMSCI_myo_rec.complete.all$mean_UEMS_blw_NLI_w52 <- NA
for (i in c(1:dim(raw_data_EMSCI_myo_rec.complete.all)[1])){
  ID = raw_data_EMSCI_myo_rec.complete.all$Patientennummer[i]
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
    temp <- raw_data_EMSCI_myo_rec.complete.all %>% dplyr::filter(Patientennummer == ID) %>%
      dplyr::select(blwNLI)
    raw_data_EMSCI_myo_rec.complete.all$mean_UEMS_blw_NLI_w52[raw_data_EMSCI_myo_rec.complete.all$Patientennummer == ID] <- rowMeans(temp, na.rm=T)[[1]]
  }
}

raw_data_EMSCI_myo_rec.complete.all$diff_ms_52_NLI <- (1 - (raw_data_EMSCI_myo_rec.complete.all$mean_UEMS_blw_NLI_w52/
                                                              raw_data_EMSCI_myo_rec.complete.all$mean_LEMS_w52))*100
raw_data_EMSCI_myo_rec.complete.all$diff_ms_52_bin_NLI <- ifelse(raw_data_EMSCI_myo_rec.complete.all$diff_ms_52_NLI>0.1, 1, 0)
table(factor(raw_data_EMSCI_myo_rec.complete.all$diff_ms_52_bin_NLI))

raw_data_EMSCI_myo_rec.complete.all %>% filter(diff_ms_52_bin_NLI ==1)

