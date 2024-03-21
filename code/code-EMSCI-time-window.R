################################################################################
# SCI - Phenomenal recovery (PR) - EMSCI description
# L. Bourguignon
# First version : 23.02.2023
# Last update : 08.11.2023
################################################################################

raw_data_EMSCI <- read.csv('/Volumes/green_groups_bds_public/Data/EMSCI/emsci_data_2020.csv')

#raw_data_EMSCI
# A_EMSCI_very_acute <- raw_data_EMSCI %>% filter(AIS == 'A', ExamStage == 'very acute')
# A_EMSCI_acute1 <- raw_data_EMSCI %>% filter(AIS == 'A', ExamStage == 'acute I')
# ID_trueA_EMSCI <- intersect(A_EMSCI_very_acute[['Patientennummer']], A_EMSCI_acute1[['Patientennummer']])
# true_A_EMSCI <- subset(raw_data_EMSCI, Patientennummer %in% ID_trueA_EMSCI)
# dim(true_A_EMSCI %>% filter(AIS == 'NT', ExamStage == 'chronic'))
# true_A_EMSCI_subset <- true_A_EMSCI %>% select(Patientennummer, AIS, ExamStage_weeks)
# true_A_EMSCI_subset$AIS_week <- paste0(true_A_EMSCI_subset$AIS, '_', true_A_EMSCI_subset$ExamStage_weeks)
# true_A_EMSCI_subset_wide <- reshape(true_A_EMSCI_subset, idvar = "Patientennummer", timevar = "ExamStage_weeks", direction = "wide")

define_df_window_EMSCI <- function(raw_data_EMSCI, window){
  # Print time window requested
  print('Time window:')
  print(window)
  
  # Define subset based on AIS grades only first
  if (length(window) == 1){
    A_EMSCI_time1 <- raw_data_EMSCI %>% filter(AIS == 'A', ExamStage_weeks == window[1])
    ID_trueA_EMSCI <- A_EMSCI_time1$Patientennummer
  } else if (length(window) == 2){
    A_EMSCI_time1 <- raw_data_EMSCI %>% dplyr::filter(AIS == 'A', ExamStage_weeks == as.numeric(window[1]))
    A_EMSCI_time2 <- raw_data_EMSCI %>% dplyr::filter(AIS == 'A', ExamStage_weeks == as.numeric(window[2]))
    ID_trueA_EMSCI <- intersect(A_EMSCI_time1[['Patientennummer']], A_EMSCI_time2[['Patientennummer']])
  } else {
    print('Time window should include 1 or 2 values only.')
  }
  
  true_A <- raw_data_EMSCI %>% dplyr::filter(Patientennummer %in% ID_trueA_EMSCI)
  true_A_time1 <- true_A %>% dplyr::filter(ExamStage_weeks == window[1])
  
  true_A[true_A == ''] <- NA
  true_A[true_A == 'NT'] <- NA
  
  print('Number of patients with AIS A in the time window requested')
  print(length(ID_trueA_EMSCI))
  
  # #Find the maximum value of LEMS recorded in the first weeks
  # values_lems <- paste0(rep('Lower0', length(window)), window)
  # subset_lems <- true_A[values_lems]
  # subset_lems <- subset_lems %>% rowwise() %>%
  #   mutate(max_early_lems = max(c_across(where(is.numeric)), na.rm = TRUE))
  # subset_lems[subset_lems == -Inf] <- NA
  # true_A$max_early_lems <- subset_lems$max_early_lems
  
  # #Find the maximum value of UEMS recorded in the first weeks
  # values_uems <- paste0(rep('Upper0', length(window)), window)
  # subset_uems <- true_A[values_uems]
  # subset_uems <- subset_uems %>% rowwise() %>%
  #   mutate(max_early_uems = max(c_across(where(is.numeric)), na.rm = TRUE))
  # subset_uems[subset_uems == -Inf] <- NA
  # true_A$max_early_uems <- subset_uems$max_early_uems
  
  levels <- c('L1', 'L2', 'L3', 'L4', 'L5', 'S1', 'S2', 'S3', 'S45')
  levels_myo <- c('L2', 'L3', 'L4', 'L5', 'S1')
  sides <- c('L', 'R')
  tests <- c('PP', 'LT')
  temp_test_ss <- paste0(rep(sides,2), rep(tests, each = 2))
  temp_test_myo <- paste0(rep(sides,1), rep('MS',2))
  col_sensory <- c(paste0(rep(temp_test_ss, each = length(levels)), '_',
    rep(levels, length(temp_test_ss))),
    paste0(rep(temp_test_myo, each = length(levels_myo)), '_',
           rep(levels_myo, length(temp_test_myo))))
  
  df_temp <- true_A_time1['Patientennummer']
  df_temp$early_sensation_left <- 0
  df_temp$bin_lems_improv <- NA
  df_temp$bin_uems_improv <- NA
  
  for (i in c(1:length(ID_trueA_EMSCI))){ # for all patients with AIS A at week 1 and 4
    #print(i)
    
    true_A_sub <- true_A %>% dplyr::filter(Patientennummer == ID_trueA_EMSCI[i])
    id = ID_trueA_EMSCI[i]
    
    true_A_sub_52 <- true_A_sub %>% dplyr::filter(ExamStage_weeks == 52)
    true_A_sub_26 <- true_A_sub %>% dplyr::filter(ExamStage_weeks == 26)
    
    if (length(window) == 1){
      true_A_sub_time1 <- true_A_sub %>% dplyr::filter(ExamStage_weeks == window[1])
      max_early_lems = true_A_sub_time1$LEMS
      max_early_uems = true_A_sub_time1$UEMS
    } else if (length(window) == 2){
      true_A_sub_time1 <- true_A_sub %>% dplyr::filter(ExamStage_weeks == window[1])
      true_A_sub_time2 <- true_A_sub %>% dplyr::filter(ExamStage_weeks == window[2])
      max_early_lems = max(true_A_sub_time1$LEMS, true_A_sub_time2$LEMS)
      max_early_uems = max(true_A_sub_time1$UEMS, true_A_sub_time2$UEMS)
    }

    # Define LEMS improvement
    
    if (is.na(true_A_sub_52$LEMS)){ # if LEMS at week 52 is missing
      value_lems = true_A_sub_26$LEMS # consider LEMS at week 26
    } else {
      value_lems = true_A_sub_52$LEMS
    }
    # calculate the difference in maximum value of LEMS recorded in the first weeks
    # and the LEMS at recovery as defined above
    if (is.na(value_lems) | is.na(max_early_lems)){ # if both week 52 and 26 LEMS are missing
      df_temp[df_temp$Patientennummer ==id, ]$bin_lems_improv <- NA # consider that the improvement info is missing for LEMS
    } else if (as.numeric(value_lems) - as.numeric(max_early_lems) > 5){ # if this difference is more than 5 points
      df_temp[df_temp$Patientennummer ==id, ]$bin_lems_improv <- 1 # consider that the patient improved in terms of LEMS
    } else {
      df_temp[df_temp$Patientennummer ==id, ]$bin_lems_improv <- 0 # otherwise consider that the patient did not improve in terms of LEMS
    }
    
    # Define UEMS improvement
    
    if (is.na(true_A_sub_52$UEMS)){ # if UEMS at week 52 is missing
      value_uems = true_A_sub_26$UEMS # consider UEMS at week 26
    } else {
      value_uems = true_A_sub_52$UEMS # otherwise consider LEMS at week 52
    }
    # calculate the difference in maximum value of UEMS recorded in the first weeks
    # and the UEMS at recovery as defined above
    if (is.na(value_uems) | is.na(max_early_uems)){ # if both week 52 and 26 UEMS are missing
      df_temp[df_temp$Patientennummer ==id, ]$bin_uems_improv <- NA # consider that the improvement info is missing for UEMS
    } else if (as.numeric(value_uems) - as.numeric(max_early_uems) > 5){ # if this difference is more than 5 points
      df_temp[df_temp$Patientennummer ==id, ]$bin_uems_improv <- 1 # consider that the patient improved in terms of UEMS
    } else {
      df_temp[df_temp$Patientennummer ==id, ]$bin_uems_improv <- 0 # otherwise consider that the patient did not improve in terms of UEMS
    }
    
    # Define if sensory function was present in the lower limbs at week 1 or 4
    
    for (col in col_sensory){
      if (is.na(true_A_sub_time1[col]) | is.na(true_A_sub_time2[col])){
        df_temp[df_temp$Patientennummer ==id,]$early_sensation_left <- NA
      } else if (!(true_A_sub_time1[col] == 0) | !(true_A_sub_time2[col] == 0)){
        df_temp[df_temp$Patientennummer ==id,]$early_sensation_left <- 1
      }
    }
  }
  return (df_temp)
}

# Define data frame of AIS A patients within the specified time window
# This dataframe contains information about UEMS and LEMS improvement (>5)
# (columns bin_uems_improv and bin_ems_improv, respectively)
# 0 = no improvement; 1 = improvement >5 points
# and whether sensory or motor function was present within the specified time
# column early_sensation_left
# 0 = no sensation/motor function left; 1 = at least one of sensory or motor function left

EMSCI_recovery <- define_df_window_EMSCI(raw_data_EMSCI, c('2','4'))

# Define AIS grade at recovery:
# based on week 52 if available, week 26 otherwise (last observation carried forward)
temp52 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'chronic',
                          Patientennummer %in% EMSCI_recovery$Patientennummer) %>%
  dplyr::select(Patientennummer, AIS, LEMS, UEMS)

temp26 <- raw_data_EMSCI %>% dplyr::filter(ExamStage == 'acute III',
                                    Patientennummer %in% EMSCI_recovery$Patientennummer) %>%
  dplyr::select(Patientennummer, AIS, LEMS, UEMS)

EMSCI_recovery <- merge(EMSCI_recovery, temp52, by = "Patientennummer")
names(EMSCI_recovery)[names(EMSCI_recovery) == "AIS"] <- "AIS.52"
names(EMSCI_recovery)[names(EMSCI_recovery) == "LEMS"] <- "LEMS.52"
names(EMSCI_recovery)[names(EMSCI_recovery) == "UEMS"] <- "UEMS.52"
EMSCI_recovery <- merge(EMSCI_recovery, temp26, by = "Patientennummer")
names(EMSCI_recovery)[names(EMSCI_recovery) == "AIS"] <- "AIS.26"
names(EMSCI_recovery)[names(EMSCI_recovery) == "LEMS"] <- "LEMS.26"
names(EMSCI_recovery)[names(EMSCI_recovery) == "UEMS"] <- "UEMS.26"
EMSCI_recovery[EMSCI_recovery == ''] <- NA
EMSCI_recovery[EMSCI_recovery == 'NT'] <- NA

for (i in c(1:dim(EMSCI_recovery)[1])){
  if (!(is.na(EMSCI_recovery$AIS.52[i]))){
    EMSCI_recovery$AIS.recovery[i] <- EMSCI_recovery$AIS.52[i]
  } else {
    EMSCI_recovery$AIS.recovery[i] <- EMSCI_recovery$AIS.26[i]
  }
  
  if (!(is.na(EMSCI_recovery$LEMS.52[i]))){
    EMSCI_recovery$LEMS.recovery[i] <- EMSCI_recovery$LEMS.52[i]
  } else {
    EMSCI_recovery$LEMS.recovery[i] <- EMSCI_recovery$LEMS.26[i]
  }
  
  if (!(is.na(EMSCI_recovery$UEMS.52[i]))){
    EMSCI_recovery$UEMS.recovery[i] <- EMSCI_recovery$UEMS.52[i]
  } else {
    EMSCI_recovery$UEMS.recovery[i] <- EMSCI_recovery$UEMS.26[i]
  }
}

# Among AIS A patients, only the ones qualified as "true severe"
# (no sensory or motor function left below L1)
# qualify for our analysis and PR definition

# Subset only the true severes
true_severe_EMSCI <- EMSCI_recovery %>% dplyr::filter(early_sensation_left == 0)
# # Subset PR patients -- old definition
# true_severe_EMSCI_recover <- true_severe_EMSCI %>%
#   dplyr::filter((AIS.recovery %in% c('C', 'D') | bin_lems_improv == 1) & 
#            (bin_uems_improv == 1 | bin_lems_improv == 1))
# Subset PR patients -- new definition
true_severe_EMSCI_recover <- true_severe_EMSCI %>%
  dplyr::filter(bin_lems_improv == 1)
# Extract IDs from PR
id_recoverers_EMSCI <- true_severe_EMSCI_recover$Patientennummer
# Add PR status label in the dataframe
true_severe_EMSCI$recover <- 'No recovery'
true_severe_EMSCI$recover[true_severe_EMSCI$Patientennummer %in% id_recoverers_EMSCI] <- "Recovery"

id_norecovery_EMSCI <- true_severe_EMSCI$Patientennummer[true_severe_EMSCI$recover == 'No recovery']

true_severe_EMSCI_complete <- true_severe_EMSCI[complete.cases(true_severe_EMSCI[c('LEMS.recovery', "recover", "bin_lems_improv")]),]
id_recoverers_EMSCI <- true_severe_EMSCI_complete %>% filter(recover == "Recovery") %>% select(Patientennummer)
id_norecovery_EMSCI <- true_severe_EMSCI_complete %>% filter(recover == "No recovery") %>% select(Patientennummer)

recover_acute1 <- raw_data_EMSCI %>% filter(Patientennummer %in% id_recoverers_EMSCI[[1]], ExamStage == 'acute I')
norecover_acute1 <- raw_data_EMSCI %>% filter(Patientennummer %in% id_norecovery_EMSCI[[1]], ExamStage == 'acute I')

table1(~ Sex + AgeAtDOI + AIS + NLI_level + LEMS + as.numeric(UEMS), data=norecover_acute1)
table1(~ Sex + AgeAtDOI + AIS + NLI_level + LEMS + as.numeric(UEMS), data=recover_acute1)
table1(~ as.numeric(LEMS.recovery) + as.numeric(UEMS.recovery) | recover, data=true_severe_EMSCI_complete)


# table(true_severe_EMSCI$recover, true_severe_EMSCI$AIS.recovery, useNA='always')
# 
# 
# temp_YOI_recover_EMSCI <- raw_data_EMSCI %>% 
#   filter(Patientennummer %in% id_recoverers_EMSCI,
#          ExamStage == 'very acute') %>%
#   select(Patientennummer, YEARDOI)
# temp_YOI_norecover_EMSCI <- raw_data_EMSCI %>% 
#   filter(!Patientennummer %in% id_recoverers_EMSCI,
#          ExamStage == 'very acute') %>%
#   select(Patientennummer, YEARDOI)
# 
# p_YOI_recover_EMSCI <- ggplot(temp_YOI_recover_EMSCI, aes(x=YEARDOI)) + 
#   geom_density()
# p_YOI_recover_EMSCI
# 
# p_YOI_norecover_EMSCI <- ggplot(temp_YOI_norecover_EMSCI, aes(x=YEARDOI)) + 
#   geom_density()
# p_YOI_norecover_EMSCI
# 
# temp_baseline <- raw_data_EMSCI %>% 
#   filter(Patientennummer %in% id_recoverers_EMSCI,
#          ExamStage == 'very acute') %>%
#   select(Patientennummer, NLI_level, AIS, AgeAtDOI, Sex)
# table(temp_baseline$NLI_level)
# table(temp_baseline$AIS)
# table(temp_baseline$Sex)
# mean(temp_baseline$AgeAtDOI, na.rm = T)
# sd(temp_baseline$AgeAtDOI, na.rm = T)
# 
# temp_EMSCI_level_recover <- true_severe_EMSCI %>% 
#   filter(Patientennummer %in% id_recoverers_EMSCI) %>% 
#   select(AIS.recovery, bin_lems_improv)
# 
# 
# temp_baseline_all <- raw_data_EMSCI %>% 
#   filter(ExamStage == 'very acute')
# 
# # Examine patients with LEMS > 0 at week 4 but currently classified as PR
# raw_data_EMSCI %>% filter(Patientennummer %in% 
#                             c('90008', '60202', '340010', '40030', '90026', '20057'), 
#                           ExamStage_weeks == 4) %>% 
#   select(Patientennummer, ExamStage_weeks, AIS, VAC, DAP, LEMS, TPP, TLT)

