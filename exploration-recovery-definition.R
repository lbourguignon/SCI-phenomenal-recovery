################################################################################
# SCI - Outstanding recovery - Changing time window for recovery definition
# L. Bourguignon
# First version : 21.03.2023
# Last update : 31.03.2023
################################################################################
## Load libraries

library(rlang)
library(table1)
library(yardstick)
library(plyr)
library(epitools)
library(lubridate)
library(dplyr)
library(stringr)
library(ggplot2)
library(SASdates)

################################################################################

max_early_value_myo <- function(col, window, df_window){
  values <- paste0(rep(col, length(window)), window)
  subset <- df_window[values]
  if (!(col %in% c('Lower', 'Upper'))){
    subset[subset == 9] <- NA
  }
  subset <- subset %>% rowwise() %>%
    dplyr::mutate(max_early_value = max(c_across(where(is.numeric)), na.rm = TRUE))
  subset[subset == -Inf] <- NA
  return(subset$max_early_value)
}

motor_improvement <- function(df_window, col, i, threshold){
  col52 = paste0(col, '52')
  col26 = paste0(col, '26')
  max_col = paste0('max_early_', col)
  if (is.na(df_window[col52][[1]][i])){ # if LEMS at week 52 is missing
    value_lems = df_window[col26][[1]][i] # consider LEMS at week 26
  } else {
    value_lems = df_window[col52][[1]][i] # otherwise consider LEMS at week 52
  }
  # calculate the difference in maximum value of LEMS recorded in the first weeks
  # and the LEMS at recovery as defined above
  if (is.na(value_lems) | is.na(df_window[max_col][[1]][i])){ # if both week 52 and 26 LEMS are missing
    out <- NA # consider that the improvement info is missing for LEMS
  } else if (value_lems - df_window[max_col][[1]][i] > threshold){ # if this difference is more than 5 points
    out <- 1 # consider that the patient improved in terms of LEMS
  } else {
    out <- 0 # otherwise consider that the patient did not improve in terms of LEMS
  }
  
  return(out)
}

true_severe <- function(raw_sygen_with_surgery, window){
  
  # Print time window requested
  print('Time window:')
  print(window)
  
  # Define subset based on AIS grades only first
  if (length(window) == 1){
    if ("0" %in% window[1]){
      true_A <- raw_sygen_with_surgery %>%
        filter(SEVCD == 1)
    } else if (!("0" %in% window)) {
      true_A <- raw_sygen_with_surgery %>%
        filter(!!as.symbol(paste0('AIS.', window[1])) == 'AIS A')
    } else {
      print('If 0 is present, it should be in the first position.')
    }
  } else if (length(window) == 2){
    if ("0" %in% window[1]){
      true_A <- raw_sygen_with_surgery %>%
        filter(SEVCD == 1,
               !!as.symbol(paste0('AIS.', window[2])) == 'AIS A')
    } else if (!("0" %in% window)) {
      true_A <- raw_sygen_with_surgery %>%
        filter(!!as.symbol(paste0('AIS.', window[1])) == 'AIS A',
               !!as.symbol(paste0('AIS.', window[2])) == 'AIS A')
    } else {
      print('If 0 is present, it should be in the first position.')
    }
  } else if (length(window) == 3){
    true_A <- raw_sygen_with_surgery %>%
      filter(SEVCD == 1,
             !!as.symbol(paste0('AIS.', window[2])) == 'AIS A',
             !!as.symbol(paste0('AIS.', window[3])) == 'AIS A')
  } else {
    print('Time window should include 1, 2 or 3 values only.')
  }
  
  print('Number of patients with AIS A in the time window requested')
  print(dim(true_A)[1])
  
  levels <- c('L1', 'L2', 'L3', 'L4', 'L5', 'S1', 'S2', 'S3', 'S45')
  sides <- c('L', 'R')
  tests <- c('PP', 'LT')
  earlyweeks <- c(paste0(rep('0', length(window)), window))
  earlyweeks <- earlyweeks[!earlyweeks %in% c('00')]
  true_A$early_sensation_left <- 0
  
  for (i in c(1:dim(true_A)[1])){ # for all patients with AIS A at week 1 and 4
    
    # Define if sensory function was present in the lower limbs at week 1 or 4
    
    if (length(earlyweeks) > 0){
      for (test in tests){
        for (lvl in levels){
          for (side in sides){
            for (week in earlyweeks){
              column = paste0(lvl, test, side, week)
              if (is.na(true_A[i, column])){
                print(paste('NA', column, 'for patient', true_A$new_PTID[i]))
              } else if (true_A[i, column] != 0) {
                true_A$early_sensation_left[i] <- 1
              }
            }
          }
        }
      }
    }
  }
  
  # Define a subset based on a stricter definition of a severe definition
  # Including not having any sensory function below L1
  
  if (length(earlyweeks) > 0){
    stricter_true_A <- true_A %>%
      filter(early_sensation_left == 0)
  } else {
    stricter_true_A <- true_A
  }
  
  stricter_true_A <- stricter_true_A %>% 
    mutate(level = if_else(str_detect(stricter_true_A$SPLVL, "T"), "T", "C"))
  
  print('Number of patients with strict severe injury in the time window requested')
  print(dim(stricter_true_A)[1])
  
  return(stricter_true_A)
  
}

make_table_1 <- function(df, window, col_recover){
  df_copy <- df
  
  df_copy$SEXCD <- factor(df_copy$SEXCD, levels=c(1,2), labels=c("Female", "Male"))
  df_copy$bin_uems_improv <- factor(df_copy$bin_uems_improv, levels=c(0, 1), labels=c("No", "Yes"))
  df_copy$bin_lems_improv <- factor(df_copy$bin_lems_improv, levels=c(0, 1), labels=c("No", "Yes"))
  df_copy$level <- factor(df_copy$level, levels=c('C', 'T'), labels=c("cervical", "thoracic"))
  df_copy$TX1_R <- factor(df_copy$TX1_R, levels=c('D1', 'D2', 'P'), labels=c("First dosage scheme", "Second dosage scheme", "Placebo"))
  df_copy$ANYANA52 <- factor(df_copy$ANYANA52, levels=c('0', '1'), labels=c("Absent", "Present"))
  df_copy$VACCD52 <- factor(df_copy$VACCD52, levels=c('0', '1'), labels=c("Absent", "Present"))
  df_copy$DECOMC01 <- factor(df_copy$DECOMC01, levels=c(0, 1), labels=c("No", "Yes"))
  df_copy$DECOCD01 <- factor(df_copy$DECOCD01, levels=c(0, 1), labels=c("No", "Yes"))
  df_copy$time_decompression_surgery <- as.numeric(df_copy$time_decompression_surgery)
  
  label(df_copy$SEXCD) <- "Sex"
  label(df_copy$AGE) <- "Age"
  label(df_copy$AIS.52) <- "AIS grade at week 52"
  label(df_copy$level) <- "Level of injury"
  label(df_copy$bin_lems_improv) <- "Improved >5pts in LEMS?ᵃ"
  label(df_copy$bin_uems_improv) <- "Improved >5pts in UEMS?ᵃ"
  label(df_copy$TX1_R) <- "Treatment group"
  label(df_copy$ANYANA52) <- "DAP at week 52"
  label(df_copy$VACCD52) <- "VAC at week 52"
  label(df_copy$DECOMC01) <- "Was decompression attempted?"
  label(df_copy$DECOCD01) <- "Was decompression obtained?"
  label(df_copy$time_decompression_surgery) <- "Time to first decompression surgery"
  
  units(df_copy$AGE) <- "years"
  units(df_copy$time_decompression_surgery) <- "minutes"
  
  caption <- paste('Time window: ', paste(window, collapse=' '))
  
  formula <- paste0('~ TX1_R + SEXCD + AGE + level + AIS.52 + ANYANA52 + VACCD52 + 
                  bin_uems_improv + bin_lems_improv + DECOMC01 + DECOCD01 + time_decompression_surgery | ', col_recover)
  
  return(table1(as.formula(formula), 
                caption = caption, data = df_copy))
  
}

df_window.1.4 <- true_severe(raw_sygen_with_surgery, c('1', '4'))

################################################################################

# How many patients had missing AIS grade at week 1 but AIS A at week 4 and 8?
test_NA_AIS1 <- raw_sygen_with_surgery %>%
  filter(AIS.1 == ' ', AIS.4 == 'AIS A', AIS.8 == 'AIS A')

dim(test_NA_AIS1)

test_NA_AIS1 %>% 
  dplyr::select(new_PTID, 
                AIS.1, AIS.4, AIS.8, AIS.16, AIS.26, AIS.52, 
                ANYANA01, ANYANA04, ANYANA08, ANYANA16, ANYANA26, ANYANA52,
                VACCD01, VACCD04, VACCD08, VACCD16, VACCD26, VACCD52, 
                SPLVL, Lower01, Lower52, Upper01, Upper52)

################################################################################

##### Investigate the recovery definition

df_window.1.4 <- true_severe(raw_sygen_with_surgery, c('1', '4'))

df_window.1.4$max_early_Lower <- max_early_value_myo('Lower0', c('1', '4'), df_window.1.4)
df_window.1.4$max_early_Upper <- max_early_value_myo('Upper0', c('1', '4'), df_window.1.4)

for (i in c(1:dim(df_window.1.4)[1])){
  # Define LEMS improvement
  df_window.1.4$bin_lems_improv_5[i] <- motor_improvement(df_window.1.4, 'Lower', i, 5)
  df_window.1.4$bin_lems_improv_2[i] <- motor_improvement(df_window.1.4, 'Lower', i, 2)
  
  # Define LEMS improvement
  df_window.1.4$bin_uems_improv_5[i] <- motor_improvement(df_window.1.4, 'Upper', i, 5)
  df_window.1.4$bin_uems_improv_2[i] <- motor_improvement(df_window.1.4, 'Upper', i, 2)
}


# Def1: AIS conversion + motor improvement OR LEMS improvement by more than 5 points

true_A_recover_new <- df_window.1.4 %>%
  filter((AIS.52 %in% c('AIS C', 'AIS D') | bin_lems_improv_5 == 1) &
           (bin_uems_improv_5 == 1 | bin_lems_improv_5 == 1))

id_recoverers <- true_A_recover_new$new_PTID
df_window.1.4$rec_AIS_motor_imp <- 'No recovery'
df_window.1.4$rec_AIS_motor_imp[df_window.1.4$new_PTID %in% id_recoverers] <- "Recovery"

# Def2: AIS conversion + motor improvement OR LEMS improvement by more than 2 points

true_A_recover_new <- df_window.1.4 %>%
  filter((AIS.52 %in% c('AIS C', 'AIS D') | bin_lems_improv_2 == 1) &
           (bin_uems_improv_2 == 1 | bin_lems_improv_2 == 1))

id_recoverers <- true_A_recover_new$new_PTID
df_window.1.4$rec_AIS_motor_imp2 <- 'No recovery'
df_window.1.4$rec_AIS_motor_imp2[df_window.1.4$new_PTID %in% id_recoverers] <- "Recovery"

# Def3: LEMS improvement alone by more than 5 points

df_window.1.4 <- df_window.1.4 %>% 
  dplyr::mutate(rec_LEMS_imp5 =
                  case_when(bin_lems_improv_5 == 1 ~ "Recovery", 
                     bin_lems_improv_5 == 0 ~ "No recovery"))

# Def4: LEMS improvement alone by more than 2 points

df_window.1.4 <- df_window.1.4 %>% 
  dplyr::mutate(rec_LEMS_imp2 =
                  case_when(bin_lems_improv_2 == 1 ~ "Recovery", 
                            bin_lems_improv_2 == 0 ~ "No recovery"))

# Def5: LEMS improvement in at least 2 myotomes

for (col in c('HIPFLL', 'KNEEXL', 'ANKDOL', 'GRETOL', 'ANKPLL', 
              'HIPFLR', 'KNEETR', 'ANKDOR', 'GRETOR', 'ANKPLR')){
  name_col = paste0('max_early_', col)
  df_window.1.4[name_col] <- max_early_value_myo(paste0(col,'0'), c('1', '4'), df_window.1.4)
}


for (i in c(1:dim(df_window.1.4)[1])){
  for (col in c('HIPFLL', 'KNEEXL', 'ANKDOL', 'GRETOL', 'ANKPLL', 
                'HIPFLR', 'KNEETR', 'ANKDOR', 'GRETOR', 'ANKPLR')){
    for (threshold in c(1, 2)){
      name_col = paste0('bin_', col, '_improv_', threshold)
      df_window.1.4[name_col][[1]][i] <- motor_improvement(df_window.1.4, col, i, threshold)
    }
  }
}

df_window.1.4 <- df_window.1.4 %>% 
  mutate(myo_recover_1 = rowSums(.[grep("_improv_1", names(.))], na.rm = T))

# df_window.1.4 %>% select('bin_HIPFLL_improv_1', 'bin_KNEEXL_improv_1', 'bin_ANKDOL_improv_1',
#                          'bin_GRETOL_improv_1', 'bin_ANKPLL_improv_1', 'bin_HIPFLR_improv_1',
#                          'bin_KNEETR_improv_1', 'bin_ANKDOR_improv_1', 'bin_GRETOR_improv_1',
#                          'bin_ANKPLR_improv_1')

df_window.1.4 <- df_window.1.4 %>% 
  dplyr::mutate(rec_myo_imp2 =
                  case_when(myo_recover_1 > 1 ~ "Recovery", 
                            myo_recover_1 < 2 ~ "No recovery"))

# Def6: LEMS improvement in at least 5 myotomes

df_window.1.4 <- df_window.1.4 %>% 
  dplyr::mutate(rec_myo_imp5 =
                  case_when(myo_recover_1 > 4 ~ "Recovery", 
                            myo_recover_1 < 5 ~ "No recovery"))

write.csv(df_window.1.4, '/Volumes/INANL/45_Lucie_Bourguignon/1_Data/1_Sygen/df_window.1.4.csv', row.names=FALSE)

# Examine characteristics by groups

table_rec_AIS_motor_imp <- make_table_1(df_window.1.4, c('1', '4'), 'rec_AIS_motor_imp') # initial definition
table_rec_AIS_motor_imp2 <- make_table_1(df_window.1.4, c('1', '4'), 'rec_AIS_motor_imp2') # initial definition but improvement reduced to 2 points instead of 5

# consider missing LEMS improvement as no recovery (n=16)
df_window.1.4$rec_LEMS_imp5[is.na(df_window.1.4$rec_LEMS_imp5)] <- 'Missing'
table_rec_LEMS_imp5 <- make_table_1(df_window.1.4, c('1', '4'), 'rec_LEMS_imp5') # LEMS with at least 5 points improvement

# consider missing LEMS improvement as no recovery (n=16)
df_window.1.4$rec_LEMS_imp2[is.na(df_window.1.4$rec_LEMS_imp2)] <- 'Missing'
table_rec_LEMS_imp2 <- make_table_1(df_window.1.4, c('1', '4'), 'rec_LEMS_imp2') # LEMS with at least 2 points improvement

table_rec_myo_imp2 <- make_table_1(df_window.1.4, c('1', '4'), 'rec_myo_imp2') # at least 2 myotomes with 1 point improvement
table_rec_myo_imp5 <-make_table_1(df_window.1.4, c('1', '4'), 'rec_myo_imp5') # at least 5 myotomes with 1 point improvement

# Look at overall of patients with the different definitions

ID_rec_AIS_motor_imp <- df_window.1.4$new_PTID[df_window.1.4$rec_AIS_motor_imp == "Recovery"]
ID_rec_AIS_motor_imp2 <- df_window.1.4$new_PTID[df_window.1.4$rec_AIS_motor_imp2 == "Recovery"]
ID_rec_LEMS_imp5 <- df_window.1.4$new_PTID[df_window.1.4$rec_LEMS_imp5 == "Recovery"]
ID_rec_LEMS_imp2 <- df_window.1.4$new_PTID[df_window.1.4$rec_LEMS_imp2 == "Recovery"]
ID_rec_myo_imp2 <- df_window.1.4$new_PTID[df_window.1.4$rec_myo_imp2 == "Recovery"]
ID_rec_myo_imp5 <- df_window.1.4$new_PTID[df_window.1.4$rec_myo_imp5 == "Recovery"]

tst <- c(unique(ID_rec_AIS_motor_imp),
         unique(ID_rec_AIS_motor_imp2), 
         unique(ID_rec_LEMS_imp5),
         unique(ID_rec_LEMS_imp2),
         unique(ID_rec_myo_imp2),
         unique(ID_rec_myo_imp5))
tst2 <- tst[duplicated(tst, nmax=2)]
tst[duplicated(tst)]

length(unique(tst[table(tst) == 6]))


df_window.1.4 %>% filter(rec_AIS_motor_imp == 'Recovery') %>%
  dplyr::select(AIS.1, AIS.52, level,
                Lower01, Lower26, Lower52,
                Upper01, Upper26, Upper52,
                LTSCOR01, LTSCOR26, LTSCOR52,
                PPSCOR01, PPSCOR26, PPSCOR52, center)

 
