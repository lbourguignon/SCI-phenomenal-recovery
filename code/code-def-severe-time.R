################################################################################
# SCI - Outstanding recovery - Changing time window for severe definition
# L. Bourguignon
# First version : 17.02.2023
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
## Load data
# Note that the Sygen file on the BMDS group share does not include date and time of surgery

raw_sygen_with_surgery <- read.csv('/Volumes/INANL/45_Lucie_Bourguignon/1_Data/1_Sygen/sygen_with_surgery.csv')

## Create variable with time to first closed treatment
raw_sygen_with_surgery$CTSTTM01_minutes <- as.numeric((raw_sygen_with_surgery$CTSTTM01%%1)*100)
raw_sygen_with_surgery$CTSTTM01_hours <- as.numeric(trunc(raw_sygen_with_surgery$CTSTTM01))
raw_sygen_with_surgery$CTSTTM01_time2 <- sprintf("%s:%s:00", 
                                                 raw_sygen_with_surgery$CTSTTM01_hours, 
                                                 as.integer(raw_sygen_with_surgery$CTSTTM01_minutes))
raw_sygen_with_surgery["CTSTTM01_time2"][raw_sygen_with_surgery["CTSTTM01_time2"] == "NA:NA:00"] <- NA
raw_sygen_with_surgery$CTSTTM01_time_combined2 <- as.POSIXct(paste(raw_sygen_with_surgery$CTSTDT01_date,
                                                                   raw_sygen_with_surgery$CTSTTM01_time2), 
                                                             format="%m/%d/%Y %H:%M:%S")

## Create variable with time to first operative treatment
raw_sygen_with_surgery$OTSTDT01_date <- num_date(raw_sygen_with_surgery$OTSTDT01)

raw_sygen_with_surgery$OTSTTM01_minutes <- as.numeric((raw_sygen_with_surgery$OTSTTM01%%1)*100)
raw_sygen_with_surgery$OTSTTM01_hours <- as.numeric(trunc(raw_sygen_with_surgery$OTSTTM01))
raw_sygen_with_surgery$OTSTTM01_time2 <- sprintf("%s:%s:00", 
                                                 raw_sygen_with_surgery$OTSTTM01_hours, 
                                                 as.integer(raw_sygen_with_surgery$OTSTTM01_minutes))
raw_sygen_with_surgery["OTSTTM01_time2"][raw_sygen_with_surgery["OTSTTM01_time2"] == "NA:NA:00"] <- NA
raw_sygen_with_surgery$OTSTTM01_time_combined2 <- as.POSIXct(paste(raw_sygen_with_surgery$OTSTDT01_date,
                                                                   raw_sygen_with_surgery$OTSTTM01_time2), 
                                                             format="%Y-%m-%d %H:%M:%S")

## Create variable with time of injury
raw_sygen_with_surgery$INJTM_minutes <- as.numeric((raw_sygen_with_surgery$INJTM%%1)*100)
raw_sygen_with_surgery$INJTM_hours <- as.numeric(trunc(raw_sygen_with_surgery$INJTM))
raw_sygen_with_surgery$INJTM_time2 <- sprintf("%s:%s:00", 
                                              raw_sygen_with_surgery$INJTM_hours, 
                                              as.integer(raw_sygen_with_surgery$INJTM_minutes))
raw_sygen_with_surgery["INJTM_time2"][raw_sygen_with_surgery["INJTM_time2"] == "NA:NA:00"] <- NA
raw_sygen_with_surgery$INJTM_time_combined2 <- as.POSIXct(paste(raw_sygen_with_surgery$INJDT_date,
                                                                   raw_sygen_with_surgery$INJTM_time2), 
                                                             format="%m/%d/%Y %H:%M:%S")

raw_sygen_with_surgery$CTST01_time_combined <- as.POSIXct(paste(raw_sygen_with_surgery$CTSTDT01_date,
                                                             raw_sygen_with_surgery$CTSTTM01_time),
                                                       format="%m/%d/%Y %H:%M:%S")

raw_sygen_with_surgery$INJ_time_combined <- as.POSIXct(paste(raw_sygen_with_surgery$INJDT_date,
                                                             raw_sygen_with_surgery$INJTM_time),
                                                       format="%m/%d/%Y %H:%M:%S")

raw_sygen_with_surgery$time_decompression_surgery_old <- difftime(raw_sygen_with_surgery$CTST01_time_combined,
                                                              raw_sygen_with_surgery$INJ_time_combined, units="mins")

raw_sygen_with_surgery$time_decompression_surgery <- difftime(raw_sygen_with_surgery$OTSTTM01_time_combined2,
                                                              raw_sygen_with_surgery$INJTM_time_combined2, units="mins")

################################################################################
## Functions

define_df_window <- function(raw_sygen_with_surgery, window){
  
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
  
  #Find the maximum value of LEMS recorded in the first weeks
  values_lems <- paste0(rep('Lower0', length(window)), window)
  subset_lems <- true_A[values_lems]
  subset_lems <- subset_lems %>% rowwise() %>%
    dplyr::mutate(max_early_lems = max(c_across(where(is.numeric)), na.rm = TRUE))
  subset_lems[subset_lems == -Inf] <- NA
  true_A$max_early_lems <- subset_lems$max_early_lems
  
  #Find the maximum value of UEMS recorded in the first weeks
  values_uems <- paste0(rep('Upper0', length(window)), window)
  subset_uems <- true_A[values_uems]
  subset_uems <- subset_uems %>% rowwise() %>%
    dplyr::mutate(max_early_uems = max(c_across(where(is.numeric)), na.rm = TRUE))
  subset_uems[subset_uems == -Inf] <- NA
  true_A$max_early_uems <- subset_uems$max_early_uems
  
  levels <- c('L1', 'L2', 'L3', 'L4', 'L5', 'S1', 'S2', 'S3', 'S45')
  sides <- c('L', 'R')
  tests <- c('PP', 'LT')
  earlyweeks <- c(paste0(rep('0', length(window)), window))
  earlyweeks <- earlyweeks[!earlyweeks %in% c('00')]
  true_A$early_sensation_left <- 0
  
  for (i in c(1:dim(true_A)[1])){ # for all patients with AIS A at week 1 and 4
    # Define LEMS improvement
    
    if (is.na(true_A$Lower52[i])){ # if LEMS at week 52 is missing
      value_lems = true_A$Lower26[i] # consider LEMS at week 26
    } else {
      value_lems = true_A$Lower52[i] # otherwise consider LEMS at week 52
    }
    # calculate the difference in maximum value of LEMS recorded in the first weeks
    # and the LEMS at recovery as defined above
    if (is.na(value_lems) | is.na(true_A$max_early_lems[i])){ # if both week 52 and 26 LEMS are missing
      true_A$bin_lems_improv[i] <- NA # consider that the improvement info is missing for LEMS
    } else if (value_lems - true_A$max_early_lems[i] > 5){ # if this difference is more than 5 points
      true_A$bin_lems_improv[i] <- 1 # consider that the patient improved in terms of LEMS
    } else {
      true_A$bin_lems_improv[i] <- 0 # otherwise consider that the patient did not improve in terms of LEMS
    }
    
    # Define UEMS improvement
    
    if (is.na(true_A$Upper52[i])){ # if UEMS at week 52 is missing
      value_uems = true_A$Upper26[i] # consider UEMS at week 26
    } else {
      value_uems = true_A$Upper52[i] # otherwise consider LEMS at week 52
    }
    # calculate the difference in maximum value of UEMS recorded in the first weeks
    # and the UEMS at recovery as defined above
    if (is.na(value_uems) | is.na(true_A$max_early_uems[i])){ # if both week 52 and 26 UEMS are missing
      true_A$bin_uems_improv[i] <- NA # consider that the improvement info is missing for UEMS
    } else if (value_uems - true_A$max_early_uems[i] > 5){ # if this difference is more than 5 points
      true_A$bin_uems_improv[i] <- 1 # consider that the patient improved in terms of UEMS
    } else {
      true_A$bin_uems_improv[i] <- 0 # otherwise consider that the patient did not improve in terms of UEMS
    }
    
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
  
  print('Number of patients with strict severe injury in the time window requested')
  print(dim(stricter_true_A)[1])
  
  # Define the subset that recovers
  true_A_recover_new <- stricter_true_A %>%
    filter((AIS.52 %in% c('AIS C', 'AIS D') | bin_lems_improv == 1) &
             (bin_uems_improv == 1 | bin_lems_improv == 1))
  
  for (i in c(1:dim(true_A_recover_new)[1])){
    for (time in c("8", "16", "26", "52")){
      col = paste0('AIS.', time)
      if (!(true_A_recover_new[i, col] == 'AIS A' | true_A_recover_new[i, col] == 'ND')){
        true_A_recover_new[i, 'time_conversion'] = time
        break
      }
    }
  }
  
  # print('Some characteristics of the recoverers')
  # print(true_A_recover_new %>%
  #         dplyr::select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, time_conversion))
  
  true_A_copy_new <- stricter_true_A
  true_A_copy_new <- true_A_copy_new %>% 
    mutate(level = if_else(str_detect(true_A_copy_new$SPLVL, "T"), "T", "C"))
  
  id_recoverers <- true_A_recover_new$new_PTID
  true_A_copy_new$recover <- 'No recovery'
  true_A_copy_new$recover[true_A_copy_new$new_PTID %in% id_recoverers] <- "Recovery"
  
  return(true_A_copy_new)
  
}

make_table_1 <- function(df, window){
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
  
  return(table1(~ TX1_R + SEXCD + AGE + level + AIS.52 + ANYANA52 + VACCD52 + 
           bin_uems_improv + bin_lems_improv + DECOMC01 + DECOCD01 + time_decompression_surgery | recover, 
           caption = caption, data = df_copy))
  
}

plot_surgery <- function(true_A_copy_new){
  sub <- true_A_copy_new %>% filter(DECOCD01 %in% c(0,1))
  sub <- sub %>% mutate(DECOCD01_label = case_when(DECOCD01 == 0 ~ 'Decompression not obtained',
                                                   DECOCD01 == 1 ~ 'Decompression obtained'))
  
  Summary.data <- sub %>% 
    dplyr::group_by(DECOCD01_label, recover) %>% 
    summarise(n = sum(!(is.na(time_decompression_surgery))), 
              Q1 = quantile(time_decompression_surgery, na.rm=TRUE)[2],
              median = quantile(time_decompression_surgery, na.rm=TRUE)[3],
              Q3 = quantile(time_decompression_surgery, na.rm=TRUE)[4])
  
  print(Summary.data)

  p <- ggplot(sub, aes(x=recover, y=as.numeric(time_decompression_surgery))) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.05) +
    geom_jitter(shape=16, position=position_jitter(0.05)) +
    theme_classic() +
    facet_grid(rows = vars(DECOCD01_label)) +
    geom_text(data=Summary.data, aes(x = recover, y = 4000, label=paste0('n=', n)), color="red", fontface =2, size = 5) + 
    labs(x = 'Recovery status', y = 'Minutes to decompression') +
    theme(axis.title.x = element_text(size = 14, face = 'bold'),
          axis.title.y = element_text(size = 14, face = 'bold'),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          strip.text = element_text(size = 14, face = 'bold'))
  return(p)
}

################################################################################
## Perform the analysis

df_window.0 <- define_df_window(raw_sygen_with_surgery, c('0'))
df_window.1 <- define_df_window(raw_sygen_with_surgery, c('1'))
df_window.4 <- define_df_window(raw_sygen_with_surgery, c('4'))
df_window.0.1 <- define_df_window(raw_sygen_with_surgery, c('0', '1'))
df_window.1.4 <- define_df_window(raw_sygen_with_surgery, c('1', '4'))
df_window.0.1.4 <- define_df_window(raw_sygen_with_surgery, c('0', '1', '4'))
df_window.1.4.8 <- define_df_window(raw_sygen_with_surgery, c('1', '4', '8'))

table1_window.0 <- make_table_1(df_window.0, c('0'))
table1_window.1 <- make_table_1(df_window.1, c('1'))
table1_window.4 <- make_table_1(df_window.4, c('4'))
table1_window.0.1 <- make_table_1(df_window.0.1, c('0', '1'))
table1_window.1.4 <- make_table_1(df_window.1.4, c('1', '4'))
table1_window.0.1.4 <- make_table_1(df_window.0.1.4, c('0', '1', '4'))
table1_window.1.4.8 <- make_table_1(df_window.1.4.8, c('1', '4', '8'))

table1_window.0
table1_window.1
table1_window.4
table1_window.0.1
table1_window.1.4
table1_window.0.1.4

surgery_window.0 <- plot_surgery(df_window.0)
surgery_window.0
surgery_window.1 <- plot_surgery(df_window.1)
surgery_window.1
surgery_window.4 <- plot_surgery(df_window.4)
surgery_window.4
surgery_window.0.1 <- plot_surgery(df_window.0.1)
surgery_window.0.1
surgery_window.1.4 <- plot_surgery(df_window.1.4)
surgery_window.1.4
surgery_window.0.1.4 <- plot_surgery(df_window.0.1.4)
surgery_window.0.1.4

write.csv(df_window.0.1.4, "/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.0.1.4.csv", row.names=FALSE)
write.csv(df_window.0, "/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.0.csv", row.names=FALSE)
write.csv(df_window.1, "/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.1.csv", row.names=FALSE)
write.csv(df_window.4, "/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.4.csv", row.names=FALSE)
write.csv(df_window.1.4, "/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/df_window.1.4.csv", row.names=FALSE)

################################################################################
## attempts to calculate the time to first decompression

# dim(df_window.0)
# sum(df_window.0$decocd01, na.rm=T)
# sum(df_window.0$decocd02, na.rm=T)
# sum(df_window.0$decocd03, na.rm=T)
# sum(df_window.0$decocd04, na.rm=T)
# sum(df_window.0$decocd08, na.rm=T)
# 
# sum(df_window.0$decomc01, na.rm=T)
# sum(df_window.0$decomc02, na.rm=T)
# sum(df_window.0$decomc03, na.rm=T)
# sum(df_window.0$decomc04, na.rm=T)
# sum(df_window.0$decomc08, na.rm=T)
# 
# df_window.0_copy <- df_window.0
# df_window.0_copy$time_decompression <- NA
# for (i in dim(df_window.0_copy)[1]){
#   if(df_window.0_copy$decocd01[i] == 1){
#     
#   }
# }
# 
# sum(is.na(raw_sygen_with_surgery$CTSTDT01))
# 
# raw_sygen_with_surgery$CTSTDT01_date
# raw_sygen_with_surgery$CTSTTM01_time
# raw_sygen_with_surgery$INJTM_time
# raw_sygen_with_surgery$INJDT_date

raw_sygen_with_surgery %>% filter(time_decompression_surgery < 0) %>% 
  select (time_decompression_surgery, CTST01_time_combined, CTSTDT01_date, CTSTTM01_time, 
          INJ_time_combined, INJDT_date, INJTM_time, INJCD, INJOTHE1)

df_window.1.4 %>% dplyr::filter(rec_AIS_motor_imp == 'Recovery') %>%
  dplyr::select(new_PTID, AIS.1, AIS.4, AIS.8, AIS.16, AIS.26, AIS.52, 
                bin_lems_improv, bin_uems_improv, time_decompression_surgery,
                level, ANYANA52, VACCD52)

df_window.1.4 %>% filter(recover == 'Recovery') %>%
  dplyr::select(new_PTID, AIS.1, AIS.4, AIS.8, AIS.16, AIS.26, AIS.52, 
                ANYANA01, ANYANA04, ANYANA08, ANYANA16, ANYANA26, ANYANA52,
                VACCD01, VACCD04, VACCD08, VACCD16, VACCD26, VACCD52,
                level)


df_window.1.4 %>% filter(recover == 'Recovery') %>%
  dplyr::select(new_PTID, OTSTDT01, PASTDT01, CTSTDT01, DECOMC01, DECOCD01)
                
                # pasttm01
                # 
                # OTSTDT01_date,
                # OTSTTM01, OTSTTM01_hours, OTSTTM01_minutes)
                
                # INJTM_time_combined2,
                # DECOMC01, DECOCD01,
                # level)
# 
# sum(is.na(raw_sygen_with_surgery$OTSTDT01))
# sum(is.na(raw_sygen_with_surgery$CTSTTM01))
# 
# raw_sygen_with_surgery$OTSTTM01_minutes <- as.numeric((raw_sygen_with_surgery$OTSTTM01%%1)*100)
# raw_sygen_with_surgery$OTSTTM01_hours <- as.numeric(trunc(raw_sygen_with_surgery$OTSTTM01))
# raw_sygen_with_surgery$OTSTDT01_date <- num_date(raw_sygen_with_surgery$OTSTDT01)

df_window.1.4 %>% filter(recover == 'Recovery') %>%
  dplyr::select(new_PTID)

df_window.4_nrecover <- df_window.4 %>% filter(recover == 'No recovery')
table(df_window.4_nrecover$bin_lems_improv + df_window.4_nrecover$bin_uems_improv)

df_window.1.4$year_injury <- substr(df_window.1.4$INJDT_date, 7, 10)
print(table1(formula(paste("~ factor(SEXCD) + AGE + Lower04 + Upper04 + factor(AIS.4) + SPLVL + factor(INJCD) + 
                           factor(INJDT_date) | factor(", 'rec_LEMS_imp5', ")")),
             data=df_window.1.4[!(is.na(df_window.1.4$rec_LEMS_imp5)),], overall=F, extra.col=list(`P-value`=pvalue)))

table(~ SEXCD | factor(rec_LEMS_imp5), data=df_window.1.4[!(is.na(df_window.1.4$rec_LEMS_imp5))], overall=F)

