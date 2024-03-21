################################################################################
# SCI - Worsening of patients
# L. Bourguignon
# First version : 09.02.2023
# Last update : 10.02.2023
################################################################################
# Load data

setwd('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/miraculous-recovery-SCI')
raw_data_Sygen <- read.csv('/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')
raw_data_Murnau <- read.csv('/Volumes/green_groups_bds_public/Data/Murnau/HematologicalBiomark_DATA_2021-01-07_0729.csv')
raw_data_EMSCI <- read.csv('/Volumes/green_groups_bds_public/Data/EMSCI/emsci_data_2020.csv')

################################################################################
# Libraries

library(networkD3)
library("survminer")
library("survival")

################################################################################
# Functions

transition_fct <- function(data_test, source, target){
  transition <- data_test %>%
    group_by_at(source) %>%
    count(!!!syms(target))
  transition <- as.data.frame(transition)
  names(transition) <- c('source', 'target', 'value')
  return(transition)
}

################################################################################
# Sygen

# Select for patients that were classified C or D at both week 1 and 4
# i.e. true C/D patients
true_CD <- raw_data_Sygen %>%
  filter(ais1 %in% c('AIS C', 'AIS D'), 
         ais4 %in% c('AIS C', 'AIS D'))

# Among true C/D patients, find the ones that regressed to A or B at week 52
true_CD_worsen <- true_CD %>%
  filter(ais52 %in% c('AIS A', 'AIS B'))
# Only one patient identified in the Sygen data 

true_CD_worsen %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52)
# Started with a AIS C and transition to AIS C from week 8 to 16
true_CD_worsen %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52)
true_CD_worsen %>%
  select(ptid, vaccd00, vaccd01, vaccd04, vaccd08, vaccd16, vaccd26, vaccd52)
true_CD_worsen %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, splvl)
true_CD_worsen %>%
  select(ptid, anyana01, anyana04, anyana08, anyana16, anyana26, anyana52)
# Regression driven by loss of VAC --> became motor complete

################################################################################
# EMSCI 

# Identify patients that were classified C or D at both week 1 and 4
# i.e. true C/D patients
CD_EMSCI_very_acute <- raw_data_EMSCI %>%
  filter(AIS %in% c('C', 'D'), ExamStage == 'very acute')
CD_EMSCI_acute1 <- raw_data_EMSCI %>%
  filter(AIS %in% c('C', 'D'), ExamStage == 'acute I')

ID_trueCD_EMSCI <- intersect(CD_EMSCI_very_acute[['Patientennummer']], 
                            CD_EMSCI_acute1[['Patientennummer']])

# Select for true C/D patients
true_CD_EMSCI <- subset(raw_data_EMSCI, Patientennummer %in% ID_trueCD_EMSCI)

# Among true C/D patients, identify the ones that regressed to A or B at week 52
temp_true_CD_EMSCI <- true_CD_EMSCI %>% filter(AIS %in% c('A', 'B'), ExamStage == 'chronic')
ID_trueCD_worsen_EMSCI <- temp_true_CD_EMSCI[['Patientennummer']]

# Select for true C/D patients that regressed to A or B at week 52
true_CD_worsen_EMSCI <- subset(true_CD_EMSCI, Patientennummer %in% ID_trueCD_worsen_EMSCI)
true_CD_worsen_EMSCI$ExamStage_weeks <- as.factor(true_CD_worsen_EMSCI$ExamStage_weeks)
# 8 patients in total

# Evolution of their AIS grades
true_CD_worsen_EMSCI_subset <- true_CD_worsen_EMSCI %>% select(Patientennummer, AIS, ExamStage_weeks)
true_CD_worsen_EMSCI_subset_wide <- reshape(true_CD_worsen_EMSCI_subset, 
                                            idvar = "Patientennummer", 
                                            timevar = "ExamStage_weeks", 
                                            direction = "wide")
true_CD_worsen_EMSCI_subset_wide[,c(1,3,2,4,5,6)]

# Evolution of their VAC
true_CD_worsen_EMSCI_VAC <- true_CD_worsen_EMSCI %>% select(Patientennummer, VAC, ExamStage_weeks)
true_CD_worsen_EMSCI_VAC_wide <- reshape(true_CD_worsen_EMSCI_VAC, 
                                            idvar = "Patientennummer", 
                                            timevar = "ExamStage_weeks", 
                                            direction = "wide")
true_CD_worsen_EMSCI_VAC_wide[,c(1,3,2,4,5,6)]

# Evolution of their DAP
true_CD_worsen_EMSCI_DAP <- true_CD_worsen_EMSCI %>% select(Patientennummer, DAP, ExamStage_weeks)
true_CD_worsen_EMSCI_DAP_wide <- reshape(true_CD_worsen_EMSCI_DAP, 
                                         idvar = "Patientennummer", 
                                         timevar = "ExamStage_weeks", 
                                         direction = "wide")
true_CD_worsen_EMSCI_DAP_wide[,c(1,3,2,4,5,6)]

# Evolution of their LEMS
true_CD_worsen_EMSCI_LEMS <- true_CD_worsen_EMSCI %>% select(Patientennummer, LEMS, ExamStage_weeks)
true_CD_worsen_EMSCI_LEMS_wide <- reshape(true_CD_worsen_EMSCI_LEMS, 
                                         idvar = "Patientennummer", 
                                         timevar = "ExamStage_weeks", 
                                         direction = "wide")
true_CD_worsen_EMSCI_LEMS_wide[,c(1,3,2,4,5,6)]

# Evolution of their UEMS
true_CD_worsen_EMSCI_UEMS <- true_CD_worsen_EMSCI %>% select(Patientennummer, UEMS, ExamStage_weeks)
true_CD_worsen_EMSCI_UEMS_wide <- reshape(true_CD_worsen_EMSCI_UEMS, 
                                          idvar = "Patientennummer", 
                                          timevar = "ExamStage_weeks", 
                                          direction = "wide")
true_CD_worsen_EMSCI_UEMS_wide[,c(1,3,2,4,5,6)]

# Look at factors that determine completeness of injury 
# i.e. VAC, DAP and S45 sensory scores (PP and LT)
true_CD_worsen_EMSCI_completeness <- true_CD_worsen_EMSCI %>% select(Patientennummer, ExamStage_weeks,
                                                                     VAC, DAP, 
                                                                     LLT_S45, RLT_S45,
                                                                     LPP_S45, RPP_S45, AIS)

# How do those factors evolve over time? What led to the regression in AIS grade?
true_CD_worsen_EMSCI_completeness_c <- true_CD_worsen_EMSCI_completeness %>% filter(ExamStage_weeks == 52)
true_CD_worsen_EMSCI_completeness_c
true_CD_worsen_EMSCI_completeness_aiii <- true_CD_worsen_EMSCI_completeness %>% filter(ExamStage_weeks == 26)
true_CD_worsen_EMSCI_completeness_aiii
true_CD_worsen_EMSCI_completeness_aii <- true_CD_worsen_EMSCI_completeness %>% filter(ExamStage_weeks == 12)
true_CD_worsen_EMSCI_completeness_aii
true_CD_worsen_EMSCI_completeness_ai <- true_CD_worsen_EMSCI_completeness %>% filter(ExamStage_weeks == 4)
true_CD_worsen_EMSCI_completeness_ai
true_CD_worsen_EMSCI_completeness_va <- true_CD_worsen_EMSCI_completeness %>% filter(ExamStage_weeks == 2)
true_CD_worsen_EMSCI_completeness_va

# Look at their level of injury at baseline
true_CD_worsen_EMSCI %>%
  filter(ExamStage == 'very acute') %>% select(Patientennummer, NLI)
