################################################################################
# SCI - Outstanding recovery - 
# L. Bourguignon
# First version : 31.03.2023
# Last update : 04.04.2023
################################################################################

################################################################################
# Library
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(naniar)
library(forcats)

################################################################################
# Load data
raw_sygen_with_surgery <- read.csv('/Volumes/INANL/45_Lucie_Bourguignon/1_Data/1_Sygen/sygen_with_surgery.csv')
path_figures_save <- '/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/figures/'

################################################################################
# Select 11 recovery patients

ID_recover <- df_window.1.4 %>% filter(rec_AIS_motor_imp == 'Recovery') %>% dplyr::select(new_PTID)
ID_norecover <- df_window.1.4 %>% filter(rec_AIS_motor_imp == 'No recovery') %>% dplyr::select(new_PTID)
sygen_11recover <- raw_sygen_with_surgery %>% dplyr::filter(new_PTID %in% ID_recover[['new_PTID']])
window.1.4_norecovery <- df_window.1.4 %>% dplyr::filter(new_PTID %in% ID_norecover[['new_PTID']])

write.csv(sygen_11recover, "/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/shiny/sygen_11recover.csv", row.names=FALSE)

sygen_11recover <- read.csv("/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/shiny/sygen_11recover.csv")

################################################################################
# Select myotomes

myotomes_left <- c('ELBFLL', 'WREXTL', 'ELBEXL', 'FINFLL', 'FINABL', 'HIPFLL', 
                  'KNEEXL', 'ANKDOL', 'GRETOL', 'ANKPLL')
myotomes_right <- c('ELBFLR', 'WREXTR', 'ELBEXR', 'FINFLR', 'FINABR', 'HIPFLR', 
                    'KNEETR', 'ANKDOR', 'GRETOR', 'ANKPLR')
levels_myotomes <- c("C5", "C6", "C7", 'C8', 'T1', 'L2', "L3", 'L4', 'L5', 'S1')
levels_dermatomes <- c("C2", "C3", "C4", "C5", "C6", "C7", 'C8', 'T1', 'T2', 
                       'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12',
                       'L1', 'L2', "L3", 'L4', 'L5', 'S1', 'S2', 'S3', 'S45')

df_heatmap_myo <- function(df, time, side){
  
  if (side == 'L'){
    myotomes <- myotomes_left
  } else if (side == 'R'){
    myotomes <- myotomes_right
  } else {
    print('Error: side must be R or L')
  }
  
  # Select relevant columns from the original data
  df <- df %>% 
    dplyr::select(one_of(paste0(myotomes, rep(time, length(myotomes)))))
  # Change column names to corresponding myotome level
  colnames(df) <- levels_myotomes
  # Add index column
  df <- cbind(index = rownames(df), df)
  rownames(df) <- 1:nrow(df)
  df_long <- melt(setDT(df), id.vars = c("index"), variable.name = "Myotomes")
  df_long$index <- paste0(df_long$index , side)
  df_long$value <- factor(df_long$value)
  
  return(df_long)
}

df_heatmap_ss <- function(df, time, side, test){
  
  temp_test <- paste0(test, side)
  col_sensory <- paste0(levels_dermatomes, 
                        rep(temp_test, length(levels_dermatomes)))

  # Select relevant columns from the original data
  df <- df %>% 
    dplyr::select(one_of(paste0(col_sensory, rep(time, length(col_sensory)))))
  # Change column names to corresponding dermatome level
  colnames(df) <- levels_dermatomes
  # Add index column
  df <- cbind(index = rownames(df), df)
  rownames(df) <- 1:nrow(df)
  df_long <- melt(setDT(df), id.vars = c("index"), variable.name = "Dermatomes")
  df_long$index <- paste0(df_long$index, side)
  df_long$value <- factor(df_long$value)
  
  return(df_long)
}

combined_df_heatmap_myo <- function(df, time){
  left_long <- df_heatmap_myo(df, time, side = 'L')
  right_long <- df_heatmap_myo(df, time, side = 'R')
  all_long <- rbind(left_long, right_long)
  all_long$Myotomes <- fct_rev(all_long$Myotomes)
  return(all_long)
}

combined_df_heatmap_ss <- function(df, time, test){
  left_long <- df_heatmap_ss(df, time, side = 'L', test)
  right_long <- df_heatmap_ss(df, time, side = 'R', test)
  all_long <- rbind(left_long, right_long)
  all_long$Dermatomes <- fct_rev(all_long$Dermatomes)
  return(all_long)
}

## Recovery
# Motor # replace by sygen_11recover for old definition of PR
sygen_11recover_motor_01_all_long <- combined_df_heatmap_myo(sygen_rec_lems_5, time = '01')
sygen_11recover_motor_04_all_long <- combined_df_heatmap_myo(sygen_rec_lems_5, time = '04')
sygen_11recover_motor_52_all_long <- combined_df_heatmap_myo(sygen_rec_lems_5, time = '52')

#PP
sygen_11recover_PP_01_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '01', 'PP')
sygen_11recover_PP_04_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '04', 'PP')
sygen_11recover_PP_52_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '52', 'PP')

#LT
sygen_11recover_LT_01_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '01', 'LT')
sygen_11recover_LT_04_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '04', 'LT')
sygen_11recover_LT_52_all_long <- combined_df_heatmap_ss(sygen_rec_lems_5, time = '52', 'LT')

## No recovery
# Motor # replace by window.1.4_norecovery for old definition of PR
window.1.4_norecovery_motor_01_all_long <- combined_df_heatmap_myo(sygen_norec_lems_5, time = '01')
window.1.4_norecovery_motor_04_all_long <- combined_df_heatmap_myo(sygen_norec_lems_5, time = '04')
window.1.4_norecovery_motor_52_all_long <- combined_df_heatmap_myo(sygen_norec_lems_5, time = '52')

#PP
window.1.4_norecovery_PP_01_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '01', 'PP')
window.1.4_norecovery_PP_04_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '04', 'PP')
window.1.4_norecovery_PP_52_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '52', 'PP')

#LT
window.1.4_norecovery_LT_01_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '01', 'LT')
window.1.4_norecovery_LT_04_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '04', 'LT')
window.1.4_norecovery_LT_52_all_long <- combined_df_heatmap_ss(sygen_norec_lems_5, time = '52', 'LT')

################################################################################
# Preparation plots

# Prepare colours
colours <- brewer.pal(n = 6, name = 'Blues')
group.colors <- c('0' = colours[1], '1' = colours[2], 
                  '2' = colours[3], '3' = colours[4], 
                  '4' = colours[5], '5' = colours[6], 
                  '9' = '#000000')

group.colors.diff <- c('0' = '#d8d8d8', '1' = '#EDF8E9', 
                       '2' = '#BAE4B3', '3' = '#74C476', 
                       '4' = "#31A354", '5' = "#006D2C", 
                       '-1' = "#FEE5D9", '-2' = "#FCAE91",
                       '-3' = "#FB6A4A", '-4' = "#DE2D26",
                       '-5' = "#A50F15")

# Reorder participants
## According to baseline TMS (old definition of PR)
#levels_baseline_TMS = c('10L', '10R', '9L', '9R', '2L', '2R', '5L', '5R', 
#                        '4L', '4R', '8L', '8R', '1L', '1R', '7L', '7R',
#                        '11L', '11R', '6L', '6R', '3L', '3R')
#levels_recovery_subgp <- c('3L', '3R', '7L', '7R', '8L', '8R', '6L', '6R',
#                           '1L', '1R', '10L', '10R', '2L', '2R', '5L', '5R',
#                           '9L', '9R', '4L', '4R', '11L', '11R')

levels_recovery_subgp <- c('2L', '2R', '4L', '4R', '5L', '5R', '3L', '3R',
                          '1L', '1R', '6L', '6R')

sygen_11recover_motor_01_all_long$index <- factor(
  sygen_11recover_motor_01_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_motor_04_all_long$index <- factor(
  sygen_11recover_motor_04_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_motor_52_all_long$index <- factor(
  sygen_11recover_motor_52_all_long$index,
  levels = levels_recovery_subgp)
                                                       
################################################################################
# Plots myotomes

# Week 01
heatmap_01 <- ggplot(sygen_11recover_motor_01_all_long, aes(index, Myotomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_01_PR-MS.png'), heatmap_01)

# Week 04
heatmap_04 <- ggplot(sygen_11recover_motor_04_all_long, aes(index, Myotomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors)+
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_04_PR-MS.pdf'), heatmap_04)

# Week 52
names(sygen_11recover_motor_52_all_long)[names(sygen_11recover_motor_52_all_long) == "value"] <- "value_recovery"
heatmap_52 <- ggplot(sygen_11recover_motor_52_all_long, aes(index, Myotomes, fill= value_recovery)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors)+
  geom_text(aes(label = value_recovery))
ggsave(paste0(path_figures_save, 'sygen_heatmap_52_PR-MS.pdf'), heatmap_52)


# Difference between week 1 and week 52
## Note: no LOCF needed, LEMS missing at week 52 for patients 1 and 6 due to
## 1 and 2 missing myotomes respectively
## Reason for 9 entered at RL3 for patient 1: RT KNEE CONTRACTED AT 90' - UBLE TO TEST
## Reason for NA at S1 for patient 6: not reported
sygen_11recover_motor <- merge(sygen_11recover_motor_52_all_long, 
                               sygen_11recover_motor_04_all_long)
sygen_11recover_motor <- sygen_11recover_motor %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(sygen_11recover_motor[,'value_recovery'][[1]])) - 
  as.numeric(as.character(sygen_11recover_motor[,4][[1]]))
sygen_11recover_motor = cbind(sygen_11recover_motor, difference=diff)

heatmap_diff <- ggplot(sygen_11recover_motor, aes(index, Myotomes, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))

ggsave(paste0(path_figures_save, 'sygen_heatmap_04-52_PR-MS.pdf'), heatmap_diff)

################################################################################
# Plots PP

# Old definition of PR
# levels_recovery_TPP = c('4L', '4R', '6L', '6R', '2L', '2R', '9L', '9R',
#                         '7L', '7R', '11L', '11R', '10L', '10R', '8L', '8R',
#                         '1L', '1R', '5L', '5R', '3L', '3R')

levels_recovery_TPP = c('6L', '6R', '1L', '1R', '5L', '5R', '3L', '3R',
                        '2L', '2R', '4L', '4R')

sygen_11recover_PP_01_all_long$index <- factor(
  sygen_11recover_PP_01_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_PP_04_all_long$index <- factor(
  sygen_11recover_PP_04_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_PP_52_all_long$index <- factor(
  sygen_11recover_PP_52_all_long$index,
  levels = levels_recovery_subgp)

# Week 01
heatmap_01_PP <- ggplot(sygen_11recover_PP_01_all_long, 
                        aes(index, Dermatomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_01_PR-PP.png'), heatmap_01_PP)

# Week 04
heatmap_04_PP <- ggplot(sygen_11recover_PP_04_all_long, 
                        aes(index, Dermatomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_04_PR-PP.png'), heatmap_04_PP)

# Week 52
## Note: no NA at week 52 so no LOCF needed here
names(sygen_11recover_PP_52_all_long)[names(sygen_11recover_PP_52_all_long) == "value"] <- "value_recovery"
heatmap_52_PP <- ggplot(sygen_11recover_PP_52_all_long, 
                        aes(index, Dermatomes, fill= value_recovery)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_52_PR-PP.png'), heatmap_52_PP)

# Diff plot
sygen_11recover_PP <- merge(sygen_11recover_PP_52_all_long, 
                               sygen_11recover_PP_04_all_long)
sygen_11recover_PP <- sygen_11recover_PP %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(sygen_11recover_PP[,'value_recovery'][[1]])) - 
  as.numeric(as.character(sygen_11recover_PP[,4][[1]]))
sygen_11recover_PP = cbind(sygen_11recover_PP, difference=diff)

heatmap_diff_PP <- ggplot(sygen_11recover_PP, aes(index, Dermatomes, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))

ggsave(paste0(path_figures_save, 'sygen_heatmap_04-52_PR-PP.png'), heatmap_diff_PP)

################################################################################
# Plots LT

# Old definition of PR
# levels_recovery_TLT = c('4L', '4R', '6L', '6R', '2L', '2R', '9L', '9R',
#                         '7L', '7R', '11L', '11R', '10L', '10R', '8L', '8R',
#                         '1L', '1R', '5L', '5R', '3L', '3R')

# order index based on MS for better comparability between plots
sygen_11recover_LT_01_all_long$index <- factor(
  sygen_11recover_LT_01_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_LT_04_all_long$index <- factor(
  sygen_11recover_LT_04_all_long$index,
  levels = levels_recovery_subgp)
sygen_11recover_LT_52_all_long$index <- factor(
  sygen_11recover_LT_52_all_long$index,
  levels = levels_recovery_subgp)

# Week 01
heatmap_01_LT <- ggplot(sygen_11recover_LT_01_all_long, 
                        aes(index, Dermatomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_01_PR-LT.png'), heatmap_01_LT)

# Week 04
heatmap_04_LT <- ggplot(sygen_11recover_LT_04_all_long, 
                        aes(index, Dermatomes, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_04_PR-LT.png'), heatmap_04_LT)

# Week 52
names(sygen_11recover_LT_52_all_long)[names(sygen_11recover_LT_52_all_long) == "value"] <- "value_recovery"
heatmap_52_LT <- ggplot(sygen_11recover_LT_52_all_long, 
                        aes(index, Dermatomes, fill= value_recovery)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
ggsave(paste0(path_figures_save, 'sygen_heatmap_52_PR-LT.png'), heatmap_52_LT)

# Diff plot
sygen_11recover_LT <- merge(sygen_11recover_LT_52_all_long, 
                            sygen_11recover_LT_01_all_long)
sygen_11recover_LT <- sygen_11recover_LT %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(sygen_11recover_LT[,'value_recovery'][[1]])) - 
  as.numeric(as.character(sygen_11recover_LT[,4][[1]]))
sygen_11recover_LT = cbind(sygen_11recover_LT, difference=diff)

heatmap_diff_LT <- ggplot(sygen_11recover_LT, aes(index, Dermatomes, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))

ggsave(paste0(path_figures_save, 'sygen_heatmap_01-52_PR-LT.png'), heatmap_diff_LT)

################################################################################
# Parallel coordinate plots

# Libraries
library(hrbrthemes)
library(GGally)
library(viridis)
library(panelr)

myotomes1_01_08 <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('01', length(myotomes_left))),
                paste0(myotomes_left, rep('04', length(myotomes_left))),
                paste0(myotomes_left, rep('08', length(myotomes_left)))))
temp_long <- long_panel(myotomes1_01_08, prefix = "0", begin = 1, end = 8, label_location = "end")
temp_long <- temp_long %>% 
  dplyr::filter(wave %in% c('1', '4', '8'))


df_myotomes1 <- data.frame(matrix(ncol = 7, nrow = length(myotomes_left)))
x <- c("myotome", "w1", "w4", "w8", "w16", "w26", "w52")
colnames(df_myotomes1) <- x
df_myotomes1$myotome <- paste0(levels_myotomes, rep('L', length(levels_myotomes)))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('01', length(myotomes_left)))))
df_myotomes1$w1 <- as.numeric(as.vector(temp[1,]))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('04', length(myotomes_left)))))
df_myotomes1$w4 <- as.numeric(as.vector(temp[1,]))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('08', length(myotomes_left)))))
df_myotomes1$w8 <- as.numeric(as.vector(temp[1,]))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('16', length(myotomes_left)))))
df_myotomes1$w16 <- as.numeric(as.vector(temp[1,]))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('26', length(myotomes_left)))))
df_myotomes1$w26 <- as.numeric(as.vector(temp[1,]))
temp <- raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(one_of(paste0(myotomes_left, rep('52', length(myotomes_left)))))
df_myotomes1$w52 <- as.numeric(as.vector(temp[1,]))
df_myotomes1$w52

## add jitter manually before plotting
ggparcoord(df_myotomes1,
           columns = 2:7, groupColumn = 1,
           showPoints = TRUE, 
           scale="globalminmax",
           title = "Parallel Coordinate Plot for the Iris Data",
           alphaLines = 0.3,
) + 
  scale_color_viridis(discrete=TRUE) +
  #theme_ipsum()+
  theme(
    plot.title = element_text(size=10)
  )

remotes::install_github("yaweige/yaweitest")
library(yaweitest)
myparallel(data = df_myotomes1,
           columns = 2:7, groupColumn = 1,
           showPoints = TRUE, 
           #scale="globalminmax",
           jitter = TRUE,
           jitterratio = 0.05)+ 
  scale_color_viridis(discrete=TRUE) +
  #theme_ipsum()+
  theme(
    plot.title = element_text(size=10)
  )


####
raw_sygen_with_surgery$BMI <- raw_sygen_with_surgery$WTM / 
  ((raw_sygen_with_surgery$HTM/100)**2)
raw_sygen_with_surgery %>% 
  dplyr::filter(new_PTID %in% ID_recover[['new_PTID']][1]) %>% 
  select(SEXCD, AGE, INJCD, RACECD, HTM, WTM, BMI, center,
         INJDT_date, INJTM_time2, center, TX1_R, EC_IM,
         SPLVL, CARDCD01, CPULCD01, EENTCD01, GASTCD01, GENICD01, 
         HEADCD01, MUSCCD01, NUMEP01, PULMCD01, SKINCD01,
         )

display.brewer.pal(3, 'Blues')
brewer.pal(3, 'Blues')

################################################################################
## Examine the reasons for 9 being reported

sygen_11recover %>% select(CTMRCOM1, CTMRCOM2, CTMRCOM3,CTMRCOM4, CTMRCOM5,
                           CTMRCOM6, CTMRCOM7, CTMRCOM8, CTMRCOM9)

## Patient 1: 
# Week 52 motor L3 R 
sygen_11recover$KENTRT52[1] # "others"
sygen_11recover$AMOTC152[1] # "RT KNEE CONTRACTED AT 90' - UBLE TO TEST"
# Week 52 PP S1 R
sygen_11recover$ASECO152[1] # "S1 PRESSURE SORE"

## Patient 3:
# Week 1 motor L3 L
sygen_11recover$KENTRT01[3]
sygen_11recover$AMOTC101[3] # "(L) TIB/FIB FX - IN CAST" --> Why are PP and LT present then?

## Patient 4:
# Week 1 PP/LT C2 R and L
sygen_11recover$ASECO101[4] # "PHILLY COLLAR"
# Week 4 PP/LT C2 R and L
sygen_11recover$ASECO104[4] # "CERVICAL COLLAR, PT CAN FEEL WHEN SHE HAVE BOWEL MOVEMENT"

## Patient 7:
# Week 1 PP/LT C2-C4 R/L
sygen_11recover$ASECO101[7] # "PHILADELPHIA COLLAR RESTRICTS TESTING."
# Week 52 PP/LT S45 R/L
sygen_11recover$ASECO152[7] # "S4-5 NOT TESTED D/T TIME CONSTRAINTS"

## Patient 8:
# Week 4 PP/LT T5-T8 R/L
sygen_11recover$ASECO104[8] # "HALO VEST T5-T8"

## Patient 9:
# Week 4 PP/LT T6-T8 R/L
sygen_11recover$ASECO104[9] # "9-HALO"

## Patient 10:
# Week 52 motor S1 R/L
sygen_11recover$APNTLT52[10]
sygen_11recover$AMOTC152[10]
sygen_11recover$ANKPLL52[10]

################################################################################
## Visualise the dates of injury
plot(raw_sygen_with_surgery$INJDT)

ggplot(df_window.1.4, aes(x = INJDT, fill = rec_AIS_motor_imp)) + 
  geom_histogram() +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)
