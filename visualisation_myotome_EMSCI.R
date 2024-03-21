library(data.table)
library(forcats)
library(RColorBrewer)
library(naniar)

levels_myotomes <- c("C5", "C6", "C7", 'C8', 'T1', 'L2', "L3", 'L4', 'L5', 'S1')
levels_dermatomes <- c("C2", "C3", "C4", "C5", "C6", "C7", 'C8', 'T1', 'T2', 
                       'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12',
                       'L1', 'L2', "L3", 'L4', 'L5', 'S1', 'S2', 'S3', 'S45')

df_heatmap <- function(time, side, test, ID){
  # test in "MS", "LT" and "PP"
  # side in "L" and "R"
  # time in "very acute", "acute I", "acute II", "acute III", and "chronic"
  
  temp_test <- paste0(side, test)
  if (test == 'MS'){
    cols <- c(paste0(rep(temp_test, length(levels_myotomes)), '_',
                   levels_myotomes), 'Patientennummer')
  } else {
    cols <- c(paste0(rep(temp_test, length(levels_dermatomes)), '_',
                   levels_dermatomes), 'Patientennummer')
  }
  
  # Select relevant columns from the original data
  df <- raw_data_EMSCI %>% 
    dplyr::filter(ExamStage == time,
           Patientennummer %in% ID) %>% 
    dplyr::select(one_of(cols))
  df[df==""] <- NA
  if (time == 'chronic'){
    for (i in c(1:dim(df)[1])){
      if (is.na(df[i, 1])){
        locf <- raw_data_EMSCI %>% 
          dplyr::filter(ExamStage == 'acute III',
                 Patientennummer %in% df$Patientennummer[i]) %>% 
          dplyr::select(one_of(cols))
        df[i,] <- locf[1,]
      }
    }
  }
  # Change column names to corresponding dermatome level
  #colnames(df) <- levels_dermatomes
  # Add index column
  df <- cbind(index = df$Patientennummer, df)
  rownames(df) <- 1:nrow(df)
  df_ohnePat <-  subset(df, select=-c(Patientennummer))
  df_long <- melt(setDT(df_ohnePat), id.vars = c("index"), variable.name = "Level")
  df_long$index <- paste0(df_long$index, side)
  df_long$value <- factor(df_long$value)
  
  return(df_long)
}

combined_df_heatmap <- function(time, test, ID){
  left_long <- df_heatmap(time, side = 'L', test, ID)
  right_long <- df_heatmap(time, side = 'R', test, ID)
  all_long <- rbind(left_long, right_long)
  all_long$Level <- substr(all_long$Level, 5, 
         nchar(as.vector(all_long$Level)))
  if (test == 'MS'){
    all_long$Level <- factor(all_long$Level, 
                             levels = levels_myotomes)
  } else {
    all_long$Level <- factor(all_long$Level, 
                             levels = levels_dermatomes)
  }
  all_long$Level <- fct_rev(all_long$Level)
  
  return(all_long)
}

## Recovery
# Motor
emsci_11recover_motor_01_all_long <- combined_df_heatmap(time = 'very acute', 'MS', id_recoverers_EMSCI[[1]])
emsci_11recover_motor_04_all_long <- combined_df_heatmap(time = 'acute I', 'MS', id_recoverers_EMSCI[[1]])
emsci_11recover_motor_52_all_long <- combined_df_heatmap(time = 'chronic', 'MS', id_recoverers_EMSCI[[1]])

# Motor stats def all cervical and thoracic
psa_matched_all_EMSCI_outliers <- psa_matched_true_severe %>% filter(tail95_blwnli == 1)
ID_out_all <- psa_matched_all_EMSCI_outliers %>% select(Patientennummer)
order_out_all <- levels(fct_reorder(factor(psa_matched_all_EMSCI_outliers$Patientennummer), psa_matched_all_EMSCI_outliers$UEMS_locf, min))
emsci_psmall_motor_01_all_long <- combined_df_heatmap(time = 'very acute', 'MS', ID_out_all[[1]])
emsci_psmall_motor_04_all_long <- combined_df_heatmap(time = 'acute I', 'MS', ID_out_all[[1]])
emsci_psmall_motor_52_all_long <- combined_df_heatmap(time = 'chronic', 'MS', ID_out_all[[1]])

#PP
emsci_11recover_PP_01_all_long <- combined_df_heatmap(time = 'very acute', 'PP', id_recoverers_EMSCI[[1]])
emsci_11recover_PP_04_all_long <- combined_df_heatmap(time = 'acute I', 'PP', id_recoverers_EMSCI[[1]])
emsci_11recover_PP_52_all_long <- combined_df_heatmap(time = 'chronic', 'PP', id_recoverers_EMSCI[[1]])

#LT
emsci_11recover_LT_01_all_long <- combined_df_heatmap(time = 'very acute', 'LT', id_recoverers_EMSCI[[1]])
emsci_11recover_LT_04_all_long <- combined_df_heatmap(time = 'acute I', 'LT', id_recoverers_EMSCI[[1]])
emsci_11recover_LT_52_all_long <- combined_df_heatmap(time = 'chronic', 'LT', id_recoverers_EMSCI[[1]])


#levels_recovery_subgp_emsci <- c('2L', '2R', '7L', '7R', '10L', '10R', '11L', '11R','9L', '9R',
#                           '17L', '17R', '4L', '4R', '3L', '3R', '15L', '15R',
#                           '1L', '1R', '12L', '12R', '8L', '8R',
#                           '14L', '14R', "6L", "6R", "5L", "5R",
#                           '16L', '16R', '13L', '13R', '18L', '18R',
#                           '19L', '19R', '20L', '20R', '21L', '21R', '22L', '22R')

# levels_recovery_subgp_emsci <- c('60096L', '60096R', '100038L', '100038R',
#                                  '90158L', '90158R', '20111L', '20111R',
#                                  '20150L', '20150R', '50007L', '50007R',
#                                  '90036L', '90036R', '300179L', '300179R',
#                                  '20057L', '20057R', '60421L', '60421R',
#                                  '555320L', '555320R', '555226L', '555226R',
#                                  '533049L', '533049R', '90240L', '90240R',
#                                  '90026L', '90026R', '230046L', '230046R',
#                                  '40030L', '40030R', '340010L', '340010R',
#                                  '60202L', '60202R', '90008L', '90008R',
#                                  '90260L', '90260R', '60054L', '60054R')

levels_recovery_subgp_emsci <- c('60054L', '60054R', '90260L', '90260R',
                                 '90008L', '90008R', '60202L', '60202R',
                                 '340010L', '340010R', '40030L', '40030R',
                                 '230046L', '230046R', '90026L', '90026R',
                                 '90240L', '90240R', '533049L', '533049R',
                                 '555226L', '555226R', '555320L', '555320R',
                                 '60421L', '60421R', '20057L', '20057R',
                                 '300179L', '300179R', '90036L', '90036R',
                                 '50007L', '50007R', '20150L', '20150R',
                                 '20111L', '20111R', '90158L', '90158R',
                                 '100038L', '100038R', '60096L', '60096R')


emsci_11recover_motor_01_all_long$index <- factor(
  emsci_11recover_motor_01_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_motor_04_all_long$index <- factor(
  emsci_11recover_motor_04_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_motor_52_all_long$index <- factor(
  emsci_11recover_motor_52_all_long$index,
  levels = levels_recovery_subgp_emsci)

emsci_psmall_motor_01_all_long$index <- factor(
  emsci_psmall_motor_01_all_long$index,
  levels = paste0(rep(order_out_all, each=2), c('L', 'R')))
emsci_psmall_motor_04_all_long$index <- factor(
  emsci_psmall_motor_04_all_long$index,
  levels = paste0(rep(order_out_all, each=2), c('L', 'R')))
emsci_psmall_motor_52_all_long$index <- factor(
  emsci_psmall_motor_52_all_long$index,
  levels = paste0(rep(order_out_all, each=2), c('L', 'R')))

heatmap_01_emsci <- ggplot(emsci_11recover_motor_01_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_01_emsci
ggsave(paste0(path_figures_save, 'emsci_heatmap_01_PR-MS.png'), heatmap_01_emsci)


heatmap_01_emsci_psmall <- ggplot(emsci_psmall_motor_01_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_01_emsci_psmall

heatmap_04_emsci_psmall <- ggplot(emsci_psmall_motor_04_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_04_emsci_psmall

heatmap_52_emsci_psmall <- ggplot(emsci_psmall_motor_52_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_52_emsci_psmall
#ggsave(paste0(path_figures_save, 'emsci_heatmap_01_PR-MS.png'), heatmap_01_emsci)


heatmap_04_emsci <- ggplot(emsci_11recover_motor_04_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_04_emsci
ggsave(paste0(path_figures_save, 'emsci_heatmap_04_PR-MS.pdf'), heatmap_04_emsci)


# note: 555226 and 90240 had NA for week 52 and thus LOCF was performed
heatmap_52_emsci <- ggplot(emsci_11recover_motor_52_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_52_emsci
ggsave(paste0(path_figures_save, 'emsci_heatmap_52_PR-MS.pdf'), heatmap_52_emsci)

# Difference between week 1 and week 52

names(emsci_psmall_motor_52_all_long)[names(emsci_psmall_motor_52_all_long) == 'value'] <- 'value_recovery'
emsci_psmall_motor <- merge(emsci_psmall_motor_52_all_long, 
                               emsci_psmall_motor_04_all_long)
emsci_psmall_motor <- emsci_psmall_motor %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(emsci_psmall_motor[,'value_recovery'][[1]])) - 
  as.numeric(as.character(emsci_psmall_motor[,4][[1]]))
emsci_psmall_motor = cbind(emsci_psmall_motor, difference=diff)


group.colors.diff <- c('0' = '#d8d8d8', '1' = '#EDF8E9', 
                  '2' = '#BAE4B3', '3' = '#74C476', 
                  '4' = "#31A354", '5' = "#006D2C", 
                  '-1' = "#FEE5D9", '-2' = "#FCAE91",
                  '-3' = "#FB6A4A", '-4' = "#DE2D26",
                  '-5' = "#A50F15")

heatmap_diff_emsci <- ggplot(emsci_psmall_motor, aes(index, Level, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))
heatmap_diff_emsci

#ggsave(paste0(path_figures_save, 'emsci_heatmap_04-52_PR-MS.pdf'), heatmap_diff_emsci)

# emsci_11recover_motor %>% group_by(index) %>% sum(difference)
# temp_sum_reco <- aggregate(emsci_11recover_motor$difference, by=list(Category=emsci_11recover_motor$index), FUN=sum)
# temp_sum_reco[order(temp_sum_reco$x),]


## PP
emsci_11recover_PP_01_all_long$index <- factor(
  emsci_11recover_PP_01_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_PP_04_all_long$index <- factor(
  emsci_11recover_PP_04_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_PP_52_all_long$index <- factor(
  emsci_11recover_PP_52_all_long$index,
  levels = levels_recovery_subgp_emsci)

heatmap_01_emsci_PP <- ggplot(emsci_11recover_PP_01_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_01_emsci_PP
ggsave(paste0(path_figures_save, 'emsci_heatmap_01_PR-PP.png'), heatmap_01_emsci_PP)

heatmap_04_emsci_PP <- ggplot(emsci_11recover_PP_04_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_04_emsci_PP
ggsave(paste0(path_figures_save, 'emsci_heatmap_04_PR-PP.png'), heatmap_04_emsci_PP)

heatmap_52_emsci_PP <- ggplot(emsci_11recover_PP_52_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_52_emsci_PP
ggsave(paste0(path_figures_save, 'emsci_heatmap_52_PR-PP.png'), heatmap_52_emsci_PP)

# Difference between week 1 and week 52

names(emsci_11recover_PP_52_all_long)[names(emsci_11recover_PP_52_all_long) == 'value'] <- 'value_recovery'
emsci_11recover_PP <- merge(emsci_11recover_PP_52_all_long, 
                               emsci_11recover_PP_01_all_long)
emsci_11recover_PP <- emsci_11recover_PP %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(emsci_11recover_PP[,'value_recovery'][[1]])) - 
  as.numeric(as.character(emsci_11recover_PP[,4][[1]]))
emsci_11recover_PP = cbind(emsci_11recover_PP, difference=diff)

heatmap_diff_emsci_PP <- ggplot(emsci_11recover_PP, aes(index, Level, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))
heatmap_diff_emsci_PP

ggsave(paste0(path_figures_save, 'emsci_heatmap_01-52_PR-PP.png'), heatmap_diff_emsci_PP)


## LT
emsci_11recover_LT_01_all_long$index <- factor(
  emsci_11recover_LT_01_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_LT_04_all_long$index <- factor(
  emsci_11recover_LT_04_all_long$index,
  levels = levels_recovery_subgp_emsci)
emsci_11recover_LT_52_all_long$index <- factor(
  emsci_11recover_LT_52_all_long$index,
  levels = levels_recovery_subgp_emsci)

heatmap_01_emsci_LT <- ggplot(emsci_11recover_LT_01_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_01_emsci_LT
ggsave(paste0(path_figures_save, 'emsci_heatmap_01_PR-LT.png'), heatmap_01_emsci_LT)

heatmap_04_emsci_LT <- ggplot(emsci_11recover_LT_04_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_04_emsci_LT
ggsave(paste0(path_figures_save, 'emsci_heatmap_04_PR-LT.png'), heatmap_04_emsci_LT)

heatmap_52_emsci_LT <- ggplot(emsci_11recover_LT_52_all_long, aes(index, Level, fill= value)) + 
  geom_tile() +
  scale_fill_manual(values=group.colors) +
  geom_text(aes(label = value))
heatmap_52_emsci_LT
ggsave(paste0(path_figures_save, 'emsci_heatmap_52_PR-LT.png'), heatmap_52_emsci_LT)

# Difference between week 1 and week 52

names(emsci_11recover_LT_52_all_long)[names(emsci_11recover_LT_52_all_long) == 'value'] <- 'value_recovery'
emsci_11recover_LT <- merge(emsci_11recover_LT_52_all_long, 
                            emsci_11recover_LT_01_all_long)
emsci_11recover_LT <- emsci_11recover_LT %>% 
  replace_with_na(replace = list(value = '9'))
diff = as.numeric(as.character(emsci_11recover_LT[,'value_recovery'][[1]])) - 
  as.numeric(as.character(emsci_11recover_LT[,4][[1]]))
emsci_11recover_LT = cbind(emsci_11recover_LT, difference=diff)

heatmap_diff_emsci_LT <- ggplot(emsci_11recover_LT, aes(index, Level, fill=factor(difference))) + 
  geom_tile() +
  scale_fill_manual(values=group.colors.diff) +
  geom_text(aes(label = difference))
heatmap_diff_emsci_LT

ggsave(paste0(path_figures_save, 'emsci_heatmap_01-52_PR-LT.png'), heatmap_diff_emsci_LT)

