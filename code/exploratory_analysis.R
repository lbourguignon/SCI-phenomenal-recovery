################################################################################
# SCI - Outstanding recovery - Definition
# L. Bourguignon
# First version : 07.02.2023
# Last update : 10.02.2023
################################################################################
# Load data

setwd('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/outstanding-recovery-SCI')
raw_data_Sygen <- read.csv('/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')
raw_data_Murnau <- read.csv('/Volumes/green_groups_bds_public/Data/Murnau/HematologicalBiomark_DATA_2021-01-07_0729.csv')
raw_data_EMSCI <- read.csv('/Volumes/green_groups_bds_public/Data/EMSCI/emsci_data_2020.csv')

################################################################################
# Libraries

library(networkD3)
library("survminer")
library("survival")
library(caret)

################################################################################
# Functions

transition_fct <- function(data_test, source, target){
  transition <- data_test %>%
    dplyr::group_by_at(source) %>%
    dplyr::count(!!!syms(target))
  transition <- as.data.frame(transition)
  names(transition) <- c('source', 'target', 'value')
  return(transition)
}

################################################################################

table(raw_data_Sygen$ais1)

true_A <- raw_data_Sygen %>%
  filter(ais1 == 'AIS A', 
         ais4 == 'AIS A')

true_A_copy <- true_A

# -------------------------
# Sankey plot : transition in AIS grades from week 01 to week 52 with 
# -------------------------
# Renaming the variable values such that they include information on the week when it was collected
true_A_copy$ais1 <- paste(true_A_copy$ais1, "_1", sep="")
true_A_copy$ais4 <- paste(true_A_copy$ais4, "_4", sep="")
true_A_copy$ais8 <- paste(true_A_copy$ais8, "_8", sep="")
true_A_copy$ais16 <- paste(true_A_copy$ais16, "_16", sep="")
true_A_copy$ais26 <- paste(true_A_copy$ais26, "_26", sep="")
true_A_copy$ais52 <- paste(true_A_copy$ais52, "_52", sep="")

# Create a table with :
# source (i.e. AIS grade at week 1),
# target (i.e. AIS grade at week 52)
# value (i.e. number of patients go from one category to another from week 1 to 52)
transition1_4 <- transition_fct(true_A_copy, 'ais1', 'ais4')
transition4_8 <- transition_fct(true_A_copy, 'ais4', 'ais8')
transition8_16 <- transition_fct(true_A_copy, 'ais8', 'ais16')
transition16_26 <- transition_fct(true_A_copy, 'ais16', 'ais26')
transition26_52 <- transition_fct(true_A_copy, 'ais26', 'ais52')

transitions_df <- rbind(transition1_4, transition4_8, transition8_16, transition16_26, transition26_52)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(transitions_df$source),
         as.character(transitions_df$target)) %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
transitions_df$IDsource <- match(transitions_df$source, nodes$name)-1
transitions_df$IDtarget <- match(transitions_df$target, nodes$name)-1
# Make the Network
sankey_transition_sygen <- sankeyNetwork(Links = transitions_df, Nodes = nodes,
                                        Source = "IDsource", Target = "IDtarget",
                                        Value = "value", NodeID = "name",
                                        sinksRight=T, iterations = 0, fontSize = 16)%>% 
  htmlwidgets::prependContent(htmltools::tags$text("AIS grade transitions Sygen"))
sankey_transition_sygen


################################################################################

true_A_murnau <- raw_data_Murnau %>%
  filter(va_ais == 'A', 
         ai_ais == 'A')

true_A_murnau_copy <- true_A_murnau

# -------------------------
# Sankey plot : transition in AIS grades from week 01 to week 52 with 
# -------------------------
# Renaming the variable values such that they include information on the week when it was collected
true_A_murnau_copy$va_ais <- paste(true_A_murnau_copy$va_ais, "_2", sep="")
true_A_murnau_copy$ai_ais <- paste(true_A_murnau_copy$ai_ais, "_4", sep="")
true_A_murnau_copy$aii_aiis <- paste(true_A_murnau_copy$aii_aiis, "_12", sep="")
true_A_murnau_copy$aiii_aiiis <- paste(true_A_murnau_copy$aiii_aiiis, "_26", sep="")
true_A_murnau_copy$c_ais <- paste(true_A_murnau_copy$c_ais, "_52", sep="")

# Create a table with :
# source (i.e. AIS grade at week 1),
# target (i.e. AIS grade at week 52)
# value (i.e. number of patients go from one category to another from week 1 to 52)
transition2_4_murnau <- transition_fct(true_A_murnau_copy, 'va_ais', 'ai_ais')
transition4_12_murnau <- transition_fct(true_A_murnau_copy, 'ai_ais', 'aii_aiis')
transition12_26_murnau <- transition_fct(true_A_murnau_copy, 'aii_aiis', 'aiii_aiiis')
transition26_52_murnau <- transition_fct(true_A_murnau_copy, 'aiii_aiiis', 'c_ais')

transitions_df_murnau <- rbind(transition2_4_murnau, transition4_12_murnau, 
                               transition12_26_murnau, transition26_52_murnau)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_murnau <- data.frame(
  name=c(as.character(transitions_df_murnau$source),
         as.character(transitions_df_murnau$target)) %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
transitions_df_murnau$IDsource <- match(transitions_df_murnau$source, nodes_murnau$name)-1
transitions_df_murnau$IDtarget <- match(transitions_df_murnau$target, nodes_murnau$name)-1
# Make the Network
sankey_transition_murnau <- sankeyNetwork(Links = transitions_df_murnau, Nodes = nodes_murnau,
                                         Source = "IDsource", Target = "IDtarget",
                                         Value = "value", NodeID = "name",
                                         sinksRight=T, iterations = 0, fontSize = 16)%>% 
  htmlwidgets::prependContent(htmltools::tags$text("AIS grade transitions Murnau"))
sankey_transition_murnau

################################################################################

raw_data_EMSCI

A_EMSCI_very_acute <- raw_data_EMSCI %>%
  filter(AIS == 'A',ExamStage == 'very acute')
A_EMSCI_acute1 <- raw_data_EMSCI %>%
  filter(AIS == 'A', ExamStage == 'acute I')

ID_trueA_EMSCI <- intersect(A_EMSCI_very_acute[['Patientennummer']], 
                            A_EMSCI_acute1[['Patientennummer']])

true_A_EMSCI <- subset(raw_data_EMSCI, Patientennummer %in% ID_trueA_EMSCI)
dim(true_A_EMSCI %>%
  filter(AIS == 'NT', ExamStage == 'chronic'))

true_A_EMSCI_subset <- true_A_EMSCI %>%
  select(Patientennummer, AIS, ExamStage_weeks)

true_A_EMSCI_subset$AIS_week <- paste0(true_A_EMSCI_subset$AIS, '_', true_A_EMSCI_subset$ExamStage_weeks)

true_A_EMSCI_subset_wide <- reshape(true_A_EMSCI_subset, idvar = "Patientennummer", timevar = "ExamStage_weeks", direction = "wide")

# Create a table with :
# source (i.e. AIS grade at week 1),
# target (i.e. AIS grade at week 52)
# value (i.e. number of patients go from one category to another from week 1 to 52)
transition2_4_EMSCI <- transition_fct(true_A_EMSCI_subset_wide, 'AIS_week.2', 'AIS_week.4')
transition4_12_EMSCI <- transition_fct(true_A_EMSCI_subset_wide, 'AIS_week.4', 'AIS_week.12')
transition12_26_EMSCI <- transition_fct(true_A_EMSCI_subset_wide, 'AIS_week.12', 'AIS_week.26')
transition26_52_EMSCI <- transition_fct(true_A_EMSCI_subset_wide, 'AIS_week.26', 'AIS_week.52')

transitions_df_EMSCI <- rbind(transition2_4_EMSCI, transition4_12_EMSCI, 
                               transition12_26_EMSCI, transition26_52_EMSCI)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_EMSCI <- data.frame(
  name=c(as.character(transitions_df_EMSCI$source),
         as.character(transitions_df_EMSCI$target)) %>% unique()
)
# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
transitions_df_EMSCI$IDsource <- match(transitions_df_EMSCI$source, nodes_EMSCI$name)-1
transitions_df_EMSCI$IDtarget <- match(transitions_df_EMSCI$target, nodes_EMSCI$name)-1
# Make the Network
sankey_transition_EMSCI <- sankeyNetwork(Links = transitions_df_EMSCI, Nodes = nodes_EMSCI,
                                          Source = "IDsource", Target = "IDtarget",
                                          Value = "value", NodeID = "name",
                                          sinksRight=T, iterations = 0, fontSize = 16)%>% 
  htmlwidgets::prependContent(htmltools::tags$text("AIS grade transitions EMSCI"))
sankey_transition_EMSCI

################################################################################

