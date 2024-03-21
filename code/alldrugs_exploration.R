################################################################################
# SCI - Outstanding recovery - Investigate drugs in recovery group
# L. Bourguignon
# First version : 20.06.2023
# Last update : 20.06.2023
# INFECTIONS
################################################################################

################################################################################
# Library
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(naniar)
library(forcats)
library(dplyr)
library(plyr)
library(stringr)

################################################################################
# Load data
raw_sygen <- read.csv(
  '/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')
df_window.1.4 <- read.csv(
  '/Volumes/green_groups_bds_public/Projects/Drug/df_window.1.4.csv')
ID_norecovry <- df_window.1.4$new_PTID[
  df_window.1.4$rec_AIS_motor_imp == 'No recovery']
ID_recovery <- df_window.1.4$new_PTID[
  df_window.1.4$rec_AIS_motor_imp == 'Recovery']


setwd('/Volumes/green_groups_bds_public/Projects/Drug/Drugs')
acetaminophen_master <- read.csv('./acetaminophen.csv')
vancomycin_master <- read.csv('./vancomycin.csv')
morphine_master <- read.csv('./morphine.csv')
midazolam_master <- read.csv('./midazolam.csv')
gentamicin_master <- read.csv('./gentamicin.csv')
lorazepam_master <- read.csv('./lorazepam.csv')
potassium_chloride_master <- read.csv('./potassium_chloride.csv')
heparin_master <- read.csv('./heparin.csv')
metoclopramide_master <- read.csv('./metoclopramide.csv')
ciprofloxacin_master <- read.csv('./ciprofloxacin.csv')
famotidine_master <- read.csv('./famotidine.csv')
albuterol_master <- read.csv('./albuterol.csv')

################################################################################
# Exploration of the infection information in the original sygen data

raw_sygen$aeons_infection02_1[raw_sygen$ptid %in% ID_norecovry]
raw_sygen$aeons_infection02_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection02_2[raw_sygen$ptid %in% ID_recovery]

raw_sygen$aeons_infection03_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection03_2[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection04_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection04_2[raw_sygen$ptid %in% ID_recovery]

raw_sygen$aeli1_infection02_1[raw_sygen$ptid %in% ID_norecovry]

raw_sygen$aeli1_infection02_1[raw_sygen$ptid %in% ID_norecovry]
raw_sygen$aeli1_infection02_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeli1_infection03_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeli1_infection04_1[raw_sygen$ptid %in% ID_recovery]


raw_sygen$aeons_infection02_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection03_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection04_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aeons_infection04_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aedur_infection03_1[raw_sygen$ptid %in% ID_recovery]
raw_sygen$aedur_infection04_1[raw_sygen$ptid %in% ID_recovery]

names(vancomycin_master)

data_heatmap_test <- vancomycin_master[vancomycin_master$NEW_ID %in% ID_recovery,]

data_heatmap_test_small <- data_heatmap_test[data_heatmap_test$NEW_ID == "19106-ABCDEFG",]
data_heatmap_test_small <- select(data_heatmap_test_small, c(4:which(colnames(data_heatmap_test_small)=="X30"))) %>% replace(is.na(.), 0)
data_heatmap_test_small <- as.matrix(data_heatmap_test_small)

# Default Heatmap
heatmap(data_heatmap_test_small)

################################################################################
# Combine all drugs prescribed into 1 file

# # New IDs after change of PR definition
# temp <- df_window.1.4[!is.na(df_window.1.4$bin_lems_improv_5),]
# ID_recovery_sygen_new_clinic <- temp$new_PTID[temp$bin_lems_improv_5 == 1]
# ID_norecovery_sygen_new_clinic <- temp$new_PTID[temp$bin_lems_improv_5 == 0]
# 
# # List of all drug files
# files <- list.files(path="/Volumes/green_groups_bds_public/Projects/Drug/Drugs",
#                     pattern="*.csv", full.names=TRUE, recursive=FALSE)
# 
# # Create empty dataframes for phenomenal recovery (PR)
# # and no phenomenal recovery (NPR) groups separately
# df.alldrugs.recovery <- data.frame()
# df.alldrugs.norecovery <- data.frame()
# 
# for (i in c(1:length(files))){ # go through all files
#   temp <- read.csv(files[i]) # read file i
#   temp <- temp[,!names(temp) %in% c("X")] # remove index column
#   temp_no_recovery <- temp[temp$NEW_ID %in% ID_recovery_sygen_new_clinic,] # select for NPR group
#   temp_recovery <- temp[temp$NEW_ID %in% ID_norecovery_sygen_new_clinic,] # select for PR group
#   # combine file i to the overall dataframes
#   df.alldrugs.recovery <- rbind(df.alldrugs.recovery, temp_recovery)
#   df.alldrugs.norecovery <- rbind(df.alldrugs.norecovery, temp_no_recovery)
# }
# 
# # Save dataframes
# write.csv(df.alldrugs.recovery, '/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_recovery_newdef.csv', row.names=FALSE)
# write.csv(df.alldrugs.norecovery, '/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_norecovery_newdef.csv', row.names=FALSE)
# 
# df.alldrugs.norecovery_copy <- df.alldrugs.norecovery
# df.alldrugs.recovery_copy <- df.alldrugs.recovery

################################################################################
# Preprocessing

df.alldrugs.norecovery_copy <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_norecovery_newdef.csv')
df.alldrugs.recovery_copy <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_recovery_newdef.csv')

df.alldrugs.norecovery_old <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_norecovery.csv')
df.alldrugs.recovery_old <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_recovery.csv')

# Convert NA to 0
df.alldrugs.norecovery_copy[is.na(df.alldrugs.norecovery_copy)] <- 0
df.alldrugs.recovery_copy[is.na(df.alldrugs.recovery_copy)] <- 0

# Combine values such that there is only 1 line per patient per drug
df.alldrugs.norecovery_copy <- plyr::ddply(df.alldrugs.norecovery_copy, c("NEW_ID", "generic_name"), plyr::numcolwise(sum))
df.alldrugs.recovery_copy <- plyr::ddply(df.alldrugs.recovery_copy, c("NEW_ID", "generic_name"), plyr::numcolwise(sum))


################################################################################
# Find common and specific drugs from the NPR and PR groups

# Select for first 30 days
df.alldrugs.norecovery_copy30 <- df.alldrugs.norecovery_copy %>% select('NEW_ID':'X30')
df.alldrugs.recovery_copy30 <- df.alldrugs.recovery_copy %>% select('NEW_ID':'X30')

# Find drugs prescribed in both groups
common_drugs_old <- intersect(unique(df.alldrugs.norecovery_old$generic_name),
                          unique(df.alldrugs.recovery_old$generic_name))
common_drugs <- intersect(unique(df.alldrugs.norecovery_copy30$generic_name),
          unique(df.alldrugs.recovery_copy30$generic_name))
#intersect(common_drugs_old, common_drugs)
## 129 drugs were prescribed in at least 1 patient in each group
### new def: 95

# Find drugs prescribed only in recovery group
setdiff(unique(df.alldrugs.recovery_copy30$generic_name),
        unique(df.alldrugs.norecovery_copy30$generic_name))
## 5 drugs present in the recovery group but not in the non recovery group
## hydrocortisone_acetate - corticosteroid
dim(df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% 'hydrocortisone_acetate', ])[1]
## sodium_citrate
dim(df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% 'sodium_citrate', ])[1]
## salicylic_acid
dim(df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% 'salicylic_acid', ])[1]
## cromolyn_sodium
dim(df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% 'cromolyn_sodium', ])[1]
## flumazenil
dim(df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% 'flumazenil', ])[1]
#### All of the 5 drugs where prescribed to only 1 patient each

# Find drugs drugs prescribed only in non-recovery group
specific_norecovery <- setdiff(unique(df.alldrugs.norecovery_copy30$generic_name),
        unique(df.alldrugs.recovery_copy30$generic_name))

df.alldrugs.norecovery_copy_specific_norecovery <- 
  df.alldrugs.norecovery_copy30[df.alldrugs.norecovery_copy30$generic_name %in% specific_norecovery, ]
table(table(df.alldrugs.norecovery_copy_specific_norecovery$generic_name))

tb_norecovery_specific <- table(df.alldrugs.norecovery_copy_specific_norecovery$generic_name)
tb_norecovery_specific[tb_norecovery_specific >=16]
## 261 drugs were given in the no recovery group but not in the recovery group
## among those 261 drugs, 95 drugs were given to only 1 patient each
## 9 drugs were given to at least 20 unique patients each
### ampicillin_and_sulbactam, penicillin + beta lactamase inhibitor --> broad spectrum
### atropine, 
### cimetidine, reduce stomach acids
### fluconazole, antifungal
### guaifenesin, anti congestant (clear mucus)
### phytonadione, = vit K, prevents bleeding
### triazolam, benzodiazepine (spasticity)

################################################################################
# Barplot visualing the prevalence of drugs prescribed both in the recovery and non recovery group

# Select only for drugs prescribed in both PR and NPR groups
df.commondrugs.norecovery <- df.alldrugs.norecovery_copy30[df.alldrugs.norecovery_copy30$generic_name %in% common_drugs, ]
df.commondrugs.recovery <- df.alldrugs.recovery_copy30[df.alldrugs.recovery_copy30$generic_name %in% common_drugs, ]

# Compute prevalence of each drug per group
df.table.norecovery <- data.frame(table(df.commondrugs.norecovery$generic_name)/length(unique(df.commondrugs.norecovery$NEW_ID)))
df.table.recovery <- data.frame(table(df.commondrugs.recovery$generic_name)/length(unique(df.commondrugs.recovery$NEW_ID)))

# Combine 2 groups into 1 dataframe
df.table.norecovery$group <- 'no recovery'
df.table.recovery$group <- 'recovery'
df.table.all <- rbind(df.table.norecovery, df.table.recovery)

# Find drugs more prevalent in the PR group compared to NPR
df.table.all$increase_amount <-  with(df.table.all, 
                                      as.integer(Freq > 
                                                   ave(Freq * (group == "recovery"), Var1, FUN = function(x) x[x!=0])))
# Store the name of those drugs
names_drugs_more_recovery <- df.table.all$Var1[df.table.all$increase_amount == 1]

# Split data based on prevalence higher in PR or NPR group
df.table.all.more_recovery <- df.table.all[df.table.all$Var1 %in% names_drugs_more_recovery, ]
df.table.all.more_norecovery <- df.table.all[!(df.table.all$Var1 %in% names_drugs_more_recovery), ]

# Reorder drug names based on prevlance values
df.table.all.more_recovery$Var1 <- factor(df.table.all.more_recovery$Var1)
df.table.all.more_recovery$Var1 <- fct_reorder(df.table.all.more_recovery$Var1, 
                                               df.table.all.more_recovery$Freq, 
                                               max)

df.table.all.more_norecovery$Var1 <- factor(df.table.all.more_norecovery$Var1)
df.table.all.more_norecovery$Var1 <- fct_reorder(df.table.all.more_norecovery$Var1, 
                                                 df.table.all.more_norecovery$Freq, 
                                                 max)

# Reorder group such that the most prevalent appears as top most or left most bar in plot
df.table.all.more_recovery$group <- factor(df.table.all.more_recovery$group, 
                                           levels = c("no recovery", "recovery"))
df.table.all.more_norecovery$group <- factor(df.table.all.more_norecovery$group, 
                                             levels = c("recovery", "no recovery"))

# Flag drugs prescribed to more than 1 patient in the PR group
names_drugs_maj_recovery <- df.table.all$Var1[df.table.all$group == 'recovery' &
                                                df.table.all$Freq > 0.1]
df.table.all.more_recovery_top <- df.table.all.more_recovery[df.table.all.more_recovery$Var1 
                                                             %in% names_drugs_maj_recovery, ]

# Plot
ggplot(data=df.table.all.more_recovery, aes(x=Var1, y=Freq, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() +
  coord_flip() +
  geom_col(width = 0.5, position = "dodge")+
  scale_fill_manual(values=c("#d4a7a5", "#800020"))

# ideas:
  ## split plot into drugs with increased prevalence in non recovery in 1 panel
  ## and drugs with increased prevalence in recovery group in the other panel
  ## order drugs by prevalence
  ## remove drug prescribed to only 1 patient in each group (recovery and non recovery)
  ## identify groups of drugs (e.g., all antibiotics)


################################################################################
## Antibiotics
################################################################################
# Description of antibio based on
# https://www.orthobullets.com/basic-science/9059/antibiotic-classification-and-mechanism
# https://clerkship.medicine.ufl.edu/files/2012/04/Antibiotic-Overview.pdf (source 1)
# https://www.researchgate.net/figure/Classification-of-antibiotics-based-on-spectrum-of-activity_tbl2_309592434 (source 2)



common_antibio_more_norecovery <- c('cefazolin', #1st generation cephalosporin                              --> narrow spectrum
                                    'trimethoprim', #folic acid synthesis inhibitor                         --> moderately broad spectrum
                                    'ceftazidime', #3rd generation cephalosporin                            --> broad spectrum
                                    'ticarcillin', #penicillin (antipseudomonal)                            --> moderately broad spectrum
                                    'tobramycin', #aminoglycoside protein synthesis inhibition              --> narrow spectrum
                                    'ampicillin') #penicillin (aminopenicillin)                             --> moderately broad spectrum

common_antibio_more_recovery <- c('acetic_acid', #                                                          --> 
                                  'amoxicillin', #penicillin (aminopenicillin)                              --> 
                                  'amoxicillin_and_clavulanate_potassium', #                                --> broad spectrum
                                  'bacitracin', #cell wall inhibitor (topical)                              --> broad spectrum (source 2)
                                  'cefotaxime', #3rd generation cephalosporin                               --> broad spectrum
                                  'ceftriaxone', #3rd generation cephalosporin                              --> broad spectrum (source 2)
                                  'ciprofloxacin', #2nd generation fluoroquinolone (DNA synthesis inhibitor)--> moderately broad spectrum
                                  'clindamycin', #macrolide (lincosamide,protein syntehsis inhibition)      --> broad spectrum
                                  'gentamicin', #aminoglycoside (protein synthesis inhibition)              --> narrow spectrum (source 2)
                                  'metronidazole', #DNA inhibitor                                           --> narrow spectrum
                                  'nafcillin', #penicillin (penicillinase-resistant-penicillin)             --> narrow spectrum 
                                  'nitrofurantoin', #                                                       --> 
                                  'ofloxacin', #2nd generation fluoroquinolone (DNA synthesis inhibitor)    --> 
                                  'piperacillin', #penicillin (antipseudomonal)                             --> moderately broad spectrum
                                  'vancomycin') # cell wall inhibitor                                       --> narrow spectrum

other_infec_more_norecovery <- c('nystatin')
other_infec_more_recovery <- c('acyclovir', 'amphotericin_B', 
                               'betamethasone_and_clotrimazole','clotrimazole',
                               'erythromycin', 'miconazole')


antibio_norecovery_only <- c('ampicillin_and_sulbactam', 'cefuroxime', 'penicillin',
                            'imipenem_and_cilastatin', 'imipenem', 'cephalexin',
                            'norfloxacin', 'aztreonam', 'chloramphenicol',
                            'cefotetan', 'tetracycline', 'cephradine', 'rifampin',
                            'oxacillin', 'benzoyl_peroxide', 'doxycycline',
                            'ceftizoxime', 'cephapirin', 'cefaclor', 'amikacin',
                            'metacyclin', 'meropenem', 'mupirocin', 'cefadroxil',
                            'cefamandole', 'clarithromycin', 'piperacillin_and_tazobactam',
                            'dicloxacillin', 'ticarcillin_clavulanate',
                            'sulfacetamide', 'mezlocillin', 'carbenicillin',
                            'azithromycin', 'cefpodoxime', 'mafenide', 'cefoperazone',
                            'cefoxitin')

all_antibio <- c(common_antibio_more_norecovery, 
                 common_antibio_more_recovery,
                 antibio_norecovery_only)

length(unique(sample(all_antibio)))

df.table.all.antibio <- df.table.all %>%
  dplyr::filter(Var1 %in% all_antibio)

#write.csv(as.data.frame(sample(all_antibio)),"/Users/blucie/Desktop/git-repos/outstanding-recovery-SCI/outputs/list-antibiotics-phenomenal-recovery-SCI.csv")

################################################################################
# Plot prevalence antibiotics

# PLot prevalence of antibio more prevalent in PR group
df.table.all.more_recovery %>%
  dplyr::filter(Var1 %in% all_antibio) %>%
  ggplot(aes(x=Var1, y=Freq, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() +
  coord_flip() +
  geom_col(width = 0.5, position = "dodge")+
  scale_fill_manual(values=c("#627313", "#C0C7A1"))

# PLot prevalence of antibio more prevalent in NPR group
df.table.all.more_norecovery %>%
  dplyr::filter(Var1 %in% common_antibio_more_norecovery) %>%
  ggplot(aes(x=Var1, y=Freq, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() +
  coord_flip() +
  geom_col(width = 0.5, position = "dodge")+
  scale_fill_manual(values=c("#800020", "#d4a7a5"))

# Plot prevalence of antibio prescribed only in NPR group
df.table.norecovery.specific <- data.frame(table(
  df.alldrugs.norecovery_copy_specific_norecovery$generic_name)/
    length(unique(df.alldrugs.norecovery_copy_specific_norecovery$NEW_ID)))
df.table.norecovery.specific$Var1 <- factor(df.table.norecovery.specific$Var1)
df.table.norecovery.specific$Var1 <- fct_reorder(df.table.norecovery.specific$Var1, 
                                                 df.table.norecovery.specific$Freq, 
                                               max)
df.table.norecovery.specific %>%
  filter(Var1 %in% antibio_norecovery_only) %>%
  ggplot(aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal() +
  coord_flip()

################################################################################
# Plot number of unique antibiotics prescribed per patient and per day
# in the first 30 days after injury

## Data preparation for PR group
# Select for antibiotics only
df.antibio.recovery <- df.alldrugs.recovery_copy %>% 
  filter(generic_name %in% all_antibio)
# Select for first 30 days
df.antibio.recovery.30days <- df.antibio.recovery %>% select('NEW_ID':'X30')
new_id_temp <- df.antibio.recovery.30days$NEW_ID
generic_name_temp <- df.antibio.recovery.30days$generic_name
# Binarise the entries (convert any value above 0 to 1)
df.antibio.recovery.30days[df.antibio.recovery.30days > 0] <- 1
df.antibio.recovery.30days$NEW_ID <- new_id_temp
df.antibio.recovery.30days$generic_name <- generic_name_temp
# Sum the number of antibiotics per patient per day
df.antibio.recovery.30days.sum <- df.antibio.recovery.30days %>%
  dplyr::group_by(NEW_ID) %>%
  dplyr::summarise(across(names(df.antibio.recovery.30days)[3:33], sum))
# Convert to long format
df.antibio.recovery.30days.sum.long <- melt(setDT(df.antibio.recovery.30days.sum), 
                                            id.vars = c("NEW_ID"), variable.name = "Day")

## Plot PR group
df.antibio.recovery.30days.sum.long %>%
  ggplot(aes(Day, factor(NEW_ID), fill = factor(value))) + 
  geom_tile(color = "white") +
  geom_text(aes(label=value)) +
  coord_equal() + 
  labs(y = 'Patient number', x = 'Day after injury') +
  scale_fill_manual(values = c(brewer.pal(9, "Blues")))

## Data preparation for NPR group
df.antibio.norecovery <- df.alldrugs.norecovery_copy %>% 
  filter(generic_name %in% all_antibio)
df.antibio.norecovery.30days <- df.antibio.norecovery %>% select('NEW_ID':'X30')
new_id_temp <- df.antibio.norecovery.30days$NEW_ID
generic_name_temp <- df.antibio.norecovery.30days$generic_name
df.antibio.norecovery.30days[df.antibio.norecovery.30days > 0] <- 1
df.antibio.norecovery.30days$NEW_ID <- new_id_temp
df.antibio.norecovery.30days$generic_name <- generic_name_temp
df.antibio.norecovery.30days.sum <- df.antibio.norecovery.30days %>%
  dplyr::group_by(NEW_ID) %>%
  dplyr::summarise(across(names(df.antibio.norecovery.30days)[3:33], sum))
df.antibio.norecovery.30days.sum.long <- melt(setDT(df.antibio.norecovery.30days.sum), 
                                            id.vars = c("NEW_ID"), variable.name = "Day")


## Plot NPR group
# Note: only a subset of the NPR group is plotted for better readibility
# The subset can be modified by changing the dplyr::filter [1:40]
df.antibio.norecovery.30days.sum.long %>%
  dplyr::filter(NEW_ID %in% unique(df.antibio.norecovery.30days.sum.long$NEW_ID)[1:40])  %>%
  ggplot(aes(Day, factor(NEW_ID), fill = factor(value))) + 
  geom_tile(color = "white") +
  geom_text(aes(label=value)) +
  coord_equal() + 
  #scale_colour_gradient(low = "#dedede", high = "#072b89") +
  labs(y = 'Patient number', x = 'Day after injury') +
  scale_fill_manual(values = c(brewer.pal(9, "Blues")))
#scale_fill_gradientn(colours = c("white", "#dedede", "#072b89"))

################################################################################
## Transform the dates of infection to days after injury

subset_infec <- raw_sygen %>% select(ptid, injdt, 
                                     aeons_infection02_1, aeons_infection02_2,
                                     aeons_infection03_1, aeons_infection03_2,
                                     aeons_infection04_1, aeons_infection03_2,
                                     aeons_infection04_3, aeons_infection04_4)

transformDate <- function(x) as.Date(x, format="%d/%m/%Y")
subset_infec_new <- data.frame(lapply(subset_infec, transformDate))
subset_infec_new$ptid <- subset_infec$ptid

DiffDay <- function(x) as.numeric(difftime(x, subset_infec_new$injdt, 
                                           units = "days"))
subset_infec_new_days <- data.frame(lapply(subset_infec_new[,3:9], DiffDay))
subset_infec_new_days$ptid <- subset_infec_new$ptid

subset_infec_new_days_recovery <- subset_infec_new_days[subset_infec_new_days$ptid %in% ID_recovery,]
subset_infec_new_days_norecovery <- subset_infec_new_days[subset_infec_new_days$ptid %in% ID_norecovry,]

nb_antibio_recovery <- as.data.frame(table(df.antibio.recovery.30days$NEW_ID))
nb_antibio_norecovery <- as.data.frame(table(df.antibio.norecovery.30days$NEW_ID))

summary(nb_antibio_recovery)
summary(nb_antibio_norecovery)

sd(nb_antibio_recovery$Freq)
sd(nb_antibio_norecovery$Freq)

table(nb_antibio_recovery$Freq)/length(unique(nb_antibio_recovery$Var1))
table(nb_antibio_norecovery$Freq)/length(unique(nb_antibio_norecovery$Var1))

################################################################################
## Chi square test vancomycin
# Select only for drugs prescribed in both PR and NPR groups
df.vanco.norecovery <- df.alldrugs.norecovery_copy[df.alldrugs.norecovery_copy$generic_name %in% "vancomycin", ] %>% select('X0':'X30')
df.vanco.recovery <- df.alldrugs.recovery_copy[df.alldrugs.recovery_copy$generic_name %in% "vancomycin", ] %>% select('X0':'X30')

df.vanco.norecovery = df.vanco.norecovery[rowSums(df.vanco.norecovery[])>0,]
df.vanco.recovery = df.vanco.recovery[rowSums(df.vanco.recovery[])>0,]

################################################################################
#Number of patients with antibio prescribed on Day 0
## NPR
sum(df.antibio.norecovery.30days$X0) #24
## PR
sum(df.antibio.recovery.30days$X0) #0

id_atb_X0 <- df.antibio.recovery.30days$NEW_ID[df.antibio.recovery.30days$X0==1]
id_atb_X1 <- df.antibio.recovery.30days$NEW_ID[df.antibio.recovery.30days$X1==1]
length(intersect(id_atb_X0, id_atb_X1))

#Number of patients with antibio prescribed on Day 1
## NPR
sum(df.antibio.norecovery.30days$X1) #78
## PR
sum(df.antibio.recovery.30days$X1) #4

################################################################################
#Number of cumulative antibiotics days in first 30 days
## NPR

head(df.antibio.norecovery.30days)
df.antibio.norecovery.30days <- df.antibio.norecovery.30days %>%
  mutate(rowsums = select(., -c(1:2)) %>% 
           rowSums(na.rm = TRUE))
test <- df.antibio.norecovery.30days %>%
  dplyr::group_by(NEW_ID) %>% 
  dplyr::mutate(Total=sum(rowsums))
plot(table(test$Total))
df.cumulative.days.antibio.NPR <- test %>% select(NEW_ID, Total)
df.cumulative.days.antibio.NPR <- df.cumulative.days.antibio.NPR %>% distinct(.keep_all = TRUE)
df.cumulative.days.antibio.NPR$status <- 'NPR'

## PR

head(df.antibio.recovery.30days)
df.antibio.recovery.30days <- df.antibio.recovery.30days %>%
  mutate(rowsums = select(., -c(1:2)) %>% 
           rowSums(na.rm = TRUE))
test <- df.antibio.recovery.30days %>%
  dplyr::group_by(NEW_ID) %>% 
  dplyr::mutate(Total=sum(rowsums))
plot(table(test$Total))
df.cumulative.days.antibio.PR <- test %>% select(NEW_ID, Total)
df.cumulative.days.antibio.PR <- df.cumulative.days.antibio.PR %>% distinct(.keep_all = TRUE)
df.cumulative.days.antibio.PR$status <- 'PR'

df.cumulative.days.antibio.all <- rbind(df.cumulative.days.antibio.NPR[,c('NEW_ID', 'Total', 'status')],
                                        df.cumulative.days.antibio.PR[,c('NEW_ID', 'Total', 'status')])
df.cumulative.days.antibio.all.first <- df.cumulative.days.antibio.all[match(unique(df.cumulative.days.antibio.all$NEW_ID), 
                                                                             df.cumulative.days.antibio.all$NEW_ID),]
df.cumulative.days.antibio.all.first$Total.cat <- cut(df.cumulative.days.antibio.all.first$Total,
                                                   breaks=c(-Inf,50,Inf), labels = c("<=150", ">150"))

table(df.cumulative.days.antibio.all.first$status, df.cumulative.days.antibio.all.first$Total.cat)
test <- chisq.test(table(df.cumulative.days.antibio.all.first$status, df.cumulative.days.antibio.all.first$Total.cat))
test

# None of the comparisons, regardless of the threshold used (based on the
# numbers observed in the PR group) yield statistically significant difference
# between the NPR and PR groups

################################################################################
#Number of days without antibio
df.antibio.recovery.30days.sum %>%
  select(., -c(1:2)) %>% 
  rowSums(df.antibio.recovery.30days.sum == 0)

df.antibio.recovery.30days.sum$days_no_atb <- rowSums(df.antibio.recovery.30days.sum == 0)
plot(table(rowSums(df.antibio.recovery.30days.sum == 0)))
df.antibio.norecovery.30days.sum$days_no_atb <- rowSums(df.antibio.norecovery.30days.sum == 0)
plot(table(rowSums(df.antibio.norecovery.30days.sum == 0)))

################################################################################
#Number of antibiotics in cervical VS thoracic injuries in NPR group

names(df.antibio.recovery.30days)
sygen_norecovery <- raw_sygen %>% filter(ptid %in% ID_norecovry)

sygen_norecovery.nli <- sygen_norecovery %>% select(ptid, splvl)
ID_norecovery_cervical <- sygen_norecovery.nli %>%
  filter(stringr::str_detect(splvl, "C")) %>%
  select(ptid)

names(sygen_norecovery.nli)[names(sygen_norecovery.nli) == 'ptid'] <- 'NEW_ID'
df.antibio.recovery.30days.nli <- merge(sygen_norecovery.nli, 
                                        df.antibio.norecovery.30days.sum,
                                        all.x = T,
                                        by = "NEW_ID")

df.antibio.recovery.30days.nli <- df.antibio.recovery.30days.nli %>%
  mutate(level = case_when(str_detect(splvl, "C") ~ 'cervical',
                           TRUE ~ "thoracic"))

ggplot(df.antibio.recovery.30days.nli,
       aes(x = splvl,
           y = days_no_atb,
           fill = level)) +
  geom_bar(stat = "identity",
           position = "dodge")

df.antibio.recovery.30days.nli%>%
  dplyr::count(days_no_atb, level) %>% 
  dplyr::mutate(freq = prop.table(n), .by = level) %>%
  ggplot(aes(x = days_no_atb,
             y = freq,
             fill = level)) +
  geom_bar(stat = "identity",
           position = "dodge")

df.antibio.recovery.30days.nli <- df.antibio.recovery.30days.nli%>%
  dplyr::count(days_no_atb, level) %>% 
  dplyr::mutate(freq = prop.table(n), .by = level)

ks.test(df.antibio.recovery.30days.nli[df.antibio.recovery.30days.nli$level == 'cervical',]$freq, 
        df.antibio.recovery.30days.nli[df.antibio.recovery.30days.nli$level == 'thoracic',]$freq)

# Seems like thoracic injuries lead to more days with no antibiotics
# --> Cervical injuries lead to more days with antibiotics
# But difference is NOT statistically significant

# --> What about PR? They are cervical but do they have the similar amount of days with no antibiotics?
# Seem hard to say based on the low sample size

################################################################################
# Number of antibiotics related to the overall severity of trauma
# as measured by the number of associated injuries among
# cardiac, cardiopulmonary, ENT, gastrointestinal, genitourinary, head, 
# musculoskeletal, skin and soft tissues

sygen.ass.injuries <- raw_sygen %>% select(ptid, cardcd01, cpulcd01, eentcd01,
                                           gastcd01, genicd01, headcd01,
                                           musccd01, numep01, pulmcd01, skincd01
                                           )
sygen.ass.injuries.PR <- sygen.ass.injuries %>% 
  filter(ptid %in% ID_recovery) %>%
  mutate(rowsums = select(., -c(1)) %>% 
           rowSums(na.rm = TRUE))
table(sygen.ass.injuries.PR$rowsums)/sum(table(sygen.ass.injuries.PR$rowsums))
sygen.ass.injuries.PR$status <- 'PR'

sygen.ass.injuries.NPR <- sygen.ass.injuries %>% 
  filter(ptid %in% ID_norecovry) %>%
  mutate(rowsums = select(., -c(1)) %>% 
           rowSums(na.rm = TRUE))
table(sygen.ass.injuries.NPR$rowsums)/sum(table(sygen.ass.injuries.NPR$rowsums))
sygen.ass.injuries.NPR$status <- 'NPR'

sygen.head.injuries.all <- rbind(sygen.ass.injuries.PR[, c('headcd01', 'status')],
      sygen.ass.injuries.NPR[, c('headcd01', 'status')])

table(sygen.head.injuries.all$status, sygen.head.injuries.all$headcd01)
test <- chisq.test(table(sygen.head.injuries.all$status, sygen.head.injuries.all$headcd01))
test

# Higher prevalence of head injuries in the PR group compared to the NPR group
# but NOT statistically significant
# Chi-squared stats: 4.5301
# pvalue = 0.033304
# Yate's corrected pvalue: 0.07485

# Correlation between the number of unique antibiotics prescribed in the first 
# 30 days after injury and the number of associated injuries
nb.antibio.30days.NPR <- as.data.frame(table(df.antibio.norecovery.30days$NEW_ID))
nb.antibio.30days.PR <- as.data.frame(table(df.antibio.recovery.30days$NEW_ID))
names(nb.antibio.30days.NPR)[names(nb.antibio.30days.NPR) == 'Var1'] <- 'ptid'
names(nb.antibio.30days.PR)[names(nb.antibio.30days.PR) == 'Var1'] <- 'ptid'
names(nb.antibio.30days.NPR)[names(nb.antibio.30days.NPR) == 'Freq'] <- 'nb.antibio'
names(nb.antibio.30days.PR)[names(nb.antibio.30days.PR) == 'Freq'] <- 'nb.antibio'

nb.antibio.30days.all <- rbind(nb.antibio.30days.NPR, nb.antibio.30days.PR)

nb.ass.injuries.NPR <- sygen.ass.injuries.NPR %>% select(ptid, rowsums, status)
names(nb.ass.injuries.NPR)[names(nb.ass.injuries.NPR) == 'rowsums'] <- 'nb.ass.injuries'
nb.ass.injuries.PR <- sygen.ass.injuries.PR %>% select(ptid, rowsums, status)
names(nb.ass.injuries.PR)[names(nb.ass.injuries.PR) == 'rowsums'] <- 'nb.ass.injuries'

nb.ass.injuries.all <- rbind(nb.ass.injuries.NPR, nb.ass.injuries.PR)

df.antibio.ass.injuries <- merge(nb.antibio.30days.all,
                                 nb.ass.injuries.all,
                                 all.x = T,
                                 by = "ptid")

ggplot(df.antibio.ass.injuries, aes(nb.ass.injuries, nb.antibio, col=status)) +
  geom_jitter(width = 0.25)+
  scale_color_manual(values=c("#C0C7A1", "#800020"))

cor.test(df.antibio.ass.injuries$nb.antibio, 
         df.antibio.ass.injuries$nb.ass.injuries, 
         method=c("s"))

# No correlations between the number of associated injuries and the number of
# unique antibiotics prescribed in the first 30 days

################################################################################
# Number of antibiotics in the first 30 days after SCI

# Number of unique antibiotics per patient in the first 30 days
table(df.antibio.recovery.30days$NEW_ID)
mean(table(df.antibio.recovery.30days$NEW_ID))
sd(table(df.antibio.recovery.30days$NEW_ID))
median(table(df.antibio.recovery.30days$NEW_ID))
IQR(table(df.antibio.recovery.30days$NEW_ID))

# Number of unique antibiotics per patient in the first 30 days
table(df.antibio.norecovery.30days$NEW_ID)
mean(table(df.antibio.norecovery.30days$NEW_ID))
sd(table(df.antibio.norecovery.30days$NEW_ID))
median(table(df.antibio.norecovery.30days$NEW_ID))
IQR(table(df.antibio.norecovery.30days$NEW_ID))

# Still seems like slightly more antibiotics were prescribed in the PR group

df.antibio.ass.injuries$nb.antibio.cat <- cut(df.antibio.ass.injuries$nb.antibio,
            breaks=c(-Inf,5,Inf), labels = c("<=5", ">5"))

table(df.antibio.ass.injuries$status, df.antibio.ass.injuries$nb.antibio.cat)
test <- chisq.test(table(df.antibio.ass.injuries$status, df.antibio.ass.injuries$nb.antibio.cat))
test
# No statistically significant difference between the number of antibiotics prescribed
# in the PR compared to NPR group, when considering the cut off of 4, 5 or 6 drugs

#Number of days without antibiotics
df.antibio.recovery.30days.sum$status <- 'PR'
df.antibio.norecovery.30days.sum$status <- 'NPR'
df.noantibio.all.30days.sum <- rbind(df.antibio.recovery.30days.sum[,c('NEW_ID', 'days_no_atb', 'status')],
                                     df.antibio.norecovery.30days.sum[,c('NEW_ID', 'days_no_atb', 'status')])
df.noantibio.all.30days.sum$days_no_atb.cat <- cut(df.noantibio.all.30days.sum$days_no_atb,
                                              breaks=c(-Inf,15,Inf), labels = c("<=15", ">15"))

table(df.noantibio.all.30days.sum$status, df.noantibio.all.30days.sum$days_no_atb.cat)
test <- chisq.test(table(df.noantibio.all.30days.sum$status, df.noantibio.all.30days.sum$days_no_atb.cat))
test
# No statistically significant different in the number of days with no antibio
# prescribed in the PR group compared to NPR group

################################################################################
# Identify prophylaxis prescriptions

# Seems like the files with information about the indications are only available
# on Kip's group share: contacted Catherine to see if she has access to the data
# or if I should arrange a new access to Kip's group share through Bobo
