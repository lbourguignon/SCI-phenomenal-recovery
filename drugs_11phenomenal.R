################################################################################
# SCI - Outstanding recovery - Investigate drugs in recovery group
# L. Bourguignon
# First version : 25.05.2023
# Last update : 25.05.2023
################################################################################

################################################################################
# Library
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(naniar)
library(forcats)
library(dplyr)

################################################################################
# Load data
sygen_11recover <- read.csv("/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/shiny/sygen_11recover.csv")
raw_sygen_with_surgery <- read.csv('/Volumes/INANL/45_Lucie_Bourguignon/1_Data/1_Sygen/sygen_with_surgery.csv')
df_window.1.4 <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_window.1.4.csv')
df_window.1.4$rec_AIS_motor_imp2
ID_norecovry <- df_window.1.4$new_PTID[df_window.1.4$rec_AIS_motor_imp == 'No recovery']

# MPSS information
sygen_11recover$I_MPM # number of hours between injury and first MPSS
sygen_11recover$BOLDOSE1 #bolus
sygen_11recover$BOLDOSE2 #bolus
sygen_11recover$INFDOSE1 #infusion
sygen_11recover$TOTBOLIN # total infusion dose including bolus in mg

################################################################################
## Find the drugs prescribed for each patient using 
## "grep -H -r 'new_PTID' ./ | cut -d: -f1" command
################################################################################


################################################################################
pat1 <- c("codeine", "fentanyl", "ceftriaxone", "potassium_chloride", "glycopyrrolate", # (surgery),
"acetaminophen_codeine", "acetaminophen", "lorazepam", "bismuth_subsalicylate",
"vecuronium bromide", "propofol", "lidocaine", "midazolam", "neostigmine",
"trimethobenzamide", "pancuronium_bromide", "vancomycin", "sucralfate", "ephedrine",
"diazepam", "hydroxyzine", "metoclopramide")
################################################################################


################################################################################
pat2 <- c("heparin", "warfarin", "famotidine", "ceftriaxone", "potassium_chloride", "acetaminophen",
"fludrocortisone", "potassium", "sulfamethoxazole_trimethoprim", "nitrofurantoin",
"ephedrine", "diazepam", "ciprofloxacin")
################################################################################


################################################################################
pat3 <- c("ketorolac", "acetaminophen", "gentamicin", "lorazepam", "furosemide", "cefotaxime",
"albumin", "vitamin_k", "morphine", "rbc", "cefazolin", "sulfamethoxazole_trimethoprim",
"ranitidine", "vancomycin", "aminophylline", "metoclopramide")

# morphine: 2.4 per day for first 5 days (day 0 to day 4)
################################################################################


################################################################################
pat4 <- c("acetaminophen_oxycodone", "pentobarbital", "famotidine", "ketorolac", "acetaminophen", 
"alprazolam", "midazolam", "zolpidem", "morphine", "cyproheptadine", "bisacodyl",
"cefazolin", "ibuprofen", "sucralfate", "ciprofloxacin")

# morphine: day 2 to 6, 96 per day
################################################################################


################################################################################
pat5 <- c("heparin", "ceftazidime", "nafcillin", "enoxaparin", "sodium_polystyrene_sulfonate",
"fresh_frozen_plasma", "calcium gluceptate", "potassium_chloride", "glycopyrrolate", # (pre-medication)
"insulin", "mannitol", "gentamicin", "dexamethasone", "cefotaxime", "midazolam",
"naloxone", "albumin", "vitamin_k", "morphine", "calcium_carbonate", "bisacodyl", "rbc",
"magnesium_sulfate", "ranitidine", "pancuronium_bromide", "vancomycin", "ciprofloxacin",
"metoclopramide")

# morphine: day 20 to day 24, 19.2 per day
################################################################################


################################################################################
pat6 <- c("heparin", "ceftazidime", "famotidine", "promethazine", "piperacillin", "gentamicin",
"propofol", "vitamin_b1", "midazolam", "albumin", "vitamin_k", "rbc",
"sulfamethoxazole_trimethoprim", "omeprazole", "sodium_phosphates",
"pancuronium_bromide", "vancomycin", "sucralfate", "aminophylline", "ciprofloxacin",
"methylprednisolone_sodium_succinate")
################################################################################


################################################################################
pat7 <- c("metronidazole", "enoxaparin", "clindamycin", "famotidine", "saline", "ketorolac",
"promethazine", "ceftriaxone", "chloral_hydrate", "acetaminophen", "gentamicin",
"ticarcillin", "amitriptyline", "lorazepam", "vecuronium bromide", "midazolam",
"zolpidem", "acetaminophen_hydrocodone", "cefazolin", "sulfamethoxazole_trimethoprim",
"vancomycin", "ofloxacin", "ibuprofen", "oxybutynin", "metoclopramide")
################################################################################


################################################################################
pat8 <- c("heparin", "ampicillin", "famotidine", "ceftriaxone", "acetaminophen", "gentamicin",
"amoxicillin", "amoxicillin_and_clavulanate_potassium", "morphine", "pseudoephedrine",
"vancomycin", "tinzaparin_sodium")

# morphine: day 0 to day 5 - 11, 11, 11, 16.5, 15, 9
################################################################################


################################################################################
pat9 <- c("heparin", "enoxaparin", "clindamycin", "famotidine", "potassium_chloride",
"glycopyrrolate", "piperacillin", "sodium_chloride", "acetaminophen_codeine",
"acetaminophen", "diphenhydramine", "gentamicin", "fludrocortisone", "lorazepam",
"furosemide", "midazolam", "morphine", "rbc", "magnesium_sulfate", "clonidine",
"vancomycin", "nitrofurantoin", "hydroxyzine", "epinephrine_sulfate", "ciprofloxacin",
"metoclopramide")

# morphine: day 3 to day 13, 72 per day
################################################################################


################################################################################
pat10 <- c("fentanyl", "influenza_vaccine", "heparin", "potassium_phosphate", "ceftazidime",
"metronidazole", "nafcillin", "enoxaparin", "clindamycin", "ketorolac", "erythromycin",
"ceftriaxone", "acetaminophen_codeine", "diphenhydramine", "gentamicin",
"vecuronium bromide", "vitamin_b1", "midazolam", "vitamin_k", "cefazolin", "ranitidine",
"vancomycin", "ofloxacin", "ibuprofen")
################################################################################


################################################################################
pat11 <- c("fentanyl", "heparin", "potassium_phosphate", "potassium_chloride", "acetaminophen",
"diphenhydramine", "gentamicin", "lorazepam", "cisapride", "midazolam", "morphine", "rbc",
"pseudoephedrine", "magnesium_sulfate", "atropine", "ranitidine", "pancuronium_bromide",
"vancomycin", "thiopental", "ciprofloxacin", "cefazolin")

# morphine: day 1 to day 6, 22.5 per day
################################################################################

x <- c(1:11)
y <- unique(c(pat1, pat2, pat3, pat4, pat5, pat6, pat7, pat8, pat9, pat10, pat11))
data <- expand.grid(X=x, Y=y)
data$value <- 0
list_drug_pat <- list(pat1, pat2, pat3, pat4, pat5, pat6, pat7, pat8, pat9, pat10, pat11)

for (i in x){
  for (j in y){
    list_temp <- list_drug_pat[i]
    if (j %in% list_temp[[1]]){
      data$value[data$X == i & data$Y == j] <- 1
    }
  }
}

order_drugs <- data %>% dplyr::group_by(Y) %>% 
  dplyr::summarise(freq = sum(value)) %>% 
  dplyr::arrange(desc(freq))

data$Y <- factor(data$Y, levels = order_drugs$Y)

order_patients <- data %>% dplyr::group_by(X) %>% 
  dplyr::summarise(freq_X = sum(value)) %>% 
  dplyr::arrange(desc(freq_X))

data$X <- factor(data$X, levels = order_patients$X)

####### Heatmaps #######
# ALL DRUGS PRESCRIBED TO 11 PHENOMENAL RECOVERIES
ggplot(data, aes(factor(X), Y, fill = factor(value))) + 
  geom_tile(color = "white") +
  coord_equal() + 
  scale_fill_manual(values = c("#dedede", "#072b89")) +
  labs(x = 'Patient number', y = 'Drug prescribed') +
  guides(fill=guide_legend(title="Drug prescribed\n1=yes, 0=no"))
  
# DRUGS PRESCRIBED TO MAJORITY OF 11 PHENOMENAL RECOVERIES (>=6)
data %>%
  group_by(Y) %>%
  filter(sum(value == 1) >= 6) %>% 
  ggplot(aes(factor(X), Y, fill = factor(value))) + 
  geom_tile(color = "white") +
  coord_equal() + 
  scale_fill_manual(values = c("#dedede", "#072b89")) +
  labs(x = 'Patient number', y = 'Drug prescribed') +
  guides(fill=guide_legend(title="Drug prescribed\n1=yes, 0=no"))

# DRUGS PRESCRIBED TO AT LEAST 2 OF 11 PHENOMENAL RECOVERIES (>=2)
data %>%
  group_by(Y) %>%
  filter(sum(value == 1) >= 2) %>% 
  ggplot(aes(factor(X), Y, fill = factor(value))) + 
  geom_tile(color = "white") +
  coord_equal() + 
  scale_fill_manual(values = c("#dedede", "#072b89")) +
  labs(x = 'Patient number', y = 'Drug prescribed') +
  guides(fill=guide_legend(title="Drug prescribed\n1=yes, 0=no"))

setwd('/Volumes/INANL/9_Catherine_Jutzeler/1_Studies/1_Drugs/Dose_time_per_Drug')
vancomycin_master <- read.csv('./vancomycin/vancomycin_wide_format_unique_ID.csv')
acetaminophen_master <- read.csv('./acetaminophen/acetaminophen_wide_format_unique_ID.csv')
midazolam_master <- read.csv('./midazolam/midazolam_wide_format_unique_ID.csv')
gentamicin_master <- read.csv('./gentamicin/gentamicin_wide_uniqueID.csv')
heparin_master <- read.csv('./heparin/heparin_wide_uniqueID.csv')
famotidine_master <- read.csv('./famotidine/famotidine_wide_format_unique_ID.csv')
ciprofloxacin_master <- read.csv('./ciprofloxacin/ciprofloxacin_wide_format_unique_ID.csv')
morphine_master <- read.csv('./morphine/morphine_wide_unique_ID.csv')

vancomycin_master %>% filter(NEW_ID == sygen_11recover$new_PTID[1])

list_master_files <- list(vancomycin_master, acetaminophen_master, midazolam_master,
                     gentamicin_master, heparin_master, famotidine_master,
                     ciprofloxacin_master, morphine_master)
list_names <- levels(data$Y)[1:8]

list_drug_11 <- c()
list_ID_11 <- c()
list_sum_11 <- c()

list_drug_norecovery <- c()
list_ID_norecovery <- c()
list_sum_norecovery <- c()

for (drug_count in c(1:length(list_master_files))){
  drug_name <- list_names[drug_count]
  drug <- list_master_files[drug_count][[1]]
  #drug <- dplyr::select(drug, -c(X))
  drug_30 <- select(drug, c(2:which(colnames(drug)=="X30"))) %>% replace(is.na(.), 0)
  drug_sub_11 <- drug_30 %>% dplyr::filter(NEW_ID %in% sygen_11recover$new_PTID)
  drug_sub_norecovery <- drug_30 %>% dplyr::filter(NEW_ID %in% ID_norecovry)
  
  drug_sub_11_sum <- drug_sub_11[,-c(1)] %>% 
    rowwise() %>% 
    mutate(rowsum = sum(c_across(everything()), na.rm = T))
  
  drug_sub_norecovery_sum <- drug_sub_norecovery[,-c(1)] %>% 
    rowwise() %>% 
    mutate(rowsum = sum(c_across(everything()), na.rm = T))
  
  for (i in c(1:dim(drug_sub_11_sum)[1])){
    list_drug_11 <- c(list_drug_11, drug_name)
    list_ID_11 <- c(list_ID_11, drug_sub_11$NEW_ID[i])
    list_sum_11 <- c(list_sum_11, drug_sub_11_sum$rowsum[i])
  }
  
  for (i in c(1:dim(drug_sub_norecovery_sum)[1])){
    list_drug_norecovery <- c(list_drug_norecovery, drug_name)
    list_ID_norecovery <- c(list_ID_norecovery, drug_sub_norecovery$NEW_ID[i])
    list_sum_norecovery <- c(list_sum_norecovery, drug_sub_norecovery_sum$rowsum[i])
  }
}

df_11_values <- data.frame(list_ID_11, list_drug_11, list_sum_11)
colnames(df_11_values) <- c('ID', 'drug', 'total_dose')

df_norecovery_values <- data.frame(list_ID_norecovery, list_drug_norecovery, list_sum_norecovery)
colnames(df_norecovery_values) <- c('ID', 'drug', 'total_dose')

df_11_values %>%
  ggplot(aes(factor(ID), drug, fill = total_dose)) + 
  geom_tile(color = "white") +
  coord_equal() + 
  #scale_colour_gradient(low = "#dedede", high = "#072b89") +
  labs(x = 'Patient number', y = 'Drug prescribed') +
  scale_fill_gradientn(colours = c("white", "#dedede", "#072b89"), 
                       values = c(0, 0.1, max(df_11_values$total_dose)))

df_11_values_complete_binary <- df_11_values %>% 
  tidyr::complete(ID, drug) %>% 
  replace(is.na(.), 0) %>%
  mutate(binary30 = ifelse(total_dose > 0, 1, 0))

# DRUGS PRESCRIBED TO MAJORITY OF 11 PHENOMENAL RECOVERIES (>=6)
df_11_values_complete_binary %>%
  ggplot(aes(drug, factor(ID), fill = factor(binary30))) + 
  geom_tile(color = "white") +
  coord_equal() + 
  scale_fill_manual(values = c("#dedede", "#072b89")) +
  labs(y = 'Patient number', x = 'Drug prescribed') +
  guides(fill=guide_legend(title="Drug prescribed\n1=yes, 0=no"))

df_11_values_complete_binary$ID <- plyr::mapvalues(df_11_values_complete_binary$ID, 
                                                   from=c("109702-ABCD", "129701-ABC", 
                                                          "19106-ABCDEFG", "19702-AB",
                                                          "199702-ABCDEFGHI",
                                                          "199709-ABCDEFGHI",
                                                          "239107-ABCDEF", "7104-ABC",
                                                          "7111-ABCDEF", "99109-ABCDEFGH",
                                                          "99704-ABCDEF"), 
                                                   to=c(1:11))

library(ggsci)
library(ggnewscale)
ggplot() + 
  geom_tile(
    data = df_11_values %>% filter(drug=="vancomycin") %>% droplevels, 
    aes(factor(ID), drug, fill = total_dose)
  ) + 
  scale_fill_distiller(palette = "Greens") + 
  new_scale_fill() + 
  geom_tile(
    data = df_11_values %>% filter(drug=="morphine") %>% droplevels, 
    aes(factor(ID), drug, fill = total_dose)
  ) + 
  scale_fill_distiller(palette ="Oranges") +
  new_scale_fill() + 
  geom_tile(
    data = df_11_values %>% filter(drug=="midazolam") %>% droplevels, 
    aes(factor(ID), drug, fill = total_dose)
  ) + 
  scale_fill_distiller(palette ="BuPu") +
  new_scale_fill() + 
  geom_tile(
    data = df_11_values %>% filter(drug=="heparin") %>% droplevels, 
    aes(factor(ID), drug, fill = total_dose)
  ) + 
  scale_colour_brewer(palette ="Blues") +
  theme(legend.position="bottom")


# library("heatmaply")
data_wide <- data %>%
  group_by(Y) %>%
  filter(sum(value == 1) >= 6) %>% 
  tidyr::spread(Y, value)
heatmaply(select(data_wide, -c(X)))
# 
# data_wide <- data %>% tidyr::spread(Y, value)
# data_wide.dist <- dist(data_wide, method='euclidean')
# data_wide.hc <- hclust(as.dist(data_wide.dist))
# 
# data_wide.dist.long <- data_wide.dist %>% as.matrix %>% melt %>%
#   mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
# 
# library(vegan)
# data(varespec)
# vare.dist <- vegdist(varespec)

ggplot(df_window.1.4) +
  geom_histogram(aes(x = I_MPM, y = after_stat(density), group=rec_AIS_motor_imp), 
                 binwidth = 0.25, colour= "black", fill = "white") +
  geom_density(fill="blue", alpha = .2, group = rec_AIS_motor_imp) 

p <- ggplot(df_window.1.4, aes(x=rec_AIS_motor_imp, y=I_MPM)) + 
  geom_violin(trim=FALSE) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth = 0.2)+
  labs(x = 'Recovery group', y = 'Time between injury and first MP dose (hours)')


test_network <- df_11_values_complete_binary %>% 
  group_by(drug) %>%
  mutate(nb.pat = sum(binary30)) %>% 
  select(-c("total_dose", "binary30", "ID")) 

nb.patient.per.drug <- head(test_network, 8)

nb.patient.per.drug %>% 
  distinct(Source) %>%
  dplyr::rename(label = Source)

