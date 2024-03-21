################################################################################
# SCI - Outstanding recovery - Sygen description
# L. Bourguignon
# First version : 10.02.2023
# Last update : 23.02.2023
################################################################################

library(table1)
library(yardstick)
library(plyr)
library(epitools)
library(lubridate)

setwd('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/outstanding-recovery-SCI')
raw_data_Sygen <- read.csv('/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')

table(raw_data_Sygen$ais1)

true_A <- raw_data_Sygen %>%
  filter(ais1 == 'AIS A', 
         ais4 == 'AIS A')

true_A_copy <- true_A

# Description of the miraculous recoveries in Sygen

true_A_recover <- true_A %>%
  filter(ais52 %in% c('AIS C', 'AIS D'))

# Look at the characteristics of the true recoveries
true_A_recover %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52)
true_A_recover %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52)
true_A_recover %>%
  select(ptid, vaccd00, vaccd01, vaccd04, vaccd08, vaccd16, vaccd26, vaccd52)
true_A_recover %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52)
true_A_recover %>%
  select(ptid, anyana01, anyana04, anyana08, anyana16, anyana26, anyana52)
true_A_recover %>%
  select(ptid, vaccd00, vaccd01, vaccd04, vaccd08, vaccd16, vaccd26, vaccd52)

# Time-to-event analysis, with time of conversion being the event of interest
# Question: When is the recovery happening?
head(true_A_recover)
for (i in c(1:dim(true_A_recover)[1])){
  for (time in c("8", "16", "26", "52")){
    col = paste0('ais', time)
    if (!(true_A_recover[i, col] == 'AIS A' | true_A_recover[i, col] == 'ND')){
      true_A_recover[i, 'time_conversion'] = time
      break
    }
  }
}

true_A_recover$time_conversion <- as.numeric(true_A_recover$time_conversion)
true_A_recover %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, time_conversion)
true_A_recover$status <- 1
fit_recovery_sygen <- survfit(Surv(time_conversion, status) ~ 1, data = true_A_recover)
ggsurvplot(fit_recovery_sygen, data = true_A_recover)
# Answer: the recovery is most often happening either between week 4 and 8 (n=9)
# or between week 26 and 52 (n=6)

# Is there a difference between treatment groups in terms of the timing of recovery?
fit_recovery_sygen_diff <- survdiff(Surv(time_conversion, status) ~ tx1_r, data = true_A_recover)
# P = 0.9

# Is there a difference between the treatment groups in terms of the likelihood of the event?
true_A$time_conversion <- 52
true_A$status <- 1
for (i in c(1:dim(true_A)[1])){
  if ((true_A[i, 'ais52'] == 'AIS C' | true_A[i, 'ais52'] == 'AIS D')){
    for (time in c("8", "16", "26", "52")){
      col = paste0('ais', time)
      if (!(true_A[i, col] == 'AIS A' | true_A[i, col] == 'ND')){
        true_A[i, 'time_conversion'] = time
        true_A[i, 'status'] = 2
        break
      }
    }
  }
}
true_A$time_conversion <- as.numeric(true_A$time_conversion)
fit_recovery_sygen_allA <- survfit(Surv(time_conversion, status) ~ tx1_r, data = true_A)
ggsurvplot(fit_recovery_sygen_allA, data = true_A)
fit_recovery_sygen_diff_allA <- survdiff(Surv(time_conversion, status) ~ tx1_r, data = true_A)
fit_recovery_sygen_diff_allA
#p=0.2, could not fit a significant difference in terms of the conversion based on the treatment group
#caution: that might also be due to the small sample size in the recovery group (n=19 in total,
# when defining recovery as conversion from AIS A at both week 2 and 4 to AIS C or D at week 52)

# Plot the individual LEMS recovery trajectories
# Are there different groups that form?
true_A_recover_lems <- true_A_recover %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52)
true_A_recover_lems_long <- gather(true_A_recover_lems, time, lems, lower00:lower52, factor_key=TRUE)
ggplot(data = true_A_recover_lems_long, aes(x = time, y = lems, group = ptid, colour= ptid)) + 
  geom_point(position=position_dodge(width=0.2))+
  geom_line(position=position_dodge(width=0.2))
# No clear groups

# Defining subgroups considering other variables than AIS grades
# Among the AIS conversions, which ones had a significant improvement in motor scores?
# "significant improvement" defined as a 5-point increase according to PMID: 23486305
# https://pubmed.ncbi.nlm.nih.gov/23486305/

true_A_recover$max_early_lems <- with(true_A_recover, pmax(lower00, lower01, lower04, na.rm = TRUE))
true_A_recover$max_early_uems <- with(true_A_recover, pmax(upper00, upper01, upper04, na.rm = TRUE))
true_A_recover$max_late_lems <- with(true_A_recover, pmax(lower26, lower52, na.rm = TRUE))
true_A_recover$max_late_uems <- with(true_A_recover, pmax(upper26, upper52, na.rm = TRUE))
for (i in c(1:dim(true_A_recover)[1])){
  if (is.na(true_A_recover$lower52[i])){
    value_lems = true_A_recover$lower26[i]
  } else {
    value_lems = true_A_recover$lower52[i]
  }
  if (value_lems - true_A_recover$max_early_lems[i] > 5){
    true_A_recover$bin_lems_improv[i] <- 1
  } else {
    true_A_recover$bin_lems_improv[i] <- 0
  }
  if (is.na(true_A_recover$upper52[i])){
    value_uems = true_A_recover$upper26[i]
  } else {
    value_uems = true_A_recover$upper52[i]
  }
  if (value_uems - true_A_recover$max_early_uems[i] > 5){
    true_A_recover$bin_uems_improv[i] <- 1
  } else {
    true_A_recover$bin_uems_improv[i] <- 0
  }
}

# Explore the relationship between AIS grade conversion, LEMS/UEMS improvement and level of injury
true_A_recover %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52, bin_lems_improv)
true_A_recover %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv)
true_A_recover %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl, anyana52, vaccd52)
true_A_recover %>%
  select(ptid, l2ltl01, l2ltl04, l2ltl52, 
         l2ltr01, l2ltr04, l2ltr52,
         l2ppl01, l2ppl04, l2ppl52, 
         l2ppr01, l2ppr04, l2ppr52)

# Thoracic injuries only
true_A_recover_thoracic <- true_A_recover %>% filter(str_detect(true_A_recover$splvl, "T"))
true_A_recover_thoracic %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)
true_A_recover_thoracic %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv, splvl)
# none of the thoracic injuries had an improvement in UEMS (expected due to ceiling effect)
true_A_recover_thoracic %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52, bin_lems_improv, splvl)
# Among thoracic injuries
# 3 patients did not improve their LEMS, let's explore that further:
true_A_recover_thoracic %>% filter(bin_lems_improv == 0) %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl, time_conversion)
# all converted from week 26 to week 52 from A to C
true_A_recover_thoracic %>% filter(bin_lems_improv == 0) %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52, bin_lems_improv, splvl)
# all had LEMS of 0 throughout
true_A_recover_thoracic %>% filter(bin_lems_improv == 0) %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv, splvl)
# all had UEMS of 50 throughout
true_A_recover_thoracic %>% filter(bin_lems_improv == 0) %>%
  select(ptid, anyana01, anyana04, anyana08, anyana16, anyana26, anyana52, time_conversion, splvl)
# all had no DAP throughout
true_A_recover_thoracic %>% filter(bin_lems_improv == 0) %>%
  select(ptid, vaccd00, vaccd01, vaccd04, vaccd08, vaccd16, vaccd26, vaccd52, time_conversion, splvl)
# all had no VAC until week 26 (included) and recovered VAC from week 26 to 52 --> drives the AIS grade conversion?
true_A_recover%>% filter(ptid == '109101-ABCDE') %>%
  select(ptid, s45ltl01, s45ltl04, s45ltl08, s45ltl16, s45ltl26, s45ltl52, 
         s45ltr01, s45ltr04, s45ltr08, s45ltr16, s45ltr26, s45ltr52,
         s45ppl01, s45ppl04, s45ppl08, s45ppl16, s45ppl26, s45ppl52, 
         s45ppr01, s45ppr04, s45ppr08, s45ppr16, s45ppr26, s45ppr52)
# all had sensory scores of 0 (LT and PP, both sides) throughout


# Cervical injuries only
true_A_recover_cervical <- true_A_recover %>% filter(str_detect(true_A_recover$splvl, "C"))
true_A_recover_cervical %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)
true_A_recover_cervical %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv, splvl)

true_A_recover_cervical$bin_lems_improv <- factor(true_A_recover_cervical$bin_lems_improv)
true_A_recover_cervical$bin_uems_improv <- factor(true_A_recover_cervical$bin_uems_improv)

# Interestingly, all cervical injuries which converted in AIS grade also improvement their motor score
# n = 5 for UEMS improved only
# n = 6 for UEMS and LEMS improved
# n = 3 for LEMS improved only
cm <- confusionMatrix(true_A_recover_cervical$bin_lems_improv, true_A_recover_cervical$bin_uems_improv)
cm <- conf_mat(true_A_recover_cervical, bin_lems_improv, bin_uems_improv)
autoplot(cm, type = "heatmap") +
  xlab('LEMS improvement') +
  ylab("UEMS improvement") +
  scale_x_discrete(labels=c("0" = "No", "1" = "Yes")) +
  scale_y_discrete(labels=c("0" = "No", "1" = "Yes"))

# Cervical injuries with no UEMS recovery (only LEMS recovery)
true_A_recover_cervical %>% filter(bin_uems_improv == 0) %>%
  select(ptid, lower00, lower01, lower04, lower08, lower16, lower26, lower52, bin_lems_improv, splvl)
true_A_recover_cervical %>% filter(bin_uems_improv == 0) %>%
  select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv, splvl)

################################################################################
#Table comparison between the patients that recover versus no recovery among the true As

true_A$ais52_locf <- ifelse(true_A$ais52 == 'ND', true_A$ais26, true_A$ais52)

# # 2 patients added to the recovery group if you perform LOCF for the missing AIS grade at week 52
# true_A <- true_A %>% 
#   mutate(recovery = if_else(ais52_locf %in% c('AIS C', 'AIS D'), "recovery", 'no recovery'))
# table1(~ factor(sexcd) + age + factor(ais52) | recovery, data=true_A)

true_A <- true_A %>% 
  mutate(level = if_else(str_detect(true_A$splvl, "T"), "T", "C"))
true_A$max_early_lems <- with(true_A, pmax(lower00, lower01, lower04, na.rm = TRUE))
true_A$max_early_uems <- with(true_A, pmax(upper00, upper01, upper04, na.rm = TRUE))
true_A$max_late_lems <- with(true_A, pmax(lower26, lower52, na.rm = TRUE))
true_A$max_late_uems <- with(true_A, pmax(upper26, upper52, na.rm = TRUE))
for (i in c(1:dim(true_A)[1])){
  if (is.na(true_A$lower52[i])){
    value_lems = true_A$lower26[i]
  } else {
    value_lems = true_A$lower52[i]
  }
  if (is.na(value_lems) | is.na(true_A$max_early_lems[i])){
    true_A$bin_lems_improv[i] <- NA
  } else if (value_lems - true_A$max_early_lems[i] > 5){
    true_A$bin_lems_improv[i] <- 1
  } else {
    true_A$bin_lems_improv[i] <- 0
  }
  if (is.na(true_A$upper52[i])){
    value_uems = true_A$upper26[i]
  } else {
    value_uems = true_A$upper52[i]
  }
  if (is.na(value_uems) | is.na(true_A$max_early_uems[i])){
    true_A$bin_uems_improv[i] <- NA
  } else if (value_uems - true_A$max_early_uems[i] > 5){
    true_A$bin_uems_improv[i] <- 1
  } else {
    true_A$bin_uems_improv[i] <- 0
  }
}

true_A <- true_A %>%
  mutate(recovery = if_else(ais52 %in% c('AIS C', 'AIS D'), "recovery", 'no recovery'))

true_A$sexcd <- factor(true_A$sexcd, levels=c(1,2), labels=c("Female", "Male"))
true_A$bin_uems_improv <- factor(true_A$bin_uems_improv, levels=c(0, 1), labels=c("No", "Yes"))
true_A$bin_lems_improv <- factor(true_A$bin_lems_improv, levels=c(0, 1), labels=c("No", "Yes"))
true_A$level <- factor(true_A$level, levels=c('C', 'T'), labels=c("cervical", "thoracic"))
true_A$tx1_r <- factor(true_A$tx1_r, levels=c('D1', 'D2', 'P'), labels=c("First dosage scheme", "Second dosage scheme", "Placebo"))
true_A$anyana52 <- factor(true_A$anyana52, levels=c('0', '1'), labels=c("Absent", "Present"))
true_A$vaccd52 <- factor(true_A$vaccd52, levels=c('0', '1'), labels=c("Absent", "Present"))

label(true_A$sexcd) <- "Sex"
label(true_A$age) <- "Age"
label(true_A$ais52) <- "AIS grade at week 52"
label(true_A$level) <- "Level of injury"
label(true_A$bin_lems_improv) <- "Improved >5pts in LEMS?ᵃ"
label(true_A$bin_uems_improv) <- "Improved >5pts in UEMS?ᵃ"
label(true_A$tx1_r) <- "Treatment group"
label(true_A$anyana52) <- "DAP at week 52"
label(true_A$vaccd52) <- "VAC at week 52"

units(true_A$age) <- "years"

footnote <- "ᵃ Comparing the highest scores from week 0, 2 or 4 and the score at week 52 if available or at week 26 otherwise"

table1(~ sexcd + age + ais52 + level + bin_uems_improv + bin_lems_improv| recovery, 
       data = true_A, footnote = footnote)

table1(~ sexcd + age + level + bin_uems_improv + bin_lems_improv + 
         tx1_r + anyana52 + vaccd52 | level*recovery, 
       data = true_A, footnote = footnote)

# Look at patients that don't convert in their AIS grades but have an improvement in their LEMS or UEMS by at least 5 points
true_A %>% 
  filter(recovery == 'no recovery',
         bin_lems_improv == 'Yes') %>% 
  select(ptid, bin_lems_improv, lower00, lower01, lower04, lower08, lower16, lower26, lower52, splvl)

true_A %>% filter(ptid == '271101-ABC') %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)
# true_A%>% filter(ptid == '271101-ABC') %>%
#   select(ptid, upper00, upper01, upper04, upper08, upper16, upper26, upper52, bin_uems_improv, splvl)
true_A%>% filter(ptid == '199709-ABCDEFGHI') %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)
true_A%>% filter(ptid == '1101-AB') %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)


true_A$ais52_locf <- ifelse(true_A$ais52 == 'ND', true_A$ais26, true_A$ais52)
true_A_recover_locf <- true_A %>%
  filter(ais52_locf %in% c('AIS C', 'AIS D'))

true_A_recover_locf %>% select(ptid, ais1, ais4, ais8, ais16, ais26, ais52)


true_A_recover %>%
  select(ptid, l2ltl01, l2ltl04, l2ltl52, 
         l2ltr01, l2ltr04, l2ltr52,
         l2ppl01, l2ppl04, l2ppl52, 
         l2ppr01, l2ppr04, l2ppr52)

levels <- c('l1', 'l2', 'l3', 'l4', 'l5', 's1', 's2', 's3', 's45')
sides <- c('l', 'r')
tests <- c('pp', 'lt')
earlyweeks <- c('01', '04')

patients_to_exclude <- c()
for (pat in c(1:dim(true_A_recover)[1])){
  for (test in tests){
    for (lvl in levels){
      for (side in sides){
        for (week in earlyweeks){
          column = paste0(lvl, test, side, week)
          #print(column)
          if (is.na(true_A_recover[pat, column])){
            print(paste('NA', column, 'for patient', pat))
          } else if (true_A_recover[pat, column] != 0){
            print(paste('Non-0 value:', column, 'for patient', true_A_recover[pat, 'ptid']))
            patients_to_exclude <- c(patients_to_exclude, true_A_recover[pat, 'ptid'])
          }
        }
      }
    }
  }
}
unique(patients_to_exclude)

true_A_recover %>% filter(ptid %in% patients_to_exclude) %>%
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)

patients_withsensory <- c()
for (pat in c(1:dim(true_A)[1])){
  for (test in tests){
    for (lvl in levels){
      for (side in sides){
        for (week in earlyweeks){
          column = paste0(lvl, test, side, week)
          #print(column)
          if (is.na(true_A[pat, column])){
            print(paste('NA', column, 'for patient', pat))
          } else if (true_A[pat, column] != 0){
            print(paste('Non-0 value:', column, 'for patient', true_A[pat, 'ptid']))
            patients_withsensory <- c(patients_withsensory, true_A[pat, 'ptid'])
          }
        }
      }
    }
  }
}
length(unique(patients_withsensory))

true_A_recover %>%
  filter(!ptid %in% patients_withsensory) %>% 
  filter(anyana52 != 0) %>% 
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl, anyana52, vaccd52)

true_A %>%
  filter(recovery == 'no recovery', bin_lems_improv == 'Yes') %>%
  filter(!ptid %in% patients_withsensory) %>% 
  select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl, anyana52, vaccd52)



################################################################################

# NEW DEFINITION OF RECOVERY
# To be defined as a true recovery, a case must:
# - not have any sensory functions in the lower limbs at week 1 and week 4
# AND, EITHER be AIS A at week 1 and 4 AND be AIS C or D at week 52 (or week 26 if missing AIS grade at week 52)
# OR, be AIS A at week 1 and 4 AND have an improvement of at least 5 points in LEMS 
# between the max LEMS in week 1 and 4 and LEMS 52 (or week 26 if week 52 is missing)

true_A <- raw_data_Sygen %>%
  filter(ais1 == 'AIS A', 
         ais4 == 'AIS A')

#Find the maximum value of LEMS recorded in the first weeks
true_A$max_early_lems <- with(true_A, pmax(lower00, lower01, na.rm = TRUE))
#Find the maximum value of UEMS recorded in the first weeks
true_A$max_early_uems <- with(true_A, pmax(upper00, upper01, na.rm = TRUE))
#Find the maximum value of LEMS recorded in the last weeks
true_A$max_late_lems <- with(true_A, pmax(lower26, lower52, na.rm = TRUE))
#Find the maximum value of UEMS recorded in the last weeks
true_A$max_late_uems <- with(true_A, pmax(upper26, upper52, na.rm = TRUE))

levels <- c('l1', 'l2', 'l3', 'l4', 'l5', 's1', 's2', 's3', 's45')
sides <- c('l', 'r')
tests <- c('pp', 'lt')
earlyweeks <- c('01', '04')
true_A$early_sensation_left <- 0

for (i in c(1:dim(true_A)[1])){ # for all patients with AIS A at week 1 and 4
  
  # Define LEMS improvement
  
  if (is.na(true_A$lower52[i])){ # if LEMS at week 52 is missing
    value_lems = true_A$lower26[i] # consider LEMS at week 26
  } else {
    value_lems = true_A$lower52[i] # otherwise consider LEMS at week 52
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
  
  if (is.na(true_A$upper52[i])){ # if UEMS at week 52 is missing
    value_uems = true_A$upper26[i] # consider UEMS at week 26
  } else {
    value_uems = true_A$upper52[i] # otherwise consider LEMS at week 52
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
  
  for (test in tests){
    for (lvl in levels){
      for (side in sides){
        for (week in earlyweeks){
          column = paste0(lvl, test, side, week)
          if (is.na(true_A[i, column])){
            print(paste('NA', column, 'for patient', true_A$ptid[i]))
          } else if (true_A[i, column] != 0) {
            true_A$early_sensation_left[i] <- 1
          }
        }
      }
    }
  }
}

stricter_true_A <- true_A %>%
  filter(early_sensation_left == 0)

dim(stricter_true_A %>% filter(bin_lems_improv == 1))

true_A_recover_new <- true_A %>%
  filter(early_sensation_left == 0 & 
           (ais52 %in% c('AIS C', 'AIS D') | bin_lems_improv == 1) &
           ( bin_uems_improv == 1 | bin_lems_improv == 1))

for (i in c(1:dim(true_A_recover_new)[1])){
  for (time in c("8", "16", "26", "52")){
    col = paste0('ais', time)
    if (!(true_A_recover_new[i, col] == 'AIS A' | true_A_recover_new[i, col] == 'ND')){
      true_A_recover_new[i, 'time_conversion'] = time
      break
    }
  }
}

true_A_recover_new %>%
  dplyr::select(ptid, ais1, ais4, ais8, ais16, ais26, ais52, bin_lems_improv, bin_uems_improv, splvl)

# true_A %>%
#   filter(early_sensation_left == 0 & (bin_lems_improv == 1)) %>%
#   select(ptid, ais1, ais52, vaccd52, anyana52)

dim(true_A_recover_new)
## Ok so we end up with 14 patients, among those I suspect 3 of them to not be true recoveries
## but if we apply the same criteria
## (i.e. absence of DAP and improvement suspicious because only occurring from week 26 to 52)
## then we also exclude the patient that improved LEMS and UEMS but did not convert in AIS grade
## So I will keep them for now but to be reconsidered

true_A_recover_new %>%
  select(ptid, ais1, ais4, ais52, bin_lems_improv, early_sensation_left, splvl, anyana52)
true_A_recover_new %>% 
  select(ptid, ais1, ais52, modben04, modben08, modben16, modben26, modben52, bin_lems_improv, splvl)

true_A_recover_new %>%
  filter(anyana52 == 0) %>%
  select(ptid, ais1, ais4, ais52, bin_lems_improv, bin_uems_improv, splvl,
         upper00, upper01, upper04, upper08, upper16, upper26, upper52, vaccd52)

## Make summary table ##

true_A_copy_new <- stricter_true_A
true_A_copy_new <- true_A_copy_new %>% 
  mutate(level = if_else(str_detect(true_A_copy_new$splvl, "T"), "T", "C"))

id_recoverers <- true_A_recover_new$ptid
true_A_copy_new$recover <- 'No recovery'
true_A_copy_new$recover[true_A_copy_new$ptid %in% id_recoverers] <- "Recovery"

table1(~ factor(sexcd) + age + factor(ais52) + factor(level) + factor(bin_uems_improv) |  factor(bin_lems_improv), 
       data = true_A_copy_new)

## Make stats comparisons ##

# Age comparison
cdat <- plyr::ddply(true_A_copy_new, "recover", summarise, age.median=median(age))
ggplot(true_A_copy_new, aes(x=age, fill=recover)) + geom_density(alpha=.3)+
  geom_vline(data=cdat, aes(xintercept=age.median,  colour=recover),
             linetype="dashed", size=1)

wilcox.test(age ~ recover, data = true_A_copy_new)

# LEMS comparison
wilcox.test(lower01 ~ recover, data = true_A_copy_new)
wilcox.test(lower04 ~ recover, data = true_A_copy_new)
wilcox.test(lower08 ~ recover, data = true_A_copy_new)
wilcox.test(lower16 ~ recover, data = true_A_copy_new)
wilcox.test(lower26 ~ recover, data = true_A_copy_new)
wilcox.test(lower52 ~ recover, data = true_A_copy_new)

true_A_recover_lems_new <- true_A_copy_new %>%
  select(ptid, recover, early_sensation_left, lower00, lower01, lower04, lower08, lower16, lower26, lower52)
true_A_recover_lems_new_long <- gather(true_A_recover_lems_new, time, lems, lower00:lower52, factor_key=TRUE)
ggplot(data = true_A_recover_lems_new_long, aes(x = time, y = lems, group = ptid, colour = factor(early_sensation_left))) + 
  geom_point(position=position_dodge(width=0.2))+
  geom_line(position=position_dodge(width=0.2))
ggplot(data = true_A_recover_lems_new_long, aes(x = time, y = lems, group = ptid, colour = factor(recover))) + 
  geom_point(position=position_dodge(width=0.2))+
  geom_line(position=position_dodge(width=0.2))

# UEMS comparison
wilcox.test(upper01 ~ recover, data = true_A_copy_new)
wilcox.test(upper04 ~ recover, data = true_A_copy_new)
wilcox.test(upper08 ~ recover, data = true_A_copy_new)
wilcox.test(upper16 ~ recover, data = true_A_copy_new)
wilcox.test(upper26 ~ recover, data = true_A_copy_new)
wilcox.test(upper52 ~ recover, data = true_A_copy_new)

# true_A_recover_uems_new <- true_A_copy_new %>%
#   filter (level == 'C') %>%
#   select(ptid, recover, early_sensation_left, upper00, upper01, upper04, upper08, upper16, upper26, upper52)
# true_A_recover_uems_new_long <- gather(true_A_recover_uems_new, time, uems, upper00:upper52, factor_key=TRUE)
# ggplot(data = true_A_recover_uems_new_long, aes(x = time, y = uems, group = ptid, colour = factor(early_sensation_left))) + 
#   geom_point(position=position_dodge(width=0.2))+
#   geom_line(position=position_dodge(width=0.2))
# ggplot(data = true_A_recover_uems_new_long, aes(x = time, y = uems, group = ptid, colour = factor(recover))) + 
#   geom_point(position=position_dodge(width=0.2))+
#   geom_line(position=position_dodge(width=0.2))

# Sex comparison
table(true_A_copy_new$recover, true_A_copy_new$sexcd)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$sexcd))

# AIS grade at week 52
table(true_A_copy_new$recover, true_A_copy_new$ais52)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$ais52))

# LEMS improvement
table(true_A_copy_new$recover, true_A_copy_new$bin_lems_improv)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$bin_lems_improv))
# statistically significant association between the improvement in LEMS and recovery 

# LEMS improvement and sensory at baseline
table(true_A_copy_new$early_sensation_left, true_A_copy_new$bin_lems_improv)
fisher.test(table(true_A_copy_new$early_sensation_left, true_A_copy_new$bin_lems_improv))
# statistically significant association between the improvement in LEMS and presence of sensation at early time point

# LEMS improvement and sensory at baseline
table(true_A_copy_new$early_sensation_left, true_A_copy_new$bin_lems_improv)
fisher.test(table(true_A_copy_new$early_sensation_left, true_A_copy_new$bin_lems_improv))
# statistically significant association between the improvement in LEMS and presence of sensation at early time point

# recovery and sensory at baseline
table(true_A_copy_new$recover, true_A_copy_new$early_sensation_left)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$early_sensation_left))
# statistically significant association between the improvement in LEMS and presence of sensation at early time point

# Associated injuries
true_A_copy_new$nb_associated_injuries <- true_A_copy_new$cardcd01 +
  true_A_copy_new$cpulcd01 +
  true_A_copy_new$eentcd01 +
  true_A_copy_new$gastcd01 +
  true_A_copy_new$genicd01 +
  true_A_copy_new$headcd01 + 
  true_A_copy_new$musccd01

true_A_copy_new <- true_A_copy_new %>% 
  mutate(nb_associated_injuries_bin = if_else(nb_associated_injuries > 0, 1, 0))

table(true_A_copy_new$recover, true_A_copy_new$nb_associated_injuries)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$nb_associated_injuries_bin))

wilcox.test(nb_associated_injuries ~ recover, data = true_A_copy_new)

inj_surg
wilcox.test(surgno ~ recover, data = true_A_copy_new)

# Treatment group
table(true_A_copy_new$recover, true_A_copy_new$tx1_r)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$tx1_r))

#Cause of injury
table(true_A_copy_new$recover, true_A_copy_new$injcd)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$injcd))

# Date of injury
x <- true_A_copy_new$injdt
date <- dmy(x)
days <- yday(date)
total_days <- days
true_A_copy_new$dayofyear <- total_days
wilcox.test(dayofyear ~ recover, data = true_A_copy_new)
ggplot(true_A_copy_new, aes(x=dayofyear, fill=recover)) + geom_density(alpha=.3)

true_A_copy_new %>% 
  filter(recover == 'Recovery') %>% 
  select(dayofyear, injcd)

# Number of days to first spinal injury
true_A_copy_new$diff_days_firstsurgery <- as.Date(true_A_copy_new$surgdt1, format="%d/%m/%Y") - 
  as.Date(true_A_copy_new$injdt, format="%d/%m/%Y")
cdat_surgery <- plyr::ddply(true_A_copy_new, "recover", summarise, 
                    days.surgery.median = median(as.numeric(diff_days_firstsurgery), na.rm = TRUE))
ggplot(true_A_copy_new, aes(x = as.numeric(diff_days_firstsurgery), fill = recover)) + geom_density(alpha = .3) +
  geom_vline(data = cdat_surgery, aes(xintercept = days.surgery.median,  colour = recover),
             linetype = "dashed", size = 1)

true_A_copy_new %>% 
  filter(recover == 'Recovery') %>% 
  select(injdt, surgdt1, diff_days_firstsurgery, age)
wilcox.test(as.numeric(diff_days_firstsurgery) ~ recover, data = true_A_copy_new)

# Infection adverse event
true_A_copy_new %>% 
  #filter(recover == 'Recovery') %>% 
  select(aeons_immune02_1, aeons_immune03_1, aeons_immune04_1, aeons_immune08_1)

# Level of injury
table(true_A_copy_new$recover, true_A_copy_new$level)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$level))

# DAP
table(true_A_copy_new$recover, true_A_copy_new$anyana52)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$anyana52))

# VAC
table(true_A_copy_new$recover, true_A_copy_new$vaccd52)
fisher.test(table(true_A_copy_new$recover, true_A_copy_new$vaccd52))

# Look at subset that doesn't recover
temp_norecovery <- true_A_copy_new %>% 
  filter(recover == 'No recovery') 
temp_norecovery 
# Level of injury
table(temp_norecovery$level)

# Number of patients truly severe with LEMS improvement
dim(true_A_copy_new %>% filter(bin_lems_improv == 1))

# Look at the time to first spinal surgery again + make plot
true_A_copy_new$diff_days_firstsurgery <- as.Date(true_A_copy_new$surgdt1, format="%d/%m/%Y") - 
  as.Date(true_A_copy_new$injdt, format="%d/%m/%Y")
cdat_surgery <- plyr::ddply(true_A_copy_new, "recover", summarise, 
                            days.surgery.median = median(as.numeric(diff_days_firstsurgery), na.rm = TRUE))
ggplot(true_A_copy_new, aes(x = as.numeric(diff_days_firstsurgery), fill = recover)) + geom_density(alpha = .3) +
  geom_vline(data = cdat_surgery, aes(xintercept = days.surgery.median,  colour = recover),
             linetype = "dashed", size = 1)

true_A_copy_new %>% filter(recover == "Recovery") %>%
  select(diff_days_firstsurgery)

ggplot(true_A_copy_new, aes(x = as.numeric(diff_days_firstsurgery), y = lower52, colour = recover)) + 
  geom_point() +
  geom_vline(data = cdat_surgery, aes(xintercept = days.surgery.median,  colour = recover),
             linetype = "dashed", size = 1)

true_A_copy_new %>%  dplyr::group_by(recover) %>% summarise(n=n(), avg = median(diff_days_firstsurgery)) -> Summary.data

p <- ggplot(true_A_copy_new, aes(x=recover, y=as.numeric(diff_days_firstsurgery))) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.05) +
  geom_jitter(shape=16, position=position_jitter(0.05)) +
  theme_classic() +
  geom_text(data=Summary.data ,aes(x = recover, y = 380, label=paste0('n=',n)), color="red", fontface =2, size = 5) + 
  labs(x = 'Recovery status', y = 'Number of days between injury and first spinal surgery')
p

# Look at quantiles of age for no recovery
temp_norecover <- true_A_copy_new %>% filter(recover == 'No recovery') %>%
  select(age)
quantile(temp_norecover$age)

