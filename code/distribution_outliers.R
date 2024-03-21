################################################################################
# SCI - Outstanding recovery - New definition PR only based on LEMS improvement
# L. Bourguignon
# First version : 20.11.2023
# Last update : 20.03.2024
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
library(reshape2)
library(MatchIt)
library(cobalt)
library(FedData)

################################################################################

setwd('/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/miraculous-recovery-SCI')

raw_sygen <- read.csv('/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')
#write.csv(raw_sygen, '/Volumes/blucie/PhD/1_SCI/10_miraculous-recovery/modified_sygen.csv')
raw_data_EMSCI <- read.csv('/Volumes/green_groups_bds_public/Data/EMSCI/emsci_data_2020.csv')

check_NLI2 = function(dataset, patient_identifier){
  # assumes columns are names level + type + side
  # e.g. C5ppl for pin prick score on C5, left hand side
  # C5ltl for light touch score on C5, left hand side
  # C5mol for motor score on C5, left hand side
  # type in c('lt','pp','mo') and side in c('l','r')
  
  # checks consecutive levels for being full score
  
  # dataset in wide format, one entry per patient
  
  # NLI definition used
  # * most caudal (highest) level for which pinprick, light touch are intact
  # AND motor score >= 3 given that all motor scores above are 5 (maximal)
  
  # known bugs
  # 1. Check if this is true still!
  #     if you have sequence TRUE NaN FALSE for full motor scores (=5)
  #     then will use NaN level as motor NLI
  #     note that for TRUE NaN NaN FALSE, will return NaN if FALSE has >= 3 motor score
  #     and TRUE level as motor NLI if FALSE has < 3 motor score
  # 2. If C5 has motor scores between 3 and 5 then use C5 as motor NLI. Not sure if incorrect
  #     as no key motor levels bove that need to be score 5 
  
  ### sensory ----
  sensory_vars_ordered_r_pp = c("c2ppr","c3ppr","c4ppr","c5ppr","c6ppr","c7ppr","c8ppr", 
                                "t1ppr", "t2ppr", "t3ppr", "t4ppr", "t5ppr", "t6ppr", "t7ppr", "t8ppr", "t9ppr", "t10ppr","t11ppr","t12ppr",
                                "l1ppr","l2ppr","l3ppr","l4ppr","l5ppr",
                                "s1ppr","s2ppr","s3ppr",
                                "s45ppr")
  sensory_levels = c('c1', sub('ppr','',sensory_vars_ordered_r_pp))
  sensory_vars_ordered_r_pp = paste0(sensory_vars_ordered_r_pp, '01')
  
  # get per patient and per dermatome, right pinprick scores
  data_sensory_pp_r = dataset %>%
    dplyr::select(all_of(c(sensory_vars_ordered_r_pp, patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'ppr') %>%
    dplyr::mutate(variable = gsub('.{5}$', '', variable)) # remove ppr01
  
  data_sensory_pp_l = dataset %>%
    dplyr::select(all_of(c(sub('ppr','ppl',sensory_vars_ordered_r_pp), patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'ppl') %>%
    dplyr::mutate(variable = gsub('.{5}$', '', variable))
  
  data_sensory_lt_r = dataset %>%
    dplyr::select(all_of(c(sub('ppr','ltr',sensory_vars_ordered_r_pp), patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'ltr') %>%
    dplyr::mutate(variable = gsub('.{5}$', '', variable))
  
  data_sensory_lt_l = dataset %>%
    dplyr::select(all_of(c(sub('ppr','ltl',sensory_vars_ordered_r_pp), patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'ltl') %>%
    dplyr::mutate(variable = gsub('.{5}$', '', variable))
  
  # merge all sensory long formats
  data_sensory = base::merge(data_sensory_pp_r, data_sensory_pp_l, 
                             by = c(patient_identifier, 'variable'))
  data_sensory = base::merge(data_sensory, data_sensory_lt_r, 
                             by = c(patient_identifier, 'variable'))
  data_sensory = base::merge(data_sensory, data_sensory_lt_l, 
                             by = c(patient_identifier, 'variable'))
  
  # check dermatomes that have intact sensation (split by left and right)
  data_sensory = data_sensory %>%
    dplyr::mutate(left_intact = (ppl==2 & ltl==2), # here some might be missing
                  right_intact = (ppr==2 & ltr==2))%>%
    dplyr::mutate(variable = factor(variable, levels = rev(sensory_levels), ordered = TRUE))
  
  # makes dermatomes ordered factors from s45 to c2
  # then checks per patient that no sensory info is missing
  # if both left is intact and right is intact then those dermatomes are kept per patient
  # the lowest (min) of these dermatomes is then the sensory NLI
  
  # find highest level with missing value
  data_sensory_missing = data_sensory %>%
    dplyr::group_by(!!rlang::sym(patient_identifier)) %>%
    dplyr::filter(is.na(left_intact) | is.na(right_intact)) %>% # per patient keep levels with missing value
    dplyr::summarise(highest_missing = max(variable))
  
  # find lowest level with intact sensation given that all above are full
  # EXCLUDING NaN, will later fix this problem
  data_sensory_NLI = data_sensory %>%
    dplyr::group_by(!!rlang::sym(patient_identifier)) %>%
    dplyr::filter(!is.na(left_intact) & !is.na(right_intact)) %>% # should give these an NLI of NaN
    dplyr::mutate(sens_intact = left_intact & right_intact) %>%
    dplyr::arrange(desc(variable)) %>% 
    # https://stackoverflow.com/questions/29273012/find-first-occurence-of-value-in-group-using-dplyr-mutate#29274880
    # first non intact sensory level, sensory level is last intact one (so one higher)
    dplyr::summarise(NLI_sensory_after_last_consec_intact = # ignoring NaN, will consider later 
                       dplyr::case_when( 
                         # exceptions: no FALSE, FALSE in first level (then no previous one)
                         sum(sens_intact == FALSE) > 0  ~
                           max(variable[sens_intact == FALSE]), # highest level with non intact sensation
                         .default = min(variable) ), # if no non-intact then take lowest non missing 
                     NLI_sensory_is_intact = 
                       dplyr::case_when(
                         sum(sens_intact == FALSE) > 0 ~ FALSE,
                         .default = TRUE)) %>%# if all are intact, then take lowest non NaN level
    dplyr::mutate(NLI_sensory_missing = 
                    dplyr::case_when(
                      NLI_sensory_is_intact == TRUE ~ NLI_sensory_after_last_consec_intact, # if was intact
                      .default = rev(sensory_levels)[as.numeric(NLI_sensory_after_last_consec_intact)+1])) %>% # one level higher
    dplyr::mutate(NLI_sensory_missing = factor(NLI_sensory_missing, levels = rev(sensory_levels), ordered = TRUE))
  
  
  data_sensory_NLI = merge(data_sensory_NLI, data_sensory_missing, 
                           by = patient_identifier, all = TRUE) %>%
    dplyr::mutate(NLI_sensory = 
                    dplyr::case_when(
                      is.na(highest_missing) ~ NLI_sensory_missing, # no missing values
                      as.numeric(highest_missing)>as.numeric(NLI_sensory_missing) ~ NA_character_, 
                      # missing values above lowest intact level
                      .default = NLI_sensory_missing), # missing value not above sensory missing
                  missing_higher_NLI = as.numeric(highest_missing)>as.numeric(NLI_sensory_missing)) %>%
    dplyr::mutate(NLI_sensory = factor(NLI_sensory, levels = rev(sensory_levels), ordered = TRUE))
  
  NLI_sens_NaN_ids = unique( (data_sensory %>% 
                                dplyr::filter(is.na(left_intact) | is.na(right_intact)))$ptid)  
  
  
  
  ### motor ----
  
  
  
  motor_vars_l = c('elbfll','wrextl', 'elbexl','finfll','finabl',
                   'hipfll','kneexl','ankdol', 'gretol','ankpll')
  motor_vars_r = c('elbflr','wrextr', 'elbexr','finflr','finabr',
                   'hipflr','kneetr','ankdor', 'gretor','ankplr') # different kneetr from kneexl
  motor_vars_l = paste0(motor_vars_l, '01')
  motor_vars_r = paste0(motor_vars_r, '01')
  
  data_motor_l = dataset %>%
    dplyr::select(all_of(c(motor_vars_l, patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'l') %>%
    dplyr::mutate(variable = gsub('.{3}$', '', variable))
  
  data_motor_r = dataset %>%
    dplyr::select(all_of(c(motor_vars_r, patient_identifier))) %>%
    reshape2::melt(id = patient_identifier, value.name = 'r') %>%
    dplyr::mutate(variable = gsub('.{3}$', '', variable)) %>%
    dplyr::mutate(variable = dplyr::case_when(grepl('kne', variable) ~ 'kneex', # name L and R are different, homogenize
                                       .default = variable))
  
  data_motor = base::merge(data_motor_l, data_motor_r, 
                           by = c(patient_identifier, 'variable'))
  
  data_motor = data_motor %>%
    dplyr::mutate(intact = (l==5 & r==5),
                  motor_3 = ( (l<5 | r<5) & l>=3 & r>=3)) %>%
    dplyr::mutate(variable = factor(variable, 
                                    levels = c('s45','s3','s2','ankpl', 'greto','ankdo','kneex','hipfl',
                                               'l1',paste0('t', seq(12, 2,by=-1)),
                                               'finab','finfl','elbex','wrext','elbfl',
                                               'c4','c3','c2', 'c1'), #rev(gsub('.{3}$', '', motor_vars_l)), 
                                    ordered = TRUE,
                                    labels = rev(sensory_levels))) # as sensory levels have all spine levels
  
  # find lowest level with intact motor score given that all above are full
  # EXCLUDING NaN, will later fix this problem
  data_motor_NLI_intact = data_motor %>%
    dplyr::group_by(!!rlang::sym(patient_identifier)) %>%
    dplyr::filter(!is.na(intact)) %>% # if intact is not NaN then motor_3 is also not, as evaluate same columns
    dplyr::arrange(desc(variable)) %>% 
    # https://stackoverflow.com/questions/29273012/find-first-occurence-of-value-in-group-using-dplyr-mutate#29274880
    # first non intact motor level, motor level is last intact one (so one higher than first non-intact)
    dplyr::summarise(NLI_motor_after_last_consec_intact = # find first FALSE, one level higher is intact
                       # ignoring NaN, will consider later 
                       dplyr::case_when( 
                         # exceptions: no FALSE, FALSE in first level (then no previous one)
                         sum(intact == FALSE) > 0  ~
                           max(variable[intact == FALSE]), # highest level with non intact motor function
                         .default = min(variable) ), # if all non NaN are intact then take lowest level
                     NLI_motor_is_intact = 
                       dplyr::case_when(
                         sum(intact == FALSE) > 0 ~ FALSE,
                         .default = TRUE),# if all are intact, then TRUE and don't take higher level in next step
                     motor_atleast3 = # find first TRUE
                       dplyr::case_when(
                         # exceptions: no FALSE, FALSE in first level (then no previous one)
                         sum(motor_3 == TRUE) > 0  ~
                           rev(sensory_levels)[as.numeric(max(variable[motor_3 == TRUE]))], # highest level with at least score 3 but not 5 on both sides
                         .default = NA_character_ )) %>% # missing if no such level exists
    dplyr::mutate(motor_intact = 
                    dplyr::case_when(
                      NLI_motor_is_intact == TRUE ~ 
                        rev(sensory_levels)[as.numeric(NLI_motor_after_last_consec_intact)], # if all non NaN were intact
                      .default = rev(sensory_levels)[as.numeric(NLI_motor_after_last_consec_intact)+1]),# one level higher
                  # so if L2 is not complete, uses L1, this follows the Note in ISNSCI exam that when full UEMS but LEMS L2 not full then use sensory evaluation to find NLI
                  # last intact is actually t1
                  motor_atleast3 = factor(motor_atleast3, levels = rev(sensory_levels), ordered = TRUE)) %>%
    
    dplyr::mutate(motor_intact = factor(motor_intact, levels = rev(sensory_levels), ordered = TRUE))
  
  print(head(data_motor_NLI_intact))
  
  # print(unique(data_motor_NLI$motor_atleast3))
  # print(unique(data_motor_NLI$motor_intact))
  
  # find highest level with missing value
  data_motor_missing = data_motor %>%
    dplyr::group_by(!!rlang::sym(patient_identifier)) %>%
    dplyr::filter(is.na(intact)) %>% # per patient keep levels with missing value
    dplyr::summarise(highest_missing = max(variable))
  
  # motor NLI is the lowest level that has at least motor score 3 given that all levels above have score 5 (full score)
  # and no missing info on full motor score or not above motor_intact
  data_motor_NLI = merge(data_motor_NLI_intact, data_motor_missing, 
                         by=patient_identifier, all=TRUE) %>%
    
    dplyr::mutate(NLI_motor = dplyr::case_when(
      as.numeric(highest_missing)>as.numeric(motor_intact) ~ NA_character_, # if have a missing value in motor level above intact lowest level 
      as.numeric(motor_intact) == (as.numeric(motor_atleast3) + 1) ~ motor_atleast3,
      .default = motor_intact)) %>%
    
    dplyr::mutate(NLI_motor = factor(NLI_motor, 
                                     levels = rev(sensory_levels), #rev(gsub('.{3}$', '', motor_vars_l)), 
                                     ordered = TRUE,
                                     labels = rev(sensory_levels)))
  
  NLI_mot_NaN_ids = unique( (data_motor %>% 
                               dplyr::filter(is.na(intact), !is.na(motor_3)))$ptid)  
  
  data_NLI_calc = merge(data_motor_NLI, data_sensory_NLI, 
                        by = patient_identifier,
                        all = TRUE) %>% 
    dplyr::mutate(NLI_motor_numeric = as.numeric(NLI_motor), # motor level is one below last intact one and we go from s45 to c2 so 
                  NLI_sensory_numeric = as.numeric(NLI_sensory)) %>%
    dplyr::mutate(NLI_numeric = pmax(NLI_motor_numeric, NLI_sensory_numeric),
                  NLI_calc = rev(sensory_levels)[NLI_numeric]) 
  return(data_NLI_calc)
}

preprocess_more_NLI = function(dataset){
  # function to impute missing values based on logic so that more NLI can be calculated
  # and less patients are discarded
  
  # * c2 sensory info missing but c3 not, then impute c2 with lowest c3, 
  #     this contributed most to not being able to find splvl
  
  out =  dataset %>% 
    dplyr::mutate(c3_sens_not_missing = !is.na(c3ppr01) & !is.na(c3ppl01) & 
                    !is.na(c3ltr01) & !is.na(c3ltl01),
                  min_c3_sens = pmin(c3ppr01, c3ppl01, c3ltr01, c3ltl01)) %>%
    dplyr::mutate(
      across(
        all_of(c('c2ppr01','c2ppl01','c2ltl01','c2ltr01')), 
        ~dplyr::case_when(is.na(.x) & c3_sens_not_missing ~ min_c3_sens, 
                   !is.na(.x) ~ .x,
                   !c3_sens_not_missing & is.na(.x) ~ NA_real_)))
  return(out)
}

sensory_vars = c("c2ppr","c3ppr","c4ppr","c5ppr","c6ppr","c7ppr","c8ppr", 
                 "c2ppl","c3ppl","c4ppl","c5ppl","c6ppl","c7ppl","c8ppl",
                 "c2ltl","c3ltl","c4ltl","c5ltl","c6ltl","c7ltl","c8ltl",
                 "c2ltr","c3ltr","c4ltr","c5ltr","c6ltr","c7ltr","c8ltr",
                 "t1ppr", "t2ppr", "t3ppr", "t4ppr", "t5ppr", "t6ppr", "t7ppr", "t8ppr", "t9ppr", "t10ppr","t11ppr","t12ppr",
                 "t1ppl", "t2ppl", "t3ppl", "t4ppl", "t5ppl", "t6ppl", "t7ppl", "t8ppl", "t9ppl", "t10ppl","t11ppl","t12ppl",
                 "t1ltr", "t2ltr", "t3ltr", "t4ltr", "t5ltr", "t6ltr", "t7ltr", "t8ltr", "t9ltr", "t10ltr","t11ltr","t12ltr",
                 "t1ltl", "t2ltl", "t3ltl", "t4ltl", "t5ltl", "t6ltl", "t7ltl", "t8ltl", "t9ltl", "t10ltl","t11ltl","t12ltl",
                 "l1ppr","l2ppr","l3ppr","l4ppr","l5ppr",
                 "l1ppl","l2ppl","l3ppl","l4ppl","l5ppl",
                 "l1ltr","l2ltr","l3ltr","l4ltr","l5ltr",
                 "l1ltl","l2ltl","l3ltl","l4ltl","l5ltl",
                 "s1ppr","s2ppr","s3ppr","s1ppl","s2ppl","s3ppl",
                 "s1ltr","s2ltr","s3ltr","s1ltl","s2ltl","s3ltl",
                 "s45ppr","s45ppl","s45ltl","s45ltr")

sensory_vars_ordered_r_pp = c("c2ppr","c3ppr","c4ppr","c5ppr","c6ppr","c7ppr","c8ppr", 
                              "t1ppr", "t2ppr", "t3ppr", "t4ppr", "t5ppr", "t6ppr", "t7ppr", "t8ppr", "t9ppr", "t10ppr","t11ppr","t12ppr",
                              "l1ppr","l2ppr","l3ppr","l4ppr","l5ppr",
                              "s1ppr","s2ppr","s3ppr",
                              "s45ppr")
sensory_vars_ordered_lr_pp = c(rbind(sub('r','l',sensory_vars_ordered_r_pp), sensory_vars_ordered_r_pp))
sensory_levels = sub('ppr','',sensory_vars_ordered_r_pp)

motor_vars = c('elbfll','elbflr','wrextl','wrextr',
               'elbexl','elbexr','finfll','finflr','finabl','finabr',
               'hipfll','hipflr','kneexl','kneetr','ankdol','ankdor',
               'gretol','gretor','ankpll','ankplr')

all_sensor_time_vars = apply(
  expand.grid(sensory_vars, c('01','04','08','16','26','52')), 1, paste, collapse = '')


splvl_notation = c('C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'T01', 'T02', 'T03', 'T04', 'T05', 'T06', 'T07', 'T08', 'T09', 'T10', 'T11', 'T12',
                   'L01','L02','L03','L04','L05','S01','S02','S03','S45')
names(splvl_notation) = c('c1', sensory_levels)

Benzel_scores = c(paste0('modben0', c(4, 8)), paste0('modben', c(16, 26, 52, 54)))

Bladder_control_scores = c(paste0('bladcd0', c(0, 1, 4, 8)), paste0('bladcd', c(16, 26, 52, 54)))
Bowel_control_scores = c(paste0('bowelc0', c(0, 1, 4, 8)), paste0('bowelc', c(16, 26, 52, 54)))

Sygen_data = raw_sygen %>%
  dplyr::mutate(
    visdt = as.Date(visdt, "%d.%m.%y"),
    yeardt = as.Date(format(as.Date(visdt, format="%d. %m. y"),"%y"), "%y"),
    injdt = as.Date(injdt, "%d.%m.%y"),
    injt_year = format(injdt,"%y"),
    injt_month = month(injdt),
    injt_month_factor = factor(injt_month),
    injt_day = day(injdt),
    injt_day_month = lubridate::yday(injdt), #paste(day(injdt), month(injdt), sep = '-'),
    injt_dayofweek = factor(weekdays(injdt),
                            levels = c('Monday', 'Tuesday', 'Wednesday',
                                       'Thursday','Friday', 'Saturday',
                                       'Sunday'),
                            ordered = FALSE),
    injt_weekend = factor(dplyr::case_when(injt_dayofweek %in% c('Saturday', 'Sunday') ~ 'weekend',
                                    is.na(injt_dayofweek) ~ NA_character_,
                                    .default = 'weekday')),
    #injt_holiday = factor(is.holiday(injdt)),
    #inj_season = factor(getSeason(injdt)),
    aracutdt = as.Date(aracutdt, "%d.%m.%y"),
    diacutdt = as.Date(diacutdt, "%d.%m.%y"),
    injtm = dplyr::case_when(injtm == 24 ~ 0, # noticed 24:00 and 00:00 being used, set default to 00:00
                      .default = injtm),
    emctm = dplyr::case_when(emctm == 24 ~ 0, # noticed 24:00 and 00:00 being used, set default to 00:00
                      .default = emctm),
    vittm2 = dplyr::case_when(vittm2 == 24 ~ 0, # noticed 24:00 and 00:00 being used, set default to 00:00
                       .default = vittm2),
    aracuttm = dplyr::case_when(aracuttm == 24 ~ 0, # noticed 24:00 and 00:00 being used, set default to 00:00
                         .default = aracuttm),
    diacuttm = dplyr::case_when(diacuttm == 24 ~ 0, # noticed 24:00 and 00:00 being used, set default to 00:00
                         .default = diacuttm),
    # combine date and time
    inj_dt_tm = as.POSIXct(as.character( paste( injdt, floor(injtm), 60*(injtm - floor(injtm)) ) ), 
                           format="%Y-%m-%d %H %M"),
    aracu_dt_tm = as.POSIXct(as.character( paste( aracutdt, floor(aracuttm), 60*(aracuttm - floor(aracuttm)) ) ), 
                             format="%Y-%m-%d %H %M"),
    diacu_dt_tm = as.POSIXct(as.character( paste( diacutdt, floor(diacuttm), 60*(diacuttm - floor(diacuttm)) ) ), 
                             format="%Y-%m-%d %H %M"),
    time_Anderson2014 = factor(dplyr::case_when( 0<=injtm & injtm <= 6 ~ 'Early AM',
                                          6<injtm & injtm <= 18 ~ 'Day',
                                          18<injtm & injtm <= (23 + 59/60) ~ 'Night',
                                          .default = NA_character_)),
    #https://now.aapmr.org/neurological-examination-and-classification-of-sci/
    COMPL_motor01 = factor(dplyr::case_when(vaccd01 == 1 ~ 'Incomplete', # based on VAC
                                     vaccd01 == 0 ~ 'Complete',
                                     .default = NA_character_)), # 16 have missing VAC data+
    COMPL_sensory01 = factor(dplyr::case_when(anyana01 == 1 | s45ltl01 %in% c(1,2) | s45ltr01  %in% c(1,2) | 
                                         s45ppl01  %in% c(1,2) | s45ppr01 %in% c(1,2) ~ 'Incomplete',
                                       anyana01 == 0 ~ 'Complete',
                                       .default = NA_character_)),
    COMPL01 = factor(dplyr::case_when(COMPL_motor01 == 'Complete' & COMPL_sensory01 == 'Complete' ~ 'Complete',
                               COMPL_motor01 == 'Incomplete' | COMPL_sensory01 == 'Incomplete' ~ 'Incomplete',
                               .default = NA_character_)),
    # plegia = factor(dplyr::case_when(splvl == 'T01' | grepl('C', splvl) ~ 'tetra', # EMSCI uses T2 as cutoff >tetra, <= para
    #                    is.na(splvl) ~ NA_character_,
    #                    .default = 'para')),
    splvl = factor(splvl, 
                   levels = c('C01','C02','C03',"C04", 'C05',
                              "C06", "C07", "C08", 'T01','T02',
                              'T03', "T04", "T05", 'T06', 'T07',
                              'T08', 'T09', "T10", "T11", 'T12'),
                   ordered = TRUE),
    ppscor52 = dplyr::case_when(ppscor52 == 999 ~ NaN, .default = ppscor52),
    ltscor52 = dplyr::case_when(ltscor52 == 999 ~ NaN, .default = ltscor52),
    ppscor01 = dplyr::case_when(ppscor01 == 999 ~ NaN, .default = ppscor01),
    ltscor01 = dplyr::case_when(ltscor01 == 999 ~ NaN, .default = ltscor01),
    AIS1 = factor(dplyr::case_when(ais1 == 'ND' | ais1 == 'NT' | ais1 == '' ~ NA_character_,
                            .default = ais1)),
    AIS4 = factor(dplyr::case_when(ais4 == 'ND' | ais4 == 'NT' | ais4 == '' ~ NA_character_,
                            .default = ais4)),
    AIS8 = factor(dplyr::case_when(ais8 == 'ND' | ais8 == 'NT' | ais8 == '' ~ NA_character_,
                            .default = ais8)),
    AIS16 = factor(dplyr::case_when(ais16 == 'ND' | ais16 == 'NT' | ais16 == '' ~ NA_character_,
                             .default = ais16)),
    AIS26 = factor(dplyr::case_when(ais26 == 'ND' | ais26 == 'NT' | ais26 == '' ~ NA_character_,
                             .default = ais26)),
    AIS52 = factor(dplyr::case_when(ais52 == 'ND' | ais52 == 'NT' | ais52 == '' ~ NA_character_,
                             .default = ais52)),
    ais1_NaN = factor(dplyr::case_when(ais1 == 'ND' | ais1 == 'NT' | ais1 == '' ~ NA_character_,
                                .default = ais1), 
                      levels = 
                        c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                      ordered = TRUE),
    ais4_NaN = factor(dplyr::case_when(ais4 == 'ND' | ais4 == 'NT' | ais4 == '' ~ NA_character_,
                                .default = ais4), 
                      levels = 
                        c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                      ordered = TRUE),
    ais8_NaN = factor(dplyr::case_when(ais8 == 'ND' | ais8 == 'NT' | ais8 == '' ~ NA_character_,
                                .default = ais8), 
                      levels = 
                        c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                      ordered = TRUE),
    ais16_NaN = factor(dplyr::case_when(ais16 == 'ND' | ais16 == 'NT' | ais16 == '' ~ NA_character_,
                                 .default = ais16), 
                       levels = 
                         c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                       ordered = TRUE),
    ais26_NaN = factor(dplyr::case_when(ais26 == 'ND' | ais26 == 'NT' | ais26 == '' ~ NA_character_,
                                 .default = ais26), 
                       levels = 
                         c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                       ordered = TRUE),
    ais52_NaN = factor(dplyr::case_when(ais52 == 'ND' | ais52 == 'NT' | ais52 == '' ~ NA_character_,
                                 .default = ais52), 
                       levels = 
                         c('AIS A','AIS B', 'AIS C', 'AIS D', 'AIS E'), 
                       ordered = TRUE),
    injury_cause = factor(dplyr::case_when(injcd == 1 ~ 'automobile',
                                    injcd == 2 ~ 'motor cycle',
                                    injcd == 3 ~ 'pedestrian',
                                    injcd == 4 ~ 'fall',
                                    injcd == 5 ~ 'water related (diving or surfing)',
                                    injcd == 6 ~ 'other sports',
                                    injcd == 7 ~ 'blunt trauma (assault)',
                                    injcd == 8 ~ 'gun shot wound',
                                    injcd == 9 ~ 'other',
                                    .default = NA_character_)),
    inf_w2 = factor((!aedur_infection02_1 %in% c('', 'UNK', NA) | #infection during week 2 or not
                       !aedur_infection02_2 %in% c('', 'UNK', NA)| #250+ patients did
                       !aedur_infection02_3 %in% c('', 'UNK', NA) |
                       !aedur_infection02_4 %in% c('', 'UNK', NA) |
                       !aedur_infection02_5 %in% c('', 'UNK', NA) |
                       !aedur_infection02_6 %in% c('', 'UNK', NA) |
                       !aedur_infection02_7 %in% c('', 'UNK', NA))),
    inf_w2_sev = factor(pmax(aesev_infection02_1, aesev_infection02_2,
                             aesev_infection02_3, aesev_infection02_4,
                             aesev_infection02_5, aesev_infection02_6,
                             aesev_infection02_7, na.rm = TRUE)),
    inf_w4 = factor((!aedur_infection04_1 %in% c('', 'UNK', NA) | #infection during week 2 or not
                       !aedur_infection04_2 %in% c('', 'UNK', NA)| #250+ patients did
                       !aedur_infection04_3 %in% c('', 'UNK', NA) |
                       !aedur_infection04_4 %in% c('', 'UNK', NA) )),
    inf_w4_sev = factor(pmax(aesev_infection04_1, aesev_infection04_2,
                             aesev_infection04_3, aesev_infection04_4, 
                             na.rm = TRUE)), 
    # parent 2011
    child = factor(dplyr::case_when(age<18 ~ 'child', age>=18 ~ 'adult', 
                             .default = NA_character_)),
    sexcd_raw = sexcd,
    sexcd = relevel(factor(sexcd), ref = 2), # set male which is larger group as reference
    racecd = factor(racecd),
    tx1_r = factor(tx1_r)
  )


locf <- function(dataset, var){
  var52 <- paste0(var, '52')
  var26 <- paste0(var, '26')
  var_recovery <- paste0(var, '_locf')
  dataset[var_recovery] <- NA
  for (i in c(1:dim(dataset)[1])){
    if (is.na(dataset[var52][[1]][i])){
      dataset[var_recovery][[1]][i] <- dataset[var26][[1]][i]
    } else {
      dataset[var_recovery][[1]][i] <- dataset[var52][[1]][i]
    }
  }
  return(dataset)
}

################################################################################
#################################### SYGEN #####################################
################################################################################

df_window.1.4 <- read.csv(
  '/Volumes/green_groups_bds_public/Projects/Drug/df_window.1.4.csv')

df_window.1.4$bin_lems_improv_5

## Redefine the PR group to exclude patients with AIS grade + improvement in UEMS 
## only and no or only marginal improvement in LEMS
sygen_rec_lems_5 <- df_window.1.4 %>% dplyr::filter(bin_lems_improv_5 == 1)
sygen_norec_lems_5 <- df_window.1.4 %>% dplyr::filter(is.na(bin_lems_improv_5) | 
                                                        bin_lems_improv_5 == 0)

## Check count UEMS < LEMS in both groups at recovery
## Sygen
sygen_rec_lems_5 %>% dplyr::select(PTNUM, Upper01, Upper04, Upper52,
                                   Lower01, Lower04, Lower26, Lower52)
# In PR group: 3 have UEMS < LEMS at recovery and one has UEMS == LEMS

sygen_norec_lems_5 %>% dplyr::select(PTNUM, Upper01, Upper04, Upper52,
                                   Lower01, Lower04, Lower26, Lower52)
# IN NPR group: none have UEMS < LEMS 
# (LEMS is very often 0 based on how we defined PR)

## EMSCI
raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('acute III', 'chronic'),
                Patientennummer %in% id_recoverers_EMSCI) %>% 
  dplyr::select(Patientennummer, ExamStage, UEMS, LEMS)
# In PR group: 2 (90260, 230046) have UEMS < LEMS at recovery 
# and one has only one point difference at acute 3 and NA at chronic

temp_norecov_emsci <- raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('acute III'),
                Patientennummer %in% id_norecovery_EMSCI) %>% 
  dplyr::select(Patientennummer, ExamStage, UEMS, LEMS)
# IN NPR group: none have UEMS < LEMS 
# (LEMS is very often 0 based on how we defined PR)

## Plot myotomes week 4 --> week 52
## Including clear sign of LOCF
## Note: see visualisation_myotome.R script

sygen_rec_lems_5 %>% dplyr::select(PTNUM, PPSCOR01, PPSCOR04, PPSCOR26, PPSCOR52,
                                   LTSCOR01, LTSCOR04, LTSCOR26, LTSCOR52)

sygen_rec_lems_5 %>% dplyr::select(PTNUM, AMOTC152, APNTLT52)

## Look at functional outcome for PR group
## EMSCI

test_recov_function_emsci <- raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('chronic'),
                Patientennummer %in% id_recoverers_EMSCI) %>% 
  #dplyr::select(Patientennummer, ExamStage, WISCI, X6min, X10m, TUG, SCIM2_TotalScore, SCIM3_TotalScore)
  dplyr::select(Patientennummer, ExamStage, names(raw_data_EMSCI)[188:233])

# 2 patients could walk (6-min 10MWT and TUG performed): 90260 and 230046
# note that they are the same patients for which UEMS < LEMS
# interestingly they also have low scores in SCIM, even for mobility
# SCIM score ranges from 18 to 66
# Need to check how the scores were distributed along the assessment 
# --> overall gives the impression that they are not so independent


################################################################################
## Look at tail distribution of delta LEMS to define PR not restricted to AIS A
################################################################################
## Sygen

# LOCF for LEMS
raw_sygen <- locf(raw_sygen, 'lower')
raw_sygen <- locf(raw_sygen, 'upper')

# Compute quantities of interest
# Delta in LEMS between week 1 and recovery
raw_sygen$diff_lems01 <- raw_sygen$lower_locf - raw_sygen$lower01
# Delta in LEMS between week 4 and recovery
raw_sygen$diff_lems04 <- raw_sygen$lower_locf - raw_sygen$lower04
# Mean improvement below NLI between week 1 and recovery
# recalculate NLI
input_sygen <- preprocess_more_NLI(Sygen_data)
calculated_NLI = check_NLI2(input_sygen, 'ptid')
calculated_NLI_subset <- calculated_NLI %>% dplyr::select('ptid', 'NLI_calc')
raw_sygen <- merge(raw_sygen, calculated_NLI_subset, by = "ptid")
raw_sygen$NLI_calc <- factor(raw_sygen$NLI_calc, 
                             levels = c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", 
                                        "t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"))

# compute the delta at the myotome level
motor_vars = c('elbfll52','elbflr','wrextl','wrextr',
               'elbexl','elbexr','finfll','finflr','finabl','finabr',
               'hipfll','hipflr','kneexl','kneetr','ankdol','ankdor',
               'gretol','gretor','ankpll','ankplr')

for (var in motor_vars){
  var_diff <- paste0('delta_', var, '_04_recovery')
  raw_sygen <- locf(raw_sygen, var)
  var04 <- paste0(var, '04')
  var_recovery <- paste0(var, '_locf')
  raw_sygen[var_diff] <- raw_sygen[var_recovery] - raw_sygen[var04]
}

# compute mean delta below NLI
motor_vars_delta <- motor_vars
for (var in motor_vars_delta){
  motor_vars_delta[motor_vars_delta == var] <- paste0('delta_', var, '_04_recovery')
}

raw_sygen$mean_delta_blwNLI <- NA
for (i in c(1:dim(raw_sygen)[1])){
  NLI = raw_sygen$NLI_calc[i]
  if (is.na(NLI)) {
    raw_sygen$mean_delta_blwNLI[i] <- NA
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) < 5){
    blwNLI <- motor_vars_delta
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) %in% c(5:8)) {
    blwNLI <- motor_vars_delta[-c(1:((match(NLI, levels(factor(raw_sygen$NLI_calc)))-4)*2))]
  } else {
    blwNLI <- motor_vars_delta[-c(1:10)]
  }
  if (!(is.na(NLI))) {
    temp <- raw_sygen %>% dplyr::filter(ptid == raw_sygen$ptid[i]) %>%
      dplyr::select(blwNLI)
    raw_sygen$mean_delta_blwNLI[i] <- rowMeans(temp, na.rm=T)
  }
}


# Mean improvement below NLI between week 4 and recovery

# Quick visualisation of the distribution of delta LEMS
ggplot(raw_sygen, aes(x=mean_delta_blwNLI)) +
  geom_histogram(binwidth=0.1, alpha=.5, position="identity") +
  geom_density()

library(statpsych)
mean_delta_blwNLI_complete <- raw_sygen$mean_delta_blwNLI[complete.cases(raw_sygen$mean_delta_blwNLI)]
test.skew(mean_delta_blwNLI_complete) # skewness=1.4382 pvalue=0 rejected
p95 <- quantile(mean_delta_blwNLI_complete, prob=0.95, na.rm=T)
mean_delta_blwNLI_complete_95 <- mean_delta_blwNLI_complete[mean_delta_blwNLI_complete<p95]
test.skew(mean_delta_blwNLI_complete_95) # skewness=1.2139 pvalue=0 rejected
p70 <- quantile(mean_delta_blwNLI_complete, prob=0.70, na.rm=T)
mean_delta_blwNLI_complete_70 <- mean_delta_blwNLI_complete[mean_delta_blwNLI_complete<p70]
test.skew(mean_delta_blwNLI_complete_70) # skewness=0.3979 pvalue=9e-04 rejected
p60 <- quantile(mean_delta_blwNLI_complete, prob=0.60, na.rm=T)
mean_delta_blwNLI_complete_60 <- mean_delta_blwNLI_complete[mean_delta_blwNLI_complete<p60]
test.skew(mean_delta_blwNLI_complete_60) # skewness=0.2143 pvalue=0.093 NOT rejected

#library(moments)
#kurtosis(raw_sygen$mean_delta_blwNLI, na.rm = T) #1.438226
#skewness(raw_sygen$mean_delta_blwNLI, na.rm = T) #4.476946

################################################################################
# Identify outliers...
################################################################################
# ...95% percentile in terms of delta LEMS week 04 to recovery (week 52 with LOCF when week 52 is missing)
id_top5perc_difflems <- as.vector(names(table(raw_sygen$ptid[raw_sygen$diff_lems04 > quantile(raw_sygen$diff_lems04, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_lems <- ifelse(raw_sygen$ptid %in% id_top5perc_difflems, 1, 0)

# ...95% percentile in terms of mean delta LEMS week 04 to recovery (week 52 with LOCF when week 52 is missing) below NLI
id_top5perc_diff_blwnli <- as.vector(names(table(raw_sygen$ptid[raw_sygen$mean_delta_blwNLI > quantile(raw_sygen$mean_delta_blwNLI, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_blwnli <- ifelse(raw_sygen$ptid %in% id_top5perc_diff_blwnli, 1, 0)

# ...95% percentile in delta LEMS week 04 to recovery within motor complete injuries at baseline
raw_sygen_baseline_complete <- raw_sygen %>% dplyr::filter(ais4 %in% c('AIS A', 'AIS B'))
id_top5perc_difflems_complete <- as.vector(names(table(raw_sygen_baseline_complete$ptid[
  raw_sygen_baseline_complete$diff_lems04 > 
    quantile(raw_sygen_baseline_complete$diff_lems04, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_lems_complete <- ifelse(raw_sygen$ptid %in% id_top5perc_difflems_complete, 1, 0)

# ...95% percentile in mean delta LEMS week 04 to recovery within motor complete injuries at baseline below NLI
id_top5perc_difflems_blwnli_complete <- as.vector(names(table(raw_sygen_baseline_complete$ptid[
  raw_sygen_baseline_complete$mean_delta_blwNLI > 
    quantile(raw_sygen_baseline_complete$mean_delta_blwNLI, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_blwnli_complete <- ifelse(raw_sygen$ptid %in% id_top5perc_difflems_blwnli_complete, 1, 0)

# ...95% percentile in delta LEMS week 04 to recovery within motor incomplete injuries at baseline
raw_sygen_baseline_incomplete <- raw_sygen %>% dplyr::filter(ais4 %in% c('AIS C', 'AIS D', 'AIS E'))
id_top5perc_difflems_incomplete <- as.vector(names(table(raw_sygen_baseline_incomplete$ptid[
  raw_sygen_baseline_incomplete$diff_lems04 > 
    quantile(raw_sygen_baseline_incomplete$diff_lems04, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_lems_incomplete <- ifelse(raw_sygen$ptid %in% id_top5perc_difflems_incomplete, 1, 0)

# ...95% percentile in mean delta LEMS week 04 to recovery within motor incomplete injuries at baseline below NLI
id_top5perc_difflems_blwnli_incomplete <- as.vector(names(table(raw_sygen_baseline_incomplete$ptid[
  raw_sygen_baseline_incomplete$mean_delta_blwNLI > 
    quantile(raw_sygen_baseline_incomplete$mean_delta_blwNLI, prob=0.95, na.rm=T)])))
#... and identify in the raw data in the raw data
raw_sygen$tail95_blwnli_incomplete <- ifelse(raw_sygen$ptid %in% id_top5perc_difflems_blwnli_incomplete, 1, 0)

# Quick check of the 95% percentile in delta LEMS in terms of aggregate motor scores
raw_sygen %>% dplyr::filter(tail95_blwnli == 1)  %>% 
  dplyr::select(ais1, ais4, ais52, splvl, NLI_calc,
                lower01, lower04, lower08, lower16, lower26, lower52,
                upper01, upper04, upper08, upper16, upper26, upper52)

# Compute the difference between UEMS and LEMS across time points
raw_sygen$diff_ms_01 <- raw_sygen$upper01 - raw_sygen$lower01
raw_sygen$diff_ms_04 <- raw_sygen$upper04 - raw_sygen$lower04
raw_sygen$diff_ms_08 <- raw_sygen$upper08 - raw_sygen$lower08
raw_sygen$diff_ms_16 <- raw_sygen$upper16 - raw_sygen$lower16
raw_sygen$diff_ms_26 <- raw_sygen$upper26 - raw_sygen$lower26
raw_sygen$diff_ms_52 <- raw_sygen$upper52 - raw_sygen$lower52
raw_sygen$diff_ms_locf <- raw_sygen$upper_locf - raw_sygen$lower_locf

# Binarize --> Is UEMS > LEMS?, 1=yes, 0=no
## Definition 1 CCS: UEMS - LEMS < 0
raw_sygen$diff_ms_locf_bin <- ifelse(raw_sygen$diff_ms_locf<0, 1, 0)
raw_sygen$diff_ms_52_bin <- ifelse(raw_sygen$diff_ms_52<0, 1, 0)
raw_sygen$diff_ms_26_bin <- ifelse(raw_sygen$diff_ms_26<0, 1, 0)
raw_sygen$diff_ms_16_bin <- ifelse(raw_sygen$diff_ms_16<0, 1, 0)
raw_sygen$diff_ms_08_bin <- ifelse(raw_sygen$diff_ms_08<0, 1, 0)
raw_sygen$diff_ms_04_bin <- ifelse(raw_sygen$diff_ms_04<0, 1, 0)
raw_sygen$diff_ms_01_bin <- ifelse(raw_sygen$diff_ms_01<0, 1, 0)

## Definition 2 CCS: UEMS - LEMS < 4
raw_sygen$diff_ms_52_bin5 <- ifelse(raw_sygen$diff_ms_52<(-4), 1, 0)
raw_sygen$diff_ms_locf_bin5 <- ifelse(raw_sygen$diff_ms_locf<(-4), 1, 0)

## Definition 3 CCS: UEMS - LEMS < 9
raw_sygen$diff_ms_52_bin10 <- ifelse(raw_sygen$diff_ms_52<(-9), 1, 0)
raw_sygen$diff_ms_locf_bin10 <- ifelse(raw_sygen$diff_ms_locf<(-9), 1, 0)

## Definition 4 CCS: UEMS - LEMS < 18
raw_sygen$diff_ms_52_bin19 <- ifelse(raw_sygen$diff_ms_52<(-18), 1, 0)
raw_sygen$diff_ms_locf_bin19 <- ifelse(raw_sygen$diff_ms_locf<(-18), 1, 0)

## Definition 5 CCS: based on NLI --> (1- (mean_UEMS_blw_NLI/mean_LEMS))*100 > 0.1

motor_vars_UEMS_52 <- c('elbfll52','elbflr52','wrextl52','wrextr52',
                        'elbexl52','elbexr52','finfll52','finflr52',
                        'finabl52','finabr52')
raw_sygen$mean_UEMS_blw_NLI_w52 <- NA
for (i in c(1:dim(raw_sygen)[1])){
  NLI = raw_sygen$NLI_calc[i]
  if (is.na(NLI)) {
    raw_sygen$mean_UEMS_blw_NLI_w52[i] <- NA
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) < 5){
    blwNLI <- motor_vars_UEMS_52
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) %in% c(5:8)) {
    blwNLI <- motor_vars_UEMS_52[-c(1:((match(NLI, levels(factor(raw_sygen$NLI_calc)))-4)*2))]
  } else {
    blwNLI <- motor_vars_UEMS_52[-c(1:10)]
  }
  if (!(is.na(NLI))) {
    temp <- raw_sygen %>% dplyr::filter(ptid == raw_sygen$ptid[i]) %>%
      dplyr::select(blwNLI)
    raw_sygen$mean_UEMS_blw_NLI_w52[i] <- rowMeans(temp, na.rm=T)
  }
}

raw_sygen$mean_LEMS_w52 <- rowMeans(raw_sygen[ , c('hipfll52','hipflr52',
                                                   'kneexl52','kneetr52',
                                                   'ankdol52','ankdor52',
                                                   'gretol52','gretor52',
                                                   'ankpll52','ankplr52')], 
                                    na.rm=TRUE)
raw_sygen$diff_ms_52_NLI <- (1 - (raw_sygen$mean_UEMS_blw_NLI_w52/
                                    raw_sygen$mean_LEMS_w52))*100
raw_sygen$diff_ms_52_bin_NLI <- ifelse(raw_sygen$diff_ms_52_NLI>0.1, 1, 0)

## Definition 5 CCS: based on NLI --> (1- (mean_UEMS_blw_NLI/mean_LEMS))*100 > 0.1
### INCLUDING LOCF

for (var in c('elbfll','elbflr','wrextl','wrextr',
              'elbexl','elbexr','finfll','finflr',
              'finabl','finabr', 'hipfll','hipflr',
              'kneexl','kneetr',
              'ankdol','ankdor',
              'gretol','gretor',
              'ankpll','ankplr')){
  raw_sygen <- locf(raw_sygen, var)
}

motor_vars_UEMS_locf <- c('elbfll_locf','elbflr_locf','wrextl_locf','wrextr_locf',
                        'elbexl_locf','elbexr_locf','finfll_locf','finflr_locf',
                        'finabl_locf','finabr_locf')
raw_sygen$mean_UEMS_blw_NLI_w_locf <- NA
for (i in c(1:dim(raw_sygen)[1])){
  NLI = raw_sygen$NLI_calc[i]
  if (is.na(NLI)) {
    raw_sygen$mean_UEMS_blw_NLI_w_locf[i] <- NA
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) < 5){
    blwNLI <- motor_vars_UEMS_locf
  } else if (match(NLI, levels(factor(raw_sygen$NLI_calc))) %in% c(5:8)) {
    blwNLI <- motor_vars_UEMS_locf[-c(1:((match(NLI, levels(factor(raw_sygen$NLI_calc)))-4)*2))]
  } else {
    blwNLI <- motor_vars_UEMS_locf[-c(1:10)]
  }
  if (!(is.na(NLI))) {
    temp <- raw_sygen %>% dplyr::filter(ptid == raw_sygen$ptid[i]) %>%
      dplyr::select(blwNLI)
    raw_sygen$mean_UEMS_blw_NLI_w_locf[i] <- rowMeans(temp, na.rm=T)
  }
}

raw_sygen$mean_LEMS_w_locf <- rowMeans(raw_sygen[ , c('hipfll_locf','hipflr_locf',
                                                   'kneexl_locf','kneetr_locf',
                                                   'ankdol_locf','ankdor_locf',
                                                   'gretol_locf','gretor_locf',
                                                   'ankpll_locf','ankplr_locf')], 
                                    na.rm=TRUE)
raw_sygen$diff_ms_locf_NLI <- (1 - (raw_sygen$mean_UEMS_blw_NLI_w_locf/
                                    raw_sygen$mean_LEMS_w_locf))*100
raw_sygen$diff_ms_locf_bin_NLI <- ifelse(raw_sygen$diff_ms_locf_NLI>0.1, 1, 0)

################################################################################


diff_ms_top5 <- raw_sygen %>% dplyr::filter(tail95_blwnli == 1)  %>% 
  dplyr::select(ptid, ais1, ais52, splvl, diff_ms_01, diff_ms_04,
                diff_ms_08, diff_ms_16, diff_ms_26, diff_ms_52, lower52, upper52, 
                diff_ms_52_bin, diff_ms_52_bin5, diff_ms_52_bin10, diff_ms_52_bin19,diff_ms_52_bin_NLI)

diff_ms_notop5<- raw_sygen %>% dplyr::filter(tail95_blwnli == 0)  %>% 
  dplyr::select(ptid, ais1, ais52, splvl, diff_ms_01, diff_ms_04,
                diff_ms_08, diff_ms_16, diff_ms_26, diff_ms_52, lower52, upper52, 
                diff_ms_52_bin, diff_ms_52_bin5, diff_ms_52_bin10, diff_ms_52_bin19, mean_delta_blwNLI)

diff_ms_notop5_but1bin_NLI <- raw_sygen %>% dplyr::filter(tail95_blwnli == 0, diff_ms_52_bin_NLI == 1)  %>% 
  dplyr::select(ptid, ais1, ais52, splvl, diff_ms_01, diff_ms_04,
                diff_ms_08, diff_ms_16, diff_ms_26, diff_ms_52, lower52, upper52, 
                diff_ms_52_bin, diff_ms_52_bin5, diff_ms_52_bin10, diff_ms_52_bin19, mean_delta_blwNLI)

# Quick visualisation of the distribution of delta LEMS
mu <- ddply(raw_sygen, "diff_ms_52_bin_NLI", summarise, grp.mean=mean(mean_delta_blwNLI, na.rm=T))
ggplot(raw_sygen, aes(x=mean_delta_blwNLI, color=factor(diff_ms_52_bin_NLI))) +
  geom_histogram(binwidth=0.1, alpha=.5, position="identity", fill='white')+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=factor(diff_ms_52_bin_NLI)),
             linetype="dashed")

# # Check characteristics of the PR group
# temp <- raw_sygen %>% dplyr::filter(!(ptid %in% id_top5perc_difflems))
# table(temp$tx1_r)
# table(diff_ms_top5$diff_ms_52)

# Convert 95% data with diff UEMS-LEMS to long format for plotting
diff_ms_top5_long <- melt(diff_ms_top5, id=c('ptid', 'ais1', 'ais52', 'splvl'))

diff_ms_top5_long[grepl("C", diff_ms_top5_long$splvl), "level"] <- "C"
diff_ms_top5_long[grepl("T", diff_ms_top5_long$splvl), "level"] <- "T"

# Plot UEMS-LEMS overtime at patient level
ggplot(diff_ms_top5_long, aes(x = variable, y = value)) +   
  geom_line(aes(group = ptid))

table(raw_sygen$diff_ms_01_bin, raw_sygen$tail95_lems)
chisq.test(table(raw_sygen$diff_ms_16_bin, raw_sygen$tail95_lems), simulate.p.value = T)

################################################################################
#Quick look at some medications in the PR vs NPR groups
vancomycin_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/vancomycin.csv')
morphine_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/morphine.csv')
acetaminophen_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/acetaminophen.csv')

length(unique(vancomycin_master$NEW_ID))

vancomycin_tail_LEMS <- vancomycin_master %>% 
  dplyr::filter(NEW_ID %in% unique(diff_ms_top5$ptid)) %>% 
  dplyr::select('NEW_ID':'X30')
length(unique(vancomycin_tail_LEMS$NEW_ID))

vancomycin_nottail_LEMS <- vancomycin_master %>% 
  dplyr::filter(!(NEW_ID %in% unique(diff_ms_top5$ptid))) %>% 
  dplyr::select('NEW_ID':'X30')

vancomycin_tail_LEMS_30days <- vancomycin_tail_LEMS[rowSums(is.na(vancomycin_tail_LEMS)) != ncol(vancomycin_tail_LEMS)-2, ]
vancomycin_nottail_LEMS_30days <- vancomycin_nottail_LEMS[rowSums(is.na(vancomycin_nottail_LEMS)) != ncol(vancomycin_nottail_LEMS)-2, ]
################################################################################

recover_ESMCI_va <- raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('very acute'),
                Patientennummer %in% id_recoverers_EMSCI)

notrecover_ESMCI_va <- raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('very acute'),
                Patientennummer %in% id_norecovery_EMSCI)
notrecover_ESMCI_c <- raw_data_EMSCI %>% 
  dplyr::filter(ExamStage %in% c('chronic'),
                Patientennummer %in% id_norecovery_EMSCI)

################################################################################
# Comparison of baseline characteristics of the top 5% in terms of mean delta MS below NLI
# compared to the other 95%

raw_sygen <- raw_sygen %>% 
  mutate(NLI_calc_bin = if_else(str_detect(raw_sygen$NLI_calc, "t"), "t", "c"))

raw_sygen[raw_sygen == ''] <- NA

# Perform propensity score matching to obtain comparable groups in terms of baseline characteristics

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

perform_psm <- function(raw_sygen, criteria){
  # Select for complete cases as PSM only works with complete data
  sygen_subset_psm <- raw_sygen %>% select(c('ptid', 'sexcd', 'age', 'lower04', 
                                             'upper04', 'ais4', 'NLI_calc_bin', 
                                             'lower_locf',
                                             criteria))
  #sygen_subset_psm[criteria] <- as.factor(sygen_subset_psm[criteria])
  sygen_complete_psm <- sygen_subset_psm[complete.cases(sygen_subset_psm), ]
  
  print('Characteristics before PSM:')
  print(table1(formula(paste("~ factor(sexcd) + age + lower04 + upper04 + factor(ais4) + NLI_calc_bin | factor(", criteria, ")")), 
         data=sygen_complete_psm, overall=F, extra.col=list(`P-value`=pvalue)))
  
  mod_test1 <- glm(formula(paste('factor(', criteria, ") ~ factor(sexcd) + age + lower04 + upper04 + factor(ais4) + NLI_calc_bin")),
                   data=sygen_complete_psm, family='binomial')
  print(summary(mod_test1))
  
  psa_nn <- matchit(formula(paste0('factor(', criteria, ') ~ factor(ais4) + upper04 + NLI_calc_bin + factor(sexcd) + age + lower04')),
                    data=sygen_complete_psm, distance='glm', method = 'nearest', m.order = 'largest', replace = F, 
                    ratio=4)
  psa_nn
  print(summary(psa_nn))
  
  love.plot(bal.tab(psa_nn),
            stat=c('m', 'v'),
            grid=T,
            thresholds = c(m=.25,v=1.25))
  
  psa_matched <- match.data(psa_nn)
  
  print('Characteristics after PSM:')
  print(table1(formula(paste("~ factor(sexcd) + age + lower04 + upper04 + factor(ais4) + NLI_calc_bin | factor(", criteria, ")")),
         data=psa_matched, overall=F, extra.col=list(`P-value`=pvalue)))
  
  raw_sygen_matched_temp <- raw_sygen %>% filter(ptid %in% psa_matched$ptid)
  
  print('Other characteristics after PSM, (not matched):')
  
  print(table1(formula(paste('~ lower_locf + upper_locf + factor(ais52) + 
                       factor(anyana04) + factor(vaccd04) + 
                       factor(anyana52) + factor(vaccd52) + factor(tx1_r) | factor(', criteria, ')')), 
         data=raw_sygen_matched_temp, overall=F, extra.col=list(`P-value`= pvalue)))
  
  print('Outcomes after PSM:')
  print(table1(formula(paste('~ factor(diff_ms_52_bin) + factor(diff_ms_52_bin5) + 
           factor(diff_ms_52_bin10) + factor(diff_ms_52_bin19) +
           factor(diff_ms_52_bin_NLI) | factor(', criteria, ')')), 
         data=raw_sygen_matched_temp, overall=F, extra.col=list(`P-value`=pvalue)))
  
  return (raw_sygen_matched_temp)
}

sygen_tail95_blwnli_psm <- perform_psm(raw_sygen, 'tail95_blwnli')
sygen_tail95_blwnli_psm$year_injury <- substr(sygen_tail95_blwnli_psm$injdt, 7, 10)
# added cause of injury and year of injury --> no differences between groups
print(table1(formula(paste("~ factor(sexcd) + age + lower04 + upper04 + factor(ais4) + NLI_calc_bin + factor(injcd) + 
                           factor(year_injury) + as.numeric(year_injury) | factor(", 'tail95_blwnli', ")")),
             data=sygen_tail95_blwnli_psm, overall=F, extra.col=list(`P-value`=pvalue)))

print(table1(formula(paste('~ factor(diff_ms_locf_bin) + factor(diff_ms_locf_bin5) + 
           factor(diff_ms_locf_bin10) + factor(diff_ms_locf_bin19) +
           factor(diff_ms_locf_bin_NLI) + lower_locf + upper_locf | factor(', 'tail95_blwnli', ')')), 
             data=sygen_tail95_blwnli_psm, overall=F, extra.col=list(`P-value`= pvalue)))

sygen_tail95_blwnli_complete <- perform_psm(raw_sygen %>% dplyr::filter(ais4 %in% c('AIS A', 'AIS B')), 'tail95_blwnli_complete')
sygen_tail95_blwnli_incomplete <- perform_psm(raw_sygen %>% dplyr::filter(ais4 %in% c('AIS C', 'AIS D', 'AIS E')), 'tail95_blwnli_incomplete')

# Combine drug files for the corresponding groups obtained

files <- list.files(path="/Volumes/green_groups_bds_public/Projects/Drug/Drugs",
                    pattern="*.csv", full.names=TRUE, recursive=FALSE)

extract_drugs_psm <- function(ID, name){
  df.alldrugs.psm.temp <- data.frame()
  for (i in c(1:length(files))){ # go through all files
      temp <- read.csv(files[i]) # read file i
      temp <- temp[,!names(temp) %in% c("X")] # remove index column
      temp_psm <- temp[temp$NEW_ID %in% ID,] # select for NPR group
      # combine file i to the overall dataframes
      df.alldrugs.psm.temp <- rbind(df.alldrugs.psm.temp, temp_psm)
  }
  # Convert NA to 0
  df.alldrugs.psm[is.na(df.alldrugs.psm)] <- 0
  # Combine values such that there is only 1 line per patient per drug
  df.alldrugs.psm <- plyr::ddply(df.alldrugs.psm, c("NEW_ID", "generic_name"), plyr::numcolwise(sum))
  # Save file
  write.csv(df.alldrugs.psm.temp, paste0('/Volumes/green_groups_bds_public/Projects/Drug/', name, '.csv'), row.names=FALSE)
  return(df.alldrugs.psm.temp)
}

# df.alldrugs.psm.all <- extract_drugs_psm(psa_matched$ptid, 'df_alldrugs_psm_new')
# df.alldrugs.psm.all_complete <- extract_drugs_psm(sygen_tail95_blwnli_complete$ptid, 'df_alldrugs_psm_complete')
# df.alldrugs.psm.all_incomplete <- extract_drugs_psm(sygen_tail95_blwnli_incomplete$ptid, 'df_alldrugs_psm_incomplete')
  
df.alldrugs.psm.all <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_psm_new.csv')
df.alldrugs.psm.all_complete <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_psm_complete.csv')
df.alldrugs.psm.all_incomplete <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/df_alldrugs_psm_incomplete.csv')

# Preprocessing drug files
# Select for the first 30 days
select.30days.drugs <- function(drug_data){
  drug_data.30days <- drug_data %>% select('NEW_ID':'X30')
  drug_data.30days <- drug_data.30days %>%
    mutate(rowsums = select(., -c(1:2)) %>% 
             rowSums(na.rm = TRUE))
  # Remove rows where sum = 0, ie. this drug was not prescribed in the first 30 days
  drug_data.30days <- drug_data.30days[drug_data.30days$rowsums > 0, ]
  return(drug_data.30days)
}
df.alldrugs.psm.all.30days <- select.30days.drugs(df.alldrugs.psm.all)
df.alldrugs.psm.all_complete.30days <- select.30days.drugs(df.alldrugs.psm.all_complete)
df.alldrugs.psm.all_incomplete.30days <- select.30days.drugs(df.alldrugs.psm.all_incomplete)

#Add classification to drug dataframe and transform to wide
labels_binarise_wide <- function(data){
  label_tail95 <- raw_sygen %>% select(ptid, tail95_blwnli, 
                                       tail95_blwnli_complete, 
                                       tail95_blwnli_incomplete, diff_ms_52_bin, 
                                       diff_ms_52_bin5,  diff_ms_52_bin10, 
                                       diff_ms_52_bin19, diff_ms_52_bin_NLI, 
                                       lower_locf, upper_locf)
  data.labels <- merge(data, label_tail95, by.x = "NEW_ID", 
                       by.y = "ptid", all.x = TRUE, all.y = FALSE)
  binary_data_wide <- data.labels[,1:2] %>%
    mutate(yesno = 1) %>%
    distinct %>%
    tidyr::spread(generic_name, yesno, fill = 0)
  binary_data_wide.labels <- merge(binary_data_wide, label_tail95, by.x = "NEW_ID", 
                                   by.y = "ptid", all.x = TRUE, all.y = FALSE)
  binary_data_wide.labels[sapply(binary_data_wide.labels, is.numeric)] <- 
    lapply(binary_data_wide.labels[sapply(binary_data_wide.labels, is.numeric)],
           as.factor)
  return(list(data.labels, binary_data_wide.labels))
}
drugs.labels_all_list <- labels_binarise_wide(df.alldrugs.psm.all.30days)
drugs.30days.long.labels_all <- drugs.labels_all_list[1][[1]]
drugs.30days.wide.bin.labels_all <- drugs.labels_all_list[2][[1]]

drugs.labels_all_complete_list <- labels_binarise_wide(df.alldrugs.psm.all_complete.30days)
drugs.30days.long.labels_all_complete <- drugs.labels_all_complete_list[1][[1]]
drugs.30days.wide.bin.labels_complete <- drugs.labels_all_complete_list[2][[1]]

drugs.labels_all_incomplete_list <- labels_binarise_wide(df.alldrugs.psm.all_incomplete.30days)
drugs.30days.long.labels_all_incomplete <- drugs.labels_all_incomplete_list[1][[1]]
drugs.30days.wide.bin.labels_incomplete <- drugs.labels_all_incomplete_list[2][[1]]

#Analysis antibiotics for each population
# 0. Check if the list of antibiotics from the PR/NPR groups cover all antibiotics observe in the newly defined groups from PSM
unique_drugs_all_psm <- unique(unique(df.alldrugs.psm.all.30days$generic_name, 
                                      df.alldrugs.psm.all_complete.30days$generic_name),
                               df.alldrugs.psm.all_incomplete.30days$generic_name)
unique_drugs_all_N.PR <- unique(df.alldrugs.recovery_copy$generic_name, df.alldrugs.norecovery_copy$generic_name)

drugs_to_classify <- setdiff(unique_drugs_all_psm, unique_drugs_all_N.PR)

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

antibio_PSM_groups <- c('neomycin', 'ampicillin_and_sulbactam', 'cefoxitin',
                        'imipenem_and_cilastatin', 'penicillin', 'cephradine',
                        'framycetin', 'cefuroxime', 'norfloxacin', 'cloxacillin',
                        'ceftizoxime', 'dicloxacillin', 'cefotetan', 'rifampin',
                        'cephapirin', 'imipenem', 'mupirocin', 'lactobacillus',
                        'cefoperazone', 'amikacin', 'cephalexin', 'tetracycline',
                        'piperacillin_and_tazobactam', 'sulfacetamide', 'aztreonam',
                        'mezlocillin', 'clarithromycin')

other_parasites_PSM <- c('fluconazole', 'ketoconazole', 
                     'neomycin_sulfate_polymyxin_B_sulfate_bacitracin_zinc_and_hydrocortisone',
                     'butoconazole_nitrate', 'sodium_hypochlorite', 'terconazole',
                     'povidone_iodine')

combined_antibio <- c(all_antibio, antibio_PSM_groups)

# 1. Comparing prevalence of common antibiotics
prevalence_common_antibio <- function(data.long, data.wide, criteria){
  common_drugs_psm <- intersect(unique(data.long %>% filter((!!as.symbol(criteria)) == 1) %>% select(generic_name))[[1]],
                                unique(data.long %>% filter((!!as.symbol(criteria)) == 0) %>% select(generic_name))[[1]])
  common_antibio_psm <- intersect(combined_antibio, common_drugs_psm)
  data.wide.antibio <- data.wide %>% dplyr::select(any_of(c(common_antibio_psm, criteria)))
  table <- table1(formula(paste('~ . | factor(', criteria, ')')), 
         data=data.wide.antibio, overall=F, extra.col=list(`P-value`= pvalue))
  return(table)
}

prevalence_common_antibio_all <- prevalence_common_antibio(drugs.30days.long.labels_all,
                                                           drugs.30days.wide.bin.labels_all,
                                                           'tail95_blwnli')
prevalence_common_antibio_complete <- prevalence_common_antibio(drugs.30days.long.labels_all_complete,
                                                                drugs.30days.wide.bin.labels_complete,
                                                           'tail95_blwnli_complete')
prevalence_common_antibio_incomplete <- prevalence_common_antibio(drugs.30days.long.labels_all_incomplete,
                                                                  drugs.30days.wide.bin.labels_incomplete,
                                                                'tail95_blwnli_incomplete')

# 2. Examine prevalence of antibiotics given in only 1 of the 2 groups
prevalence_1gp_antibio <- function(data.long, criteria){
  all_drugs_psm <- unique(as.vector(data.long %>% select(generic_name)))[[1]]
  all_antibio_psm <- intersect(as.vector(all_drugs_psm), as.vector(combined_antibio))
  common_drugs_psm <- intersect(unique(data.long %>% filter((!!as.symbol(criteria)) == 1) %>% select(generic_name))[[1]],
                                unique(data.long %>% filter((!!as.symbol(criteria)) == 0) %>% select(generic_name))[[1]])
  common_antibio_psm <- intersect(combined_antibio, common_drugs_psm)
  
  unique_antibio_1gp <- setdiff(all_antibio_psm, common_antibio_psm)
  table <- table(data.long %>% filter(generic_name %in% unique_antibio_1gp) %>% select(generic_name))
  return(table)
}

prevalence_1gp_antibio_all <- prevalence_1gp_antibio(drugs.30days.long.labels_all, 'tail95_blwnli')
prevalence_1gp_antibio_all
drugs.30days.long.labels_all %>% filter(generic_name == 'ceftizoxime') %>% select(tail95_blwnli)

prevalence_1gp_antibio_complete <- prevalence_1gp_antibio(drugs.30days.long.labels_all_complete, 'tail95_blwnli_complete')
prevalence_1gp_antibio_complete
drugs.30days.long.labels_all_complete %>% filter(generic_name == 'piperacillin') %>% select(tail95_blwnli_complete, NEW_ID)
chisq.test(table(drugs.30days.wide.bin.labels_complete$tail95_blwnli_complete, 
                 drugs.30days.wide.bin.labels_complete$piperacillin))$p.value

prevalence_1gp_antibio_incomplete <- prevalence_1gp_antibio(drugs.30days.long.labels_all_incomplete, 'tail95_blwnli_incomplete')
prevalence_1gp_antibio_incomplete
drugs.30days.long.labels_all_incomplete %>% filter(generic_name == 'piperacillin') %>% select(tail95_blwnli_incomplete, NEW_ID)
chisq.test(table(drugs.30days.wide.bin.labels_incomplete$tail95_blwnli_incomplete, 
                 drugs.30days.wide.bin.labels_incomplete$piperacillin))$p.value

# 3. Number of patients with antibiotics on day 0 or 1
#Number of patients with antibio prescribed on Day 0
temp <- unique(drugs.30days.long.labels_all %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X0, tail95_blwnli))
temp[is.na(temp)] <- 0
ID_X0 <- temp[temp$X0>0,]$NEW_ID
table(temp$X0>0, temp$tail95_blwnli)
chisq.test(table(temp$X0>0, temp$tail95_blwnli))
# When considering all Sygen: 7 (6%) from control group and 3 (10%) from outlier group had antibios at day 0

temp <- unique(drugs.30days.long.labels_all_complete %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X0, tail95_blwnli_complete))
temp[is.na(temp)] <- 0
table(temp$X0>0, temp$tail95_blwnli_complete)
chisq.test(table(temp$X0>0, temp$tail95_blwnli_complete))
# When considering complete injuries: 4 (5%) from control group and 2 (10%) from outlier group had antibios at day 0

temp <- unique(drugs.30days.long.labels_all_incomplete %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X0, tail95_blwnli_incomplete))
temp[is.na(temp)] <- 0
table(temp$X0>0, temp$tail95_blwnli_incomplete)
chisq.test(table(temp$X0>0, temp$tail95_blwnli_incomplete))
# When considering incomplete injuries: 6 (15%) from control group and 0 (0%) from outlier group had antibios at day 0

#Number of patients with antibio prescribed on Day 1
temp <- unique(drugs.30days.long.labels_all %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X1, tail95_blwnli))
temp[is.na(temp)] <- 0
table(temp$X1>0, temp$tail95_blwnli)
ID_X1 <- temp[temp$X1>0,]$NEW_ID
chisq.test(table(temp$X1>0, temp$tail95_blwnli))
# When considering all Sygen: 16 (13%) from control group and 8 (24%) from outlier group had antibios at day 0

temp <- unique(drugs.30days.long.labels_all_complete %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X1, tail95_blwnli_complete))
temp[is.na(temp)] <- 0
table(temp$X1>0, temp$tail95_blwnli_complete)
chisq.test(table(temp$X1>0, temp$tail95_blwnli_complete))
# When considering complete injuries: 12 (16%) from control group and 5 (26%) from outlier group had antibios at day 0

temp <- unique(drugs.30days.long.labels_all_incomplete %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X1, tail95_blwnli_incomplete))
temp[is.na(temp)] <- 0
table(temp$X1>0, temp$tail95_blwnli_incomplete)
chisq.test(table(temp$X1>0, temp$tail95_blwnli_incomplete))
# When considering incomplete injuries: 6 (31%) from control group and 6 (67%) from outlier group had antibios at day 0

temp <- unique(drugs.30days.long.labels_all %>% filter(generic_name %in% combined_antibio) %>% select(NEW_ID, X0, X1, tail95_blwnli))
temp[is.na(temp)] <- 0
temp <- temp %>% mutate(Max = pmax(X0, X1))
temp_ordered <- temp[order(temp$NEW_ID, -abs(temp$Max) ), ]
temp_ordered_unique <- temp_ordered[ !duplicated(temp_ordered$NEW_ID), ]

table(temp_ordered_unique$Max>0, temp_ordered_unique$tail95_blwnli)
chisq.test(table(temp$Max>0, temp$tail95_blwnli))


# 4. Number of unique antibiotics in the first 30 days

# Number of unique antibiotics per patient in the first 30 days - all
table(drugs.30days.long.labels_all_antibio$NEW_ID)
mean(table(drugs.30days.long.labels_all_antibio$NEW_ID))
sd(table(drugs.30days.long.labels_all_antibio$NEW_ID))
median(table(drugs.30days.long.labels_all_antibio$NEW_ID))
IQR(table(drugs.30days.long.labels_all_antibio$NEW_ID))
# by subgroups
drugs.30days.long.labels_all_antibio <- drugs.30days.long.labels_all %>% dplyr::filter(generic_name %in% combined_antibio)
drugs.30days.long.labels_all_antibio_control <- drugs.30days.long.labels_all_antibio %>% dplyr::filter(tail95_blwnli == 0)
summary(as.data.frame(table(drugs.30days.long.labels_all_antibio_control$NEW_ID)))
#min: 1, max:8, median[q1-q3] = 3[2-4], mean=3.575
drugs.30days.long.labels_all_antibio_recov <- drugs.30days.long.labels_all_antibio %>% dplyr::filter(tail95_blwnli == 1)
summary(as.data.frame(table(drugs.30days.long.labels_all_antibio_recov$NEW_ID)))
#min: 1, max:7, median[q1-q3] = 3[2-5.5], mean=3.926

ks.test(as.data.frame(table(drugs.30days.long.labels_all_antibio_recov$NEW_ID))$Freq,
        as.data.frame(table(drugs.30days.long.labels_all_antibio_control$NEW_ID))$Freq)

drugs.30days.long.labels_all_complete_antibio <- drugs.30days.long.labels_all_complete %>% dplyr::filter(generic_name %in% combined_antibio)
drugs.30days.long.labels_all_complete_antibio_control <- drugs.30days.long.labels_all_complete_antibio %>% dplyr::filter(tail95_blwnli_complete == 0)
summary(as.data.frame(table(drugs.30days.long.labels_all_complete_antibio_control$NEW_ID)))
#min: 1, max:23, median[q1-q3] = 5[3-8], mean=6.105
drugs.30days.long.labels_all_complete_antibio_recov <- drugs.30days.long.labels_all_complete_antibio %>% dplyr::filter(tail95_blwnli_complete == 1)
summary(as.data.frame(table(drugs.30days.long.labels_all_complete_antibio_recov$NEW_ID)))
#min: 1, max:16, median[q1-q3] = 5[4-8], mean=6.368

drugs.30days.long.labels_all_incomplete_antibio <- drugs.30days.long.labels_all_incomplete %>% dplyr::filter(generic_name %in% combined_antibio)
drugs.30days.long.labels_all_incomplete_antibio_control <- drugs.30days.long.labels_all_incomplete_antibio %>% dplyr::filter(tail95_blwnli_incomplete == 0)
summary(as.data.frame(table(drugs.30days.long.labels_all_incomplete_antibio_control$NEW_ID)))
#min: 2, max:17, median[q1-q3] = 5.5[3-9.75], mean=6.912
drugs.30days.long.labels_all_incomplete_antibio_recov <- drugs.30days.long.labels_all_incomplete_antibio %>% dplyr::filter(tail95_blwnli_incomplete == 1)
summary(as.data.frame(table(drugs.30days.long.labels_all_incomplete_antibio_recov$NEW_ID)))
#min: 2, max:20, median[q1-q3] = 5[2.75-6], mean=6.125

# 5. Number of cumulative antibiotics-days in the first 30 days
# Number of antibiotics in the first 30 days after SCI
drugs.30days.long.labels_all_antibio_bin <- drugs.30days.long.labels_all_antibio %>%
  dplyr::mutate(across(X0:X30, ~case_when(.x > 0 ~ 1, TRUE ~ 0))) 

drugs.30days.long._all_antibio_bin <- drugs.30days.long.labels_all_antibio_bin %>%
  dplyr::select(NEW_ID:X30)
drugs.30days.long._all_antibio_bin_sum <- plyr::ddply(drugs.30days.long._all_antibio_bin, c("NEW_ID"), plyr::numcolwise(sum))

drugs.30days.long._all_antibio_bin_sum <- drugs.30days.long._all_antibio_bin_sum %>%
  mutate(rowsums = select(., -c(1)) %>% 
           rowSums(na.rm = TRUE))
## RE ADD THE LABEL OF INTEREST
drugs.30days.long._all_antibio_bin_sum_label <- merge(drugs.30days.long._all_antibio_bin_sum, drugs.30days.long.labels_all_antibio %>% dplyr::select(NEW_ID, tail95_blwnli), by = "NEW_ID")
drugs.30days.long._all_antibio_bin_sum_label <- drugs.30days.long._all_antibio_bin_sum_label[!duplicated(drugs.30days.long._all_antibio_bin_sum_label), ]
table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$rowsums)
drugs.30days.long._all_antibio_bin_sum_label$cum_atb <- cut(drugs.30days.long._all_antibio_bin_sum_label$rowsums,
                                                                    breaks=c(-Inf,17,Inf), labels = c("<=25", ">25"))
table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$cum_atb)
test <- chisq.test(table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$cum_atb))
test

temp1 <- drugs.30days.long._all_antibio_bin_sum_label %>% filter(tail95_blwnli == 1) %>%
  select(rowsums)
temp2 <- drugs.30days.long._all_antibio_bin_sum_label %>% filter(tail95_blwnli == 0) %>%
  select(rowsums)

ks.test(temp1[[1]], temp2[[1]])

# 6. Number of days without antibiotics
drugs.30days.long._all_antibio_bin_sum$days_no_atb <- rowSums(drugs.30days.long._all_antibio_bin_sum == 0)

## RE ADD THE LABEL OF INTEREST
drugs.30days.long._all_antibio_bin_sum_label <- merge(drugs.30days.long._all_antibio_bin_sum, drugs.30days.long.labels_all_antibio %>% dplyr::select(NEW_ID, tail95_blwnli), by = "NEW_ID")

drugs.30days.long._all_antibio_bin_sum_label$days_no_atb.cat <- cut(drugs.30days.long._all_antibio_bin_sum_label$days_no_atb,
                                                   breaks=c(-Inf,13,Inf), labels = c("<=15", ">15"))

drugs.30days.long._all_antibio_bin_sum_label <- drugs.30days.long._all_antibio_bin_sum_label[!duplicated(drugs.30days.long._all_antibio_bin_sum_label), ]

table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$days_no_atb)
table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$days_no_atb.cat)
test <- chisq.test(table(drugs.30days.long._all_antibio_bin_sum_label$tail95_blwnli, drugs.30days.long._all_antibio_bin_sum_label$days_no_atb.cat))
test



df_window.1.4