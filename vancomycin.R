################################################################################
# SCI - Outstanding recovery - Investigate drugs in recovery group
# L. Bourguignon
# First version : 20.07.2023
# Last update : 20.07.2023
# Vancomycin
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

################################################################################
# Load data

raw_sygen <- read.csv('/Volumes/green_groups_bds_public/Data/Sygen/JohnKramersProject_DATA_2019-10-07_0111.csv')
vancomycin_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/vancomycin.csv')
acetaminophen_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/acetaminophen.csv')
acetaminophen_oxycodone_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/acetaminophen_and_oxycodone.csv')
morphine_master <- read.csv('/Volumes/green_groups_bds_public/Projects/Drug/Drugs/morphine.csv')

################################################################################
#Function preprocessing 

prepro_fct <- function(data){
  long_ID <- unique(data$NEW_ID)
  print(paste0('Number of patients which received the drug in the 1st year: ', length(long_ID)))
  
  data.30days <- data %>% select('NEW_ID':'X30')
  new_id_temp <- data.30days$NEW_ID
  generic_name_temp <- data.30days$generic_name
  # Binarise the entries (convert any value above 0 to 1)
  data.30days[data.30days > 0] <- 1
  data.30days$NEW_ID <- new_id_temp
  data.30days$generic_name <- generic_name_temp
  # Sum the number of days of data per patient across the first 30 days after SCI
  data.30days <- data.30days %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(sum = sum(c_across(3:33), na.rm = T))
  
  data.30days.received <- data.30days[data.30days$sum>0,]
  
  id_received <- unique(data.30days.received$NEW_ID)
  print(paste0('Number of patients which received the drug in the first 30 days: ', length(id_received)))
  
  sygen_copy <- raw_sygen
  sygen_copy <- sygen_copy %>% 
    mutate(received = case_when(ptid %in% id_received ~ 1, 
                             TRUE ~ 0))
  
  subset_received_sygen <- raw_sygen[raw_sygen$ptid %in% id_received, ]
  subset_notreceived_sygen <- raw_sygen[!(raw_sygen$ptid %in% id_received), ]
  
  print('--Baseline--')
  print('AIS distribution for patients with treatment')
  print('Number')
  print(table(subset_received_sygen$ais1))
  print('Proportions')
  print(table(subset_received_sygen$ais1)/dim(subset_received_sygen)[1])
  print('AIS distribution for patients with no treatment')
  print('Number')
  print(table(subset_notreceived_sygen$ais1))
  print('Proportions')
  print(table(subset_notreceived_sygen$ais1)/dim(subset_notreceived_sygen)[1])
  
  print('--Recovery--')
  print('AIS distribution for patients with treatment')
  print('Number')
  print(table(subset_received_sygen$ais52))
  print('Proportions')
  print(table(subset_received_sygen$ais52)/dim(subset_received_sygen)[1])
  print('AIS distribution for patients with no treatment')
  print('Number')
  print(table(subset_notreceived_sygen$ais52))
  print('Proportions')
  print(table(subset_notreceived_sygen$ais52)/dim(subset_notreceived_sygen)[1])
  
  return(list(subset_received_sygen, subset_notreceived_sygen))
}

################################################################################
# Vancomycin
output.vanco <- prepro_fct(vancomycin_master)
received.vanco <- output.vanco[1][[1]]
notreceived.vanco <- output.vanco[2][[1]]

prop.test(x= as.vector(table(notreceived.vanco$ais1))[2:5], 
          n= as.vector(table(raw_sygen$ais1))[2:5])
# The proportion of patients receiving vancomycin is significantly different 
# in the different AIS grade categories
# p-value = 6.336e-05

chisq.test(as.vector(table(notreceived.vanco$ais1))[2:5], 
           p = as.vector(table(received.vanco$ais1))[2:5], 
           rescale.p = TRUE)
# Statistically significant difference is proportion of AIS grades between 
# receiver and none-receiver of vancomycin
# p-value < 2.2e-16

################################################################################
# Acetaminophen

output.acetaminophen <- prepro_fct(acetaminophen_master)
received.acetaminophen <- output.acetaminophen[1][[1]]
notreceived.acetaminophen <- output.acetaminophen[2][[1]]

prop.test(x= as.vector(table(notreceived.acetaminophen$ais1))[2:5], 
          n= as.vector(table(raw_sygen$ais1))[2:5])
# The proportion of patients receiving acetaminophen is NOT significantly different 
# in the different AIS grade categories
# p-value = 0.0516

chisq.test(as.vector(table(notreceived.acetaminophen$ais1))[2:5], 
           p = as.vector(table(received.acetaminophen$ais1))[2:5], 
           rescale.p = TRUE)
# Statistically significant difference is proportion of AIS grades between 
# receiver and none-receiver of acetaminophen
# p-value = 0.003445

################################################################################
# Acetaminophen + oxycodone

output.acetaminophen.oxycodone <- prepro_fct(acetaminophen_oxycodone_master)
received.acetaminophen.oxycodone <- output.acetaminophen.oxycodone[1][[1]]
notreceived.acetaminophen.oxycodone <- output.acetaminophen.oxycodone[2][[1]]

prop.test(x= as.vector(table(notreceived.acetaminophen.oxycodone$ais1))[2:5], 
          n= as.vector(table(raw_sygen$ais1))[2:5])
# The proportion of patients receiving acetaminophen + oxycodone is significantly different 
# in the different AIS grade categories (trend towards being significant)
# p-value = 0.04786

chisq.test(as.vector(table(notreceived.acetaminophen.oxycodone$ais1))[2:5], 
           p = as.vector(table(received.acetaminophen.oxycodone$ais1))[2:5], 
           rescale.p = TRUE)
# Statistically significant difference is proportion of AIS grades between 
# receiver and none-receiver of acetaminophen + oxycodone
# p-value = 0.0001635

################################################################################
# Morphine

output.morphine <- prepro_fct(morphine_master)
received.morphine <- output.morphine[1][[1]]
notreceived.morphine <- output.morphine[2][[1]]

prop.test(x= as.vector(table(notreceived.morphine$ais1))[2:5], 
          n= as.vector(table(raw_sygen$ais1))[2:5])
# The proportion of patients receiving morphine is significantly different 
# in the different AIS grade categories 
# --> increased proportion of patients with no morphine in the first 30 days
# from AIS A to D
# p-value = 0.005107

chisq.test(as.vector(table(notreceived.morphine$ais1))[1:4], 
           p = as.vector(table(notreceived.morphine$ais1))[1:4], 
           rescale.p = TRUE)
# NO statistically significant difference is proportion of AIS grades between 
# receiver and none-receiver of morphine
# p-value = 1

################################################################################
################################################################################
lems52_sygen_raw <- ggplot(rbind(received.vanco, notreceived.vanco), 
                           aes(x=lower52, group=received, colour=received)) + 
  geom_density(alpha=.2) + xlim(0,50)

sygen_copy_lems <- sygen_copy %>% select(ptid,
                                         received,
                                         ais1,
                                         lower01,
                                         lower04,
                                         lower08,
                                         lower16,
                                         lower26,
                                         lower52)

sygen_copy_lems_long <- tidyr::gather(sygen_copy_lems, time, value, 
                                      lower01:lower52, factor_key=TRUE)

plot <- ggplot(sygen_copy_lems_long, aes(x=time, y=value, group=ptid, colour=factor(vanco))) + 
  facet_wrap(~ais1, nrow=2) +
  geom_jitter() +
  geom_line(size=1) +
  stat_summary(aes(y=value, group=vanco), fun.data=mean_se, geom = "ribbon", alpha = 0.6) +
  stat_summary(aes(y=value, group=vanco), fun=mean, geom="line", size = 1) +
  ylim(0, 55) +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold", margin = margin(r=18)),
        plot.background = element_rect(fill = "#FFFFFF")) +
  xlab('Exam Stage')
plot


slope <- function(df)
{
  df = df[complete.cases(df),] #remove rows (= 1 time point for 1 individual for 1 blood marker) with missing values
  remove_AIS = c('E', '', 'AIS E', 'ND')
  df = df[!df$ais1 %in% remove_AIS,] #filter out patients with AIS E or missing grades at baseline (Sygen) or acute 1 stage (Murnau)
  
  #df$variable = as.numeric(as.character(df$variable))
  
  test <- lmer(value ~ received*ais1 + variable + age + level + sexcd + (1|ptid), data = df, na.action = na.omit) # fit ANOVA test
  m.lst <- lstrends(test, "received", var = "ais1")
  
  return(summary(pairs(m.lst)))
}
library(lme4)
library(lsmeans)
df <- sygen_copy %>% select(ptid,
                            sexcd,
                            age,
                            splvl,
                            received,
                            ais1,
                            lower01,
                            lower04,
                            lower08,
                            lower16,
                            lower26,
                            lower52)

df <- df%>% 
  mutate(level = if_else(str_detect(df$splvl, "T"), "T", "C"))

df_long <- melt(df, id=c('ptid', 'sexcd', 'age', 'splvl', 'level',
                         'received', 'ais1'))

df_long$variable<- as.numeric(str_remove(df_long$variable, "lower"))

temp <- slope(df_long)

library(car)

test <- Anova(lm1 <- lmer(value ~ received*ais1 + variable + age + level + sexcd + (1|ptid), 
                          data = df, na.action = na.omit), 
              type="III")
test
