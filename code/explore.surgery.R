
library(table1)

## Create variable with time to first closed treatment
df_window.1.4$CTSTTM01_minutes <- as.numeric((df_window.1.4$CTSTTM01%%1)*100)
df_window.1.4$CTSTTM01_hours <- as.numeric(trunc(df_window.1.4$CTSTTM01))
df_window.1.4$CTSTTM01_time2 <- sprintf("%s:%s:00", 
                                                 df_window.1.4$CTSTTM01_hours, 
                                                 as.integer(df_window.1.4$CTSTTM01_minutes))
df_window.1.4["CTSTTM01_time2"][df_window.1.4["CTSTTM01_time2"] == "NA:NA:00"] <- NA
df_window.1.4$CTSTTM01_time_combined2 <- as.POSIXct(paste(df_window.1.4$CTSTDT01_date,
                                                                   df_window.1.4$CTSTTM01_time2), 
                                                             format="%m/%d/%Y %H:%M:%S")

## Create variable with time to first operative treatment
df_window.1.4$OTSTDT01_date <- num_date(df_window.1.4$OTSTDT01)

df_window.1.4$OTSTTM01_minutes <- as.numeric((df_window.1.4$OTSTTM01%%1)*100)
df_window.1.4$OTSTTM01_hours <- as.numeric(trunc(df_window.1.4$OTSTTM01))
df_window.1.4$OTSTTM01_time2 <- sprintf("%s:%s:00", 
                                                 df_window.1.4$OTSTTM01_hours, 
                                                 as.integer(df_window.1.4$OTSTTM01_minutes))
df_window.1.4["OTSTTM01_time2"][df_window.1.4["OTSTTM01_time2"] == "NA:NA:00"] <- NA
df_window.1.4$OTSTTM01_time_combined2 <- as.POSIXct(paste(df_window.1.4$OTSTDT01_date,
                                                                   df_window.1.4$OTSTTM01_time2), 
                                                             format="%Y-%m-%d %H:%M:%S")

## Create variable with time of injury
df_window.1.4$INJTM_minutes <- as.numeric((df_window.1.4$INJTM%%1)*100)
df_window.1.4$INJTM_hours <- as.numeric(trunc(df_window.1.4$INJTM))
df_window.1.4$INJTM_time2 <- sprintf("%s:%s:00", 
                                              df_window.1.4$INJTM_hours, 
                                              as.integer(df_window.1.4$INJTM_minutes))
df_window.1.4["INJTM_time2"][df_window.1.4["INJTM_time2"] == "NA:NA:00"] <- NA
df_window.1.4$INJTM_time_combined2 <- as.POSIXct(paste(df_window.1.4$INJDT_date,
                                                                df_window.1.4$INJTM_time2), 
                                                          format="%m/%d/%Y %H:%M:%S")

df_window.1.4$CTST01_time_combined <- as.POSIXct(paste(df_window.1.4$CTSTDT01_date,
                                                                df_window.1.4$CTSTTM01_time),
                                                          format="%m/%d/%Y %H:%M:%S")

df_window.1.4$INJ_time_combined <- as.POSIXct(paste(df_window.1.4$INJDT_date,
                                                             df_window.1.4$INJTM_time),
                                                       format="%m/%d/%Y %H:%M:%S")

df_window.1.4$time_decompression_surgery_old <- difftime(df_window.1.4$CTST01_time_combined,
                                                                  df_window.1.4$INJ_time_combined, units="mins")

df_window.1.4$time_decompression_surgery <- difftime(df_window.1.4$OTSTTM01_time_combined2,
                                                              df_window.1.4$INJTM_time_combined2, units="mins")

df_window.1.4$CORDCD01 <- factor(df_window.1.4$CORDCD01)
df_window.1.4$CORDOC01 <- factor(df_window.1.4$CORDOC01)
df_window.1.4$DECOCD01 <- factor(df_window.1.4$DECOCD01)
df_window.1.4$DECOMC01 <- factor(df_window.1.4$DECOMC01)
df_window.1.4$DURACD01 <- factor(df_window.1.4$DURACD01)
df_window.1.4$DUROCD01 <- factor(df_window.1.4$DUROCD01)
df_window.1.4$EXCDCD01 <- factor(df_window.1.4$EXCDCD01)
df_window.1.4$EXCICD01 <- factor(df_window.1.4$EXCICD01)
df_window.1.4$FUSCD01 <- factor(df_window.1.4$FUSCD01)
df_window.1.4$FUSICD01 <- factor(df_window.1.4$FUSICD01)
df_window.1.4$IMMOBC01 <- factor(df_window.1.4$IMMOBC01)
df_window.1.4$INFXCD01 <- factor(df_window.1.4$INFXCD01)
df_window.1.4$INOPCD01 <- factor(df_window.1.4$INOPCD01)
df_window.1.4$INTFCD01 <- factor(df_window.1.4$INTFCD01)
df_window.1.4$INTRAC01 <- factor(df_window.1.4$INTRAC01)
df_window.1.4$LAMCD101 <- factor(df_window.1.4$LAMCD101)
df_window.1.4$time_decompression_surgery_24h <- as.numeric(df_window.1.4$time_decompression_surgery/24)
df_window.1.4 %>% 
  dplyr::mutate(time_decompression_surgery_24h_binary = case_when(time_decompression_surgery_24h <  24 ~ 0,
                             time_decompression_surgery_24h >= 24 ~ 1))
df_window.1.4$time_decompression_surgery_24h_binary <- factor(df_window.1.4$time_decompression_surgery_24h_binary)

formula <- paste0('~ CORDCD01 + CORDOC01 + DECOCD01 + DECOMC01 + DURACD01 + DUROCD01 +  EXCDCD01 + EXCICD01 + EXDICD01 + FUSCD01 + FUSICD01 + IMMOBC01 + 
                  INFXCD01 + INOPCD01 + INTFCD01 + INTRAC01 + LAMCD101 + time_decompression_surgery_24h| ', "rec_AIS_motor_imp")

table1(as.formula(formula), data = df_window.1.4)

c("CORDCD01","CORDOC01", "DECOCD01", "DECOMC01", "DURACD01", "DUROCD01", 
  "EXCDCD01", "EXCICD01", "EXDICD01", "FUSCD01", "FUSICD01", "IMMOBC01", 
  "INFXCD01", "INOPCD01", "INTFCD01", "INTRAC01", "LAMCD101", 
  "time_decompression_surgery")

df_window.1.4$time_decompression_surgery


table(raw_sygen$sevcd)
table(raw_sygen$tx1_r)
