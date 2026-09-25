# --------------------------
# Part 1 : SAS Raw data set Clean 
# Author: Congyu Zhang
# --------------------------
#1 . Load the Environment 
getRversion()

rm(list=ls())
packages <- c("dplyr","haven","tidyr","here")
lapply(packages,library,character.only = TRUE)

#2. Data exploration and clean
# --------------------------
# subject characteristic
# --------------------------
adsl2 <- read_sas(here("SAS raw data", data_file = 'adsl2.sas7bdat')) 

demo <- adsl2 %>%
  select(UNI_TRNC,AGE,SEX,RACE,ETHNIC,REGIONDI,ACTARM,AGEGRP,AGE65,BWT,BHT,BBMI,BECOG,PROFL,MUTSTAT,
         MSTAGE,PRIORADJ,SCRNLDH)

# --------------------------
# PRO data clean
# --------------------------
qs <- read_sas(here("SAS raw data","qs.sas7bdat")) # Key PRO data

pro <- qs %>% 
  filter(QSCAT == "EORTC QLQ-C30 V3", QSTESTCD != "QSALL") %>%  #Remove QSALL(Entire Questionnaire was not done)
  select(UNI_TRNC,QSCAT,QSTESTCD,QSTEST,QSSTRESN,QSBLFL,VISITDY,QSDY,VISIT,UNI_ID) %>%
  # remove follow-up period - treatment changed
  filter(!VISIT %in% c("FOLLOW-UP 4 WEEKS","FOLLOW-UP 12 WEEKS")) %>%
  # remove mislabeled VISIT
  filter(!VISIT %in% c("Cycle greater than 18 \u0096 Day 1")) %>%
  # remove missing observation time (one patient with PRO but no QADY data)
  filter(QSDY != "NA") %>%
  #remove duplicate rows: with same study time and QSTEST
  distinct(UNI_TRNC,QSTESTCD,QSDY,.keep_all = TRUE)  

# left_join demographics and pro data
pro_final <- pro %>%
          left_join(demo,by = "UNI_TRNC")

# Export the PRO dataset
#write.csv(pro_final, here("data","all_pro.csv"),row.names = FALSE)

# --------------------------
# Tumor response data clean
# --------------------------
adrs <- read_sas(here("SAS raw data", 'adrs.sas7bdat')) # tumor RECIST response

# Tumor response data
tumor <- adrs %>%
  # remove crossover patients and patients received actual treatment
  filter(CROSSFL == "N", ACTARM != "") %>%
  # unit for WT (kg) HT (cm) BBMI (kg/m^2)
  select(UNI_TRNC,TRT01A,PARAM,PARAMCD,AVAL,AVALC,ADY,AVISIT) 

# Best Confirmed Overall Response - Investigator without missing data
ORR <- tumor %>%
  filter(PARAMCD == "BESRSPI", AVAL != "NA") 

ORR %>%
  group_by(UNI_TRNC,PARAM,AVAL) %>%
  filter(n() >1)  #no duplicate checked

# Export the ORR dataset
# write.csv(ORR, here("data","ORR.csv"),row.names = FALSE)

# --------------------------
# Treatment duration data clean
# --------------------------
# Treatment Duration of each Drug and Reason for Discontinuation (199 patients)
treatment_time <- adrs %>%
  # remove crossover patients and patients without actual treatment
  filter(CROSSFL == "N", ACTARM != "") %>%
  # unit for WT (kg) HT (cm) BBMI (kg/m^2)
  select(UNI_TRNC,ACTARM,TRTSDTM,TRTEDTM,TRTDUR,DISCSTUD,STDSSDT,STDDRS,DISCAE,
         VEMSDTM,VEMEDTM,VEMDUR,DISCVEM,VEMDRS,VEMSSDT,MEKSDTM,MEKEDTM,MEKDUR,
         DISCMEK,MEKDRS,MEKSSDT) %>%
  distinct(UNI_TRNC,.keep_all = TRUE)

colSums(is.na(treatment_time))
# 2 patients missing cobimetinib discontinuation date

# Export the Treatment time dataset
# write.csv(treatment_time, here("data","treatment_time.csv"),row.names = FALSE)


# --------------------------
# TTE data clean
# --------------------------
adte <- read_sas(here("SAS raw data", "adte.sas7bdat")) # TTE data

length(unique(adte$UNI_ID))

TTE <- adte %>%
  # remove crossover patients and patients received actual treatment
  filter(CROSSFL == "N", ACTARM != "") %>%
  # unit for WT (kg) HT (cm) BBMI (kg/m^2)
  select(UNI_TRNC,TRT01A,PARAM,PARAMCD,EVNTDESC,CNSR,ADY) 

# survival (OS) data - two BECOG missing
survival <- TTE %>%
  filter(PARAMCD == "OS")

# PFS data - Earliest Contributing Event to Investigator PFS - two BECOG missing
PFS <- TTE %>%
  filter(PARAMCD == "PFSINV") 

# Export the TTE data
# write.csv(survival, here("data","survival.csv"),row.names = FALSE)
# write.csv(PFS, here("data","PFS.csv"),row.names = FALSE)

# --------------------------
# Tumor size data
# --------------------------
tr <- read_sas(here("SAS raw data", data_file = 'tr.sas7bdat')) # tumor results

tr_clean <- tr %>%
        select(UNI_TRNC,everything(),-STUDYID,-DOMAIN,-UNI_ID,-TRSPID)

# Export the tumor size data
#write.csv(tr_clean, here("data","tumor_size.csv"),row.names = FALSE)
