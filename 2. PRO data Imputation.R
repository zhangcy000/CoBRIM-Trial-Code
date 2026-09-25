# --------------------------
# Part 2 : PRO data imputation for POPPK model
# Author: Congyu Zhang
# --------------------------
#1 . Load the Environment and import the raw data
getRversion()

rm(list=ls())
packages <- c("dplyr","haven","tidyr","here","gtsummary","stringr","ggplot2", "purrr")
lapply(packages,library,character.only = TRUE)

pro <- read.csv(here("data","all_pro.csv")) %>%
      select(-UNI_ID, -QSCAT)

# ----- demographics information
demo_table <- pro %>%
  group_by(UNI_TRNC) %>%
  slice(1) %>%
  ungroup() %>%
  select(ACTARM, AGE, SEX, RACE, ETHNIC, REGIONDI, AGEGRP, AGE65, BWT, BHT, BECOG, MUTSTAT,
         MSTAGE,PRIORADJ,SCRNLDH) %>%
  tbl_summary(
    by = ACTARM,
    type = list(
      c(AGE, BWT, BHT) ~ "continuous",
      c(SEX, RACE, ETHNIC, REGIONDI, AGEGRP, AGE65,BECOG, MUTSTAT,
        MSTAGE,PRIORADJ,SCRNLDH) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "median = {median} ({min}, {max})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits =all_continuous() ~ 2
  )

print(demo_table)

demo_table2 <- pro %>%
  group_by(UNI_TRNC) %>%
  slice(1) %>%
  ungroup() %>%
  select(ACTARM, AGE, SEX, RACE, ETHNIC, REGIONDI, AGEGRP, AGE65, BWT, BHT, BECOG, MUTSTAT,
         MSTAGE,PRIORADJ,SCRNLDH) %>%
  tbl_summary(
    type = list(
      c(AGE, BWT, BHT) ~ "continuous",
      c(SEX, RACE, ETHNIC, REGIONDI, AGEGRP, AGE65,BECOG, MUTSTAT,
        MSTAGE,PRIORADJ,SCRNLDH) ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "median = {median} ({min}, {max})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits =all_continuous() ~ 2
  )
print(demo_table2)
 
# check the organization of the questionnaire (complete datasets (30 questions at each time))
count_id_time <- pro %>%
  group_by(UNI_TRNC,QSDY) %>%
  summarise(
    count = n()
  )

pro %>%
  distinct(QSTESTCD,QSTEST) %>%
  arrange(QSTESTCD)

sum(is.na(pro$QSSTRESN)) # 0 missing

pro %>%
  group_by(UNI_TRNC,QSTESTCD,QSDY) %>%
  filter(n() >1)  #no duplicate checked

# check whether VISIT are coded same for all questions
visit_check <- pro %>%
  group_by(UNI_TRNC,QSDY) %>%
  summarise(
    N_visit = n_distinct(VISIT),
    .groups = "drop",
  ) %>%
  filter(N_visit > 1) 

#2. Calculate each sub-group item scores
# transform the dataset to wider format
demo_col <- c("UNI_TRNC", "VISIT", "ACTARM", "AGE", "SEX", "RACE", "ETHNIC", "REGIONDI", "AGEGRP", "AGE65", 
              "BWT", "BHT", "BECOG", "MUTSTAT", "MSTAGE", "PRIORADJ", "SCRNLDH")

pro_wide <- pro %>%
          filter(!(UNI_TRNC == "8f1f25a0dd" & QSDY == 155)) %>%
          # change QSDY = -1 to 1
          mutate(QSDY = ifelse(QSDY == -1, 1, QSDY)) %>%
          mutate(item_name = paste0("Q", as.integer(str_extract(QSTESTCD, "\\d{1,2}$")))) %>%
          select(all_of(demo_col), QSDY, item_name, QSSTRESN) %>%
          pivot_wider(
            names_from = item_name,
            values_from = QSSTRESN
          )

# sum-up the item scores
pro_score <- pro_wide %>%
      mutate(
        # QoL
        QL2 = rowSums(cbind(Q29,Q30)),
        
        # Functional scales
        PF2 = rowSums(cbind(Q1,Q2,Q3,Q4,Q5)),
        RF2 = rowSums(cbind(Q6,Q7)),
        EF = rowSums(cbind(Q21,Q22,Q23,Q24)),
        CF = rowSums(cbind(Q20,Q25)),
        SF = rowSums(cbind(Q26,Q27)),
        
        # Symptom scales/items
        FA = rowSums(cbind(Q10,Q12,Q18)),
        NV = rowSums(cbind(Q14,Q15)),
        PA = rowSums(cbind(Q9,Q19)),
        DY = Q8,
        SL = Q11,
        AP = Q13,
        CO = Q16,
        DI = Q17,
        FI = Q28,
        
        # total score
        FUN_SUM = rowSums(cbind(PF2, RF2, EF, CF, SF)),
        SYMP_SUM = rowSums(cbind(FA, NV, PA, DY, SL, AP, CO, DI, FI))      
        ) %>%
      arrange(UNI_TRNC,QSDY) %>%
      group_by(UNI_TRNC) %>%
      mutate(T0 = min(QSDY),TAFR = QSDY - T0, QOL2 = QL2 - 2,
             FUN_SUM2 = FUN_SUM - 15, SYMP_SUM2 = SYMP_SUM - 13) %>%
      select(UNI_TRNC, QSDY, TAFR, QL2, QOL2, FUN_SUM, FUN_SUM2, SYMP_SUM,SYMP_SUM2,
             PF2, RF2, EF, CF, SF, FA, NV, PA, DY, SL, AP, CO, DI, FI,
             everything())

# ===============================================================
# 1. QOL data only - generate the CSV file for POPPK model
# ===============================================================
# filter the global health status (QoL) PRO and transform to pivot_longer format
pro_ghs_longer <- pro_score %>%
      rename(QoL = QL2) %>%
      pivot_longer(
            cols = c(Q29,Q30),
            names_to = "QNAME",
            values_to = "QSCORE"
          ) %>%
      #mutate(QSCORE = if_else(QSCORE == 1, 2, QSCORE)) %>%
      select(UNI_TRNC, QSDY, TAFR, QNAME, QSCORE, all_of(demo_col))

ghs_filter <- pro_ghs_longer %>%
        filter(VISIT %in% c("CYCLE 1 DAY 1", "CYCLE 1 DAY 15", "CYCLE 2 DAY 1", "CYCLE 2 DAY 15",
                      "CYCLE 4 DAY 1", "CYCLE 6 DAY 1", "CYCLE 8 DAY 1"))

# Baseline dataset
pro_ghs_base <- pro_ghs_longer %>%
              filter(TAFR == 0)

ggplot(pro_ghs_base, aes(x = factor(QSCORE))) +
    geom_bar() +
    facet_wrap(~QNAME) +
    theme_bw()

#write.csv(pro_ghs_base, here("data","QoL_base_6grp.csv"),row.names = FALSE)
#write.csv(pro_ghs_longer, here("data","QoL_full.csv"),row.names = FALSE)
#write.csv(ghs_filter, here("data","QoL_full_filter.csv"),row.names = FALSE)

# check the VISIT and QSDY
pro_ghs_longer %>%
  group_by(VISIT) %>%
  summarise(day = paste(unique(QSDY), collapse = ", "))

# ===============================================================
# 2. QOL data only - EDA 
# ===============================================================
# linear transformation
pro_qol_combine <- pro_score %>%
              select(UNI_TRNC, QSDY, TAFR, QL2, all_of(demo_col)) %>%
              mutate(QoL_linear = round(((QL2/2 - 1)/6)*100, 2)) %>%
          group_by(UNI_TRNC) %>%
          mutate(
            base = QL2[which.min(QSDY)],
            base_change = QL2 - base,
            base_linear = QoL_linear[which.min(QSDY)],
            base_change_linear = QoL_linear - base_linear,
          ) %>%
          ungroup %>%
          select(UNI_TRNC, QSDY, TAFR, QL2, base, base_change, QoL_linear, base_linear, base_change_linear, everything()) %>%
          # rename the columns
          rename(ID = UNI_TRNC,TIME = QSDY) %>%
          rename_with(toupper)

# ============================================================================
# generation id, time file for Simulx
qol_export <- pro_qol_combine %>%
          select(ID, TIME, VISIT) %>%
          rename(id = ID, time = TIME)

qol_time <- qol_export %>% select(-VISIT)

#write.csv(qol_export, here("data","qol_original.csv"), row.names = FALSE)
#write.csv(qol_time, here("data","qol_time.txt"), row.names = FALSE)
# ============================================================================

# check the number of observation of each patient
count_obs <- pro_qol_combine %>%
          group_by(ID) %>%
          summarise(
            count = n()
          ) %>%
          ungroup()

table(count_obs$count)  
  
# ===============================================
# ---- Generate Plot
# Plot by Visit - from cycle 1 to cycle 8
# ===============================================
table(pro_qol_combine$VISIT)

qol_visit_plot <- pro_qol_combine %>%
          filter(VISIT %in% c("CYCLE 1 DAY 1", "CYCLE 1 DAY 15", "CYCLE 2 DAY 1", "CYCLE 2 DAY 15",
                              "CYCLE 4 DAY 1", "CYCLE 6 DAY 1", "CYCLE 8 DAY 1")) %>%
          select(-TIME, -TAFR) %>%
          mutate(Day = case_when(
                  grepl("^CYCLE 1 DAY 1$",VISIT,ignore.case = TRUE) ~ "0",
                  grepl("^CYCLE 1 DAY 15$",VISIT,ignore.case = TRUE) ~ "15",
                  grepl("^CYCLE 2 DAY 1$",VISIT,ignore.case = TRUE) ~ "29",
                  grepl("^CYCLE 2 DAY 15$",VISIT,ignore.case = TRUE) ~ "43",
                  grepl("^CYCLE 4 DAY 1$",VISIT,ignore.case = TRUE) ~ "85",
                  grepl("^CYCLE 6 DAY 1$",VISIT,ignore.case = TRUE) ~ "141",
                  grepl("^CYCLE 8 DAY 1$",VISIT,ignore.case = TRUE) ~ "197",
                  TRUE ~ "missing"
                ),Day = as.numeric(Day)) %>%
          select(ID, Day, everything())

# --- statistical analysis comparing two arms (CFB only)
qol_plot_mean <- qol_visit_plot %>%
  group_by(ACTARM, VISIT, Day) %>%
  summarise(
    mean_change = mean(BASE_CHANGE_LINEAR, na.rm = TRUE),
    se_change = sd(BASE_CHANGE_LINEAR, na.rm = TRUE)/sqrt(n()),
    sd = sd(BASE_CHANGE_LINEAR, na.rm = TRUE),
    pt_num = n(),
    .groups = "drop"
  ) %>%
  mutate(VISIT = str_replace_all(VISIT, "\n", ""))

qol_pvalue <- qol_visit_plot %>%
  filter(VISIT != "CYCLE 1 DAY 1") %>%
  group_by(VISIT, Day) %>%
  summarise(
    p_value_t = t.test(BASE_CHANGE_LINEAR ~ ACTARM)$p.value,
    p_value_wil = wilcox.test(BASE_CHANGE_LINEAR ~ ACTARM)$p.value,
    .groups = "drop"
  )

stat_qol <- qol_plot_mean %>%
          filter(Day != "0") %>%
          select(ACTARM, VISIT, Day, mean_change, sd) %>%
          pivot_wider(names_from = ACTARM, values_from = c(mean_change, sd)) %>%
          rename(mean_treatment = 'mean_change_VEMURAFENIB + GDC0973',
                 sd_treatment = 'sd_VEMURAFENIB + GDC0973',
                 mean_placebo = 'mean_change_VEMURAFENIB + PLACEBO',
                 sd_placebo = 'sd_VEMURAFENIB + PLACEBO') %>%
          left_join(qol_pvalue, by = c("VISIT","Day")) %>%
          mutate(p.adj_t = p.adjust(p_value_t, method = "BH"),
                 p.adj_wil = p.adjust(p_value_wil, method = "BH"))

#write.csv(stat_qol, here("data","visit_comparison.csv"), row.names = FALSE)

# ggplot
qol_plot_mean <- qol_plot_mean %>%
          mutate(VISIT = str_replace(VISIT, " DAY", "\nDay"))  %>%
          mutate(Treatment = dplyr::recode(ACTARM,
          "VEMURAFENIB + GDC0973" = "Cobimetinib and vemurafenib",
          "VEMURAFENIB + PLACEBO" = "Placebo and vemurafenib"))

p1 <- ggplot(qol_plot_mean, aes(x = Day, y = mean_change, color = Treatment, 
                                group = Treatment, shape = Treatment)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", alpha= 0.8) +
  geom_hline(yintercept = -10, linetype = "dashed", color = "black", alpha= 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = mean_change - se_change,
                    ymax = mean_change + se_change),
                width = 6, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  #geom_text(data = qol_plot_mean %>%
  #          filter(ACTARM == "VEMURAFENIB + GDC0973"),
  #          aes(label = VISIT),
  #          vjust = -1, size = 4, color = "black") +
  scale_x_continuous(breaks = seq(0,220,30),limits = c(0,220),name = "Time (days)") +
  scale_y_continuous(breaks = seq(-20,20,10),limits = c(-20,20),name = "Change from Baseline, mean (SE)") +
  scale_color_manual(values = c("#0099CC","#E74D3D")) +
  scale_shape_manual(values = c(15,16)) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank()
  )

print(p1)

#ggsave(here("plot", "1.2 QoL by Visit_se.tiff"), p1, width = 12, height = 6, unit = "in", dpi = 300)

# ===============================================
# --- statistical analysis comparing two arms (raw data and CFB)
# ===============================================
cycle8_qol <- qol_visit_plot %>%
            filter(VISIT == "CYCLE 8 DAY 1")

# ----- baseline comparison  
cycle8_qol %>%
  group_by(ACTARM) %>%
  summarise(
    mean = mean(BASE_LINEAR),
    median = median(BASE_LINEAR),
    sd = sd(BASE_LINEAR),
    se = sd/sqrt(n())
  ) 

#ACTARM                 mean median    sd    se
#1 VEMURAFENIB + GDC0973  69.3   66.7  18.6  2.47
#2 VEMURAFENIB + PLACEBO  74.8   75    20.4  3.18

t_base <- t.test(BASE_LINEAR ~ ACTARM, data = cycle8_qol, conf.level = 0.95)

diff(t_base$estimate)
#5.498759

t_base$conf.int
#[1] -13.510542, 2.513023

t_base$p.value
#[1] 0.1758673

cycle8_qol %>%
  group_by(ACTARM) %>%
  summarise(
    median = median(BASE_LINEAR),
    q25 = quantile(BASE_LINEAR, 0.25),
    q75 = quantile(BASE_LINEAR, 0.75),
    IQR = q75 - q25
  ) 

#ACTARM                median   q25   q75   IQR
#1 VEMURAFENIB + GDC0973   66.7  58.3  83.3    25
#2 VEMURAFENIB + PLACEBO   75    66.7  91.7    25

wilcox_base <- wilcox.test(BASE_LINEAR ~ ACTARM, data = cycle8_qol, conf.int = TRUE,
                           conf.level = 0.95, exact = FALSE)

wilcox_base$estimate
# -8.329927 
wilcox_base$conf.int
#-1.666007e+01  3.790392e-05
wilcox_base$p.value
#[1] 0.1195014

# ----- Cycle 8 day 1 comparison  
cycle8_qol %>%
  group_by(ACTARM) %>%
  summarise(
    mean_qol = mean(QOL_LINEAR),
    sd_qol = sd(QOL_LINEAR),
    se_qol = sd_qol/sqrt(n()),
    
    mean_cfb = mean(BASE_CHANGE_LINEAR),
    sd_cfb = sd(BASE_CHANGE_LINEAR),
    se_cfb = sd_cfb/sqrt(n()),
  ) 

#ACTARM                   mean_qol sd_qol se_qol mean_cfb sd_cfb se_cfb
#1 VEMURAFENIB + GDC0973     68.9   18.7   2.47   -0.438   17.3   2.29
#2 VEMURAFENIB + PLACEBO     67.1   19.2   3.00   -7.72    23.0   3.59

# ----- Observed: QOL_linear
t_qol <- t.test(QOL_LINEAR ~ ACTARM, data = cycle8_qol, conf.level = 0.95)

# t_qol$estimate
# #mean in group VEMURAFENIB + GDC0973 mean in group VEMURAFENIB + PLACEBO 
# #68.85982                            67.07366 
# 
# diff(t_qol$estimate)
# #-1.786166

t_qol$conf.int
#[1] -5.936650  9.508982
t_qol$p.value
#[1] 0.6467911

cycle8_qol %>%
  group_by(ACTARM) %>%
  summarise(
    median_qol = median(QOL_LINEAR),
    q25_qol = quantile(QOL_LINEAR, 0.25),
    q75_qol = quantile(QOL_LINEAR, 0.75),
    IQR_qol = q75_qol - q25_qol
  ) 

#ACTARM                median_qol q25_qol q75_qol IQR_qol
#1 VEMURAFENIB + GDC0973       66.7      50    83.3    33.3
#2 VEMURAFENIB + PLACEBO       66.7      50    83.3    33.3

wilcox_qol <- wilcox.test(QOL_LINEAR ~ ACTARM, data = cycle8_qol, conf.int = TRUE,
                           conf.level = 0.95, exact = FALSE)

wilcox_qol$estimate
# 5.625924e-06

wilcox_qol$conf.int
#[1] -8.329979  8.330090

wilcox_qol$p.value
#[1] 0.7016522

# ----- Observed: BASE_CHANGE_LINEAR
t_cfb <- t.test(BASE_CHANGE_LINEAR ~ ACTARM, data = cycle8_qol, conf.level = 0.95)

t_cfb$estimate
#mean in group VEMURAFENIB + GDC0973 mean in group VEMURAFENIB + PLACEBO 
#-0.4382456                          -7.723170

diff(t_cfb$estimate)
#-7.284925

t_cfb$conf.int
#[1] -1.209228 15.779078

t_cfb$p.value
#[1] 0.09161293

cycle8_qol %>%
  group_by(ACTARM) %>%
  summarise(
    median_cfb = median(BASE_CHANGE_LINEAR),
    q25_cfb = quantile(BASE_CHANGE_LINEAR, 0.25),
    q75_cfb = quantile(BASE_CHANGE_LINEAR, 0.75),
    IQR_cfb = q75_cfb - q25_cfb
  ) 

wilcox_cfb <- wilcox.test(BASE_CHANGE_LINEAR ~ ACTARM, data = cycle8_qol, conf.int = TRUE,
                          conf.level = 0.95, exact = FALSE)

wilcox_cfb$estimate
# 8.330023 

wilcox_cfb$conf.int
#[1] -6.601669e-05  1.666998e+01

wilcox_cfb$p.value
#[1] 0.07130368








