# --------------------------
# Part 3b : Statistical analysis duplicate for second candidate model
# Author: Congyu Zhang
# --------------------------
#1 . Load the Environment and import the raw data
getRversion()

rm(list=ls())
packages <- c("ggforce","haven","ggsci","here","patchwork","ggplot2","dplyr", "rstatix", "survival",
              "survminer","ggpubr", "tidyr", "car", "broom")
lapply(packages,library,character.only = TRUE)

select <- dplyr::select

# ---------------------------------------------------------
#2.statistical analysis - parameters with clinical endpoint
# ---------------------------------------------------------
# extract individual parameters
# --- Mean PRO total score prm
prm <- read.csv(here("data", "poppk", "irt", "prm_test4.txt")) %>%
  select(id, BASE_mode, EMAX_mode, KD_mode, AGE, BHT, BWT, SEX, RACE, ACTARM,
         ETHNIC, MSTAGE, MUTSTAT, PRIORADJ, REGIONDI, SCRNLDH, tBECOG) %>%
  rename(Base = BASE_mode, Emax = EMAX_mode, Kd = KD_mode) %>%
  mutate(ACTARM = factor(case_when(
    ACTARM == "VEMURAFENIB + GDC0973" ~ "Cobimetinib and vemurafenib",
    ACTARM == "VEMURAFENIB + PLACEBO" ~ "Placebo and vemurafenib",
    .default = NA_character_
  ), levels = c("Cobimetinib and vemurafenib", "Placebo and vemurafenib")))

# boxplot to compare the parameter distribution
prm_long <- prm %>%
  select(id, ACTARM, Emax) %>%
  pivot_longer(cols = c(Emax),
               names_to = "parameter",
               values_to = "value") %>%
  mutate(Treatment = dplyr::recode(ACTARM,
          "Treatment" = "Cobimetinib and vemurafenib",
          "Placebo" = "Placebo and vemurafenib"))

p1 <- ggplot(prm_long, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ parameter) +
  scale_y_continuous(limits = c(-3,3)) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 5.5,
                     label.x = 0.7, label.y = 2.6) +
  labs(y = "Predicited Value") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16, color = "black")
  ) +
  scale_fill_jco()
print(p1)

#ggsave("plot/Emax.tiff", plot = p1, width = 5, height = 7, unit ="in", dpi = 300)

prm %>%
  group_by(ACTARM) %>%
  summarise(median_emax = median(Emax))

vars <- c("Base", "Emax", "Kd")
# ===============================================
#----- 1. ORR Best Confirmed Overall Response - Investigator(Visit time different)
# ===============================================
orr_dat <- read.csv(here("data", "ORR.csv")) %>% 
  rename(id = UNI_TRNC, ORR = AVALC) %>%
  inner_join(prm, by = "id") %>%
  select(-TRT01A, -PARAM, -PARAMCD, -AVAL, -AVISIT) %>%
  # define the responder vs. non-repsonder
  mutate(resp = ifelse(ORR %in% c("CR", "PR"), 1, 0))

orr_dat %>%
  group_by(ACTARM, ORR) %>%
  summarise(n = n(), .groups = "drop")

# -- logistic regression 
res <- do.call(rbind, lapply(vars, function(v){
  fit <- glm(as.formula(paste("resp ~", v)),
             family = binomial, data = orr_dat)
  ci <-confint(fit)[v,]
  data.frame(
    parameter  = v,
    OR = exp(coef(fit)[v]),
    OR_l = exp(ci[1]),
    OR_h = exp(ci[2]),
    p_value = coef(summary(fit))[v, "Pr(>|z|)"]
  )
}))

res_adj <- do.call(rbind, lapply(vars, function(v){
  fit <- glm(as.formula(paste("resp ~", v, "+ ACTARM")),
             family = binomial, data = orr_dat)
  ci <-confint(fit)[v,]
  data.frame(
    parameter  = v,
    OR = exp(coef(fit)[v]),
    OR_l = exp(ci[1]),
    OR_h = exp(ci[2]),
    p_value = coef(summary(fit))[v, "Pr(>|z|)"]
  )
}))

print(res, digits = 3)
print(res_adj, digits = 3)

# ===============================================
#----- 2. Survival analysis
# ===============================================
os <- read.csv(here("data", "survival.csv")) %>%
  #filter(ADY <= 600) %>%
  rename(id = UNI_TRNC, time = ADY) %>%
  mutate(event = ifelse(CNSR == 0, 1, 0), time_month = time/30,
         TRT01A = relevel(factor(TRT01A), ref = "VEMURAFENIB + PLACEBO")) %>%
  select(id, time, time_month, event, TRT01A)

# Cox - compare Baseline with OS
base_qol <- read.csv(here("data", "QoL_base.csv")) %>%
              rename(id = UNI_TRNC) %>%
              mutate(tBECOG = ifelse(is.na(BECOG), 0, BECOG))

q29 <- base_qol %>% filter(QNAME == "Q29")
q30 <- base_qol %>% filter(QNAME == "Q30")

os_q29 <- os %>% inner_join(q29, by = "id")
os_q30 <- os %>% inner_join(q30, by = "id")

# Cox model with Q29 and Q30
cox_q29 <- coxph(Surv(time_month, event) ~ QSCORE + ACTARM + BECOG, data = os_q29)
summary(cox_q29)

cox_q30 <- coxph(Surv(time_month, event) ~ QSCORE + ACTARM + BECOG, data = os_q30)
summary(cox_q30)

#=======================================================
# Cox - compare parameters with OS
#=======================================================
# keep the shared id 
os_dat <- os %>% 
        inner_join(prm, by = "id") %>%
        mutate(
          logKd_z = as.numeric(scale(log(Kd))),
          Base_z = as.numeric(scale(Base)),
          Emax_z = as.numeric(scale(Emax)),
          tBECOG = factor(tBECOG),
        )

# check the colinearility of the parameters
cor(os_dat[, c("logKd_z", "Base_z", "Emax_z")], use = "complete.obs")

#1. Multivariable analysis
# --- unadjustment analysis
cox_multi_unadj <- coxph(
  Surv(time, event) ~ Base_z + Emax_z + logKd_z, data = os_dat
)
summary(cox_multi_unadj)

# comment:
#coef exp(coef) se(coef)      z Pr(>|z|)  
#Base_z  -0.24409   0.78342  0.11733 -2.080   0.0375 *
#Emax_z  -0.12448   0.88295  0.11277 -1.104   0.2697  
#logKd_z -0.03014   0.97031  0.08868 -0.340   0.7339

cox.zph(cox_multi_unadj)

# --- adjustment analysis
cox_multi_adj <- coxph(
  Surv(time, event) ~ Base_z + Emax_z + logKd_z + ACTARM + strata(tBECOG), data = os_dat
)
summary(cox_multi_adj)

# comment:
#coef exp(coef) se(coef)      z Pr(>|z|)  
#Base_z          -0.19546   0.82245  0.12192 -1.603    0.109
#Emax_z          -0.07937   0.92370  0.11191 -0.709    0.478
#logKd_z         -0.03563   0.96500  0.08977 -0.397    0.691
#ACTARMTreatment -0.23060   0.79406  0.18541 -1.244    0.214

cox.zph(cox_multi_adj)

# --- adjustment analysis - 2 (adjust both ECOG and ACTARM)
cox_multi_adj_2 <- coxph(
  Surv(time, event) ~ Base_z + Emax_z + logKd_z + ACTARM + tBECOG, data = os_dat
)
summary(cox_multi_adj_2)

cox.zph(cox_multi_adj_2)

#2. Forest Plot
vars_cox <- c("Base_z", "Emax_z", "logKd_z")

unadj_results <- tidy(cox_multi_unadj, exponentiate = TRUE, conf.int = TRUE) %>%
          filter(term %in% vars_cox) %>%
          mutate(model = "Unadjusted")

adj_results <- tidy(cox_multi_adj, exponentiate = TRUE, conf.int = TRUE) %>%
          filter(term %in% vars_cox) %>%
          mutate(model = "Adjusted by ECOG and Treatment")

combine_results <- bind_rows(unadj_results, adj_results) %>%
                      mutate(
                        term = dplyr::recode(term, 
                                      "Base_z" = "Base",
                                      "Emax_z" = "Emax",
                                      "logKd_z" = "Kd"),
                        model = factor(model, levels = c("Unadjusted",
                                                         "Adjusted by ECOG and Treatment")))
# generate forest plot
os_forest <- ggplot(combine_results, aes(x = estimate, y = term, color = model)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                position = position_dodge(width = 0.4), width = 0.2, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  #scale_x_log10() +
  theme_bw() +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
  ) +
  scale_color_jco()

print(os_forest)

#ggsave("plot/os_forest.tiff", plot = os_forest, width = 10, height = 5, unit ="in", dpi = 300)

  
# -- km plot for inner_combine data
km_fit <- survfit(Surv(time_month, event) ~ ACTARM, data = os_dat)

km_plot <- ggsurvplot(
  km_fit,
  data = os_dat,
  conf.int = FALSE,
  legend = "none",
  legend.labs = c("Cobimetinib and vemurafenib", "Placebo and vemurafenib"),
  palette = "jco",
  pval = TRUE,
  pval.size = 5,
  
  xlab = "Time (months)",
  xlim = c(0, 72),
  break.time.by = 12,
  
  ylab = "Overall Survival (%)",
  ylim = c(0,1),
  yticklabs = c("0%", "25%", "50%", "75%", "100%"),
  surv.scale = "percent",
  
  size = 1.2,
  censor = TRUE,
  censor.shape = "+",
  censor.size = 6,
  
  font.x = c(12, "bold", "black"),
  font.y = c(12, "bold", "black"),
  font.tickslab = c(12, "plain", "black"),
  
  risk.table = FALSE,
  risk.table.y.text = FALSE,
  risk.table.height = 0.22,
  risk.table.fontsize = 5,
  tables.theme = theme_cleantable() +
    theme(axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color ="black", face = "bold")),
  
  ggtheme = theme_classic()
)

km_plot$plot <- km_plot$plot +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30", linewidth = 0.8)

print(km_plot)

tiff("plot/km_os_clean.tiff", width = 6.5, height = 4, unit ="in", res = 300)
print(km_plot)
dev.off()


#----- 3. PFS analysis
pfs <- read.csv(here("data", "PFS.csv")) %>%
  filter(!EVNTDESC == "Randomization") %>%
  rename(id = UNI_TRNC, time = ADY) %>%
  mutate(event = ifelse(CNSR ==0, 1, 0), time_month = time/30,
         TRT01A = relevel(factor(TRT01A), ref = "VEMURAFENIB + PLACEBO")) %>%
  select(id, time, time_month, event, TRT01A)

#=======================================================
# Cox - compare Baseline with PFS
#=======================================================
pfs_q29 <- pfs %>% inner_join(q29, by = "id")
pfs_q30 <- pfs %>% inner_join(q30, by = "id")

# Cox model with Q29 and Q30
cox_q29_pfs <- coxph(Surv(time_month, event) ~ QSCORE + ACTARM, data = pfs_q29)
summary(cox_q29_pfs)

cox_q30_pfs <- coxph(Surv(time_month, event) ~ QSCORE + ACTARM, data = pfs_q30)
summary(cox_q30_pfs)

#=======================================================
# Cox - compare parameters with PFS
#=======================================================
# keep the shared id 
pfs_dat <- pfs %>% 
      inner_join(prm, by = "id") %>%
      mutate(
        logKd_z = as.numeric(scale(log(Kd))),
        Base_z = as.numeric(scale(Base)),
        Emax_z = as.numeric(scale(Emax)),
        tBECOG = factor(tBECOG),
      )

# check the colinearility of the parameters
cor(pfs_dat[, c("logKd_z", "Base_z", "Emax_z")], use = "complete.obs")

#1. Multivariable analysis
# --- unadjustment analysis
cox_multi_unadj_pfs <- coxph(
  Surv(time, event) ~ Base_z + Emax_z + logKd_z, data = pfs_dat
)
summary(cox_multi_unadj_pfs)
cox.zph(cox_multi_unadj_pfs)

# --- adjustment analysis
cox_multi_adj_pfs <- coxph(
  Surv(time, event) ~ Base_z + Emax_z + logKd_z + ACTARM + strata(tBECOG), data = pfs_dat
)
summary(cox_multi_adj_pfs)
cox.zph(cox_multi_adj_pfs)

#2. Forest Plot
unadj_results_pfs <- tidy(cox_multi_unadj_pfs, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% vars_cox) %>%
  mutate(model = "Unadjusted")

adj_results_pfs <- tidy(cox_multi_adj_pfs, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term %in% vars_cox) %>%
  mutate(model = "Adjusted by ECOG and Treatment")

combine_results_pfs <- bind_rows(unadj_results_pfs, adj_results_pfs) %>%
  mutate(
    term = dplyr::recode(term, 
                         "Base_z" = "Base",
                         "Emax_z" = "Emax",
                         "logKd_z" = "Kd"),
    model = factor(model, levels = c("Unadjusted",
                                     "Adjusted by ECOG and Treatment")))
# generate forest plot
pfs_forest <- ggplot(combine_results_pfs, aes(x = estimate, y = term, color = model)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                position = position_dodge(width = 0.4), width = 0.2, linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  #scale_x_log10() +
  theme_bw() +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, color = NULL) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, color = "black"),
    strip.text = element_text(size = 15, face = "bold", color = "black")
  )

print(pfs_forest)

#ggsave("plot/os_forest.tiff", plot = os_forest, width = 10, height = 5, unit ="in", dpi = 300)

