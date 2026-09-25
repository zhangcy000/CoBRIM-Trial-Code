# --------------------------
# Part : Simulation 
# Author: Congyu Zhang
# --------------------------
# Load the Environment 
getRversion()
#[1] ‘4.4.0’
rm(list=ls())

packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "here", "ggh4x")
lapply(packages,library,character.only = TRUE)

# ============================================================
# Part I - Simulation from Cycle 1 day 1 to Cycle 8 day 1 (original time point)
# ============================================================
#1. Import the original data and simulated data
qol <- read.csv(here("data","QoL_full.csv"))

qol_sim <- read.csv(here("data","sim_qol.csv")) %>%
            select(-ID, -UID) 

#3. Data Interpretation 
# First step is to filter out patient ID at each time point
visit_time_map <- tibble(
  VISIT = c("CYCLE 1 DAY 1", "CYCLE 1 DAY 15", "CYCLE 2 DAY 1", "CYCLE 2 DAY 15",
            "CYCLE 4 DAY 1", "CYCLE 6 DAY 1", "CYCLE 8 DAY 1")
)

id_by_visit <- qol %>%
        filter(VISIT %in% visit_time_map$VISIT) %>%
        distinct(UNI_TRNC, VISIT, QSDY) %>%
        rename(TIME = QSDY)

qol_filter <- qol_sim %>%
            inner_join(id_by_visit, by = c("original_id" = "UNI_TRNC", "TIME" = "TIME")) 

# calculate the average qol score of each patient
qol_avg <- qol_filter %>%
            group_by(original_id, TIME, VISIT, obsid, ACTARM, BECOG) %>%
            summarise(obs_mean = mean(obs, na.rm = TRUE), .groups = "drop") %>%
            mutate(QNAME = dplyr::recode(obsid, "y1" = "Q29", "y2" = "Q30"),
                   Treatment = dplyr::recode(ACTARM, "VEMURAFENIB + PLACEBO" = "Placebo and vemurafenib", 
                                              "VEMURAFENIB + GDC0973" = "Cobimetinib and vemurafenib")) %>%
            select(-obsid, -ACTARM) %>%
    # transform to wider formation
            pivot_wider(names_from = QNAME, values_from = obs_mean) %>%
    # linear transformation of QoL 
            mutate(raw_score = (Q29 + Q30)/2,
                   QoL = round((((raw_score - 1)/6) * 100),2)) %>%
            group_by(original_id) %>%
            mutate(
              base = QoL[which.min(TIME)],
              base_change = QoL - base
            ) %>%
            rename_with(toupper)


# Generate the Plot as raw data (cycle 1 day 1 to cycle 8 day 1)
qol_plot <- qol_avg %>%
          select(-TIME) %>%
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
          select(ORIGINAL_ID, Day, everything())

# --- statistical analysis comparing two arms (CFB only)
qol_plot_mean <- qol_plot %>%
  group_by(Day, TREATMENT, VISIT) %>%
  summarise(
    mean_change = mean(BASE_CHANGE, na.rm = TRUE),
    se = sd(BASE_CHANGE, na.rm = TRUE)/sqrt(n()),
    sd = sd(BASE_CHANGE, na.rm = TRUE),
    pt_num = n(),
    .groups = "drop"
  )

qol_pvalue <- qol_plot %>%
  filter(Day != "0") %>%
  group_by(VISIT) %>%
  summarise(
    p_value_t = t.test(BASE_CHANGE ~ TREATMENT)$p.value,
    p_value_wil = wilcox.test(BASE_CHANGE ~ TREATMENT)$p.value,
    .groups = "drop"
  )

stat_qol <- qol_plot_mean %>%
  filter(Day != "0") %>%
  select(TREATMENT, VISIT, Day, mean_change, sd) %>%
  pivot_wider(names_from = TREATMENT, values_from = c(mean_change, sd)) %>%
  left_join(qol_pvalue, by = c("VISIT")) %>%
  mutate(p.adj_t = p.adjust(p_value_t, method = "BH"),
         p.adj_wil = p.adjust(p_value_wil, method = "BH"))

#write.csv(stat_qol, here("data","visit_comparison_sim_rawtime.csv"), row.names = FALSE)

p1 <- ggplot(qol_plot_mean, aes(x = Day, y = mean_change, color = TREATMENT, 
                                group = TREATMENT, shape = TREATMENT)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", alpha= 0.8) +
  geom_hline(yintercept = -10, linetype = "dashed", color = "black", alpha= 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = mean_change - se,
                    ymax = mean_change + se),
                width = 6, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0,220,30),limits = c(0,220),name = "Time (days)") +
  scale_y_continuous(breaks = seq(-10,10,10),limits = c(-10,10),name = "Change from Baseline, mean (SE)") +
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

#ggsave(here("plot", "QoL by Visit_se_sim_rawtime.tiff"), p1, width = 12, height = 6, unit = "in", dpi = 300)

# --- statistical analysis comparing two arms (raw data and CFB)
cycle8_qol <- qol_plot %>%
  filter(VISIT == "CYCLE 8 DAY 1")

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(shapiro_p = shapiro.test(BASE_CHANGE)$p.value)

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(shapiro_p = shapiro.test(QOL)$p.value)

# ----- baseline comparison  
cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    mean = mean(BASE),
    sd = sd(BASE),
    se = sd/sqrt(n())
  ) 

#TREATMENT                    mean    sd    se
#1 Cobimetinib and vemurafenib  68.7  15.7  2.07
#2 Placebo and vemurafenib      74.8  15.9  2.48

t_base <- t.test(BASE ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_base$estimate)
#6.05472

t_base$conf.int
#[1] -12.4786109   0.3691714

t_base$p.value
#[1] 0.06436778

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median = median(BASE),
    q25 = quantile(BASE, 0.25),
    q75 = quantile(BASE, 0.75),
    IQR = q75 - q25
  ) 

#TREATMENT                   median   q25   q75   IQR
#1 Cobimetinib and vemurafenib   69.4  59.4  80.7  21.3
#2 Placebo and vemurafenib       74.1  62.9  88.9  26.1

wilcox_base <- wilcox.test(BASE ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                           conf.level = 0.95, exact = FALSE)

wilcox_base$estimate
# -6.219969

wilcox_base$conf.int
#-12.8799444   0.5999383

wilcox_base$p.value
#[1] 0.07704929

# ----- Cycle 8 day 1 comparison  
cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    mean_qol = mean(QOL),
    sd_qol = sd(QOL),
    se_qol = sd_qol/sqrt(n()),
    
    mean_cfb = median(BASE_CHANGE),
    sd_cfb = sd(BASE_CHANGE),
    se_cfb = sd_cfb/sqrt(n()),
  ) 

#TREATMENT                   mean_qol sd_qol se_qol mean_cfb sd_cfb se_cfb
#1 Cobimetinib and vemurafenib     68.7   17.1   2.27   -0.470   15.6   2.06
#2 Placebo and vemurafenib         70.2   15.0   2.34   -6.24    16.0   2.50

# ----- Observed: QOL
t_qol <- t.test(QOL ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_qol$estimate)
#1.537843

t_qol$conf.int
#[1] -8.001962  4.926275

t_qol$p.value
#[1] 0.637702

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median_qol = median(QOL),
    q25_qol = quantile(QOL, 0.25),
    q75_qol = quantile(QOL, 0.75),
    IQR_qol = q75_qol - q25_qol
  ) 

#TREATMENT                   median_qol q25_qol q75_qol IQR_qol
#1 Cobimetinib and vemurafenib       68.3    58.8    79.9    21.1
#2 Placebo and vemurafenib           70.7    62.4    80.9    18.5

wilcox_qol <- wilcox.test(QOL ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                          conf.level = 0.95, exact = FALSE)

wilcox_qol$estimate
# -1.550039 

wilcox_qol$conf.int
#[1] -7.589974  4.789969

wilcox_qol$p.value
#[1] 0.629433

# ----- Observed: BASE_CHANGE
t_cfb <- t.test(BASE_CHANGE ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_cfb$estimate)
#-4.516876 

t_cfb$conf.int
#[1]-1.930554 10.964306

t_cfb$p.value
#[1] 0.1672709

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median_cfb = median(BASE_CHANGE),
    q25_cfb = quantile(BASE_CHANGE, 0.25),
    q75_cfb = quantile(BASE_CHANGE, 0.75),
    IQR_cfb = q75_cfb - q25_cfb
  ) 

#TREATMENT                   median_cfb q25_cfb q75_cfb IQR_cfb
#1 Cobimetinib and vemurafenib     -0.470   -8.11    6.30    14.4
#2 Placebo and vemurafenib         -6.24   -13.2     1.43    14.7

wilcox_cfb <- wilcox.test(BASE_CHANGE ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                          conf.level = 0.95, exact = FALSE)

wilcox_cfb$estimate
# 5.070022  

wilcox_cfb$conf.int
#[1] -0.3399857 10.6299697

wilcox_cfb$p.value
#[1]0.06628817

# ============================================================
# Part II - Simulation for QoL longitudinal Change (Treatment + BECOG)
# ============================================================
rm(list=ls())
qol_subgroup <- read.csv(here("data","subgroup.csv")) 

# === 1. Baseline probability for each score (ECOG 0 vs ECOG 1)
qol_base <- qol_subgroup %>%
          filter(TIME == 0) %>%
          select(rep, ID, obs, obsid, BECOG) %>%
          mutate(QNAME = dplyr::recode(obsid, "y1" = "Overall Health", "y2" = "Quality of Life"))

qol_base_mean <- qol_base %>%
              group_by(BECOG, QNAME) %>%
              summarise(mean_score = mean(obs), .groups = "drop") %>%
              mutate(ECOG = factor(BECOG, levels = c(0,1), labels = c("ECOG 0", "ECOG 1"))) 

prob_base <- qol_base %>%
              group_by(BECOG, QNAME, obs) %>%
              summarise(n = n(), .groups = "drop") %>%
              group_by(BECOG, QNAME) %>%
              mutate(prob = n/sum(n)) %>%
              ungroup() %>%
              mutate(ECOG = factor(BECOG, levels = c(0,1), labels = c("ECOG 0", "ECOG 1"))) %>%
            # complete the probability for score 1 - 7
              complete(ECOG, QNAME, obs = 1:7, fill = list(prob = 0, n = 0))

bubble1 <- ggplot(prob_base, aes(x = obs, y = QNAME, size = prob, color = ECOG)) +
          geom_point(shape = 21, fill = NA, stroke = 1) +
          geom_point(data = qol_base_mean, aes(x = mean_score, y = QNAME, color = ECOG),
                     shape = 4, size = 2, stroke = 1.3) +
          labs(y = NULL, x = "Score") +
          scale_size_continuous(name = "Probability", 
                                range = c(0,12), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
          scale_color_manual(name = "ECOG Status", values = c("ECOG 0" = "blue", "ECOG 1" = "red"),
                             guide = guide_legend(override.aes = list(size = 3, stroke = 1.3))) +
          scale_x_continuous(breaks = 1:7) +
          theme_bw() +
          theme(
              legend.position = "top",
              legend.box = "vertical",
              legend.spacing.y = unit(0.3, "cm"),
              axis.text = element_text(color = "black", size = 16),
              axis.title = element_text(color = "black", size = 16),
              legend.title = element_text(size = 16),
              legend.text = element_text(size = 16),
          )
          
print(bubble1)

#ggsave(here("plot", "base_score.tiff"), bubble1, width = 12, height = 6, unit = "in", dpi = 300)

# === 2. Steady State probability for each score (Treatment Vs. Placebo)
qol_ss <- qol_subgroup %>%
  filter(TIME == 360 & BECOG == 0) %>%
  select(rep, ID, obs, obsid, ACTARM) %>%
  mutate(QNAME = dplyr::recode(obsid, "y1" = "Overall Health", "y2" = "Quality of Life"))

qol_ss_mean <- qol_ss %>%
  group_by(ACTARM, QNAME) %>%
  summarise(mean_score = mean(obs), .groups = "drop") %>%
  mutate(Treatment = factor(ACTARM, levels = c("VEMURAFENIB + GDC0973","VEMURAFENIB + PLACEBO"), 
                            labels = c("Cobimetinib and vemurafenib", "Placebo and vemurafenib"))) 

prob_ss <- qol_ss %>%
  group_by(ACTARM, QNAME, obs) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ACTARM, QNAME) %>%
  mutate(prob = n/sum(n)) %>%
  ungroup() %>%
  mutate(Treatment = factor(ACTARM, levels = c("VEMURAFENIB + GDC0973","VEMURAFENIB + PLACEBO"), 
                            labels = c("Cobimetinib and vemurafenib", "Placebo and vemurafenib"))) %>%
  # complete the probability for score 1 - 7
  complete(Treatment, QNAME, obs = 1:7, fill = list(prob = 0, n = 0))

bubble2 <- ggplot(prob_ss, aes(x = obs, y = QNAME, size = prob, color = Treatment)) +
  geom_point(shape = 21, fill = NA, stroke = 1) +
  geom_point(data = qol_ss_mean, aes(x = mean_score, y = QNAME, color = Treatment),
             shape = 3, size = 2, stroke = 1.3) +
  labs(y = NULL, x = "Score") +
  scale_size_continuous(name = "Probability", 
                        range = c(0,12), breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1)) +
  scale_color_manual(name = "Treatment", values = c("Cobimetinib and vemurafenib" = "blue", 
                                                    "Placebo and vemurafenib" = "red"),
                     guide = guide_legend(override.aes = list(size = 3, stroke = 1.3))) +
  scale_x_continuous(breaks = 1:7) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.box = "vertical",
    legend.spacing.y = unit(0.3, "cm"),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(color = "black", size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
  )

print(bubble2)

#ggsave(here("plot", "ss_score.tiff"), bubble2, width = 12, height = 6, unit = "in", dpi = 300)

# === 3. Longitudinal Change (treatment + BECOG)
qol_wide_grp <- qol_subgroup %>%
              select(rep, TIME, obs, obsid, group) %>%
              pivot_wider(names_from = obsid, values_from = obs) %>%
              rename(Q29 = y1, Q30 = y2) %>%
              mutate(
                raw_score = (Q29 + Q30)/2,
                QoL = round((((raw_score - 1)/6) * 100),2)
              )

qol_summary <- qol_wide_grp %>%
          group_by(group, TIME) %>%
          summarise(
            mean_QoL = mean(QoL),
            sd_QoL = sd(QoL),
            lower_CI = mean_QoL - 1.96* sd_QoL,
            upper_CI = mean_QoL + 1.96* sd_QoL,
            .groups = "drop"
          )

qol_summary2 <- qol_summary %>%
        mutate(Treatment = dplyr::recode(group,
                                   "TRTECOG0" = "Cobimetinib and vemurafenib",
                                   "TRTECOG1" = "Cobimetinib and vemurafenib",
                                   "PlaceboECOG0" = "Placebo and vemurafenib",
                                   "PlaceboECOG1" = "Placebo and vemurafenib"),
               ECOG = dplyr::recode(group,
                              "TRTECOG0" = "ECOG0",
                              "TRTECOG1" = "ECOG1",
                              "PlaceboECOG0" = "ECOG0",
                              "PlaceboECOG1" = "ECOG1"))

p2 <- ggplot(qol_summary2, aes(x = TIME, y = mean_QoL, color = Treatment, fill = Treatment)) +
      geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.35, color = NA) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~ECOG) +
      labs(x = "Time (Days)", y = "QoL Score, mean (95% CI)") +
      scale_y_continuous(breaks = seq(20,100,20),limits = c(20,100)) +
      scale_x_continuous(breaks = seq(0,360,60),limits = c(0,360)) +
      scale_color_manual(values = c("Cobimetinib and vemurafenib" = "#2166AC",
                                    "Placebo and vemurafenib" = "#B2182B")) +
      scale_fill_manual(values = c("Cobimetinib and vemurafenib" = "#2166AC",
                                   "Placebo and vemurafenib" = "#B2182B")) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 16, face = "bold"),
        legend.position = "right",
        axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
      )

print(p2)

qol_210 <- qol_summary %>%
          filter(TIME == 210)
print(qol_210)
#group         TIME mean_QoL sd_QoL lower_CI upper_CI
#1 PlaceboECOG0   210     64.9   8.16     48.9     80.9
#2 PlaceboECOG1   210     45.7   9.80     26.4     64.9
#3 TRTECOG0       210     77.3   6.40     64.8     89.9
#4 TRTECOG1       210     60.5   8.41     44.0     77.0

#ggsave(here("plot", "Mean QoL_grp.tiff"), p2, width = 13, height = 5, unit = "in", dpi = 300)


# compare the QoL difference between arms (subgroup by ECOG)
qol_arm_wide <- qol_wide_grp %>%
            select(rep, TIME, group, QoL) %>%
            pivot_wider(names_from = group, values_from = QoL)

diff_arm <- qol_arm_wide %>%
        mutate(diff_ECOG0 = TRTECOG0 - PlaceboECOG0,
               diff_ECOG1 = TRTECOG1 - PlaceboECOG1) %>%
        select(rep, TIME, diff_ECOG0, diff_ECOG1)

summary_ECOG0 <- diff_arm %>%
          group_by(TIME) %>%
          summarise(
            mean_diff = mean(diff_ECOG0),
            sd_diff = sd(diff_ECOG0),
            lower_CI = mean_diff - 1.96* sd_diff,
            upper_CI = mean_diff + 1.96* sd_diff,
            .groups = "drop"
          )

summary_ECOG1 <- diff_arm %>%
          group_by(TIME) %>%
          summarise(
            mean_diff = mean(diff_ECOG1),
            sd_diff = sd(diff_ECOG1),
            lower_CI = mean_diff - 1.96* sd_diff,
            upper_CI = mean_diff + 1.96* sd_diff,
            .groups = "drop"
          )

# compare treatment difference (no ECOG subgroup)
qol_arm2 <- qol_wide_grp %>%
          mutate(arm = case_when(
            grepl("^TRT", group) ~ "Treatment",
            grepl("^Placebo", group) ~ "Placebo",
            TRUE ~ NA_character_
          ))

arm_diff <- qol_arm2 %>%
          group_by(rep, TIME, arm) %>%
          summarise(QoL_mean = mean(QoL), .groups = "drop") %>%
          pivot_wider(names_from = arm, values_from = QoL_mean) %>%
          mutate(diff = Treatment - Placebo) %>%
          select(rep, TIME, diff)

summary_arm <- arm_diff %>%
          group_by(TIME) %>%
          summarise(
            mean_diff = mean(diff),
            sd_diff = sd(diff),
            lower_CI = mean_diff - 1.96* sd_diff,
            upper_CI = mean_diff + 1.96* sd_diff,
            .groups = "drop"
          )






