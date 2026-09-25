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
# Part I - Simulation from Cycle 1 day 1 to Cycle 8 day 1 (default simulation)
# ============================================================
#1. Import the original data and simulated data
qol <- read.csv(here("data","QoL_full.csv"))

qol_sim <- read.csv(here("data","cycle8_full.csv")) %>%
            select(-ID, -UID)

#3. Data Interpretation 
# First step is to filter out patient ID at each time point
visit_time_map <- tibble(
  VISIT = c("CYCLE 1 DAY 1", "CYCLE 1 DAY 15", "CYCLE 2 DAY 1", "CYCLE 2 DAY 15",
            "CYCLE 4 DAY 1", "CYCLE 6 DAY 1", "CYCLE 8 DAY 1"),
  TIME = c(0, 15, 30, 45, 90, 150, 210)
)

id_by_visit <- qol %>%
        filter(VISIT %in% visit_time_map$VISIT) %>%
        distinct(UNI_TRNC, VISIT) %>%
        left_join(visit_time_map, by = "VISIT")

qol_filter <- qol_sim %>%
            inner_join(id_by_visit, by = c("original_id" = "UNI_TRNC", "TIME" = "TIME")) 

# calculate the average qol score of each patient
qol_avg <- qol_filter %>%
            group_by(original_id, TIME, obsid, ACTARM, BECOG) %>%
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
            filter(TIME %in% c(0,15,30,45,90,150,210)) %>%
            mutate(VISIT = case_when(
              TIME == 0 ~ "Cycle 1 Day 1",
              TIME == 15 ~ "Cycle 1 Day 15",
              TIME == 30 ~ "Cycle 2 Day 1",
              TIME == 45 ~ "Cycle 2 Day 15",
              TIME == 90 ~ "Cycle 4 Day 1",
              TIME == 150 ~ "Cycle 6 Day 1",
              TIME == 210 ~ "Cycle 8 Day 1",
              TRUE ~ NA_character_
            ))

# --- statistical analysis comparing two arms (CFB only)
qol_plot_mean <- qol_plot %>%
  group_by(TIME, TREATMENT, VISIT) %>%
  summarise(
    mean_change = mean(BASE_CHANGE, na.rm = TRUE),
    se = sd(BASE_CHANGE, na.rm = TRUE)/sqrt(n()),
    sd = sd(BASE_CHANGE, na.rm = TRUE),
    pt_num = n(),
    .groups = "drop"
  )

qol_pvalue <- qol_plot %>%
  filter(TIME != "0") %>%
  group_by(VISIT) %>%
  summarise(
    p_value_t = t.test(BASE_CHANGE ~ TREATMENT)$p.value,
    p_value_wil = wilcox.test(BASE_CHANGE ~ TREATMENT)$p.value,
    .groups = "drop"
  )

stat_qol <- qol_plot_mean %>%
  filter(TIME != "0") %>%
  select(TREATMENT, VISIT, TIME, mean_change, sd) %>%
  pivot_wider(names_from = TREATMENT, values_from = c(mean_change, sd)) %>%
  left_join(qol_pvalue, by = c("VISIT")) %>%
  mutate(p.adj_t = p.adjust(p_value_t, method = "BH"),
         p.adj_wil = p.adjust(p_value_wil, method = "BH"))

#write.csv(stat_qol, here("data","visit_comparison_sim.csv"), row.names = FALSE)

p1 <- ggplot(qol_plot_mean, aes(x = TIME, y = mean_change, color = TREATMENT, 
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

#ggsave(here("plot", "QoL by Visit_se_sim_2.tiff"), p1, width = 12, height = 6, unit = "in", dpi = 300)


## --- statistical analysis comparing two arms (raw data and CFB)
cycle8_qol <- qol_plot %>%
  filter(VISIT == "Cycle 8 Day 1")

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
#1 Cobimetinib and vemurafenib  69.0  16.0  2.12
#2 Placebo and vemurafenib      75.3  16.6  2.59

t_base <- t.test(BASE ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_base$estimate)
# 6.311955

t_base$conf.int
#[1] -12.972425   0.348514

t_base$p.value
#[1] 0.06295053

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median = median(BASE),
    q25 = quantile(BASE, 0.25),
    q75 = quantile(BASE, 0.75),
    IQR = q75 - q25
  ) 

#TREATMENT                   median   q25   q75   IQR
#1 Cobimetinib and vemurafenib   69.2  59.4  81.4  21.9
#2 Placebo and vemurafenib       74.9  62.6  89.9  27.3

wilcox_base <- wilcox.test(BASE ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                           conf.level = 0.95, exact = FALSE)

wilcox_base$estimate
# -6.520033 
wilcox_base$conf.int
#-13.3299876   0.6199786
wilcox_base$p.value
#[1] 0.08073234

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
#1 Cobimetinib and vemurafenib     68.7   17.1   2.27    0.320   17.0   2.25
#2 Placebo and vemurafenib         70.1   15.2   2.38   -6.42    17.2   2.68

# ----- Observed: QOL
t_qol <- t.test(QOL ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_qol$estimate)
#1.398703 

t_qol$conf.int
#[1] -7.926594  5.129187

t_qol$p.value
#[1] 0.671421

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median_qol = median(QOL),
    q25_qol = quantile(QOL, 0.25),
    q75_qol = quantile(QOL, 0.75),
    IQR_qol = q75_qol - q25_qol
  ) 

#TREATMENT                   median_qol q25_qol q75_qol IQR_qol
#1 Cobimetinib and vemurafenib       68.4    58.8    79.6    20.8
#2 Placebo and vemurafenib           71.8    61.9    81.3    19.4

wilcox_qol <- wilcox.test(QOL ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                          conf.level = 0.95, exact = FALSE)

wilcox_qol$estimate
# -1.510025

wilcox_qol$conf.int
#[1] -7.489982  4.929973

wilcox_qol$p.value
#[1] 0.6474421

# ----- Observed: BASE_CHANGE
t_cfb <- t.test(BASE_CHANGE ~ TREATMENT, data = cycle8_qol, conf.level = 0.95)

diff(t_cfb$estimate)
#-4.913252

t_cfb$conf.int
#[1] -2.047249 11.873753

t_cfb$p.value
#[1] 0.1641396

cycle8_qol %>%
  group_by(TREATMENT) %>%
  summarise(
    median_cfb = median(BASE_CHANGE),
    q25_cfb = quantile(BASE_CHANGE, 0.25),
    q75_cfb = quantile(BASE_CHANGE, 0.75),
    IQR_cfb = q75_cfb - q25_cfb
  ) 

#TREATMENT                   median_cfb q25_cfb q75_cfb IQR_cfb
#1 Cobimetinib and vemurafenib      0.320   -8.48    7.03    15.5
#2 Placebo and vemurafenib         -6.42   -16.1     1.38    17.4

wilcox_cfb <- wilcox.test(BASE_CHANGE ~ TREATMENT, data = cycle8_qol, conf.int = TRUE,
                          conf.level = 0.95, exact = FALSE)

wilcox_cfb$estimate
# 6.010032  

wilcox_cfb$conf.int
#[1] 0.03009355 11.66000084

wilcox_cfb$p.value
#[1] 0.04928596


# ============================================================
# Part III - NLME model demonstration
# ============================================================
rm(list=ls())
#-----------------------------------------
#-- Demonstration of IIV
#-----------------------------------------
#1. Import the population prm and individual prm
pop_prm <- read.csv(here("data","poppk","irt","pop_prm.txt")) %>%
            filter(parameter %in% c("BASE_pop", "beta_EMAX_ACTARM_VEMURAFENIB___PLACEBO", "KD_pop")) %>%
            select(parameter, value) %>%
            pivot_wider(names_from = parameter, values_from = value) %>%
            rename(BASE = BASE_pop, EMAX = beta_EMAX_ACTARM_VEMURAFENIB___PLACEBO, KD = KD_pop)

idv_prm <- read.csv(here("data","poppk","irt","prm_test4.txt")) %>%
              filter(ACTARM == "VEMURAFENIB + PLACEBO") %>%
              select(id, BASE_mode, EMAX_mode, KD_mode, tBECOG, ACTARM)

# randomly pick three patients
set.seed(1)

idv_sample <- idv_prm %>%
          distinct(id, .keep_all = TRUE) %>%
          slice_sample(n = 3) %>%
          select(id, BASE_mode, EMAX_mode, KD_mode) %>%
          mutate(plot_id = row_number())

# generate the plot
time_seq <- seq(0, 240, by = 1)

pop_curve <- tibble(
  time = time_seq,
  value = pop_prm$BASE + pop_prm$EMAX * (1 - exp(-pop_prm$KD * time_seq)),
  type = "Population"
)

idv_curve <- idv_sample %>%
          rowwise()%>%
          reframe(
            time = time_seq,
            value = BASE_mode + EMAX_mode * (1 - exp(-KD_mode * time_seq)),
            type = paste0("Individual ", plot_id)
          )
            
plot_data <- bind_rows(pop_curve, idv_curve) %>%
              mutate(type = factor(type, levels = c("Population", "Individual 1",
                                                    "Individual 2", "Individual 3")))

IIV_plot <- ggplot(plot_data, aes(x = time, y = value, group = type, color = type)) +
          geom_line(data = filter(plot_data, type == "Population"), linewidth = 1.2) +
          geom_line(data = filter(plot_data, type != "Population"), 
                    linewidth = 0.8) +
          scale_color_manual(
                name ="",
                values = c(
                    "Population" = "black",
                    "Individual 1" = "#7CAE00",
                    "Individual 2" = "#00BFC4",
                    "Individual 3" = "#C77CFF"
                    )
          ) +
          scale_x_continuous(breaks = seq(0,240,60),limits = c(0,240)) +
          labs(x = "Time (day)", y = "QOL Latent") +
          theme_bw() +
          theme(
            legend.position = "right",
            axis.text = element_text(color = "black", size = 15),
            axis.title = element_text(color = "black", size = 15),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            panel.grid = element_blank()
          ) 

print(IIV_plot)

#ggsave(here("plot", "NLME1.tiff"), IIV_plot, width = 8, height = 4, unit = "in", dpi = 300)

#-----------------------------------------
#-- Demonstration of epison
#-----------------------------------------
rm(list=ls())

#1. Import the original data and simulated data
qol <- read.csv(here("data","QoL_full.csv"))

qol_sim <- read.csv(here("data","cycle8_full.csv")) %>%
            select(-ID, -UID)

pop_prm <- read.csv(here("data","poppk","irt","pop_prm.txt")) 
idv_prm <- read.csv(here("data","poppk","irt","prm_test4.txt"))


q29_thershold <- pop_prm %>%
          filter(parameter %in% c("Q1th1_pop", "Q1th2_pop", "Q1th3_pop", 
                                  "Q1th4_pop", "Q1th5_pop", "Q1th6_pop")) %>%
          arrange(parameter) %>%
          pull(value)


one_id <- sample(unique(qol_sim$original_id),1)

obs_one <- qol_sim %>%
          filter(original_id == one_id, obsid == "y1", rep == 1) %>%
          select(TIME, obs) %>%
          rename(time = TIME, obs_score= obs) %>%
          filter(time %in% c(0,15,30,45,60,90,120,150,180,210,240))

one_prm <- idv_prm %>% filter(id == one_id)

time_interval <- seq(0, 240, by = 1)

theta_curve <- tibble(
  time = time_interval,
  theta = one_prm$BASE_mode + one_prm$EMAX_mode * (1 - exp(-one_prm$KD_mode * time_interval))
)

expected_score <- function(theta, thresholds){
  cum_probs <- sapply(thresholds, function(b) 1/(1+exp(-(theta-b))))
  cum_probs <-c(1, cum_probs, 0)
  cat_probs <- -diff(cum_probs)
  sum(cat_probs * (1:7))
}

theta_curve <- theta_curve %>%
          mutate(expected_score = sapply(theta, expected_score, thresholds = q29_thershold))


match_data <- obs_one %>%
          left_join(theta_curve, by = "time")


theta_plot <- ggplot() +
          geom_line(data = theta_curve, aes(x = time, y = expected_score), 
                    color = "steelblue", linewidth = 1) +
          geom_point(data = obs_one, aes(x = time, y = obs_score),
                     color = "black", size = 2) +
          scale_x_continuous(breaks = seq(0,240,60),limits = c(0,240),name = "Time (day)") +
          scale_y_continuous(breaks = seq(1,7,1),limits = c(1,7),name = "Q29:Overall Health Score") +
          theme_bw() +
          theme(
            axis.text = element_text(color = "black", size = 15),
            axis.title = element_text(color = "black", size = 15),
            panel.grid = element_blank()
          )

print(theta_plot)

#ggsave(here("plot", "NLME2.tiff"), theta_plot, width = 4, height = 4, unit = "in", dpi = 300)







