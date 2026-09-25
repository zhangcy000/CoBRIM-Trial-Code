# --------------------------
# Part : VPC Plot for IRT - QoL Model
# Author: Congyu Zhang
# --------------------------
#1 . Load the Environment and import the raw data
getRversion()

rm(list=ls())
packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "here", "ggh4x")
lapply(packages,library,character.only = TRUE)

#2. VPC generation
#=========================================================
# Function to generate the data for VPC
#=========================================================
draw_vpc_categorical <- function(project.name, item.number, simdata.list){
    # Read empirical distribution data
    empirical <- read.csv(here("data", "poppk", "irt", 
                    paste0(project.name, "/ChartsData/VisualPredictiveCheck/y", item.number, "_distribution.txt")),
                    header = TRUE)
    
    # Get simulation data for this item
    simdata <- simdata.list[[item.number]]
    
    # Extract bin edges from empirical data
    bins <- unique(c(empirical$binsTimeBefore, empirical$binsTimeAfter))
    bins <- sort(bins)
    binsmiddle <- (bins[2:length(bins)] + bins[1:(length(bins)-1)])/2
    
    # Assign simulated data to bins
    simdata$binindex <- findInterval(simdata$time, bins, all.inside = TRUE)
    
    # Get the column name for this item's simulated values
    y_col <- paste0("sim_y", item.number)
    
    # Calculate proportions of each level in each bin for each replicate
    prop_list <- list()

    for (cat in 1:7) {
      prop_cat <- aggregate(
        simdata[[y_col]] == cat,
        by = list (bin = simdata$binindex, rep = simdata$rep),
        FUN = mean
      )
      names(prop_cat)[3] <- "proportion"
      prop_cat$category <- cat
      prop_list[[cat]] <- prop_cat
    }
    
    # Combine all categories
    all_props <- do.call(rbind, prop_list)
    
    # Calculate median, 2.5% and 97.%% percentiles across replicates
    total_score_rep <- all_props %>%
      group_by(bin,rep) %>%
      summarise(total_score = sum(proportion * category), .groups = "drop")
    
    total_score_intervals <- total_score_rep %>%
      group_by(bin) %>%
      summarise(
        lower = quantile(total_score, 0.025),
        median = quantile(total_score, 0.50),
        upper = quantile(total_score, 0.975)
      ) %>%
      mutate(bin_middle = binsmiddle[bin])
    
    # Calculate total score for empirical data
    empirical_total <- empirical %>%
         mutate(
           category_value = case_when(
             category == "[1-1]" ~ 1,
             category == "[2-2]" ~ 2,
             category == "[3-3]" ~ 3,
             category == "[4-4]" ~ 4,
             category == "[5-5]" ~ 5,
             category == "[6-6]" ~ 6,
             category == "[7-7]" ~ 7,
           ),
           bin_middle = (binsTimeBefore + binsTimeAfter)/2
         ) %>%
      group_by(bin_middle) %>%
      summarise(
        empirical_score = sum(propCategory_empirical * category_value)
      )
    
    # Return data for plotting
    vpc_data <- data.frame(
      TIME = total_score_intervals$bin_middle,
      exp_l = total_score_intervals$lower,
      exp_h = total_score_intervals$upper,
      exp_obs = empirical_total$empirical_score[match(total_score_intervals$bin_middle,empirical_total$bin_middle)],
      ITEM = paste0("Item_", item.number)
    )
    
    return(vpc_data)
}

# Function to create combined data for QoL items
vpc_data_all <- function(project.name, simdata.list) {
  
  # Create list to store data
  vpc_data_list <- list()
  
  # Generate VPC data for both QoL items
  for (i in 1:2){
    vpc_data_list[[i]] <- draw_vpc_categorical(project.name, i, simdata.list)
  }
  
  # Combine all data
  vpc_exp <- do.call(rbind, vpc_data_list)
 
  return(vpc_exp) 
}


project.name <- "test4"

# Import all simulation data for two QoL items
simdata.list <- list()

for (i in 1:2){
  sim_yi <- read.csv(here("data", "poppk", "irt", 
                    paste0(project.name, "/ChartsData/VisualPredictiveCheck/y", i, "_simulations.txt")),
                     header = TRUE) %>%
            select(-split, -color, -filter)
  
  simdata.list[[i]] <- sim_yi
}

vpc_basic <- vpc_data_all(project.name, simdata.list)

vpc_plot <- vpc_basic %>%
          mutate(Question = case_when(
            ITEM == "Item_1" ~ "Q29: Overall Health",
            ITEM == "Item_2" ~ "Q30: Overll QoL",
          ))

p1 <- ggplot(vpc_plot, aes(x = TIME)) +
    geom_ribbon(aes(ymin = exp_l, ymax = exp_h), fill = "grey60", alpha = 0.5) +
    geom_line(aes(y = exp_obs), color = "royalblue4", linewidth = 1.2) +
    labs(x = "Time (Days)", y = "Item Score") +
    scale_x_continuous(limits = c(0,650), breaks = seq(0, 650, 150)) +
    scale_y_continuous(limits = c(1,7), breaks = seq(1,7,1)) +
    facet_wrap(~Question) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 15, color = "black"),
      strip.text = element_text(size = 15, face = "bold", color = "black")
    )

print(p1)

#ggsave(here("plot", "vpc.tiff"), p1, width = 8, height = 4, unit = "in", dpi = 300)

#3. Distribution of parameters and random effects
rm(list=ls())

# ==== 3.1 random effects
re_dis <- read.csv(here("data", "poppk", "irt", "test4", "ChartsData",
                     "DistributionOfTheStandardizedRandomEffects", "pdf.txt")) %>%
          select(-split)

re_long <- re_dis %>%
            pivot_longer(
              cols = everything(),
              names_to = c("parameter", ".value"),
              names_pattern = "standEta_(.*)_(abscissa|pdf)"
            ) %>%
            group_by(parameter) %>%
            arrange(abscissa, .by_group = TRUE) %>%
            mutate(bin_width = c(diff(abscissa), diff(abscissa)[1])) %>%
            ungroup()
            
re_hist_plot <- ggplot(re_long, aes(x = abscissa, y = pdf, color = parameter)) + 
                geom_col(aes(width = bin_width), fill = "steelblue", color = "white", alpha =0.9) +
                stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "black", linewidth = 0.8,
                              inherit.aes = FALSE, data = data.frame(abscissa = seq(-4,4,0.01)), 
                              aes(x = abscissa)) +
                labs(y = "Probability Density Function", x = "Random Effect \u03b7") +
                facet_wrap(~parameter, scales = "free_x") +
                theme_bw() +
                theme(
                    strip.text = element_text(size = 15, face = "bold"),
                    axis.title = element_text(size = 15, face = "bold"),
                    axis.text = element_text(size = 15, color = "black")
              )

print(re_hist_plot)

#ggsave(here("plot", "RE.tiff"), re_hist_plot, width = 10, height = 5, unit = "in", dpi = 300)

# ==== 3.2 NPDE
npde_q29 <- read.csv(here("data", "poppk", "irt", "test4", "ChartsData",
                        "DistributionOfTheResiduals", "y1_pdf.txt")) %>%
            select(-split) %>%
            mutate(QNAME = "Q29: Overall Health")

npde_q30 <- read.csv(here("data", "poppk", "irt", "test4", "ChartsData",
                          "DistributionOfTheResiduals", "y2_pdf.txt")) %>%
            select(-split) %>%
            mutate(QNAME = "Q30: Overall QoL")

npde <- rbind(npde_q29, npde_q30)

npde_plot <- npde %>%
  group_by(QNAME) %>%
  arrange(npde_abscissa, .by_group = TRUE) %>%
  mutate(bin_width = c(diff(npde_abscissa), diff(npde_abscissa)[1])) %>%
  ungroup()

npde_hist_plot <- ggplot(npde_plot, aes(x = npde_abscissa, y = npde_pdf)) + 
  geom_col(aes(width = bin_width), fill = "steelblue", color = "white", alpha =0.9) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "black", linewidth = 0.8,
                inherit.aes = FALSE, data = data.frame(npde_abscissa = seq(-2,2,0.01)), 
                aes(x = npde_abscissa)) +
  labs(y = "Probability Density Function", x = "NPDE") +
  facet_wrap(~QNAME) +
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1)) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 15, color = "black")
  )

print(npde_hist_plot)

#ggsave(here("plot", "NPDE.tiff"), npde_hist_plot, width = 10, height = 5, unit = "in", dpi = 300)







