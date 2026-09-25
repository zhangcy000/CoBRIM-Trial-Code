# --------------------------
# Part : ICC Curve for Question 29 and Question 30
# Author: Congyu Zhang
# --------------------------
#1 . Load the Environment and import the raw data
getRversion()

rm(list=ls())
packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "here", "ggh4x")
lapply(packages,library,character.only = TRUE)

#2. Create Graded Response Equation
item_names <- c("Overall Health", "Quality of Life")

gr_model <- function(item_list, item_names){
  
  all_data <- data.frame()
  theta <- seq(-4, 3, length.out = 100)
  
  for(item_idx in 1:length(item_list)){
    item <- item_list[[item_idx]]
    a <- item$a
    b <- item$b
    n_cat <- length(b) + 1
    
    name <- item_names[item_idx]
    
    for (k in 0:(n_cat - 1)){
      prob <- numeric(length(theta))
      
      for (i in 1:length(theta)){
        if(k == 0){
          prob[i] <- 1 - 1/(1 + exp(-a * (theta[i] - b[1])))
        } else if(k == n_cat-1){
          prob[i] <- 1/(1+exp(-a * (theta[i] - b[k])))
        } else {
          prob[i] <- 1/(1+exp(-a * (theta[i] - b[k]))) - 1/(1+exp(-a * (theta[i] - b[k+1])))
        }
      }
      
      item_data <- data.frame(
        theta = theta,
        probability = prob,
        score = factor(k+1, levels = 1:n_cat),
        item = name
      )
      
      all_data <- rbind(all_data, item_data)
    }
  }
  
  return(all_data)
}


#3. Plot function
gr_plot <- function(item_list, item_names){
  
  # Generate the data
  plot_data <- gr_model(item_list, item_names)
  plot_data$item <- factor(plot_data$item, levels = item_names)
  
  levels_color <- c("red3", "darkorange2", "orange", "gold", "yellowgreen", "green3", "darkgreen")
  
  # Generate the plot
  prob <- ggplot(plot_data, aes(x = theta, y = probability, color = score)) +
          geom_line(size = 1) +
          facet_wrap(~item, ncol = 2) +
          labs(x = "Latent Variable", y = "Probability") +
          theme_bw() +
          theme(
            strip.text = element_text(size = 15, face = "bold"),
            legend.position = "top",
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16),
            axis.title = element_text(size = 15, face = "bold"),
            axis.text = element_text(size = 15, color = "black")
          ) +
          scale_color_manual(values = levels_color, breaks = 1:7, labels = as.character(1:7)) +
          guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
          scale_y_continuous(breaks = seq(-4,3,1)) +
          scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1))
  
  return(prob)
}

# 4. Import the ICC parameters
ICC <- read.csv(here("data","poppk","irt","ICC.txt")) %>%
          select(parameter, value)

shared_a <- ICC$value[ICC$parameter == "DIS1_pop"]

items <- lapply(1:2, function(n){
  b <- ICC$value[ICC$parameter %in% paste0("Q", n, "th", 1:6, "_pop")]
  list (a = shared_a, b = cumsum(b))
})

names(items) <- paste0("item", 1:length(items))

# 5. Generate the plot
p1 <- gr_plot(items, item_names)
print(p1)

#ggsave("plot/ICC.tiff", plot = p1, width = 8, height = 4, unit ="in", dpi = 300)





