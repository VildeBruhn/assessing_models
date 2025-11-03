###########################################
##                                       ##
##  Plots and statistics on the results  ##
##                                       ##
###########################################

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

rm(list = ls())

library(lme4)
library(lmerTest)
library(readr)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(ggplot2)



#--------------------------------------------------------------------
# Mean support for stasis depending on type of fossil and lifestyle 
#--------------------------------------------------------------------

#Import metadata and results of the test (only models with single processes, no adequacy)
load("./model_test_uni.Rdata")
load("./ln_data_meta_uni.Rdata")

#Extract support for stasis for each time series
support_Stasis <- lapply(model_test, function(df) {
  return(df["Stasis", "Akaike.wt"])
})

tsIDs <- unique(names(ln_data_meta))
for (tsID in tsIDs) {
  ln_data_meta[[tsID]]$support_stasis <- support_Stasis[[tsID]]
}

# Separate the time series in micro and macrofossils
ln_data_meta_macro <- Filter(function(x) x$microfossil[1] == "no", ln_data_meta)
ln_data_meta_micro <- Filter(function(x) x$microfossil[1] == "yes", ln_data_meta)

# Separate the time series in benthic and planktic microfossils
ln_data_meta_micro_benthic <- Filter(function(x) x$environment[1] == "marine benthic", ln_data_meta_micro)
ln_data_meta_micro_planktic <- Filter(function(x) x$environment[1] == "marine planktic", ln_data_meta_micro)

# Separate each dataset into two types of traita
ln_data_meta_macro_size <- Filter(function(x) x$trait_type[1] %in% c("linear", "area"), ln_data_meta_macro)
ln_data_meta_macro_shape <- Filter(function(x) x$trait_type[1] %in% c("ratio", "angle"), ln_data_meta_macro)
ln_data_meta_micro_benthic_size <- Filter(function(x) x$trait_type[1] %in% c("linear", "area"), ln_data_meta_micro_benthic)
ln_data_meta_micro_benthic_shape <- Filter(function(x) x$trait_type[1] %in% c("ratio", "angle"), ln_data_meta_micro_benthic)
ln_data_meta_micro_planktic_size <- Filter(function(x) x$trait_type[1] %in% c("linear", "area"), ln_data_meta_micro_planktic)
ln_data_meta_micro_planktic_shape <- Filter(function(x) x$trait_type[1] %in% c("ratio", "angle"), ln_data_meta_micro_planktic)

# Get the sample size for each dataset
ln_data_meta_micro_benthic_size_stasis <- unlist(lapply(ln_data_meta_micro_benthic_size, function(x) x$support_stasis[1]))
ln_data_meta_micro_benthic_shape_stasis <- unlist(lapply(ln_data_meta_micro_benthic_shape, function(x) x$support_stasis[1]))
ln_data_meta_micro_planktic_size_stasis <- unlist(lapply(ln_data_meta_micro_planktic_size, function(x) x$support_stasis[1]))
ln_data_meta_micro_planktic_shape_stasis <- unlist(lapply(ln_data_meta_micro_planktic_shape, function(x) x$support_stasis[1]))
ln_data_meta_macro_size_stasis <- unlist(lapply(ln_data_meta_macro_size, function(x) x$support_stasis[1]))
ln_data_meta_macro_shape_stasis <- unlist(lapply(ln_data_meta_macro_shape, function(x) x$support_stasis[1]))

display_order <- c(1:6)

fossiltype <- c("micro benthic", "micro benthic", "micro planktic", "micro planktic", "macro", "macro")

trait <- c("size", "shape", "size", "shape", "size", "shape")

sample_size <- c(length(ln_data_meta_micro_benthic_size_stasis), length(ln_data_meta_micro_benthic_shape_stasis), length(ln_data_meta_micro_planktic_size_stasis),
                length(ln_data_meta_micro_planktic_shape_stasis), length(ln_data_meta_macro_size_stasis), length(ln_data_meta_macro_shape_stasis))
                
# Get the mean for each dataset
mean <- c(mean(ln_data_meta_micro_benthic_size_stasis), mean(ln_data_meta_micro_benthic_shape_stasis), mean(ln_data_meta_micro_planktic_size_stasis),
                mean(ln_data_meta_micro_planktic_shape_stasis), mean(ln_data_meta_macro_size_stasis), mean(ln_data_meta_macro_shape_stasis))

# Get the confidence interval for each dataset
sd <- c(sd(ln_data_meta_micro_benthic_size_stasis), sd(ln_data_meta_micro_benthic_shape_stasis), sd(ln_data_meta_micro_planktic_size_stasis),
       sd(ln_data_meta_micro_planktic_shape_stasis), sd(ln_data_meta_macro_size_stasis), sd(ln_data_meta_macro_shape_stasis))

se <- sd / sqrt(sample_size)
t <- qt(0.975, df = sample_size - 1)

# Make dataframe
ln_data_support_stasis <- data.frame(display_order, fossiltype, trait, sample_size, mean, sd, se, t)
ln_data_support_stasis$display_order <- factor(ln_data_support_stasis$display_order, levels = 1:6, ordered = TRUE)

# Create plot
boxplot_fossiltype_trait = ggplot(ln_data_support_stasis, aes(x = interaction(fossiltype, trait), y = mean, shape = trait, fill = trait)) +
  geom_point(size = 4, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  scale_shape_manual(values = c(21, 22)) +  # circle = 21, square = 22
  scale_fill_manual(values = c("black", "grey70")) +
  labs(y = "Support for Stasis Model", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "lightgray", linetype = "dashed", size = 0.5),
    plot.title = element_text(size = 38),
    axis.title = element_text(size = 38),
    axis.text = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

# save the graphs
png("./results_paleoTS_v0.6.1/plot/plot_fossiltype_trait.png", width = 2000, height = 1500)
grid.arrange(boxplot_fossiltype_trait, nrow = 1)
dev.off()
