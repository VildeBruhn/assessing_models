##############################################################
##   Assessing models of evolution (mode shift included)    ##
##                       STATISTICS                         ##
##############################################################

# install packages
#install.packages("evoTS") #version 1.0.3
#install.packages("devtools")
#devtools::install_github("klvoje/adePEM")
#install.packages("tidyverse")
#install.packages("paleoTS") #version 0.6.1
#install.packages("wesanderson")

# load packages
library(tidyverse)
library(devtools)
library(wesanderson)
library(gridExtra)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
library(scales)

rm(list = ls())

# set working directory
setwd("C:/Users/marionth/OneDrive - Universitetet i Oslo/Skrivebord/PhD/Github/assessing_models_evolution/assessing_models")
source("./assessing_models_uni_functions.R")


#--------------
# IMPORT FILES
#--------------
ln_data_meta_shift <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# remove time series with less than 14 steps
ln_data_meta_shift <- subset(ln_data_meta_shift, steps >= 14)

# remove modern time series
ln_data_meta_shift <- subset(ln_data_meta_shift, period_start != "Present")

# remove Syverson
ln_data_meta_shift <- subset(ln_data_meta_shift, URL != "https://doi.org/10.1017/pab.2024.37")

# make list based on ID
ln_data_meta_shift <- lapply(split(ln_data_meta_shift,ln_data_meta_shift$tsID), function(x) as.list(x))


# load data from analyses
load("./ln_data_shift.Rdata")
load("./aicc_shift_passed.Rdata")
load("./adeq_shift_passed.Rdata")


#--------------
# PREPARE DATA
#--------------


# get aicc model info into metadata
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW, "GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW, "URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis, "stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Strict_stasis, "strict stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Decel, "decel")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Accel, "accel")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU, "OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_mov_opt_anc, "OU mov opt anc")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_mov_opt, "OU mov opt")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_Stasis, "Stasis-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_URW, "Stasis-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_GRW, "Stasis-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, Stasis_OU, "Stasis-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_URW, "URW-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_GRW, "URW-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_OU, "URW-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_GRW, "GRW-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_OU, "GRW-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_OU, "OU-OU")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_GRW, "OU-GRW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_URW, "OU-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, OU_Stasis, "OU-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_URW, "GRW-URW")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, GRW_Stasis, "GRW-Stasis")
ln_data_meta_shift <- model_aicc(ln_data_meta_shift, URW_Stasis, "URW-Stasis")

# get the type of model (without or with shift) info into metadata
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Strict_stasis, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Decel, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Accel, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt_anc, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_mov_opt, "no shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, Stasis_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_OU, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_GRW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, OU_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_URW, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, GRW_Stasis, "mode shift")
ln_data_meta_shift <- model_type(ln_data_meta_shift, URW_Stasis, "mode shift")

# get adequate model info into metadata
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_adeq_passed, "GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_adeq_passed, "URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_adeq_passed, "stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Strict_stasis_adeq_passed, "strict stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Decel_adeq_passed, "decel")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Accel_adeq_passed, "accel")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_adeq_passed, "OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_mov_opt_adeq_passed, "OU mov opt")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_Stasis_adeq_passed, "Stasis-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_URW_adeq_passed, "Stasis-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_GRW_adeq_passed, "Stasis-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, Stasis_OU_adeq_passed, "Stasis-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_URW_adeq_passed, "URW-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_GRW_adeq_passed, "URW-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_OU_adeq_passed, "URW-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_GRW_adeq_passed, "GRW-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_OU_adeq_passed, "GRW-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_OU_adeq_passed, "OU-OU")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_GRW_adeq_passed, "OU-GRW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_URW_adeq_passed, "OU-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, OU_Stasis_adeq_passed, "OU-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_URW_adeq_passed, "GRW-URW")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, GRW_Stasis_adeq_passed, "GRW-Stasis")
ln_data_meta_shift <- model_adeq(ln_data_meta_shift, URW_Stasis_adeq_passed, "URW-Stasis")

# bind data to dataframe
unit_list <- c("popID", "total_N", "steps", "interval_MY", "trait_type",
               "lat", "lon", "model_aicc", "model_type", "model_adequate")
plot_data <- bind(ln_data_meta_shift, unit_list)

# remove time series with NA
plot_data <- plot_data %>% drop_na(model_aicc)

# ordering the model according to the number of parameters
level_order <- c("no shift", "mode shift")
plot_data$model_type <- factor(plot_data$model_type, levels = level_order)
plot_data$model_type <- relevel(plot_data$model_type, ref = "no shift")


#--------------------
# STATS RELATIVE FIT
#--------------------

###### interval MY ######
intv_my <- plot_data[c("model_type", "interval_MY")]

# LMM
lmm_model_time <- lmer(interval_MY ~ model_type + (1| popID), data = plot_data)
#lmm_model_time <- lm(interval_MY ~ model_type, data = plot_data)
summary(lmm_model_time)

# plot
pdf("./results_paleoTS_v0.6.1/plot/interval_my_shift_aicc.pdf")
intv_plot <- ggplot(intv_my, aes(x = interval_MY, y = factor(model_type, levels = level_order),
                                 )) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("Interval (MY)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                axis.title = element_text(size = 13),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot
dev.off()

###### steps ######
steps <- plot_data[c("model_type", "steps")]

# LMM
lmm_model_steps <- lmer(steps ~ model_type + (1| popID), data = plot_data)
#lmm_model_steps <- lm(steps ~ model_type, data = plot_data)
summary(lmm_model_steps)

# plot
pdf("./results_paleoTS_v0.6.1/plot/steps_shift_aicc.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = factor(model_type, levels = level_order),
)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("ln(Steps)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                     axis.title = element_text(size = 13),
                                                     axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res <- plot_data[c("model_type", "steps", "interval_MY")]
res$resolution <- res$steps/res$interval_MY

# LMM
plot_data$resolution = plot_data$steps/plot_data$interval_MY
lmm_model_resolution <- lmer(resolution ~ model_type + (1| popID), data = plot_data)
#lmm_model_resolution <- lm(resolution ~ model_type, data = plot_data)
summary(lmm_model_resolution)

# plot
pdf("./results_paleoTS_v0.6.1/plot/resolution_shift_aicc.pdf")
res_plot <- ggplot(res, aes(x = log(resolution), y = factor(model_type, levels = level_order),
)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("ln(Resolution)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                 axis.title = element_text(size = 13),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()


# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_shift_aicc.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()



# put the output of the linear mixed effects models in the same table
models <- rownames(summary(lmm_model_time)$coefficients)
interval <- as.data.frame(summary(lmm_model_time)$coefficients)
steps <- as.data.frame(summary(lmm_model_steps)$coefficients)
resolution <- as.data.frame(summary(lmm_model_resolution)$coefficients)

interval_df <- data.frame(term = models, i.Est = interval$Estimate,
                          i.SE = interval$"Std", i.p.value = interval$"Pr")
interval_df$i.Est[2:nrow(interval_df)] <- interval_df$i.Est[1] + interval_df$i.Est[2:nrow(interval_df)]

steps_df <- data.frame(s.Est = steps$Estimate,
                       s.SE = steps$"Std", s.p.value = steps$"Pr")
steps_df$s.Est[2:nrow(steps_df)] <- steps_df$s.Est[1] + steps_df$s.Est[2:nrow(steps_df)]

resolution_df <- data.frame(r.Est = resolution$Estimate,
                            r.SE = resolution$"Std", r.p.value = resolution$"Pr")
resolution_df$r.Est[2:nrow(resolution_df)] <- resolution_df$r.Est[1] + resolution_df$r.Est[2:nrow(resolution_df)]

lmm_result_table <- cbind(interval_df, steps_df, resolution_df)
lmm_result_table[, -1] <- round(lmm_result_table[, -1], 3)

write.csv(lmm_result_table, file = "./results_paleoTS_v0.6.1/table_lmm_shift.csv", row.names = FALSE)



#-------------------
# PLOT ABSOLUTE FIT
#-------------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

# ordering the model according to the number of parameters
level_order <- c("no shift", "mode shift")
plot_data$model_type <- factor(plot_data$model_type, levels = level_order)
plot_data$model_type <- relevel(plot_data$model_type, ref = "no shift")


###### interval MY ######
intv_my <- plot_data2[c("model_type", "interval_MY")]

# LMM
lmm_model_time2 <- lmer(interval_MY ~ model_type + (1| popID), data = plot_data2)
#lmm_model_time2 <- lm(interval_MY ~ model_type, data = plot_data2)
summary(lmm_model_time2)

# plot
pdf("./results_paleoTS_v0.6.1/plot/interval_my_shift_adeq.pdf")
intv_plot <- ggplot(intv_my, aes(x = interval_MY, y = factor(model_type, levels = level_order),
)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("Interval (MY)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                     axis.title = element_text(size = 13),
                                                     axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot
dev.off()

###### steps ######
steps <- plot_data2[c("model_type", "steps")]

# LMM
lmm_model_steps2 <- lmer(steps ~ model_type + (1| popID), data = plot_data2)
#lmm_model_steps2 <- lm(steps ~ model_type, data = plot_data2)
summary(lmm_model_steps2)

# plot
pdf("./results_paleoTS_v0.6.1/plot/steps_shift_adeq.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = factor(model_type, levels = level_order),
)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("ln(Steps)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                 axis.title = element_text(size = 13),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res <- plot_data2[c("model_type", "steps", "interval_MY")]
res$resolution <- res$steps/res$interval_MY

# LMM
plot_data2$resolution = plot_data2$steps/plot_data2$interval_MY
lmm_model_resolution2 <- lmer(resolution ~ model_type + (1| popID), data = plot_data2)
#lmm_model_resolution2 <- lm(resolution ~ model_type, data = plot_data2)
summary(lmm_model_resolution2)

# plot
pdf("./results_paleoTS_v0.6.1/plot/resolution_shift_adeq.pdf")
res_plot <- ggplot(res, aes(x = log(resolution), y = factor(model_type, levels = level_order),
)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("no shift", "mode shift")) + 
  xlab("ln(Resolution)") + ylab("Model type") + theme(axis.text = element_text(size = 10),
                                                      axis.title = element_text(size = 13),
                                                      axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()



# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_shift_adeq.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()

# put the output of the linear mixed effects models in the same table
models2 = rownames(summary(lmm_model_time2)$coefficients)
interval2 = as.data.frame(summary(lmm_model_time2)$coefficients)
steps2 = as.data.frame(summary(lmm_model_steps2)$coefficients)
resolution2 = as.data.frame(summary(lmm_model_resolution2)$coefficients)

interval_df2 <- data.frame(term = models2, i.Est = interval2$Estimate,
                           i.SE = interval2$"Std", i.p.value = interval2$"Pr")
interval_df2$i.Est[2:nrow(interval_df2)] <- interval_df2$i.Est[1] + interval_df2$i.Est[2:nrow(interval_df2)]

steps_df2 <- data.frame(s.Est = steps2$Estimate,
                        s.SE = steps2$"Std", s.p.value = steps2$"Pr")
steps_df2$s.Est[2:nrow(steps_df2)] <- steps_df2$s.Est[1] + steps_df2$s.Est[2:nrow(steps_df2)]

resolution_df2 <- data.frame(r.Est = resolution2$Estimate,
                             r.SE = resolution2$"Std", r.p.value = resolution2$"Pr")
resolution_df2$r.Est[2:nrow(resolution_df2)] <- resolution_df2$r.Est[1] + resolution_df2$r.Est[2:nrow(resolution_df2)]

lmm_result_table2 <- cbind(interval_df2, steps_df2, resolution_df2)
lmm_result_table2[, -1] <- round(lmm_result_table2[, -1], 3)

write.csv(lmm_result_table2, file = "./results_paleoTS_v0.6.1/table_lmm_shift_adeq.csv", row.names = FALSE)


#---------------------------------
# PLOT STATISTICS ON THE DATASET
#---------------------------------

# bind data to dataframe
unit_list <- c("popID", "taxa", "period_start", "steps", "interval_MY")
plot_dataset <- bind(ln_data_meta_shift, unit_list)
plot_dataset$resolution <- plot_dataset$steps/plot_dataset$interval_MY

###### Taxa plot ###### 
plot_dataset$taxa <- replace(plot_dataset$taxa, plot_dataset$taxa == "chondrichthyan", "fish")

unique(plot_dataset$taxa) #no cephalopod, no echinoderm, no bird

taxa_levels <- c(
  "foraminifer", "coccolith", "radiolarian", "diatom",
  "bryozoan", "bivalve", "gastropod", "ostracod", "brachiopod", "trilobite", "graptolite",
  "mammal", "conodont", "fish"
) 

taxa_cols <- c(
  # Protists
  foraminifer = "#598B8C",
  coccolith    = "#719EA0",
  radiolarian  = "#99C0C2",
  diatom       = "#BFE0E1",
  # Invertebrates
  bryozoan    = "#968F2C",
  bivalve     = "#A89F33",
  gastropod   = "#B7AF3B",
  #cephalopod  = "#C6BF43",
  ostracod    = "#DCCB4E",
  brachiopod  = "#E3D460",
  trilobite   = "#EADC72",
  #echinoderm  = "#F1E583",
  graptolite  = "#F8EC95",
  # Vertebrates
  mammal    = "#B84400",
  #bird     = "#DB6000",
  conodont = "#F07C26",
  fish      = "#FF983D"
)

plot_dataset %>%
  filter(!is.na(taxa)) %>%
  count(taxa, name = "n") %>%
  mutate(
    fraction = round(n / sum(n), 3),
    pct = percent(fraction),
    ypos = cumsum(fraction) - fraction/2
  ) -> df_taxa


df_taxa <- df_taxa %>%
  mutate(taxa = factor(taxa, levels = taxa_levels))  
used_levels <- levels(df_taxa$taxa)                   
used_cols   <- taxa_cols[used_levels]

percentages <- setNames(df_taxa$pct, df_taxa$taxa)
legend_labels <- paste0(used_levels, " (", percentages[used_levels], ")")

taxa_dataset_plot <- ggplot(df_taxa, aes(x = 1, y = fraction, fill = taxa)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  coord_polar(theta = "y", direction = -1) +
  scale_fill_manual(values = used_cols, labels = legend_labels) +
  theme_void() +
  theme(
    legend.key.spacing.y = unit(0.1, "cm"),
  ) +
  labs(title = "Taxa", fill = "taxa") +
  guides(fill=guide_legend(ncol=2, byrow=FALSE))

###### Age plot ###### 
unique(plot_dataset$period_start) #no permian, no triassic

period_levels <- c(
  "Cambrian", "Ordovician", "Silurian", "Devonian",
  "Carboniferous", "Jurassic", "Cretaceous",
  "Paleogene", "Neogene", "Quaternary"
)

period_cols <- c(
  "Cambrian"      = rgb(127, 160, 86, maxColorValue = 255),
  "Ordovician"    = rgb(0, 146, 112, maxColorValue = 255),
  "Silurian"      = rgb(179, 225, 182, maxColorValue = 255),
  "Devonian"      = rgb(203, 140, 55, maxColorValue = 255),
  "Carboniferous" = rgb(103, 165, 153, maxColorValue = 255),
  #"Permian"       = rgb(240, 64, 40, maxColorValue = 255),
  #"Triassic"      = rgb(129, 43, 146, maxColorValue = 255),
  "Jurassic"      = rgb(52, 178, 201, maxColorValue = 255),
  "Cretaceous"    = rgb(127, 198, 78, maxColorValue = 255),
  "Paleogene"     = rgb(253, 154, 82, maxColorValue = 255),
  "Neogene"       = rgb(255, 230, 25, maxColorValue = 255),
  "Quaternary"    = rgb(249, 249, 127, maxColorValue = 255)
)

plot_dataset %>%
  filter(!is.na(period_start)) %>%
  mutate(period_start = factor(period_start, levels = period_levels)) %>%
  count(period_start, name = "n") %>%
  filter(!is.na(period_start)) %>%
  mutate(
    fraction = n / sum(n),
    pct = percent(fraction),
    ypos = cumsum(fraction) - fraction / 2
  ) -> df_periods

used_levels <- levels(df_periods$period_start)
used_cols <- period_cols[used_levels]

percentages <- setNames(df_periods$pct, df_periods$period_start)
legend_labels <- paste0(used_levels, " (", percentages[used_levels], ")")

age_dataset_plot <- ggplot(df_periods, aes(x = 1, y = fraction, fill = period_start)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  coord_polar(theta = "y", direction = -1) +
  scale_fill_manual(values = used_cols, labels = legend_labels) +
  theme_void() +
  theme(
    legend.key.spacing.y = unit(0.1, "cm"),
  ) +
  labs(title = "Geological period", fill = "Period")

###### interval plot ######
intv_dataset_plot <- ggplot(plot_dataset, aes(x = interval_MY)) +
  geom_histogram(bins = 20, color = "black", fill = "#8D8680", linewidth = 0.2) +
  labs(x = "Interval (My)",
       y = "Time series count") +
  theme_classic()

###### steps plot ######
steps_dataset_plot <- ggplot(plot_dataset, aes(x = steps)) +
  geom_histogram(bins = 20, color = "black", fill = "#8D8680", linewidth = 0.2) +
  scale_x_log10() +
  labs(x = "Steps",
       y = "Time series count") +
  theme_classic()

###### resolution plot ######
res_dataset_plot <- ggplot(plot_dataset, aes(x = resolution)) +
  geom_histogram(bins = 20, color = "black", fill = "#8D8680", linewidth = 0.2) +
  scale_x_log10() +
  labs(x = "Resolution",
       y = "Time series count") +
  theme_classic()

# save the dataset figure
plot_dataset_final = list(taxa_dataset_plot, age_dataset_plot, intv_dataset_plot, steps_dataset_plot, res_dataset_plot)

plot_dataset_display = grid.arrange(
  grobs = plot_dataset_final,
  widths = c(1, 5, 5, 5, 5, 5, 5, 1),
  heights = c(1, 10, 1, 8, 1),
  layout_matrix = rbind(c(NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, 1, 1, 1, 2, 2, 2, NA),
                        c(NA, NA, NA, NA, NA, NA, NA, NA),
                        c(NA, 3, 3, 4, 4, 5, 5, NA),
                        c(NA, NA, NA, NA, NA, NA, NA, NA))
)

ggsave("./results_paleoTS_v0.6.1/plot/dataset_shift_v2.pdf", plot_dataset_display,
       width = 9, height = 6, units = "in", dpi = 300)

