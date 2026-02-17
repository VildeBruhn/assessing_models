#########################################
##   Assessing models of evolution     ##
##            STATISTICS               ##
#########################################

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
library(deeptime)

rm(list = ls())

# set working directory
setwd("/Users/vildebruhnkinneberg/Documents/GitHub/assessing_models_evolution/assessing_models/")
source("./assessing_models_uni_functions.R")


#--------------
# IMPORT FILES
#--------------
df <- read_delim("./timeseries/metadata.txt", col_names = TRUE, delim = "\t")

# remove time series with less than 10 steps
df <- subset(df, steps >= 7)

# remove modern time series
df <- subset(df, period_start != "Present")

# remove Syverson
df <- subset(df, URL != "https://doi.org/10.1017/pab.2024.37")

# make list based on ID
df <- lapply(split(df,df$tsID), function(x) as.list(x))


# load data from analyses
load("./aicc_uni_passed.Rdata")
load("./model_test_uni.Rdata")
load("./adeq_uni_passed.Rdata")


#--------------
# PREPARE DATA
#--------------


# get aicc model info into metadata
df <- model_aicc(df, GRW, "GRW")
df <- model_aicc(df, URW, "URW")
df <- model_aicc(df, stasis, "stasis")
df <- model_aicc(df, strict_stasis, "strict stasis")
df <- model_aicc(df, decel, "decel")
df <- model_aicc(df, accel, "accel")
df <- model_aicc(df, OU, "OU")
df <- model_aicc(df, OU_mov_opt_anc, "OU mov opt anc")
df <- model_aicc(df, OU_mov_opt, "OU mov opt")

# get adequate model info into metadata
df <- model_adeq(df, GRW_adeq_passed, "GRW")
df <- model_adeq(df, URW_adeq_passed, "URW")
df <- model_adeq(df, stasis_adeq_passed, "stasis")
df <- model_adeq(df, strict_stasis_adeq_passed, "strict stasis")
df <- model_adeq(df, decel_adeq_passed, "decel")
df <- model_adeq(df, accel_adeq_passed, "accel")
df <- model_adeq(df, OU_adeq_passed, "OU")
df <- model_adeq(df, OU_mov_opt_anc_adeq_passed, "OU mov opt anc")
df <- model_adeq(df, OU_mov_opt_adeq_passed, "OU mov opt")

# bind data to dataframe
unit_list <- c("stID", "popID", "species", "description", "total_N", "steps", 
               "interval_MY", "model_aicc", "model_adequate")
plot_data <- bind(df, unit_list)

# remove time series with NA
plot_data <- plot_data %>% drop_na(model_aicc)

# collapse strict stasis to stasis
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "strict stasis", "stasis")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "strict stasis", "stasis")

# collapse OU mov opt anc to OU mov opt
plot_data$model_aicc <- replace(plot_data$model_aicc, plot_data$model_aicc == "OU mov opt anc", "OU mov opt")
plot_data$model_adequate <- replace(plot_data$model_adequate, plot_data$model_adequate == "OU mov opt anc", "OU mov opt")


# add parameters column
plot_data$parameters <- plot_data$model_aicc
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "stasis", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "URW", 2)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "GRW", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "accel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "decel", 3)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU", 4)
plot_data$parameters <- replace(plot_data$parameters, plot_data$parameters == "OU mov opt", 5)

# ordering the model according to the number of parameters
level_order <- c("stasis", "URW", "GRW", "accel", "decel", "OU", "OU mov opt")
plot_data$model_aicc <- factor(plot_data$model_aicc, levels = level_order)


# set colors
col_val1 <- c("#85B7B9", "#DCCB4E")
col_val2 <- c("#F11B00","#E5A208", "#ADC397", "#3A9AB2")
col_val_extra <- c("#3A9AB2", "#85B7B9", "#ADC397", "#DCCB4E", "#E5A208", "#ED6E04", "#F11B00")


#---------------
# RELATIVE FIT
#---------------


###### interval MY ######
intv_my <- plot_data

# LMM
lmm_intv <- lmer(interval_MY ~ 1 + model_aicc + (1| popID), data = intv_my)
summary(lmm_intv)

sink(file = "./results_paleoTS_v0.6.1/lmm_interval_uni_results_aicc.txt")
summary(lmm_intv)
sink()

# stats 
intv_my_mean <- aggregate(intv_my$interval_MY, list(intv_my$model_aicc), mean)
intv_my_mean <- as.data.frame(intv_my_mean)
names(intv_my_mean) <- c("model", "mean interval (MY)")
intv_my_median <- aggregate(intv_my$interval_MY, list(intv_my$model_aicc), median)
intv_my_median <- as.data.frame(intv_my_median)
names(intv_my_median) <- c("model", "median interval (MY)")

sink(file = "./results_paleoTS_v0.6.1/interval_my_results_aicc.txt")
intv_my_mean
intv_my_median
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/interval_my_uni_aicc.pdf")
intv_plot <- ggplot(intv_my, aes(x = interval_MY, y = model_aicc,
                    fill = parameters)) +
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                axis.title = element_text(size = 13),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot
dev.off()

###### steps ######
steps <- plot_data

# stats
steps_mean <- aggregate(steps$interval_MY, list(steps$model_aicc), mean)
steps_mean <- as.data.frame(steps_mean)
names(steps_mean) <- c("model", "mean steps")
steps_median <- aggregate(steps$interval_MY, list(steps$model_aicc), median)
steps_median <- as.data.frame(steps_median)
names(steps_median) <- c("model", "median steps")

sink(file = "./results_paleoTS_v0.6.1/steps_uni_results_aicc.txt")
steps_mean
steps_median
sink()

# LMM
lmm_steps <- lmer(steps ~ model_aicc + (1| popID), data = steps)
summary(lmm_steps)

sink(file = "./results_paleoTS_v0.6.1/lmm_steps_uni_results_aicc.txt")
summary(lmm_steps)
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/steps_uni_aicc.pdf")
steps_plot <- ggplot(steps, aes(x = log(steps), y = model_aicc, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                        axis.title = element_text(size = 13),
                                        axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res <- plot_data
res$resolution <- res$steps/res$interval_MY

# LMM
lmm_res <- lmer(resolution ~ 1 + model_aicc + (1| popID), data = res)
summary(lmm_res)

sink(file = "./results_paleoTS_v0.6.1/lmm_res_uni_results_aicc.txt")
summary(lmm_res)
sink()

# stats
res_mean <- aggregate(res$resolution, list(res$model_aicc), mean)
res_mean <- as.data.frame(res_mean)
names(res_mean) <- c("model", "mean res")
res_median <- aggregate(res$resolution, list(res$model_aicc), median)
res_median <- as.data.frame(res_median)
names(res_median) <- c("model", "median res")

sink(file = "./results_paleoTS_v0.6.1/resolution_uni_results_aicc.txt")
res_mean
res_median
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/resolution_uni_aicc.pdf")
res_plot <- ggplot(res, aes(x = log(resolution), y = model_aicc, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                        axis.title = element_text(size = 13),
                                                        axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()



# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_aicc.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()



# put the output of the linear mixed effects models in the same table
models <- rownames(summary(lmm_intv)$coefficients)
interval <- as.data.frame(summary(lmm_intv)$coefficients)
steps <- as.data.frame(summary(lmm_steps)$coefficients)
resolution <- as.data.frame(summary(lmm_res)$coefficients)

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

write.csv(lmm_result_table, file = "./results_paleoTS_v0.6.1/table_lmm_uni.csv", row.names = FALSE)



#---------------
# ABSOLUTE FIT
#---------------


# remove NA data
plot_data2 <- plot_data 
plot_data2 <- plot_data2 %>% drop_na(model_adequate)

# ordering the model according to the number of parameters
plot_data2$model_adequate <- factor(plot_data2$model_adequate, levels = level_order)

###### interval MY ######
intv_my2 <- plot_data2


# stats
intv_my2_mean <- aggregate(intv_my2[, 2], list(intv_my2$model_adequate), mean)
intv_my2_mean <- as.data.frame(intv_my2_mean)
names(intv_my2_mean) <- c("model", "mean interval (MY)")
intv_my2_median <- aggregate(intv_my2[, 2], list(intv_my2$model_adequate), median)
intv_my2_median <- as.data.frame(intv_my2_median)
names(intv_my2_median) <- c("model", "median interval (MY)")

sink(file = "./results_paleoTS_v0.6.1/interval_my_results_adeq.txt")
intv_my2_mean
intv_my2_median
sink()

# LMM
lmm_intv2 <- lmer(interval_MY ~ 1+ model_adequate + (1| popID), data = intv_my2)
summary(lmm_intv2)

sink(file = "./results_paleoTS_v0.6.1/lmm_interval_adeq_uni_results.txt")
summary(lmm_intv2)
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/interval_my_uni_adeq.pdf")
intv_plot2 <- ggplot(intv_my2, aes(x = interval_MY, y = model_adequate,
                                 fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) + 
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("Interval (MY)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                axis.title = element_text(size = 13),
                                                axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
intv_plot2
dev.off()

###### steps ######
steps2 <- plot_data2

# stats
steps2_mean <- aggregate(steps2[, 2], list(steps2$model_adequate), mean)
steps2_mean <- as.data.frame(steps2_mean)
names(steps2_mean) <- c("model", "mean steps2")
steps2_median <- aggregate(steps2[, 2], list(steps2$model_adequate), median)
steps2_median <- as.data.frame(steps2_median)
names(steps2_median) <- c("model", "median steps2")

sink(file = "./results_paleoTS_v0.6.1/steps_uni_results_adeq.txt")
steps2_mean
steps2_median
sink()

# LMM
lmm_steps2 <- lmer(steps ~ model_adequate + (1| popID), data = steps2)
summary(lmm_steps2)

sink(file = "./results_paleoTS_v0.6.1/lmm_steps_adeq_uni_results.txt")
summary(lmm_steps2)
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/steps_uni_adeq.pdf")
steps_plot2 <- ggplot(steps2, aes(x = log(steps), y = model_adequate, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Steps)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                            axis.title = element_text(size = 13),
                                            axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
steps_plot
dev.off()

###### resolution ######
res2 <- plot_data2
res2$resolution <- res2$steps/res2$interval_MY

# stats
res2_mean <- aggregate(res2[, 5], list(res2$model_adequate), mean)
res2_mean <- as.data.frame(res2_mean)
names(res2_mean) <- c("model", "mean res")
res2_median <- aggregate(res2[, 5], list(res2$model_adequate), median)
res2_median <- as.data.frame(res2_median)
names(res2_median) <- c("model", "median res")

sink(file = "./results_paleoTS_v0.6.1/resolution_uni_results_adeq.txt")
res2_mean
res2_median
sink()

# LMM
lmm_res2 <- lmer(resolution ~ 1+ model_aicc + (1| popID), data = res2)
summary(lmm_res2)

sink(file = "./results_paleoTS_v0.6.1/lmm_res_adeq_uni_results.txt")
summary(lmm_res2)
sink()

# plot
pdf("./results_paleoTS_v0.6.1/plot/resolution_uni_adeq.pdf")
res_plot <- ggplot(res2, aes(x = log(resolution), y = model_adequate, fill = parameters)) + 
  geom_boxplot() + theme_classic() + scale_y_discrete(labels = c("Stasis", "URW", 
                                                                 "GRW", "Accel.", "Decel.", 
                                                                 "OU", "OU mov. opt.")) +
  scale_fill_discrete(name = "Parameters", labels = c("2", "3", "4",
                                                      "5"), palette = col_val2) +
  xlab("ln(Resolution)") + ylab("Model") + theme(axis.text = element_text(size = 10),
                                                 axis.title = element_text(size = 13),
                                                 axis.title.x = element_text(margin = margin(t = 17, r = 0, b = 0, l = 0)))
res_plot
dev.off()



# put plots in same figure
pdf(width = 15.0, height = 5.5, file = "./results_paleoTS_v0.6.1/plot/interval_steps_res_adeq.pdf")
grid.arrange(intv_plot,
             steps_plot,
             res_plot, nrow = 1)
dev.off()

# put the output of the linear mixed effects models in the same table
models2 <- rownames(summary(lmm_intv2)$coefficients)
interval2 <- as.data.frame(summary(lmm_intv2)$coefficients)
steps2 <- as.data.frame(summary(lmm_steps2)$coefficients)
resolution2 <- as.data.frame(summary(lmm_res2)$coefficients)

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

write.csv(lmm_result_table2, file = "./results_paleoTS_v0.6.1/table_lmm_uni_adeq.csv", row.names = FALSE)


###### DAICc gap ######
# Extract the AICc results
aicc <- lapply(model_test, function(x) x[(names(x) %in% c("AICc"))])

# Extract the lowest AICc results for each model
aicc_min <- lapply(aicc, function(x) {
  which.min(as.numeric(unlist(x)))
})

# Calculate DeltaAICc between the best and second best model
aicc_gap <- list()
aicc_mingap <- list()
for (i in 1:length(aicc)) {
  for (j in 1:nrow(aicc[[1]][])) {
    if (aicc[[i]][j,] != aicc[[i]][aicc_min[[i]],]) {
      aicc_gap[[j]] = abs(aicc[[i]][j,] - aicc[[i]][aicc_min[[i]],])
    } 
  }
  aicc_mingap[i] <- min(unlist(aicc_gap))
  names(aicc_mingap)[i] <- names(aicc_min)[i]
}

plot_data$deltaAICc <- unlist(aicc_mingap[plot_data$data_frame])
plot_data$adequacy <- ifelse(is.na(plot_data$model_adequate), "inadequate", "adequate")

adequacy_order <- c("adequate", "inadequate")
plot_data$adequacy <- factor(plot_data$adequacy, levels = adequacy_order)

# stats
#How many time series have a deltaAICc inferior or equal to 2
TS_deltaAICc2_c <- sum(plot_data$deltaAICc <= 2)
TS_deltaAICc2_p <- TS_deltaAICc2_c/nrow(plot_data)*100

# LMM
lmm_model_deltaAICc <- lmer(deltaAICc ~ 1+ adequacy + (1| popID), data = plot_data)
#lmm_model_deltaAICc <- lm(deltaAICc ~ adequacy, data = plot_data)
summary(lmm_model_deltaAICc)

sink(file = "./results_paleoTS_v0.6.1/deltaAICc_uni_results_adeq.txt")
round(summary(lmm_model_deltaAICc)$coefficient, 2)
paste("Total number of time series investigated:", nrow(plot_data))
paste("Time series with DAICc <= 2:   ", TS_deltaAICc2_c, "   ", TS_deltaAICc2_p, "%")
sink()

# plot
pdf(width = 5.0 , height = 4.5, file = "./results_paleoTS_v0.6.1/plot/DAICc_uni_adeq.pdf")
Daicc_plot = ggplot(plot_data, aes(x = adequacy, y = deltaAICc, fill = adequacy)) +
  geom_boxplot() +
  labs(x = expression(bold("Adequacy status"))) +
  theme_classic() + ylab(~ paste(bold(Delta), bold(" AICc gap"), " (second best model - first best model)")) +
  scale_y_continuous(limits = c(0, 8)) +
  scale_fill_manual(values = c("adequate" = "#ADC397", "inadequate" = "#E5A208"), name = "") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.text = element_text(size = 11))
Daicc_plot
dev.off()



#---------------------------------
# PLOT STATISTICS ON THE DATASET
#---------------------------------

# bind data to dataframe
unit_list <- c("popID", "taxa", "period_start", "steps", "interval_MY")
plot_dataset <- bind(df, unit_list)
plot_dataset$resolution <- plot_dataset$steps/plot_dataset$interval_MY

###### Taxa plot ######
plot_dataset$taxa <- replace(plot_dataset$taxa, plot_dataset$taxa == "chondrichthyan", "fish")

taxa_levels <- c(
  "foraminifer", "coccolith", "radiolarian", "diatom",
  "bryozoan", "bivalve", "gastropod", "cephalopod", "ostracod", "brachiopod", "trilobite", "echinoderm", "graptolite",
  "mammal", "bird", "conodont", "fish"
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
  cephalopod  = "#C6BF43",
  ostracod    = "#DCCB4E",
  brachiopod  = "#E3D460",
  trilobite   = "#EADC72",
  echinoderm  = "#F1E583",
  graptolite  = "#F8EC95",
# Vertebrates
  mammal    = "#B84400",
  bird     = "#DB6000",
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
  labs(title = "Taxa", fill = "") +
  guides(fill=guide_legend(ncol=2, byrow=FALSE))

###### Age plot ###### 
period_levels <- c(
  "Cambrian", "Ordovician", "Silurian", "Devonian",
  "Carboniferous", "Permian", "Triassic", "Jurassic", "Cretaceous",
  "Paleogene", "Neogene", "Quaternary"
)

period_cols <- c(
  "Cambrian"      = rgb(127, 160, 86, maxColorValue = 255),
  "Ordovician"    = rgb(0, 146, 112, maxColorValue = 255),
  "Silurian"      = rgb(179, 225, 182, maxColorValue = 255),
  "Devonian"      = rgb(203, 140, 55, maxColorValue = 255),
  "Carboniferous" = rgb(103, 165, 153, maxColorValue = 255),
  "Permian"       = rgb(240, 64, 40, maxColorValue = 255),
  "Triassic"      = rgb(129, 43, 146, maxColorValue = 255),
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

age_dataset_plot <- ggplot(df_periods, aes(x = 1, y = fraction, fill = period_start)) +
  geom_col(width = 1, color = "black", linewidth = 0.2) +
  coord_polar(theta = "y", direction = -1) +
  scale_fill_manual(values = used_cols, labels = legend_labels) +
  theme_void() +
  theme(
    legend.key.spacing.y = unit(0.1, "cm"),
  ) +
  labs(title = "Geological period", fill = "")

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
       y = "") +
  theme_classic()

###### resolution plot ######
res_dataset_plot <- ggplot(plot_dataset, aes(x = resolution)) +
  geom_histogram(bins = 20, color = "black", fill = "#8D8680" , linewidth = 0.2) +
  scale_x_log10() +
  labs(x = "Resolution",
       y = "") +
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

ggsave("./results_paleoTS_v0.6.1/plot/dataset_uni_v3.pdf", plot_dataset_display,
       width = 9, height = 6, units = "in", dpi = 300)
