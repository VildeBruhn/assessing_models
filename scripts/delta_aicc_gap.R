####################################################
##    Assessing models of phenotypic evolution    ##
##                                                ##
##                 DELTA AICc GAP                 ##
##                                                ##
##                Marion Thaureau                 ##
##              Created 01.06.2026                ##
####################################################

# install packages
#install.packages("tidyverse")
#install.packages("lme4")

# load packages
library(tidyverse)
library(lme4)

rm(list = ls())

# set working directory
PATH = "[PATH_TO_MAIN_FOLDER]"
setwd(PATH)

# import functions
source("./scripts/functions.R")

###### DAICc gap ######

# import files
load("./data/model_test_single_mode.Rdata")
load("./data/plot_data.Rdata")

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

# add to data frame
plot_data$deltaAICc <- unlist(aicc_mingap[plot_data$data_frame])
plot_data$adequacy <- ifelse(is.na(plot_data$model_adequate), "inadequate", "adequate")

adequacy_order <- c("adequate", "inadequate")
plot_data$adequacy <- factor(plot_data$adequacy, levels = adequacy_order)

# get number og time series with deltaAICc lower or equal to 2
TS_deltaAICc2_c <- sum(plot_data$deltaAICc <= 2)
TS_deltaAICc2_p <- TS_deltaAICc2_c/nrow(plot_data)*100

# LMM
lmm_model_deltaAICc <- lmer(deltaAICc ~ 1+ adequacy + (1| popID), data = plot_data)
summary(lmm_model_deltaAICc)

# plot
#pdf(width = 5.0 , height = 4.5, file = "[PATH_TO_OUTPUT_FOLDER]/DAICc_gap.pdf")
Daicc_plot = ggplot(plot_data, aes(x = adequacy, y = deltaAICc, fill = adequacy)) + 
  geom_boxplot() + labs(x = "") + theme_classic() + ylab(~ paste(bold(Delta), " AICc gap")) +
  scale_y_continuous(limits = c(0, 8)) + scale_x_discrete(labels = c("Adequate", "Inadequate"), ) +
  scale_fill_manual(values = c("#99C0C2", "#F1E583")) + guides(fill="none") +
  #scale_fill_manual(values = c("adequate" = "#99C0C2", "inadequate" = "#F1E583"), name = "") +
  theme(axis.title = element_text(size = 12), axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(size = 12), legend.text = element_text(size = 11))
Daicc_plot
#dev.off()
