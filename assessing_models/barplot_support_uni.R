######################
## Code for Barplot ##
######################

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2


# importing results of test without shift and TS with minimum 10 datapoints
load("./model_test_uni.Rdata")
metadata_trait <- read_delim("./timeseries/metadata_trait.txt", col_names = TRUE, delim = "\t")


#-------------------------------------
# Plot of the support for each model
#-------------------------------------

### extracting AICc weight list for each time series ###
support_Stasis <- lapply(model_test, function(df) {
  return(df["Stasis", "Akaike.wt"])
})

support_URW <- lapply(model_test, function(df) {
  return(df["URW", "Akaike.wt"])
})

support_GRW <- lapply(model_test, function(df) {
  return(df["GRW", "Akaike.wt"])
})

support_Strict_stasis <- lapply(model_test, function(df) {
  return(df["StrictStasis", "Akaike.wt"])
})

support_Decel <- lapply(model_test, function(df) {
  return(df["Decel", "Akaike.wt"])
})

support_Accel <- lapply(model_test, function(df) {
  return(df["Accel", "Akaike.wt"])
})

support_OU <- lapply(model_test, function(df) {
  return(df["OU", "Akaike.wt"])
})

support_OU_mov_opt_anc <- lapply(model_test, function(df) {
  return(df["OU model with moving optimum (ancestral state at optimum)", "Akaike.wt"])
})

support_OU_mov_opt <- lapply(model_test, function(df) {
  return(df["OU model with moving optimum", "Akaike.wt"])
})


support_Stasis_values <- unlist(support_Stasis)
df_Stasis <- data.frame(values = support_Stasis_values)

hist_Stasis = ggplot(df_Stasis, aes(x = values)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Stasis", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0))) + 
  theme_minimal()


support_URW_values <- unlist(support_URW)
df_URW <- data.frame(values = support_URW_values)

hist_URW = ggplot(df_URW, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for URW", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_GRW_values <- unlist(support_GRW)
df_GRW <- data.frame(values = support_GRW_values)

hist_GRW = ggplot(df_GRW, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for GRW", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 300), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Strict_stasis_values <- unlist(support_Strict_stasis)
df_Strict_stasis <- data.frame(values = support_Strict_stasis_values)

hist_Strict_stasis = ggplot(df_Strict_stasis, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Strict_stasis", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 400), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Accel_values <- unlist(support_Accel)
df_Accel <- data.frame(values = support_Accel_values)

hist_Accel = ggplot(df_Accel, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Accel", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 400), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Decel_values <- unlist(support_Decel)
df_Decel <- data.frame(values = support_Decel_values)

hist_Decel = ggplot(df_Decel, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Decel", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 400), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_values <- unlist(support_OU)
df_OU <- data.frame(values = support_OU_values)

hist_OU = ggplot(df_OU, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 500), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_anc_values <- unlist(support_OU_mov_opt_anc)
df_OU_mov_opt_anc <- data.frame(values = support_OU_mov_opt_anc_values)

hist_OU_mov_opt_anc = ggplot(df_OU_mov_opt_anc, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt_anc", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 500), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_values <- unlist(support_OU_mov_opt)
df_OU_mov_opt <- data.frame(values = support_OU_mov_opt_values)

hist_OU_mov_opt = ggplot(df_OU_mov_opt, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 500), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()

# Set common theme for the plots
custom_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),  # Title size and space
  axis.title.x = element_text(size = 14),  # X-axis title size
  axis.title.y = element_text(size = 14),  # Y-axis title size
  axis.text.x = element_text(size = 12),   # X-axis text size
  axis.text.y = element_text(size = 12),   # Y-axis text size
  panel.spacing = unit(2, "lines")  # Space between plots
)

# Apply the custom theme to each histogram
hist_Stasis <- hist_Stasis + custom_theme
hist_URW <- hist_URW + custom_theme
hist_GRW <- hist_GRW + custom_theme
hist_Strict_stasis <- hist_Strict_stasis + custom_theme
hist_Accel <- hist_Accel + custom_theme
hist_Decel <- hist_Decel + custom_theme
hist_OU <- hist_OU + custom_theme
hist_OU_mov_opt_anc <- hist_OU_mov_opt_anc + custom_theme
hist_OU_mov_opt <- hist_OU_mov_opt + custom_theme


# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_support_uni_barplot.png", width = 2400, height = 1800)
grid.arrange(hist_Stasis, hist_URW, hist_GRW, hist_Strict_stasis, hist_Accel, hist_Decel, hist_OU, hist_OU_mov_opt_anc, hist_OU_mov_opt, ncol = 3)
dev.off()


#----------------------------------------------------------------
# Barplot of the support for each model depending on type of trait
#----------------------------------------------------------------

metadatashort <- metadata_trait[metadata_trait$tsID %in% names(model_test), ]

# Separate TS of size and shape
metadata_size <- metadatashort[metadatashort$trait_category %in% c("size", "area"), ]
metadata_shape <- metadatashort[metadatashort$trait_category %in% c("shape", "angle"), ]

model_test_size <- model_test[names(model_test) %in% metadata_size$tsID]
model_test_shape <- model_test[names(model_test) %in% metadata_shape$tsID]


###### FOR SIZE TRAITS ######
### extracting AICc weight list for each time series ###
support_Stasis_size <- lapply(model_test_size, function(df) {
  return(df["Stasis", "Akaike.wt"])
})

support_URW_size <- lapply(model_test_size, function(df) {
  return(df["URW", "Akaike.wt"])
})

support_GRW_size <- lapply(model_test_size, function(df) {
  return(df["GRW", "Akaike.wt"])
})

support_Strict_stasis_size <- lapply(model_test_size, function(df) {
  return(df["StrictStasis", "Akaike.wt"])
})

support_Decel_size <- lapply(model_test_size, function(df) {
  return(df["Decel", "Akaike.wt"])
})

support_Accel_size <- lapply(model_test_size, function(df) {
  return(df["Accel", "Akaike.wt"])
})

support_OU_size <- lapply(model_test_size, function(df) {
  return(df["OU", "Akaike.wt"])
})

support_OU_mov_opt_anc_size <- lapply(model_test_size, function(df) {
  return(df["OU model with moving optimum (ancestral state at optimum)", "Akaike.wt"])
})

support_OU_mov_opt_size <- lapply(model_test_size, function(df) {
  return(df["OU model with moving optimum", "Akaike.wt"])
})

support_Stasis_values_size <- unlist(support_Stasis_size)
df_Stasis_size <- data.frame(values = support_Stasis_values_size)

hist_Stasis_size = ggplot(df_Stasis_size, aes(x = values)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Stasis for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) + 
  theme_minimal()


support_URW_values_size <- unlist(support_URW_size)
df_URW_size <- data.frame(values = support_URW_values_size)

hist_URW_size = ggplot(df_URW_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for URW for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_GRW_values_size <- unlist(support_GRW_size)
df_GRW_size <- data.frame(values = support_GRW_values_size)

hist_GRW_size = ggplot(df_GRW_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for GRW for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Strict_stasis_values_size <- unlist(support_Strict_stasis_size)
df_Strict_stasis_size <- data.frame(values = support_Strict_stasis_values_size)

hist_Strict_stasis_size = ggplot(df_Strict_stasis_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Strict_stasis for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Accel_values_size <- unlist(support_Accel_size)
df_Accel_size <- data.frame(values = support_Accel_values_size)

hist_Accel_size = ggplot(df_Accel_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Accel for size trait", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 200), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Decel_values_size <- unlist(support_Decel_size)
df_Decel_size <- data.frame(values = support_Decel_values_size)

hist_Decel_size = ggplot(df_Decel_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Decel for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 200), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_values_size <- unlist(support_OU_size)
df_OU_size <- data.frame(values = support_OU_values_size)

hist_OU_size = ggplot(df_OU_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_anc_values_size <- unlist(support_OU_mov_opt_anc_size)
df_OU_mov_opt_anc_size <- data.frame(values = support_OU_mov_opt_anc_values_size)

hist_OU_mov_opt_anc_size = ggplot(df_OU_mov_opt_anc_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt_anc for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_values_size <- unlist(support_OU_mov_opt_size)
df_OU_mov_opt_size <- data.frame(values = support_OU_mov_opt_values_size)

hist_OU_mov_opt_size = ggplot(df_OU_mov_opt_size, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue4", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt for size traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 250), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


###### FOR SHAPE TRAITS ######
### extracting AICc weight list for each time series ###
support_Stasis_shape <- lapply(model_test_shape, function(df) {
  return(df["Stasis", "Akaike.wt"])
})

support_URW_shape <- lapply(model_test_shape, function(df) {
  return(df["URW", "Akaike.wt"])
})

support_GRW_shape <- lapply(model_test_shape, function(df) {
  return(df["GRW", "Akaike.wt"])
})

support_Strict_stasis_shape <- lapply(model_test_shape, function(df) {
  return(df["StrictStasis", "Akaike.wt"])
})

support_Decel_shape <- lapply(model_test_shape, function(df) {
  return(df["Decel", "Akaike.wt"])
})

support_Accel_shape <- lapply(model_test_shape, function(df) {
  return(df["Accel", "Akaike.wt"])
})

support_OU_shape <- lapply(model_test_shape, function(df) {
  return(df["OU", "Akaike.wt"])
})

support_OU_mov_opt_anc_shape <- lapply(model_test_shape, function(df) {
  return(df["OU model with moving optimum (ancestral state at optimum)", "Akaike.wt"])
})

support_OU_mov_opt_shape <- lapply(model_test_shape, function(df) {
  return(df["OU model with moving optimum", "Akaike.wt"])
})


support_Stasis_values_shape <- unlist(support_Stasis_shape)
df_Stasis_shape <- data.frame(values = support_Stasis_values_shape)

hist_Stasis_shape = ggplot(df_Stasis_shape, aes(x = values)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine2", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Stasis for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) + 
  theme_minimal()


support_URW_values_shape <- unlist(support_URW_shape)
df_URW_shape <- data.frame(values = support_URW_values_shape)

hist_URW_shape = ggplot(df_URW_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine2", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for URW for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_GRW_values_shape <- unlist(support_GRW_shape)
df_GRW_shape <- data.frame(values = support_GRW_values_shape)

hist_GRW_shape = ggplot(df_GRW_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine2", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for GRW for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Strict_stasis_values_shape <- unlist(support_Strict_stasis_shape)
df_Strict_stasis_shape <- data.frame(values = support_Strict_stasis_values_shape)

hist_Strict_stasis_shape = ggplot(df_Strict_stasis_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "aquamarine2", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Strict_stasis for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Accel_values_shape <- unlist(support_Accel_shape)
df_Accel_shape <- data.frame(values = support_Accel_values_shape)

hist_Accel_shape = ggplot(df_Accel_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue3", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Accel for shape trait", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_Decel_values_shape <- unlist(support_Decel_shape)
df_Decel_shape <- data.frame(values = support_Decel_values_shape)

hist_Decel_shape = ggplot(df_Decel_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue3", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for Decel for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_values_shape <- unlist(support_OU_shape)
df_OU_shape <- data.frame(values = support_OU_values_shape)

hist_OU_shape = ggplot(df_OU_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue3", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_anc_values_shape <- unlist(support_OU_mov_opt_anc_shape)
df_OU_mov_opt_anc_shape <- data.frame(values = support_OU_mov_opt_anc_values_shape)

hist_OU_mov_opt_anc_shape = ggplot(df_OU_mov_opt_anc_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue3", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt_anc for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()


support_OU_mov_opt_values_shape <- unlist(support_OU_mov_opt_shape)
df_OU_mov_opt_shape <- data.frame(values = support_OU_mov_opt_values_shape)

hist_OU_mov_opt_shape = ggplot(df_OU_mov_opt_shape, aes(x = values)) +  
  geom_histogram(breaks = seq(0, 1, by = 0.1), fill = "dodgerblue3", color = "black", boundary = 0, closed = "right") +
  labs(title = "Histogram of Support for OU_mov_opt for shape traits", 
       x = "Akaike weight", 
       y = "Frequency") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1), expand = expansion(mult = c(0, 0))) +  
  scale_y_continuous(limits = c(0, 150), expand = expansion(mult = c(0, 0))) +  
  theme_minimal()

# Set common theme for the plots
custom_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", margin = margin(b = 10)),  # Title size and space
  axis.title.x = element_text(size = 14),  # X-axis title size
  axis.title.y = element_text(size = 14),  # Y-axis title size
  axis.text.x = element_text(size = 12),   # X-axis text size
  axis.text.y = element_text(size = 12),   # Y-axis text size
  panel.spacing = unit(2, "lines")  # Space between plots
)

# Apply the custom theme to each histogram
hist_Stasis_size <- hist_Stasis_size + custom_theme
hist_URW_size <- hist_URW_size + custom_theme
hist_GRW_size <- hist_GRW_size + custom_theme
hist_Strict_stasis_size <- hist_Strict_stasis_size + custom_theme
hist_Accel_size <- hist_Accel_size + custom_theme
hist_Decel_size <- hist_Decel_size + custom_theme
hist_OU_size <- hist_OU_size + custom_theme
hist_OU_mov_opt_anc_size <- hist_OU_mov_opt_anc_size + custom_theme
hist_OU_mov_opt_size <- hist_OU_mov_opt_size + custom_theme

# Apply the custom theme to each histogram
hist_Stasis_shape <- hist_Stasis_shape + custom_theme
hist_URW_shape <- hist_URW_shape + custom_theme
hist_GRW_shape <- hist_GRW_shape + custom_theme
hist_Strict_stasis_shape <- hist_Strict_stasis_shape + custom_theme
hist_Accel_shape <- hist_Accel_shape + custom_theme
hist_Decel_shape <- hist_Decel_shape + custom_theme
hist_OU_shape <- hist_OU_shape + custom_theme
hist_OU_mov_opt_anc_shape <- hist_OU_mov_opt_anc_shape + custom_theme
hist_OU_mov_opt_shape <- hist_OU_mov_opt_shape + custom_theme


# save the graphs
png("./results_paleoTS_v0.6.1/plot/results_support_traitcategory_barplot.png", width = 2400, height = 1800)
grid.arrange(hist_Stasis_size, hist_URW_size, hist_GRW_size, hist_Stasis_shape, hist_URW_shape, hist_GRW_shape, hist_Strict_stasis_size, hist_Accel_size, hist_Decel_size,
               hist_Strict_stasis_shape, hist_Accel_shape, hist_Decel_shape,  hist_OU_size, hist_OU_mov_opt_anc_size, hist_OU_mov_opt_size, hist_OU_shape, hist_OU_mov_opt_anc_shape, hist_OU_mov_opt_shape, ncol = 3)

dev.off()

