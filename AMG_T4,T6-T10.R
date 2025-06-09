library(dplyr)
library(ggplot2)
library(tidyr)  
library(reshape2)
library(stringr)

setwd("C:/Users/Дарина/Desktop/Everest/virom")

calculate_p_value <- function(r, n) {
  if (is.na(r) || n < 3 || abs(r) == 1) {
    return(NA)  
  }
  t_value <- (r * sqrt(n - 2)) / sqrt(1 - r^2)
  p_value <- 2 * pt(-abs(t_value), df = n - 2)  
  return(p_value)
}

 data <- read.csv("virus_C_all_removed.csv")

colnames(data) <- make.names(colnames(data))

result <- data %>%
  group_by(time.point.9, ProteinID) %>%
  summarise(
    Mean_relab = mean(Relative.abundance, na.rm = TRUE), 
    n_values = n(), 
    .groups = "drop"
  ) %>% 
  ungroup() %>% as.data.frame()

periods <- c("T4", "T6", "T7", "T8", "T9", "T10")

filtered_result <- result %>%
  filter(time.point.9 %in% periods)

filtered_result <- filtered_result %>%
  group_by(ProteinID) %>%
  filter(any(Mean_relab != 0, na.rm = TRUE)) %>%
  ungroup()

protein_time_combinations <- expand.grid(
  ProteinID=unique(filtered_result$ProteinID),
  time.point.9=periods
)

filtered_result <- protein_time_combinations %>% 
  left_join(., filtered_result, by=c("ProteinID", "time.point.9"))

dcast_filtered_result <- reshape2::dcast(filtered_result %>% select(ProteinID,time.point.9,Mean_relab),
                                         ProteinID~time.point.9, value.var="Mean_relab") %>% 
  select(c("ProteinID", all_of(periods)))

dcast_filtered_result <- dcast_filtered_result %>% 
  rowwise() %>% 
  mutate(
    Mean = mean(c(T4, T6, T7, T8, T9, T10), na.rm = TRUE),
    StandardDeviation = sd(c(T4, T6, T7, T8, T9, T10), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(
    Z_T4 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T4 - Mean) / StandardDeviation, 0),
    Z_T6 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T6 - Mean) / StandardDeviation, 0),
    Z_T7 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T7 - Mean) / StandardDeviation, 0),
    Z_T8 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T8 - Mean) / StandardDeviation, 0),
    Z_T9 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T9 - Mean) / StandardDeviation, 0),
    Z_T10 = ifelse(!is.na(StandardDeviation) & StandardDeviation != 0, (T10 - Mean) / StandardDeviation, 0) 
  ) %>% unique()

Zscore_filtered_result <- dcast_filtered_result[,c(1,grep("Z_", colnames(dcast_filtered_result)))]
save(Zscore_filtered_result, file="ZAMG_T4,T6-T10.RData")

correlation_results <- reshape2::melt(Zscore_filtered_result, id = "ProteinID") %>%
  mutate(time_index=case_when(variable=="Z_T4" ~ 1,
                              variable=="Z_T6" ~ 2,
                              variable=="Z_T7" ~ 3,
                              variable=="Z_T8" ~ 4,
                              variable=="Z_T9" ~ 5,
                              variable=="Z_T10" ~ 6,
                              TRUE ~ NA)) %>% 
  group_by(ProteinID) %>%
  summarise(
    Zscore_Mean_relab_values = list(value),
    available_periods = sum(!is.na(value)),
    correlation = ifelse(available_periods >= 2 & sd(unlist(value), na.rm = TRUE) > 0, 
                         cor(unlist(value), unlist(time_index), use = "complete.obs"),
                         NA),
    p_value = ifelse(!is.na(correlation), calculate_p_value(correlation, available_periods), NA),
    trend = case_when(
      correlation > 0 ~ "Increase",
      correlation < 0 ~ "Decrease",
      TRUE ~ "No Change"
    ),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  as.data.frame()

correlation_results <- correlation_results %>%
  mutate(
    Zscore_Mean_relab_values = sapply(Zscore_Mean_relab_values, function(x) paste(x, collapse = "; "))
  )%>%
  select(-available_periods) %>% 
  unique()


annotation_data <- data %>%
  select(ProteinID, AMG_KO_name) %>%
  filter(!is.na(AMG_KO_name) & AMG_KO_name != "" & AMG_KO_name != "-") %>%
  mutate(annotation = AMG_KO_name) %>%  
  select(-AMG_KO_name)

correlation_results <- correlation_results %>%
  left_join(annotation_data %>% select(ProteinID, annotation), by = "ProteinID") %>%
  filter(!is.na(annotation) & annotation != "" & annotation != "-") %>%
  unique()

significant_correlation_results <- correlation_results %>%
  filter(p_value < 0.05 & (correlation < -0.7 | correlation > 0.7))

significant_correlation_results <- significant_correlation_results %>%
  mutate(
    Zscore_Mean_relab_values_numeric = lapply(Zscore_Mean_relab_values, function(x) {
      vec <- unlist(strsplit(as.character(x), "; "))
      vec[vec == "NA"] <- NA
      as.numeric(vec)
    }),
    
    Zscore_Mean_relab_average = sapply(Zscore_Mean_relab_values_numeric, function(x) mean(x, na.rm = TRUE))
  )%>%
  select(-Zscore_Mean_relab_values_numeric) %>% 
  as.data.frame()

ave_Mean_relab_df <- filtered_result %>% 
  group_by(ProteinID) %>% 
  summarise(ave_Mean_relab=mean(Mean_relab, na.rm = TRUE)) %>% 
  ungroup() %>% 
  as.data.frame()

significant_correlation_results <- merge(significant_correlation_results, ave_Mean_relab_df, by="ProteinID")

top_100_significant_correlation_results <- significant_correlation_results %>%
  arrange(desc(ave_Mean_relab)) %>%
  head(100)

trim_annotation <- function(x) {
  words <- unlist(str_split(x, "\\s+"))  
  if (length(words) > 10) {
    return(paste(words[1:10], collapse = " ")) 
  } else {
    return(x) 
  }
}

top_100_significant_correlation_results <- top_100_significant_correlation_results %>%
  mutate(
    annotation = sapply(annotation, trim_annotation),  
    annotation = str_wrap(annotation, width = 40),  
    group = ifelse(row_number() <= 50, "First 50", "Second 50"),  
    significant = ifelse(p_value < 0.01, "*", "")
  )

top_100_extended <- merge(Zscore_filtered_result, top_100_significant_correlation_results, by = "ProteinID")

top_100_long <- top_100_extended %>%
  tidyr::pivot_longer(cols = starts_with("Z_T"),
                      names_to = "Timepoint", values_to = "Z_score") %>%
  mutate(Timepoint = factor(Timepoint, levels = paste0("Z_T", c(4, 6, 7, 8, 9, 10)),
                            labels = paste0("T", c(4, 6, 7, 8, 9, 10))))

min_z <- min(top_100_long$Z_score, na.rm = TRUE)
max_z <- max(top_100_long$Z_score, na.rm = TRUE)



# pdf("B_AMG_T4,T6-T10.pdf",width = 17,height = 14)
# ggplot() +
#   geom_text(data = top_100_extended, aes(x = "Significant", y = ProteinID, label = significant), size = 4, color = "black") +
#   geom_point(data = top_100_extended, aes(x = "Correlation", y = ProteinID, size = abs(correlation), color = correlation)) +
#   geom_tile(data = top_100_long, aes(x = Timepoint, y = ProteinID, fill = Z_score), color = "black") +
#   scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", midpoint = 0, limits = c(min_z, max_z), name = "Mean relab(Z-score)") +
#   scale_color_gradientn(limits=c(-1, 1), colors = c("blue", "lightblue", "white", "pink", "red"), breaks = c(-1, -0.5, 0, 0.5, 1), labels = format(c(-1, -0.5, 0, 0.5, 1))) +
#   scale_size_continuous(range = c(1, 3), name = "correlation") +
#   scale_y_discrete(breaks = top_100_extended$ProteinID, labels = top_100_extended$annotation, expand = c(0,0))+
#   scale_x_discrete(expand = c(0,0))+
#   facet_wrap(~group, ncol = 2, scales = "free_y", strip.position = "top") +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_text(size = 15, lineheight = 1),
#     axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     strip.text = element_text(size = 8, face = "bold"),
#     plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
#     legend.position = "bottom",
#     legend.title = element_text(size = 7),
#     legend.text = element_text(size = 7),
#     legend.key.size = unit(0.5, "cm"),
#     plot.margin = margin(1, 150, 1, 150)
#   ) +
#   labs(title = "Correlation of Protein Relative Abundance across T4,T6-T10 Timepoints (AMG)")
# dev.off()


# Remove all objects except the original dataset
# rm(list = setdiff(ls(), "data"))
