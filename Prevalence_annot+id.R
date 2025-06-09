library(dplyr)
library(tidyr)
library(openxlsx)
setwd("C:/Users/Дарина/Desktop/Everest/virom")
data <- read.csv("virus_B_all_removed.csv")


df_filtered <- data %>%
  filter(!is.na(AMG_KO_name),
         trimws(AMG_KO_name) != "",
         trimws(AMG_KO_name) != "-",
         !(time.point.9 %in% c("Ta", "Tb")))

unique_proteins_full <- df_filtered %>%
  group_by(ProteinID) %>%
  slice(1) %>%
  ungroup()

cat("Уникальные строки с аннотированными протеинами:\n")



df_filtered <- df_filtered %>%
  mutate(Relative.abundance = ifelse(Relative.abundance != 0, 1, 0))

df_filtered <- df_filtered %>%
  mutate(person_time = paste0(time.point.9, "_", people))

pivot_df <- df_filtered %>%
  select(person_time, ProteinID, Relative.abundance) %>%
  group_by(person_time, ProteinID) %>%
  summarise(value = max(Relative.abundance), .groups = "drop") %>%
  pivot_wider(names_from = ProteinID, values_from = value, values_fill = 0)

pivot_df <- pivot_df %>%
  mutate(
    timepoint_num = as.numeric(gsub("^T", "", sub("_.*", "", person_time))),
    person_num = as.numeric(gsub("^S", "", sub(".*_", "", person_time)))
  ) %>%
  arrange(timepoint_num, person_num) %>%
  select(-timepoint_num, -person_num)


protein_data <- pivot_df[, -1]
total_rows <- nrow(protein_data)

prevalence_vector <- colSums(protein_data == 1, na.rm = TRUE) / total_rows


prevalence_row <- c("prevalent", as.list(prevalence_vector))
names(prevalence_row) <- names(pivot_df)


pivot_df_with_prevalence <- rbind(pivot_df, prevalence_row)

# write.csv(pivot_df_with_prevalence, "pivot_df_with_prevalence.csv", row.names = FALSE)
# write.xlsx(pivot_df_with_prevalence, file = "pivot_df_with_prevalence.xlsx")

df_no_prev <- pivot_df_with_prevalence %>%
  filter(person_time != "prevalent")

df_prevalence_by_tp <- df_no_prev %>%
  mutate(timepoint = sub("_.*", "", person_time))

prevalence_by_timepoint <- df_prevalence_by_tp %>%
  group_by(timepoint) %>%
  summarise(across(-person_time, ~ mean(as.numeric(.)), .names = "{.col}"))

prevalence_by_timepoint <- prevalence_by_timepoint %>%
  mutate(timepoint_num = as.numeric(gsub("^T", "", timepoint))) %>%
  arrange(timepoint_num) %>%
  select(-timepoint_num)

protein_cols <- colnames(prevalence_by_timepoint)[-1]
nonzero_proteins <- sapply(prevalence_by_timepoint[protein_cols], function(x) any(as.numeric(x) > 0))

filtered_prevalence <- prevalence_by_timepoint[, c("timepoint", names(nonzero_proteins[nonzero_proteins]))]

num_removed <- sum(!nonzero_proteins)
num_remaining <- sum(nonzero_proteins)
cat("Отфильтровано протеинов без представленности:", num_removed, "\n")
cat("Осталось протеинов:", num_remaining, "\n")

library(ggplot2)
library(tidyr)
library(dplyr)

plot_data <- filtered_prevalence %>%
  pivot_longer(cols = -timepoint, names_to = "ProteinID", values_to = "Prevalence") %>%
  mutate(Prevalence = as.numeric(Prevalence))

protein_annotation_map <- df_filtered %>%
  select(ProteinID, AMG_KO_name) %>%
  filter(trimws(AMG_KO_name) != "", AMG_KO_name != "-", !is.na(AMG_KO_name)) %>%
  distinct() %>%
  group_by(ProteinID) %>%
  summarise(AMG_KO_name = paste(unique(AMG_KO_name), collapse = "; "), .groups = "drop")

plot_data <- plot_data %>%
  left_join(protein_annotation_map, by = "ProteinID") %>%
  mutate(ProteinLabel = ifelse(
    is.na(AMG_KO_name),
    ProteinID,
    paste0(AMG_KO_name, " (", ProteinID, ")")
  )) %>%
  mutate(
    ProteinLabel = ifelse(nchar(ProteinLabel) > 140,
                          paste0(substr(ProteinLabel, 1, 139), "..."),
                          ProteinLabel)
  )


plot_data$timepoint <- factor(plot_data$timepoint, levels = paste0("T", sort(as.numeric(gsub("T", "", unique(plot_data$timepoint))))))

plot_data$ProteinLabel <- factor(
  plot_data$ProteinLabel,
  levels = plot_data %>%
    group_by(ProteinLabel) %>%
    summarise(mean_prev = mean(Prevalence, na.rm = TRUE)) %>%
    arrange(desc(mean_prev)) %>%
    pull(ProteinLabel)
)

save_heatmaps <- TRUE       
plot_width <- 15            
plot_height <- 10           
plot_dpi <- 300            
plot_units <- "in"          

protein_labels <- levels(plot_data$ProteinLabel)
chunks <- split(protein_labels, cut(seq_along(protein_labels), breaks = 5, labels = FALSE))

for (i in seq_along(chunks)) {
  subset_labels <- chunks[[i]]
  subset_data <- plot_data %>% filter(ProteinLabel %in% subset_labels)
  
  p <- ggplot(subset_data, aes(x = timepoint, y = ProteinLabel, fill = Prevalence)) +
    geom_tile(color = "gray85") +
    scale_fill_gradient(
      low = "white",
      high = "#4b0082",
      limits = c(0, 1),
      name = "Prevalence"
    ) +
    labs(
      title = paste("Protein Prevalence Heatmap Across Timepoints (Annotation+ProteinID) ", i, "/5"),
      x = "Timepoint",
      y = "Functional Annotation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank()
    )
  
  print(p)
  
  if (save_heatmaps) {
    ggsave(
      filename = paste0("B_Prevalence_Annot_ID_", i, ".pdf"),
      plot = p,
      width = plot_width,
      height = plot_height,
      dpi = plot_dpi,
      units = plot_units
    )
    cat("Saved: Prevalence_heatmap_", i, ".pdf\n")
  }
}

