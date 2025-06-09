library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(stringr)

setwd("C:/Users/Дарина/Desktop/Everest/virom")

data <- read.csv("virus_C_all_removed.csv")

periods <- c("T4", "T6", "T7", "T8", "T9", "T10")

filtered_data <- data %>%
  filter(time.point.9 %in% periods) %>% 
  group_by(time.point.9, ProteinID) %>%
  summarise(
    Mean_relab = mean(Relative.abundance, na.rm = TRUE), 
    n_values = n(), 
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  group_by(ProteinID) %>%
  filter(any(Mean_relab != 0, na.rm = TRUE)) %>%  
  ungroup() %>% 
  as.data.frame()

filtered_data$time.point.9 <- factor(filtered_data$time.point.9, levels = periods)

data_filtered <- data %>%
  inner_join(filtered_data, by = c("ProteinID", "time.point.9")) %>% 
  arrange(ProteinID) 

trim_annotation <- function(x) {
  words <- unlist(str_split(x, "\\s+")) 
  if (length(words) > 10) {
    return(paste(words[1:10], collapse = " "))  
  } else {
    return(x) 
  }
}

data_filtered <- data_filtered %>%
  filter(!is.na(AMG_KO_name) & AMG_KO_name != "" & AMG_KO_name != "-") %>%
  mutate(annotation = str_wrap(sapply(AMG_KO_name, trim_annotation), width = 100))

annotation <- data_filtered %>% 
  select(ProteinID, annotation) %>% 
  unique()
rownames(annotation) <- annotation$ProteinID

wilcox_results <- data_filtered %>%
  group_by(ProteinID) %>%
  filter(n_distinct(time.point.9) > 1) %>%
  summarise(p_values = list({
    time_points <- unique(time.point.9)
    
  
    if (length(time_points) < 2) {
      return(tibble(t1 = character(), t2 = character(), p_value = numeric()))
    }
    
    result_list <- combn(time_points, 2, function(pair) {
      t1 <- pair[1]
      t2 <- pair[2]
      
      data_t1 <- Relative.abundance[time.point.9 == t1]
      data_t2 <- Relative.abundance[time.point.9 == t2]
      
      if (length(data_t1) == 0 | length(data_t2) == 0) {
        return(tibble(t1 = t1, t2 = t2, p_value = NA))
      }
      
      p_value <- if (length(data_t1) > 1 & length(data_t2) > 1) {
        wilcox.test(data_t1, data_t2, exact = FALSE)$p.value
      } else {
        NA
      }
      
      tibble(t1 = t1, t2 = t2, p_value = p_value)
    }, simplify = FALSE)
    
    bind_rows(result_list)
  })) %>% 
  unnest(p_values)

wilcox_results <- wilcox_results %>% 
  mutate(signif = case_when(
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

significant_proteins <- wilcox_results %>%
  filter(p_value < 0.05) %>%
  group_by(t1, t2) %>% 
  pull(ProteinID) %>% 
  unique()

heatmap_data <- filtered_data %>%
  filter(ProteinID %in% significant_proteins) %>%
  dcast(ProteinID ~ time.point.9, value.var = "Mean_relab")

max_time_data <- heatmap_data %>%
  pivot_longer(cols = starts_with("T"), names_to = "time.point", values_to = "Mean_relab") %>%
  group_by(ProteinID) %>%
  summarise(max_time = time.point[which.max(Mean_relab)],
            max_expression = max(Mean_relab))

heatmap_data <- heatmap_data %>%
  left_join(max_time_data, by = "ProteinID")

heatmap_data$max_time <- factor(heatmap_data$max_time, levels = periods)
sorted_data <- heatmap_data %>%
  group_by(max_time) %>%
  arrange(max_time, desc(max_expression), .by_group = TRUE) %>%
  ungroup()

heatmap_matrix <- as.matrix(sorted_data[, periods])
rownames(heatmap_matrix) <- annotation[sorted_data$ProteinID,]$annotation
heatmap_matrix[is.na(heatmap_matrix)] <- 0  

p <- pheatmap(heatmap_matrix,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_rownames = TRUE,
         main = "Heatmap of Significant Proteins",
         color = colorRampPalette(colors = c("blue","white","red"))(100)
)

pdf("P_valueAMG.T4_T6-T10.pdf",width = 15,height = 25)
gridExtra::grid.arrange(p$gtable)
dev.off()

openxlsx::write.xlsx(wilcox_results, file = "AMG_T4_T6-T10_pvalue.xlsx")

# Remove all objects except the original dataset
#rm(list = setdiff(ls(), "data"))
