library(openxlsx)
library(dplyr)
library(pheatmap)
library(tidyr)
library(stringr)
library(ggplot2)

setwd("C:/Users/Дарина/Desktop/Everest/virom")

periods <- c("T4", "T6", "T7", "T8", "T9", "T10")

heatmap_data <- read.xlsx("HHeatmap_data_all_T4,T6-T10.xlsx")

annotation_data <- read.csv("Groups.csv")

data <- merge(annotation_data, heatmap_data, by.x = "query", by.y = "ProteinID", all.y = T) %>% 
  mutate(COG_category = ifelse(is.na(COG_category), "-", COG_category))

bar_data <- data %>% 
  group_by(COG_category) %>% 
  summarise(COG_number = n()) %>% 
  arrange(desc(COG_number)) %>% 
  ungroup() %>% 
  as.data.frame()

bar_data$COG_category <- factor(bar_data$COG_category, 
                                levels = bar_data$COG_category[order(bar_data$COG_number, decreasing = TRUE)])


pdf("all_T4,T6-T10.COG_bar.pdf", width = 14, height = 2)
ggplot(bar_data, aes(x = COG_category, y = COG_number)) +
  geom_bar(stat = "identity", fill = "darkblue", width = 0.6) +
  labs(title = "COG Category of Protein",
       x = "COG_category",
       y = "Number of proteins") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
        axis.text.y = element_text(size = 8))+
  scale_y_continuous(expand = c(0,0))
dev.off()

filter_bar_data <- bar_data %>% 
  arrange(desc(COG_number)) %>% 
  slice_head(n = 25)

data <- data %>% 
  filter(COG_category %in% filter_bar_data$COG_category)

data$COG_category <- factor(data$COG_category, 
                            levels = filter_bar_data$COG_category[order(filter_bar_data$COG_number, decreasing = TRUE)])


max_time_data <-  data %>%
  pivot_longer(cols = all_of(periods), names_to = "time.point", values_to = "Mean_relab") %>%
  group_by(query) %>%
  summarise(max_time = time.point[which.max(Mean_relab)],
            max_expression = max(Mean_relab))

heatmap <- data %>%
  left_join(max_time_data, by = "query")


heatmap$max_time <- factor(heatmap$max_time, levels = periods)

sorted_data <- heatmap %>%
  arrange(max_time, COG_category, desc(max_expression), .by_group = TRUE)

heatmap_matrix <- sorted_data %>%
  select(all_of(periods)) %>%
  as.matrix()
rownames(heatmap_matrix) <- sorted_data$query

annotation_row <- data.frame(Pathway = sorted_data$COG_category)
rownames(annotation_row) <- rownames(heatmap_matrix)


top_heatmap <- pheatmap(heatmap_matrix,
                        annotation_row = annotation_row, 
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        show_rownames = FALSE,
                        scale = "row", 
                        color = colorRampPalette(c("blue", "white", "red"))(100))

pdf("all_T4,T6-T10.COG_heatmap.top25.pdf", width = 5, height = 6)
gridExtra::grid.arrange(top_heatmap$gtable)
dev.off()


data_list <- split(data, ~ COG_category)

heatmap_list <- lapply(data_list, function(x){

  sub_max_time_data <-  x %>%
    pivot_longer(cols = all_of(periods), names_to = "time.point", values_to = "Mean_relab") %>%
    group_by(query) %>%
    summarise(max_time = time.point[which.max(Mean_relab)],
              max_expression = max(Mean_relab))

  sub_heatmap <- x %>%
    left_join(sub_max_time_data, by = "query")

  sub_heatmap$max_time <- factor(sub_heatmap$max_time, levels = periods)

  sub_sorted_data <- sub_heatmap %>%
    group_by(max_time) %>%
    arrange(max_time, desc(max_expression), .by_group = TRUE) %>%
    ungroup()

  sub_heatmap_matrix <- sub_sorted_data %>%
    select(all_of(periods)) %>%
    as.matrix()

  rownames(sub_heatmap_matrix) <- str_wrap(sub_sorted_data$annotaion, width = 100)

  COG_ID <- unique(x$COG_category)

  p <- pheatmap(sub_heatmap_matrix,
           scale = "row",
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_rownames = TRUE,
           main = paste0("Heatmap of Significant Proteins of COG: ", COG_ID),
           color = colorRampPalette(colors = c("blue","white","red"))(100)
  )

  p$gtable

})

pdf("all_T4,T6-T10.COG_heatmap.pdf", width = 10, height = 8)
for (i in 1:length(heatmap_list)) {
  gridExtra::grid.arrange(heatmap_list[[i]])
}
dev.off()

