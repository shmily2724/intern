library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
#==================prepare data====================================================
file_path <- "/media/shmily/writable/working/result/Kinh_all_haplotype_result.xlsx"

overview <- read_excel(file_path, sheet = "Overview") %>%
  rename_with(trimws)
# wide -> long 

meta_cols <- c("Country", "Reference", "Population", "Sample size")
hap_cols <- setdiff(names(overview), meta_cols)

length(hap_cols)
head(hap_cols)

hg_long <- overview %>%
  mutate(across(all_of(hap_cols), as.numeric)) %>%
  pivot_longer(
    cols = all_of(hap_cols),
    names_to = "Haplogroup",
    values_to = "Count"
  ) %>%
  mutate(
    Haplogroup = str_trim(Haplogroup),
    present = Count > 0
  )
# giữ true haplogroup present
hg_present <- hg_long %>%
  filter(present)
#group hg 
hg_class <- hg_present %>%
  group_by(Haplogroup) %>%
  summarise(
    n_pop = n_distinct(Population),
    populations = paste(sort(unique(Population)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    class = case_when(
      n_pop == 3 ~ "shared_3pop",
      n_pop == 2 ~ "shared_2pop",
      n_pop == 1 & populations == "Central" ~ "private_Central",
      n_pop == 1 & populations == "North"   ~ "private_North",
      n_pop == 1 & populations == "South"   ~ "private_South"
    )
  )
#đếm 
hg_class %>%
  count(class)

hg_3pop <- hg_class %>% filter(class == "shared_3pop")
hg_2pop <- hg_class %>% filter(class == "shared_2pop")
hg_private <- hg_class %>% filter(grepl("^private_", class))
#===============================draw plot======================================
#create datframe 
plot_df <- hg_present %>%
  left_join(hg_class, by = "Haplogroup") %>%
  mutate(
    Frequency = Count / `Sample size` * 100,
    class_label = case_when(
      class == "shared_3pop" ~ "3 populations",
      class == "shared_2pop" ~ "2 populations",
      grepl("^private_", class) ~ "1 population"
    )
  )
#modify
hg_order <- plot_df %>%
  group_by(Haplogroup) %>%
  summarise(max_freq = max(Frequency), .groups = "drop") %>%
  arrange(desc(max_freq)) %>%
  pull(Haplogroup)

plot_df$Haplogroup <- factor

plot_df <- plot_df %>%
  mutate(
    class_label = case_when(
      class == "shared_3pop" ~ "3 populations",
      class == "shared_2pop" ~ "2 populations",
      grepl("^private_", class) ~ "1 population"
    )
  )
#draw plot 
p <- ggplot(plot_df, aes(x = Haplogroup, y = Frequency,
                         color = class_label, shape = Population)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Haplogroup),
    size = 3,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.size = 0.2
  ) +
  scale_color_manual(values = c(
    "3 populations" = "red",
    "2 populations" = "goldenrod",
    "1 population" = "springgreen4"
  )) +
  scale_shape_manual(values = c(
    "South" = 16,
    "North" = 17,
    "Central" = 4
  )) +
  labs(
    x = NULL,
    y = "Haplogroup frequency (%)",
    color = "Distribution",
    shape = "Population",
    title = "Haplogroup distribution in three Kinh populations"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
p
