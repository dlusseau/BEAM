# Plotting trials for the CD case (area, metier, year, quarter)
tot_reliability3 <- fread("results/tot_bpue_reliability_CD.csv")

tot_reliability3 <- tot_reliability3 %>%
  filter(overall_reliability == TRUE)

library(dplyr)
library(ggplot2)
library(forcats)
library(scales)

# Option A - Heatmap with tot_mean values ----

plot_data <- tot_reliability3 %>%
  filter(!is.na(tot_mean), !is.na(tot_lwr), !is.na(tot_upr)) %>% 
  mutate(year_quarter = sprintf("%d-Q%d", year, quarter))

cd_heatmap <- ggplot(plot_data, aes(x = year_quarter, y = forcats::fct_rev(factor(metierl4)),
    fill = tot_mean)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_wrap(vars(areacode), scales = "free_y") +
  scale_fill_viridis_c(option = "C", trans = scales::log10_trans(),
    name = "tot_mean") +
  labs(x = "Year and quarter", y = "Métier") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"))

# Save
ggsave(cd_heatmap, file = "results/graphs_CD/tot_CD_heatmap.png", height = 8, width = 14)

# Option A.2 with no transformed scale in the predictions
# ggplot(plot_data, aes(x = year_quarter, y = forcats::fct_rev(factor(metierl4)),
#                       fill = tot_mean)) +
#   geom_tile(colour = "white", linewidth = 0.25) +
#   facet_wrap(vars(areacode), scales = "free_y") +
#   scale_fill_viridis_c(option = "C", name = "tot_mean") +
#   labs(x = "Year and quarter", y = "Métier") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         panel.grid = element_blank(),
#         strip.background = element_rect(fill = "grey95"))

# Option B - Panel per year, area in the secondary y-axis ----
tot_reliability3_heatmap_year <- ggplot(plot_data, aes(x = factor(quarter, levels = 1:4, labels = paste0("Q", 1:4)),
                      y = forcats::fct_rev(factor(metierl4)), fill = tot_mean)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(
    rows = vars(areacode),
    cols = vars(year),
    scales = "free_y",
    space = "free_y") +
  scale_fill_viridis_c(option = "C", trans = scales::log10_trans()) +
  labs(x = "Quarter", y = "Métier level 4") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "right")

ggsave(tot_reliability3_heatmap_year, file = "results/graphs_CD/tot_CD_heatmap_year.png",
       width = 16, height = 8)


# Incertidumbre ----
plot_data_cis <- plot_data %>%
  mutate(interval_width = tot_upr - tot_lwr,
         relative_interval_width = case_when(tot_mean > 0 ~ interval_width/tot_mean,
      TRUE ~ NA_real_))


uncertainty <- ggplot(filter(plot_data_cis), aes(x = factor(quarter, levels = 1:4,
                                             labels = paste0("Q", 1:4)),
                                  y = forcats::fct_rev(factor(metierl4)),
                                  fill = relative_interval_width)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(rows = vars(areacode), cols = vars(year), scales = "free_y",
             space = "free_y") +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "Quarter", y = "Métier") +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95"))

ggsave(uncertainty, filename = "results/graphs_CD/tot_CD_heatmap_year_uncertainty.png",
       width = 16, height = 8)

# Option C - Bars format per year and area ----
data_plot <- tot_reliability3 %>%
  filter(!is.na(tot_mean), !is.na(tot_lwr), !is.na(tot_upr)) %>%
  dplyr::mutate(
    quarter = factor(quarter, levels = 1:4, labels = paste0("Q", 1:4)),
    areacode = as.character(areacode),
    metierl4 = toupper(as.character(metierl4))
  ) %>%
  droplevels()


tot_reliability3_plot_bars <- ggplot(data_plot, aes(x = forcats::fct_rev(factor((metierl4))),
    y = tot_mean, fill = quarter, group = quarter)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr),
    position = position_dodge2(width = 0.75, preserve = "single", padding = 0.10),
    width = 0.65, colour = "black", linewidth = 0.35) +
  facet_grid(rows = vars(areacode), cols = vars(year), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Q1" = "#A6CEE3", "Q2" = "#1F78B4", "Q3" = "#B2DF8A", "Q4" = "#33A02C"),
    drop = FALSE, name = "Quarter") +
  scale_y_log10(labels = scales::label_number(accuracy = 1, big.mark = " ")) +
  coord_flip() +
  labs(x = "Métier level 4", y = "Predicted total bycatch (individuals)") +
  theme_bw() +
  theme(legend.position = "bottom",
                 legend.title = element_text(face = "plain"),
    plot.caption = element_text(face = "italic", hjust = 1),
    strip.background = element_rect(fill = "grey95", colour = "grey40",
                                             linewidth = 0.4),
    strip.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 8))

ggsave(tot_reliability3_plot_bars, filename = "results/graphs_CD/tot_CD_barsquarter.png",
       width = 16, height = 8)

# Option D - Plot by area ----
data_plot <- tot_reliability3 %>%
  filter(!is.na(tot_mean), !is.na(tot_lwr), !is.na(tot_upr)) %>% 
  mutate(year_quarter = sprintf("%d-Q%d", year, quarter))

for(i in unique(data_plot$areacode)){
  plot_data_this <- data_plot %>% filter(areacode %in% i)
  plot_this <- ggplot(plot_data_this, aes(x = year_quarter, 
                                          y = forcats::fct_rev(factor(metierl4)),
                        fill = tot_mean)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    #facet_wrap(vars(areacode), scales = "free_y") +
    scale_fill_viridis_c(option = "C", trans = scales::log10_trans(),
                         name = "tot_mean") +
    labs(x = "Year and quarter", y = "Métier",
         title = paste0(i)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95"))
  #print(plot_this)
ggsave(plot_this, filename = paste0("results/graphs_CD/tot_CD_heatmap_perareacode", i, ".png"))
}


# Option E - Per metier and area, one bar per year (x-axis) & quarter in the secondary y-axis ----
for(j in unique(data_plot$metierl4)){
  data_plot_this <- data_plot %>% filter(metierl4 %in% j)
  tot_reliability3_plot_bars_year <- ggplot(data_plot_this, aes(x = year, y = tot_mean, fill = year, group = year)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr),
                position = position_dodge2(width = 0.75, preserve = "single", padding = 0.10),
                width = 0.65, colour = "black", linewidth = 0.35) +
  facet_grid(rows = vars(quarter), cols = vars(areacode), scales = "free_y", space = "free_y") +
  # scale_fill_manual(") +
  scale_y_log10(labels = scales::label_number(accuracy = 1, big.mark = " ")) +
  #labs(x = "Métier level 4", y = "Predicted total bycatch (individuals)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "plain"),
        plot.caption = element_text(face = "italic", hjust = 1),
        strip.background = element_rect(fill = "grey95", colour = "grey40",
                                        linewidth = 0.4),
        strip.text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8)) +
    ggtitle(paste0(j), subtitle = "No trends should be inferred, y-axis in log-scale")
  #print(tot_reliability3_plot_bars_year)
  ggsave(tot_reliability3_plot_bars_year, filename = paste0("results/graphs_CD/tot_CD_barspermetier_", j,
                                                  ".png"),
         width = 14, height = 10)
}

# They do not represent trends, log scale