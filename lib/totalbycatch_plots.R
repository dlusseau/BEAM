############################################################
##### Produces plot post reliability check
#####################################################

library(ggplot2)

# Taking `tot_reliability.csv`
data <- fread("results/tot_reliability.csv")

# Taking taxa from `ecoreg_species`
ecoreg_sps <- fread("data/ecoregion_metier_sps_risk_bycatch_2023.csv") # taxon column 

# Prepare the data to plot 
# Upper case scientific names
data$scientific <- paste(toupper(substr(data$species, 1, 1)), 
                         substr(data$species, 2, 
                                nchar(data$species)), sep="")
data$ecoregion <- paste(toupper(substr(data$ecoregion, 1, 1)), 
                        substr(data$ecoregion, 2, 
                               nchar(data$ecoregion)), sep="")
data$metierl4 <- toupper(data$metierl4)

# Create label
data$label <- paste(data$ecoregion, ", ",
                    data$scientific, ", ",
                    data$metierl4, sep = "")


# Total bycatch per taxa ----
# Remove no reliable tot estimates
data_total <- data %>% 
  filter(overall_reliability == TRUE)

# Order by ecoregion, species & metierl4
data_plot_tot <- data_total %>%
  arrange(desc(ecoregion), scientific, metierl4) %>%
  mutate(label = forcats::fct_inorder(label))

mammals <-
  ggplot(subset(data_plot_tot, taxon %in% "mammals"), 
         aes(x = label, 
             y = tot_mean)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr), width = 0.5, 
                fill = "lightblue",fatten=1.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  coord_flip() +
  xlab("Ecoregion, Species, Metier level 4") + 
  ylab("Total Bycatch (individuals)") +
  theme_minimal() +
  labs(caption = "Note that x-axis is on a logarithmic scale") +
  theme(plot.caption = element_text(face = "italic"))
mammals

turtles <-
  ggplot(subset(data_plot_tot, taxon %in% "turtles"), 
         aes(x = label, 
             y = tot_mean)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr), width = 0.5, 
                fill = "lightblue",fatten=1.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  coord_flip() +
  xlab("Ecoregion, Species, Metier level 4") + 
  ylab("Total Bycatch (individuals)") +
  theme_minimal() +
  labs(caption = "Note that x-axis is on a logarithmic scale") +
  theme(plot.caption = element_text(face = "italic"))
turtles

elasmobranchs <-
  ggplot(subset(data_plot_tot, taxon %in% "elasmobranchs"), 
         aes(x = label, 
             y = tot_mean)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr), width = 0.5, 
                fill = "lightblue",fatten=1.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  coord_flip() +
  xlab("Ecoregion, Species, Metier level 4") + 
  ylab("Total Bycatch (individuals)") +
  theme_minimal() +
  labs(caption = "Note that x-axis is on a logarithmic scale") +
  theme(plot.caption = element_text(face = "italic"))
elasmobranchs

seabirds <-
  ggplot(subset(data_plot_tot, taxon %in% "seabirds"), 
         aes(x = label, 
             y = tot_mean)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr), width = 0.5, 
                fill = "lightblue",fatten=1.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  coord_flip() +
  xlab("Ecoregion, Species, Metier level 4") + 
  ylab("Total Bycatch (individuals)") +
  theme_minimal() +
  labs(caption = "Note that x-axis is on a logarithmic scale") +
  theme(plot.caption = element_text(face = "italic"))
seabirds

fish <-
  ggplot(subset(data_plot_tot, taxon %in% "fish"), 
         aes(x = label, 
             y = tot_mean)) +
  geom_crossbar(aes(ymin = tot_lwr, ymax = tot_upr), width = 0.5, 
                fill = "lightblue",fatten=1.5) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  coord_flip() +
  xlab("Ecoregion, Species, Metier level 4") + 
  ylab("Total Bycatch (individuals)") +
  theme_minimal() +
  labs(caption = "Note that x-axis is on a logarithmic scale") +
  theme(plot.caption = element_text(face = "italic"))
fish

# Save the results
ggsave(filename = "results/figures/tot_bycatch_mammals_certains.png", plot = mammals,
       width = 12, height = 8, dpi = 300)
ggsave(filename = "results/figures/tot_bycatch_turtles_certains.png", plot = turtles,
       width = 12, height = 8, dpi = 300)
ggsave(filename = "results/figures/tot_bycatch_elasmobranchs_certains.png", plot = elasmobranchs,
       width = 15, height = 12, dpi = 300)
ggsave(filename = "results/figures/tot_bycatch_seabirds_certains.png", plot = seabirds,
       width = 12, height = 8, dpi = 300)
ggsave(filename = "results/figures/tot_bycatch_fishs_certains.png", plot = fish,
       width = 15, height = 12, dpi = 300)
