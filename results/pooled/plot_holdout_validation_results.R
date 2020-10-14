# Plot results of hold out testing
rm(list=ls())
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis/results/pooled')
require(tidyverse); require(ggplot2); require(data.table)

# Function to parse results
parse_df <- function(dat, dataset_name) {
  return(dat %>% 
           mutate(holdout_set = dataset_name) %>%
           rename(feature_space = V1) %>%
           mutate(feature_space = recode(feature_space, "Raw" = "Before decontamination", "Raw CR" = "After decontamination")) %>%
           select(feature_space, holdout_set, 
                  external_test_precision, external_test_recall,
                  external_test_AUPRC) %>%
           gather(metric, value = performance, external_test_precision, external_test_recall, external_test_AUPRC) %>%
           filter(!grepl("RA", feature_space))
         )
}

karius_df <- parse_df(fread("hold_karius_out_model_results.csv"), "Karius")
grumaz_df <- parse_df(fread("hold_grumaz_out_model_results.csv"), "Grumaz-16/19")
kapusta_df <- parse_df(fread("hold_kapusta_out_model_results.csv"), "Gosiewski-17")

plot_df <- rbind(karius_df, grumaz_df, kapusta_df)
plot_df$metric <- factor(plot_df$metric, levels = unique(plot_df$metric))
plot_df$feature_space <- factor(plot_df$feature_space, levels = unique(plot_df$feature_space))

# Generate barplot
metric_names <- c(
  "external_test_precision" = "Precision",
  "external_test_recall" = "Recall",
  "external_test_AUPRC" = "AUPRC"
)

ggplot(plot_df, aes(x = holdout_set, y = performance, fill = feature_space)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  facet_grid(metric~., labeller = as_labeller(metric_names)) +
  geom_text(aes(label = round(performance, 2)), position = position_dodge(width=0.9), vjust=0, size = 2) +
  labs(x = "Holdout set", y = "Performance", fill = "Feature Space") +
  theme(axis.title.x = element_text(vjust = -1))

ggsave("holdout_results_barplot.png", dpi=600, width = 7, height = 5)
