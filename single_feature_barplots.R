# Plot AUROC scores of single feature models compared to full model
rm(list=ls())
setwd('~/git_repos/Polymicrobial-Signature-of-Sepsis/results/drop_confirmed_features_full/')

df <- read.csv('single_feature_auroc_no_confirmed.csv', stringsAsFactors = F)
colnames(df) <- c("Model", "AUROC")
df$type <- "Single feature"

model <- read.csv("karius_drop_confirmed_features_results.csv", stringsAsFactors = F)
model <- model[c(3, 7), c("X", "test_AUROC")]
colnames(model) <- colnames(df)[1:2]

# Add Karius-Neat
to_append <- data.frame(Model = c("Karius-Neat"), AUROC = c(0.995))
model <- rbind(model, to_append)
model$type <- "Multi-feature"

plot_df <- rbind(model, df)
plot_df$Model[plot_df$Model == "Without Confirmed CR"] <- "Karius-Without"
plot_df$Model[plot_df$Model == "Only Confirmed"] <- "Karius-Only"
plot_df$Model
plot_df$Model <- factor(plot_df$Model, levels=plot_df$Model[order(plot_df$AUROC, decreasing = T)])
plot_df$AUROC <- plot_df$AUROC - 0.5

require(ggplot2)
ggplot(plot_df, aes(x = Model, y = AUROC, fill = type)) + 
  geom_bar(stat = "identity") +
  labs(x = "Feature Space", y = "AUROC") +
  ylim(c(0, 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), 
                     labels = as.character(seq(0.5, 1.0, 0.1)), 
                     limits = c(0, 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        plot.margin = unit(c(0,0,0,1), "cm")) 

ggsave("AUROC_barplot.png", dpi=600, width = 7, height = 5)
