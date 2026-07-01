library(tidyverse)
library(arrow)
library(modelr)
library(broom)
library(ggpubr)
library(ggcorrplot)
library(plotly)
library(ggfortify)

main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 24),
  axis.text = element_text(size = 20),
  legend.title = element_text(size = 24),
  legend.text = element_text(size = 20),
  strip.text = element_text(size = 24),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(color = "black", fill = NA, size = 1)
)

file <- "data/results.parquet"

df_ <- read_parquet(file) %>%
  mutate(across(everything(), ~ replace_na(., 0)))



# PCA --------------------------------------------------------------------


params <- c("S1", "a1", "b1", "S2", "a2", "b2")
cor_df <- df_ %>%
  select(all_of(params)) %>%
  cor(use = "pairwise.complete.obs")

ggcorrplot(
  cor_df, hc.order = F, type = "upper", outline.col = "white", lab = TRUE, 
  lab_size = 4, title = "Correlation Matrix", ggtheme = theme_minimal()
)


pca <- df_ %>%
  select(all_of(params)) %>%
  prcomp(scale. = TRUE)

eigenvectors <- pca$rotation %>%
  as.data.frame() %>%
  mutate(
    across(everything(), ~ .x * 3),
    variable = rownames(.)
  )

summary(pca)

kms <- kmeans(df_[params], centers = 2)

df_pca <- df_ %>%
  bind_cols(pca$x) %>%
  mutate(cluster = kms$cluster)

df_pca %>%
  ggplot(aes(x = PC1, y = PC2, color = assay)) +
  geom_point(size = 2) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = eigenvectors, linewidth = 1, color = "blue") +
  geom_label(aes(label = variable), data = eigenvectors, size = 6, color = "blue") +
  stat_ellipse(aes(group=cluster), color="black", level=0.95) +
  # scale_color_gradient(low="blue", high="red") +
  guides(
    color = guide_legend(override.aes = list(size = 6)) 
  ) +
  main_theme +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8, .95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.direction = "vertical",
  )


