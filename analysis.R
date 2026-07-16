library(tidyverse)
library(arrow)
library(modelr)
library(broom)
library(ggpubr)
library(ggcorrplot)
library(plotly)
library(ggfortify)
library(Rdimtools)

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
  mutate(
    sample_type = case_when(
      str_detect(sample_type, "oral") ~ "oral swab",
      str_detect(sample_type, "nasal") ~ "nasal swab",
      TRUE ~ sample_type
    ),
    across(everything(), ~ replace_na(., 0))
  )


# PCA --------------------------------------------------------------------


params <- c("S1", "a1", "b1", "S2", "a2", "b2")
cor_df <- df_ %>%
  select(all_of(params)) %>%
  cor(use = "pairwise.complete.obs")

ggcorrplot(
  cor_df,
  hc.order = F, type = "upper", outline.col = "white", lab = TRUE,
  lab_size = 4, title = "Correlation Matrix", ggtheme = theme_minimal()
)


pca <- df_ %>%
  select(all_of(params)) %>%
  do.spc(df_$assay, ndim = 1)

# eigenvectors <- pca$rotation %>%
#   as.data.frame() %>%
#   mutate(
#     across(everything(), ~ .x * 10),
#     variable = rownames(.)
#   )

# summary(pca)

# kms <- kmeans(df_[params], centers = 2)

df_pca <- df_ %>%
  mutate(
    Y1 = pca$Y[, 1],
    # Y2 = pca$Y[, 2],
  )
# bind_cols(pca$Y, names = "")
# mutate(cluster = kms$cluster)

df_pca %>%
  arrange(mortem) %>%
  ggplot(aes(x = mortem, y = Y1, color = assay)) +
  geom_boxplot() +
  # geom_point(size = 2, alpha=0.1) +
  # geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = eigenvectors, linewidth = 1, color = "blue") +
  # geom_text(aes(label = variable), data = eigenvectors, size = 6, color = "blue") +
  # stat_ellipse(aes(group=cluster), color="black", level=0.95) +
  # facet_grid(cols = vars(mortem), rows = vars(assay)) +
  # scale_color_gradient(low="blue", high="red") +
  # coord_fixed() +
  guides(
    color = guide_legend(override.aes = list(size = 6, alpha = 1))
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
ggsave("figures/pca.png", width = 16, height = 6)
