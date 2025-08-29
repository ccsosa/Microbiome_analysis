library(dplyr)
library(ggplot2)
library(ggrepel)
library(forcats)

# Suponiendo que df$P ya es BH-ajustada
df <- x2

# Suponemos que df ya tiene taxon, level, P (ajustado BH)
df <- df %>%
  filter(!is.na(P)) %>%
  mutate(
    neglog10P = -log10(P),
    grouping = ifelse(is.na(level) | level == "", "Unknown", level)
  )

# Mantener orden original de aparición de grouping
df$grouping <- factor(df$grouping, levels = unique(df$grouping))

# Posiciones acumuladas
df <- df %>%
  group_by(grouping) %>%
  mutate(pos_in_group = row_number()) %>%
  ungroup()

group_sizes <- df %>% count(grouping, name = "n")
group_offsets <- group_sizes %>%
  mutate(offset = lag(cumsum(n), default = 0)) %>%
  select(grouping, offset)

df <- df %>%
  left_join(group_offsets, by = "grouping") %>%
  mutate(pos = pos_in_group + offset)

# Centros para etiquetas en eje X
group_centers <- df %>%
  group_by(grouping) %>%
  summarise(center = mean(range(pos)), .groups = "drop")

# Taxa significativos
df_sig <- df %>% filter(P < 0.05)

# Selección top-N por grupo (ejemplo N=10)
N <- 10
df_top_per_group <- df_sig %>%
  group_by(grouping) %>%
  slice_min(order_by = P, n = N, with_ties = FALSE) %>%
  ungroup()

# Plot
ggplot(df, aes(x = pos, y = neglog10P, color = grouping)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_point(data = df_sig, aes(x = pos, y = neglog10P),
             shape = 21, color = "black", size = 2, inherit.aes = FALSE) +
  scale_x_continuous(breaks = group_centers$center,
                     labels = levels(df$grouping)) +
  labs(x = NULL, y = expression(-log[10](italic(P))),
       title = "Manhattan plot of Wilcoxon test results (FDR-adjusted)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggrepel::geom_text_repel(
    data = df_top_per_group,
    aes(label = taxon),
    size = 3, box.padding = 0.3, max.overlaps = 100000
  )
