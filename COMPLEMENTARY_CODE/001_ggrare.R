
ggrare_function <- function(physeq,name){
message("Saving rarefaction curves plot (using ggrare function)")
ggrare_plot <-
  ranacapa::ggrare(
    physeq,
    step = 100,
    se = TRUE,
    parallel = T,
    color = "file",
    label = "orig_id"
  )
rare_data <- ggrare_plot$data
# 
# Paso 2: median and median absolute deviation
summary_stats <- rare_data %>%
  group_by(file, Size) %>%
  summarise(
    median_richness = median(.S, na.rm = TRUE),
    mad_richness = mad(.S, na.rm = TRUE),
    .groups = "drop"
  )

# Paso 3: smooth median and mad
smooth_data <- summary_stats %>%
  group_by(file) %>%
  arrange(Size) %>%
  mutate(
    smooth_median = predict(loess(median_richness ~ Size, span = 0.4)),
    smooth_mad = predict(loess(mad_richness ~ Size, span = 0.4))
  ) %>%
  ungroup()

#Plot
ggrare_plot <- ggplot() +
  # Líneas individuales (una por muestra)
  geom_line(data = rare_data,
            aes(x = Size, y = .S, group = file, color = file),
            alpha = 0.3, size = 0.9) +
  
  # Banda: mediana ± MAD suavizada
  geom_ribbon(data = smooth_data,
              aes(x = Size,
                  ymin = smooth_median - smooth_mad,
                  ymax = smooth_median + smooth_mad,
                  fill = file),
              alpha = 0.2, color = NA) +
  
  # Línea de mediana suavizada
  geom_line(data = smooth_data,
            aes(x = Size, y = smooth_median, color = file),
            size = 3,) +
  ggsci::scale_color_jco() +
  ggsci::scale_fill_jco() +
  
  labs(x = "Number of Sequences", y = "Richness") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  )

message("Plotting rarefaction curves")
ggsave(
  paste0(graph_dir, "/", name,"_001_rarefaction_curve_ggrare.pdf"),
  ggrare_plot,
  scale = 0.9,
  width = 25,
  height = 14,
  units = "in",
  dpi = 600
)
}
