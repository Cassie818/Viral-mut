library(gprofiler2)
library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)


process_gost_data <- function(genes, title_name) {
  gostres <- gost(query = genes, 
                  organism = "hsapiens", 
                  significant = TRUE, 
                  exclude_iea = TRUE, 
                  user_threshold = 0.05, 
                  correction_method = "g_SCS", 
                  sources = c("GO:MF"))
  
  if (is.null(gostres)) return(NULL)
  
  df <- gostres$result %>%
    filter(!term_name %in% c("calcium ion transmembrane transporter activity", 
                             "glutamate-gated calcium ion channel activity")) %>%
    arrange(p_value) %>%
    head(10) %>%
    mutate(term_name = factor(term_name, levels = rev(unique(term_name[order(intersection_size, -p_value)]))))
  
  return(df)
}

high_pli_df <- process_gost_data(c('CACNA1D', 'CAMTA1', 'CYBB', 'DICER1', 'DYNC1H1', 
                                   'EXT1', 'GRIN2A', 'GRIN2B', 'IL2RG', 'INF2', 
                                   'KDM3B', 'KDM6A', 'MORC2', 'MYRF', 'NR0B1', 
                                   'OFD1', 'PKD1', 'POGZ', 'PURA', 'RB1', 
                                   'SERPING1', 'SLC16A2', 'SMC1A', 'TAF1', 'TNFRSF1A', 
                                   'WAS'), "High")

low_pli_df <- process_gost_data(c('ACAD9', 'CC2D2A', 'DEAF1', 'DOK7', 'DYNC2H1', 
                                  'EYS', 'GNPTAB', 'HBA1', 'HEXB', 'HLCS', 
                                  'IL12RB1', 'ITGB4', 'KCNJ11', 'KCNT1', 'MAN2B1', 
                                  'MSH6', 'MYO15A', 'MYOC', 'NBEAL2', 'NEU1', 
                                  'NTRK1', 'OCA2', 'OTOF', 'PC', 'PIK3R2', 
                                  'RAG1', 'SLC4A11', 'SUMF1', 'TECTA', 'TGM1', 
                                  'TYMP', 'USH2A', 'WDR62', 'WFS1'), "Low")


build_plot <- function(data, title, color_low = "#542788") {
  ggplot(data, aes(x = intersection_size, y = term_name, fill = p_value)) +
    geom_col(width = 0.75) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) + 
    scale_fill_gradient(low = color_low, high = "#F0F0F0", 
                        name = "p-value",
                        labels = function(x) sprintf("%.3f", x),
                        guide = guide_colorbar(barwidth = 0.8, barheight = 4)) +
    labs(title = title, x = "Count", y = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.y = element_text(size = 11, color = "black", lineheight = 0.8),
      axis.text.x = element_text(size = 11),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey95"),
      plot.margin = margin(5, 5, 5, 5)
    )
}


p_high <- build_plot(high_pli_df, "PLM better (high pLI)", color_low = "#4D004B")
p_low <- build_plot(low_pli_df, "PLM better (low pLI)", color_low = "#084594")


final_plot <- (p_high / p_low) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

print(final_plot)