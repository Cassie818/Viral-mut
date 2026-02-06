library(gprofiler2)
library(ggplot2)
library(stringr)
library(dplyr) 


esm_gene_list <- c('ACAD9', 'CACNA1D', 'CAMTA1', 'CC2D2A', 'CD40LG', 
                   'CFH', 'CPLANE1', 'CYBB', 'DEAF1', 'DICER1', 
                   'DOK7', 'DYNC1H1', 'DYNC2H1', 'EXT1', 'EYS', 
                   'GNPTAB', 'GRIN2A', 'GRIN2B', 'HBA1', 'HEXB',
                   'HLCS', 'IL12RB1', 'IL2RG', 'INF2', 'ITGB4', 
                   'KCNJ11', 'KCNT1', 'KDM3B', 'KDM6A', 'MAN2B1', 
                   'MECP2', 'MORC2', 'MSH6', 'MYO15A', 'MYOC', 
                   'MYRF', 'NBEAL2', 'NEU1', 'NR0B1', 'NSD2', 
                   'NTRK1', 'OCA2', 'OFD1', 'OTOF', 'PC', 
                   'PIK3R2', 'PKD1', 'POGZ', 'PURA', 'RAG1', 
                   'RB1', 'SERPING1', 'SLC16A2', 'SLC4A11', 'SMC1A', 
                   'SUMF1', 'TAF1', 'TECTA', 'TGM1', 'TNFRSF1A', 
                   'TP53', 'TYMP', 'USH2A', 'WAS', 'WDR62', 
                   'WFS1')

gostres <- gost(query = esm_gene_list, 
                organism = "hsapiens", 
                ordered_query = FALSE, 
                multi_query = FALSE, 
                significant = TRUE, 
                exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.05, 
                correction_method = "g_SCS", 
                domain_scope = "annotated", 
                custom_bg = NULL, 
                numeric_ns = "", 
                sources = c("GO:MF")) 


res_df <- gostres$result
res_df <- res_df %>%
  filter(term_name != "protein dimerization activity") %>%      
  filter(term_name != "voltage-gated channel activity") %>%
  filter(p_value < 0.05) 

res_df <- res_df[order(res_df$p_value), ]
top_terms <- head(res_df, 10) 

top_terms$term_name <- factor(top_terms$term_name, 
                              levels = top_terms$term_name[order(top_terms$intersection_size)])

ggplot(top_terms, aes(x = intersection_size, y = term_name, fill = p_value)) +
  geom_col() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_gradient(low = "purple", high = "#E0E0E0", name = "p_value")+
  labs(title = "PLM better",
       x = "Count", 
       y = NULL) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey90"), 
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 1.6,
    plot.title = element_text(
      size = 12,  
      hjust = 0.5,    
    ),
  )