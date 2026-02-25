library(gprofiler2)
library(ggplot2)
library(stringr)
library(dplyr) 

calm_gene_list <- c('ABCD1', 'AFG3L2', 'AR', 'ARX', 'ATP1A1', 
                    'CACNA1A', 'CACNA1C', 'DHX37', 'DNM1', 'DNMT1', 
                    'EYA1', 'EZH2', 'FOXG1', 'FOXP3', 'GARS1', 
                    'GATA1', 'GATA2', 'GJA3', 'GNB1', 'HGD', 
                    'KCNH1', 'MEF2C', 'PCYT1A', 'PDE4D', 'PPP2R1A', 
                    'SCN3A', 'SNRNP200', 'SOX9', 'TGFBR1')

gostres <- gost(query = calm_gene_list, 
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

if (is.null(gostres)) {
  stop("None")
}

res_df <- gostres$result %>%
  filter(term_name != "DNA-binding transcription activator activity") %>%
  arrange(p_value) 

top_terms <- head(res_df, 10)

top_terms <- top_terms[, sapply(top_terms, is.atomic)]

top_terms$term_name <- reorder(top_terms$term_name, -top_terms$p_value)

ggplot(top_terms, aes(x = intersection_size, y = term_name, fill = p_value)) +
  geom_col() +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_gradient(low = "purple", high = "#E0E0E0", name = "p_value") +
  labs(title = "CLM better",
       x = "Count", 
       y = NULL) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey90"), 
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 1.6,
    plot.title = element_text(size = 12, hjust = 0.5)
  )