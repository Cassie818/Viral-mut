library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

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

gene_entrez <- bitr(esm_gene_list, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

target_genes <- gene_entrez$ENTREZID

ego <- enrichGO(gene          = target_genes,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",     
                pAdjustMethod = "BH",      
                pvalueCutoff  = 0.05,      
                qvalueCutoff  = 0.2,      
                readable      = TRUE) 

ego_sim <- simplify(ego, 
                    cutoff = 0.5, 
                    by = "p.adjust", 
                    select_fun = min)


p <- barplot(ego_sim, 
             showCategory = 5) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_distiller(palette = "BuPu", direction = -1) + 
  labs(
    title = "PLM better",
    x = "Count",
    fill = "p.adjust"
  ) +
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

print(p)