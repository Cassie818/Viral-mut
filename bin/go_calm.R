library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)

calm_gene_list <- c('ABCD1', 'AFG3L2', 'AR', 'ARX', 'ATP1A1', 
                    'CACNA1A', 'CACNA1C', 'DHX37', 'DNM1', 'DNMT1', 
                    'EYA1', 'EZH2', 'FOXG1', 'FOXP3', 'GARS1', 
                    'GATA1', 'GATA2', 'GJA3', 'GNB1', 'HGD', 
                    'KCNH1', 'MEF2C', 'PCYT1A', 'PDE4D', 'PPP2R1A', 
                    'SCN3A', 'SNRNP200', 'SOX9', 'TGFBR1')

gene_entrez <- bitr(calm_gene_list, 
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
                    cutoff = 0.7, 
                    by = "p.adjust", 
                    select_fun = min)


p <- barplot(ego_sim, 
             showCategory = 5) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + 
  labs(
    title = "CLM better",
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
