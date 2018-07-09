#JRB reproduce combined gsea output

source("replot_functions.R")

## load and combine KEGG and GO GSEA results
gsea_KEGG = load_gsea(rnk_file = "KEGG/gsea_data.gsea_data.rnk", 
                      sets_file = "KEGG/gene_sets.gmt", 
                      results_file = "KEGG/results.edb")
gsea_GO = load_gsea(rnk_file = "GO/metric2.rnk", 
                    sets_file = "GO/gene_sets.gmt", 
                    results_file = "GO/results.edb")

gsea_full = combine_gsea(gsea_KEGG, gsea_GO)

## define and extract sets of interest
#control gene sets to plot here
set_oi = c("KEGG_LYSOSOME", "GO_PLASMA_MEMBRANE_PROTEIN_COMPLEX")

full_dt = rbindlist(lapply(set_oi, function(set_name){
  extract_set(set_name, gsea_full)
}))

## prepare plots

#control colors used
g_colors = seqsetvis::safeBrew(length(set_oi), "Dark2")
names(g_colors) = set_oi

#enforce consistent xlimits
xlim = c(0, max(full_dt$x))

p_prof = ggplot(full_dt[type == "profile"], aes(x = x, y = y, color = set)) + 
  geom_path(size = 2) +
  coord_cartesian(xlim = xlim, expand = FALSE) +
  scale_color_manual(values = g_colors) +
  labs(y = "ES", x = "", color = "") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(size = 1.5), 
        axis.ticks = element_line(size = 1.5, colour = "black"),
        axis.ticks.length = unit(.45, "cm"), 
        legend.position = "top", 
        legend.text = element_text(size = 10),
        legend.direction = "vertical",
        legend.justification = "right",
        axis.title = element_text(size = 16))

hidt = full_dt[type == "hit_index"]

p_ticks = ggplot(hidt, 
                 aes(x = x, 
                     y = set, 
                     color = set)) + 
  coord_cartesian(xlim = xlim, expand = FALSE) +
  geom_tile() + 
  labs(x = "hit index", y = "") +
  scale_color_manual(values = g_colors) +
  guides(color = "none") +
  theme(panel.background = element_blank(), axis.line.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.line.x = element_line(size = 1.5), 
        axis.ticks.x = element_line(size = 1.5, colour = "black"),
        axis.ticks.length = unit(.45, "cm"),
        axis.title = element_text(size = 16))

## assemble plots to final

g2 <- ggplotGrob(p_prof)
g3 <- ggplotGrob(p_ticks)

# force xcoords to match
g2$widths <- unit.pmax(g2$widths, g3$widths)
g3$widths <- unit.pmax(g2$widths, g3$widths)

#use heights to redistribute
grid.arrange(g2, g3, heights = c(3, 1.5))



