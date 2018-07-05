#make a volcano plot like https://bioinformatics.stackexchange.com/questions/4080/volcano-plot-in-r
#NOTE I tried with base R graphics and it got messy fast, switched to ggplot2 and now life is good.


if(!file.exists("volcano.RData")){
    system("wget http://galaxy.med.uvm.edu/static/shared/volcano.RData --no-check-certificate")
}
load("volcano.RData")

gene_list$group = "none"
gene_list[left_biased,]$group = "left"
gene_list[right_biased,]$group = "right"

gene_list$pvalue[gene_list$pvalue == 0] = min(gene_list$pvalue[gene_list$pvalue != 0])

#jrb two issues:
#jrb 1) your pvalues are in linear scale and this expects log10
#jrb 2) pvalues are rounded down to 0 around 10^-306 which breaks log10

library(ggplot2)
library(ggrepel)
library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

p = ggplot(gene_list, aes(x = log2FoldChange, y = pvalue, color = group)) +
    geom_point(alpha = .2) +
    scale_y_continuous(trans = reverselog_trans(10)) +
    scale_color_manual(values = c("none" = "black", "left" = "blue", "right" = "red")) +
    theme_classic() +
    labs(x = "Fold Change (log2)", y = "P value") + guides(color = "none")

p
#jrb add up/down annotation
de_labs = paste0(c("Down\n", "Up\n"), c(sum(left_biased), sum(right_biased)))
p = p + annotate("label", x = c(-12, 12), y = rep(10^-250, 2), label = de_labs)
p
gene_list$x = gene_list$X
gene_list$X = toupper(gene_list$x)
goi = merge(gene_list, interesting_genes, by.x = "X", by.y = "V1")

plab_simple =  p + geom_text_repel(data = goi, 
                                   aes(label = X),
                                   color = "black")

plab_allLines =  p + geom_text_repel(data          = goi,
                                     color = "black",
                                  nudge_y = 10^2,
                                  aes(label = X),
                                  segment.size  = 0.2,
                                  segment.color = "grey20",
                                  direction     = "x") 
plab_simple
plab_allLines
# interesting_genes = interesting_genes$V1
# 
# fold_changes = gene_list$log2FoldChange
# names(fold_changes) = gene_list$X
# 
# 
# pvalues = gene_list$pvalue
# pvalues[pvalues == 0] = min(pvalues[pvalues != 0])
# pvalues = log10(pvalues)
# names(pvalues) = gene_list$X
# 
# fc_interesting_genes = fold_changes[interesting_genes]
# pval_interesting_genes = pvalues[interesting_genes]
# 
# #jrb i'm adding xspan and yspan variables
# xspan = round(max(abs(fold_changes)))
# xspan = c(-xspan, xspan)
# yspan = c(0, floor(min(pvalues)))
# 
# ### plot corpus
# plot(fold_changes, pvalues,
#      pch = point_type, cex = point_size, col = palette[1],
#      xlab = "Fold Change (log2)", ylab = "P value",
#      xlim = xspan, ylim = yspan, yaxt="n")
# aty <- axTicks(2)
# labels <- sapply(aty,function(i)
#     as.expression(bquote(10^ .(i)))
# )
# labels[6] <- 1
# axis(2,at=aty,labels=labels)
# 
# ### plot deferentially expressed genes in colour + genes on interest with black
# points(fold_changes[left_biased], pvalues[left_biased],
#        pch = point_type, cex = point_size, col = palette[2])
# points(fold_changes[right_biased], pvalues[right_biased],
#        pch = point_type, cex = point_size, col = palette[3])
# points(fc_interesting_genes, pval_interesting_genes, pch = point_type, cex = point_size)
# 
# ### add text
# #jrb labelling more than 4 points with this method is a real headache
# #jrb i'm going to leave it at the first 4 genes (already annoying since this assumed certain coordinates)
# #jrb if labelling everything is real priority we can work on that more.
# 
# xpos = function(xcoord, xspan = xspan){
#     min(xspan) + diff(xspan) * rep(c(1, 5)/6, each = 2)
# }
# 
# ypos = function(xcoord, xspan = xspan){
#     min(yspan) + diff(yspan) * rep(c(1, 5)/6, each = 2)
# }
# 
# 
# labelx <- rep(xspan, each = 2) * .66
# labely <- c(-1.5, -0.5, -0.5, -1.5)
# labels <- c('Il17s', 'Casz1', 'Rorc', 'Il22')
# text(labelx, labely, labels)
# text(-4, -4, "down\n194 genes")
# text(4, -4, "up\n156 genes")
# 
# ### add lines between text and data points
# for(i in 1:4){
#     lines(c(labelx[i], fc_interesting_genes[i]), c(labely[i], pval_interesting_genes[i]))
# }
