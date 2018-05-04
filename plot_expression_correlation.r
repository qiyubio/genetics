GE=read.delim('genes_4samples_wt.fpkm_known_gene_plot_order')

#pseudocount
s=1e-3

new<-GE[which(GE$prdm9_11_FPKM >0.1 | GE$prdm9_12_FPKM >0.1 |GE$hop2_11_FPKM >0.1 | GE$hop2_12_FPKM >0.1|GE$ewt_FPKM >0.1),]
GE<-new
round(cor(GE[grep('FPKM',names(GE))]),9)
round(cor(GE[grep('FPKM',names(GE))],method='spearman'),9)

panel.cor <- function(x, y, method='spearman', digits=2, prefix="R=", cex.cor, ...)
  #panel.cor <- function(x, y, method='pearson', digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)
}


png("FPKM_scatter_cor_log10_spearman_genes.png",res=200,width = 1300, height = 1300)
pairs(log10(GE[grep('FPKM',names(GE))]+s), horInd = c(1,2,3,4,5),verInd=c(1,2,3,4,5), lower.panel=function(...) smoothScatter(..., nrpoints=0, add=TRUE,colramp = colorRampPalette(c("white","blue","cyan","green","yellow"))), upper.panel=panel.cor, labels=c(expression(paste(italic("Prdm9"),"-/- #1")),expression(paste(italic("Prdm9"),"-/- #2")),expression(paste(italic("Hop2"),'-/- #1')),expression(paste(italic("Hop2"),'-/- #2')),"WT adult",size=100))
dev.off()

