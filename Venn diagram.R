install.packages("VennDiagram")
library(VennDiagram)
limma <- read.csv("genesusinglimma.csv", stringsAsFactors = F)
nsolver <- read.csv("genesusingnsolver.csv", stringsAsFactors = F)
limmafoldchange <- subset(limma, abs(logFC) > 0.58)
nsolverfoldchange <- subset(nsolver, abs(logFC) > 0.58)
pdf("venn_diagram2.pdf")
plotgenes <- venn.plot <- venn.diagram(list(nsolverfoldchange$Probe.Name, limmafoldchange$Gene), NULL, fill=c("red", "green"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Nsolver", "Limma"))
grid.draw(plotgenes)
dev.off()
limmalist <- limmafoldchange$Gene
nsolverlist <- nsolverfoldchange$Probe.Name
intersectlist <- intersect(limmalist,nsolverlist)
intersectlist
