Nanostring <- read.csv("Nanostring.csv", stringsAsFactors = F) # reads the CSV file
Nanostring <- Nanostring[order(Nanostring$Avg.Count),] # getting an idea for the average counts
install.packages("dplyr")
library(dplyr)
Nanostring <-filter(Nanostring, Avg.Count > 18, Avg.Count > 18,X..Samples.above.Threshold > 32)
summary(Nanostring)
rownames(Nanostring) <- Nanostring$Probe.Name
Nanostring <- Nanostring[8:19]
Matrix <- data.matrix(Nanostring)
keep <-rowSums(Matrix) > 250
Nanostring <- Nanostring[keep,]
Matrix <- data.matrix(Nanostring)
match("Apoe",rownames(Matrix))
Matrix <- Matrix[-460,]
boxplot(log2(Matrix))
Matrix <- Matrix[, -c(1,2,3,7,8)]
Matrix <- Matrix[,-7]
Logmatrix <- log2(Matrix)
boxplot(Logmatrix)
barplot(colMeans(Logmatrix), ylim = c(0,10))
barplot(apply(Logmatrix, 2, sd))

design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)

install.packages("limma")
library(limma)
# Ordinary fit
fit <- lmFit(Logmatrix,design)
fit <- eBayes(fit)
listfit <- topTable(fit,coef=2, number = 40, sort.by = "logFC")
listfittotgenes <- topTable(fit,coef=2, number = 459, sort.by = "logFC")
listfitpval <- topTable(fit,coef=2, number = 40, sort.by = "p")
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)
# attemptin volcano plot
library(calibrate)
listfittotgenes <- mutate(listfittotgenes, Gene = row.names(listfittotgenes))

with(listfittotgenes, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2.5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(listfittotgenes, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(listfittotgenes, abs(logFC)>0.84), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(listfittotgenes, adj.P.Val<.05 & abs(logFC)>0.84), points(logFC, -log10(P.Value), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
with(subset(listfittotgenes, adj.P.Val<.05 & abs(logFC)>0.84), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))

with(subset(listfittotgenes, abs(logFC)>0.84), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))
