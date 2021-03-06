Nanostring <- read.csv("Nanostring.csv", stringsAsFactors = F) # reads the CSV file
Nanostring <- Nanostring[order(Nanostring$Avg.Count),] # getting an idea for the average counts
install.packages("dplyr")
library(dplyr)
# we are gonna remove genes with count less than 18 and samples above threshold 32% (above 17 in 32% of12 lanes)
Nanostring <-filter(Nanostring, Avg.Count > 18, Avg.Count > 18,X..Samples.above.Threshold > 32)
summary(Nanostring)
# we are going to give the rows of the database a name, it will help R when we will conver dataframe to Matrix to perform operation
rownames(Nanostring) <- Nanostring$Probe.Name
# getting read of all the colum we do not need and leave only the samples
Nanostring <- Nanostring[8:19]
# we are gonna tranform the databse into a Matrix 
Matrix <- data.matrix(Nanostring)
# getting read of low count genes by performing row sums function in the matrix (this function cannot be performed on sataframes, then since the Matrix and the Nanostirng dataframe have the same rownames we can only keep the  rownames generated in line 15 in the nanostirng databse in line 16, this is done with subsetting
keep <-rowSums(Matrix) > 250
Nanostring <- Nanostring[keep,]
#reconver the dataframe into a matrix again, i go back and forth because there are operations in the dataframe that you cannot perform in matrix
Matrix <- data.matrix(Nanostring)
#Apoe is the only outlier gene identified in anothe script, I am going to remove that
match("Apoe",rownames(Matrix))
Matrix <- Matrix[-460,]
# I am going to boxplot the matrix to check how the counts are cmoing by
boxplot(log2(Matrix))
# removing from the matrix samples I am not interested in
Matrix <- Matrix[, -c(1,2,3,7,8)]
Matrix <- Matrix[,-7]
Logmatrix <- log2(Matrix)
# check again how the matrix looks like without the samples I am not intersted in
boxplot(Logmatrix)
#boxplotting the mean count of each samples, and then the SD, if the samples are weel normalized they should not vary so much
barplot(colMeans(Logmatrix), ylim = c(0,10))
barplot(apply(Logmatrix, 2, sd))
# from this step on we are gonna calculate differential gene expression (DEG) between different goups using the Limma package
# Design is the contrast matrix, I have not understood it very much but it works
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)

install.packages("limma")
library(limma)
# Ordinary fit they are the operation required for limma
fit <- lmFit(Logmatrix,design)
fit <- eBayes(fit)
listfit <- topTable(fit,coef=2, number = 40, sort.by = "logFC")
# with this command I ll have the list of the fold change of all the 459 genes sorted by fold change (I put the parameter number = 459 = number of rows = number of genes)
listfittotgenes <- topTable(fit,coef=2, number = 459, sort.by = "logFC")
# sam as aboves but sorted by P value
listfitpval <- topTable(fit,coef=2, number = 40, sort.by = "p")
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)
# attemptin volcano plot
library(calibrate)
# adding a new column to the file generated by limma, this colum will be used as label for different genes in the vulcano plot
listfittotgenes <- mutate(listfittotgenes, Gene = row.names(listfittotgenes))
# this is actually the volcano plot plotting different color, read the code for more explanation
with(listfittotgenes, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2.5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(listfittotgenes, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(listfittotgenes, abs(logFC)>0.84), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(listfittotgenes, adj.P.Val<.05 & abs(logFC)>0.84), points(logFC, -log10(P.Value), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot, showing both > 1.5 fold change and adj p value < 0.05
with(subset(listfittotgenes, adj.P.Val<.05 & abs(logFC)>0.84), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))
# Label points with the textxy function from the calibrate plot, showing both > 1.5 fold change, adj p value is not taken into consideration
with(subset(listfittotgenes, abs(logFC)>0.84), textxy(logFC, -log10(P.Value), labs=Gene, cex=.8))
