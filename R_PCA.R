# 分组展示样本PCA结果
library(ggplot2)
library(ggforce)
#Data input
args <-commandArgs(trailingOnly = TRUE)
myData <-read.table("PCA_output.eigenvec", sep=' ' , header=TRUE,)
eigenval <-read.table("PCA_output.eigenval")
#Sum the eigenvalues
total <-sum(eigenval)
#Set output filename
filename <-paste("PCA_output.eigenvec", '.pdf', sep='')
#Set output format
pdf(filename, width=8, height=6)
#Set axis labels
xlabel <-paste('PC1(', round((eigenval$V1[1]/total)*100, 2), '%)', sep='')
ylabel <-paste('PC2(', round((eigenval$V1[2]/total)*100, 2), '%)', sep='')

#Plot PCA
ggplot(myData, aes(x=PC1, y=PC2, color=Category)) + 
  geom_point() + 
  xlab(xlabel) + 
  ylab(ylabel) +
  scale_size_identity()
dev.off()
