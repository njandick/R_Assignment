# R_Assignment
# Part 1
fang_et_al_genotypes <- read.table("/Users/Nicole/desktop/BCB546X-Spring2016/R_Assignment/fang_et_al_genotypes.txt", header=TRUE) 
is.list (fang_et_al_genotypes)
typeof(fang_et_al_genotypes)
length(fang_et_al_genotypes) 
nrow(fang_et_al_genotypes)
dim(fang_et_al_genotypes)
colnames(fang_et_al_genotypes)
rownames(fang_et_al_genotypes)
head(fang_et_al_genotypes)
snp_position <- read.delim("/Users/Nicole/desktop/BCB546X-Spring2016/R_Assignment/snp_position.txt")
is.list (snp_position)
typeof(snp_position)
length(snp_position) 
nrow(snp_position)
dim(snp_position)
colnames(snp_position)
head(snp_position)
teosinte_genotype <- fang_et_al_genotypes[(fang_et_al_genotypes$Group == "ZMPBA") | (fang_et_al_genotypes$Group =="ZMPIL")| (fang_et_al_genotypes$Group == "ZMPJA"), ]
maize_genotype <- fang_et_al_genotypes[(fang_et_al_genotypes$Group == "ZMMIL") | (fang_et_al_genotypes$Group =="ZMMMR")| (fang_et_al_genotypes$Group == "ZMMLR"), ]
teosinte_genotype <- t(teosinte_genotype)
maize_genotype <- t(maize_genotype)
teosinte_genotype <- t(teosinte_genotype)
snp_position <- snp_position [c(-2,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15)]
snp_position <- data.frame(snp_position)
maize_genotype <- data.frame(maize_genotype)
maize_genotype_no_header <- maize_genotype[-(1:3),]
teosinte_genotype_no_header <- teosinte_genotype[-(1:3),]
maize_genotype_no_header<-cbind(colnames(fang_et_al_genotypes[4:ncol(fang_et_al_genotypes)]),maize_genotype_no_header)
colnames(maize_genotype_no_header)[1]<-paste("SNP_ID")
total_maize_genotypes <- merge(maize_genotype_no_header, snp_position, by= "SNP_ID")
teosinte_genotype_no_header<-cbind(colnames(fang_et_al_genotypes[4:ncol(fang_et_al_genotypes)]),teosinte_genotype_no_header)
colnames(teosinte_genotype_no_header)[1]<-paste("SNP_ID")
total_teosinte_genotypes <- merge(teosinte_genotype_no_header, snp_position, by="SNP_ID")
unique(total_maize_genotypes[,1575])
for(i in 1:10){
  temp_chromosome_df<-total_maize_genotypes[total_maize_genotypes[,1575]==i,]
  
  ordered_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(ordered_temp_chromosome_df, sprintf("maize_chromosome_%d.csv",i), row.names=FALSE)
  
}
maize_chromosome_1<-read.csv("maize_chromosome_1.csv")
unique(total_teosinte_genotypes[,977])
for(i in 1:10){
  temp_t_chromosome_df<-total_teosinte_genotypes[total_teosinte_genotypes[,977]==i,]
  
  ordered_temp_t_chromosome_df<- temp_t_chromosome_df[order(as.numeric(temp_t_chromosome_df$Position)),]
  
  write.csv(ordered_temp_t_chromosome_df, sprintf("teosinte_chromosome_%d.csv",i), row.names=FALSE)
}
teosinte_chromosome_1<-read.csv("teosinte_chromosome_1.csv")

for(i in 1:10){
  temp_chromosome_df<-total_maize_genotypes[total_maize_genotypes[,1575]==i,]
 
  temp_chromosome_df[temp_chromosome_df=="?"] <-"-"
  
  descending_order_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(descending_order_temp_chromosome_df, sprintf("descending_position_maize_%d.csv",i), row.names=FALSE)
}
descending_position_maize_test<-read.csv("descending_position_maize_1.csv")


for(i in 1:10){
  temp_chromosome_df<- total_teosinte_genotypes[total_teosinte_genotypes[,977]==i,]
  
  temp_chromosome_df[temp_chromosome_df=="?"]<-"-"
  
  descending_order_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(descending_order_temp_chromosome_df, sprintf("descending_position_teosinte_%d.csv",i), row.names= FALSE) 
}
descending_position_test <- read.csv("descending_position_teosinte_1.csv")

# Part 2
library(reshape2)
library(ggplot2)
melt(snp_position)
melt(fang_et_al_genotypes)
snp_count_chromosome_df<-data.frame()
for (i in 1:10){
  temp_chromosome_df<-snp_position[snp_position[,2]==i,]
  
  chromosome_count<-nrow(temp_chromosome_df)
  
  chromosome_count<-cbind(i,chromosome_count)
  
  snp_chromosome_count_df<-rbind(snp_chromosome_count_df, chromosome_count)
}
colnames(snp_count_chromosome_df)
colnames(snp_count_chromosome_df)<-c("Chromosome","Number_of_SNPs")
colnames(snp_count_chromosome_df)
ggplot(snp_count_chromosome_df)+geom_point(aes(x=Chromosome, y=Number_of_SNPs))+ ggtitle("Number of SNPs per Chromosome")


