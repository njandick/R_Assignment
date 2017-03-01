# R_Assignment
## Part 1
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
# looked over both the fang_et_al_genotypes file and the snp_position file
# found that fang_et_al_genotypes has 986 columns and 2782 rows and is considered a list
# found that snp_position file has 15 columns and 983 rows, is also considered a list and is a delimeted file


# sorted fang_et_al_genotypes into the teosinte and maize genotypes, got rid of the columns in snp_position that were not needed 
teosinte_genotype <- fang_et_al_genotypes[(fang_et_al_genotypes$Group == "ZMPBA") | (fang_et_al_genotypes$Group =="ZMPIL")| (fang_et_al_genotypes$Group == "ZMPJA"), ]
maize_genotype <- fang_et_al_genotypes[(fang_et_al_genotypes$Group == "ZMMIL") | (fang_et_al_genotypes$Group =="ZMMMR")| (fang_et_al_genotypes$Group == "ZMMLR"), ]
teosinte_genotype <- t(teosinte_genotype)
maize_genotype <- t(maize_genotype)
snp_position <- snp_position [c(-2,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15)]

# For some reason R was not reading these files as data frames, so I switched it so they were, got rid of the header in the gnotype files and then merged the snp and genotype data
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

# Seperated the maize genotype data into files for each chromosome
unique(total_maize_genotypes[,1575])
for(i in 1:10){
  temp_chromosome_df<-total_maize_genotypes[total_maize_genotypes[,1575]==i,]
  
  ordered_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(ordered_temp_chromosome_df, sprintf("maize_chromosome_%d.csv",i), row.names=FALSE)
  
}
maize_chromosome_1<-read.csv("maize_chromosome_1.csv")

# Seperated the teosinte genotype data into files for each chromosome
unique(total_teosinte_genotypes[,977])
for(i in 1:10){
  temp_t_chromosome_df<-total_teosinte_genotypes[total_teosinte_genotypes[,977]==i,]
  
  ordered_temp_t_chromosome_df<- temp_t_chromosome_df[order(as.numeric(temp_t_chromosome_df$Position)),]
  
  write.csv(ordered_temp_t_chromosome_df, sprintf("teosinte_chromosome_%d.csv",i), row.names=FALSE)
}
teosinte_chromosome_1<-read.csv("teosinte_chromosome_1.csv")

# Seperated the maize data, put in descending order and changed ? to -
for(i in 1:10){
  temp_chromosome_df<-total_maize_genotypes[total_maize_genotypes[,1575]==i,]
 
  temp_chromosome_df[temp_chromosome_df=="?"] <-"-"
  
  descending_order_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(descending_order_temp_chromosome_df, sprintf("descending_position_maize_%d.csv",i), row.names=FALSE)
}
descending_position_maize_test<-read.csv("descending_position_maize_1.csv")

# Seperated the teosinte data, put in descending order and changed ? to -

for(i in 1:10){
  temp_chromosome_df<- total_teosinte_genotypes[total_teosinte_genotypes[,977]==i,]
  
  temp_chromosome_df[temp_chromosome_df=="?"]<-"-"
  
  descending_order_temp_chromosome_df<- temp_chromosome_df[order(as.numeric(temp_chromosome_df$Position)),]
  
  write.csv(descending_order_temp_chromosome_df, sprintf("descending_position_teosinte_%d.csv",i), row.names= FALSE) 
}
descending_position_test <- read.csv("descending_position_teosinte_1.csv")

## Part 2
# loaded the necessary libraries and created a ggplot showing the Number of SNPs per Chromosome
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

# created a ggplot to show the Percentage of Homozygous, Heterozygous and NA SNPs in the Data

total<- as.data.frame(total)
homozygous <- as.data.frame(matrix(ncol=1,nrow=nrow(total)))
na <-as.data.frame(matrix(ncol=1,nrow=nrow(total)))
heterozygous<-as.data.frame(matrix(ncol=1,nrow=nrow(total))) 
homozygous_df<-cbind(total[c(1:3)],homozygous,heterozygous,na,total[c(4:ncol(total))])
colnames(homozygous_df)[4]<-paste("Homozygous_percent")
colnames(homozygous_df)[5]<-paste("Heterozygous_percent")
colnames(homozygous_df)[6]<-paste("NA_percent")
homozygous_df[homozygous_df=="?/?"] <- NA
for (i in 1:nrow(homozygous_df))
{
  SNPs<-homozygous_df[i,c(7:ncol(homozygous_df))]
  counter<-0
  
  for(j in 1:ncol(SNPs)){
    if (is.na(SNPs[1,j])==FALSE){
      if((SNPs[1,j] == "A/A") | (SNPs[1,j] == "T/T")| (SNPs[1,j] == "C/C") | (SNPs[1,j] == "G/G")){
        counter<-counter +1;
      }
    }
  }
  homozygous_df[i,6]<-sum(is.na(SNPs))/ncol(SNPs)
  homozygous_df[i,4]<-counter/ncol(SNPs)
  homozygous_df[i,5]<-1-(homozygous_df[i,4])-(homozygous_df[i,6])
}
homozygous_sorted_group<- homozygous_df[order(homozygous_df$Group),]
homozygous_sorted_sample_ID<-homozygous_df[order(homozygous_df$Sample_ID),]
temp_df<- melt(homozygous_sorted_group[c(1:6)], id= c("Sample_ID", "JG_OTU", "Group"))
ggplot(temp_df, aes(x=temp_df$Sample_ID, y=temp_df$value)) + geom_bar(stat="identity")+ xlab("Sample_Id") + ylab("Proportion")+ ggtitle("Homozygous,Heterozygous and NA SNPs"


# looked at how many SNPs had multiple positions in the data and graphed it in a ggplot
snp_mult_pos_df<-data.frame()
for (i in 1:nrow=nrow(snp_position))
  {
  temp_mult_df<-snp_position[snp_position[,6]==i,]
  mult_pos<-nrow(temp_mult_df)
  mult_pos<-cbind(i,mult_pos)
  snp_mult_pos_df<-rbind(snp_mult_pos_df, mult_pos)
}
colnames(snp_mult_pos_df)
colnames(snp_mult_pos_df)<-c("i","mult_pos")
colnames(snp_mult_pos_df)
ggplot(snp_mult_pos_df)+geom_point(aes(x=i, y=mult_pos))+ ggtitle("Number of SNPs with Multiple Positions")


