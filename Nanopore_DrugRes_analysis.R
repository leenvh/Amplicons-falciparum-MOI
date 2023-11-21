require(dplyr)
require(scales)
require(data.table)
require(ggplot2)

setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION")
Variants_Exp51<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp51.csv")
Variants_Exp52<-read.csv("~/Google Drive/My Drive/PhD_LSHTM/Projects/AmpliconSequencing/MinION/variants_Exp52.csv")

ggplot(Variants_Exp51, aes(x = as.factor(gene), y = depth)) +
  geom_boxplot()+
  scale_y_sqrt(breaks = c(1, 50, 100, 500, 1000, 30000, 60000,90000))+
  theme_classic()

ggplot(Variants_Exp52, aes(x = as.factor(gene), y = depth)) +
  geom_boxplot()+
  scale_y_sqrt(breaks = c(1, 50, 100, 500, 1000, 30000, 60000,90000))+
  theme_classic()



