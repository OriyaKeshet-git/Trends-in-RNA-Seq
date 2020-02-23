library(ggplot2)
library(dplyr)
library(rentrez)
#----------------create graph for All RNA related articles per year divided by  All published articles pre year---------------------------

#Initialize count:
artCount = c()

#Find all published articles and all RNA related published articles per year:
for(y in 2010:2018)
{
  query = paste( y, "[PDAT]", sep = "")
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  tmp = c(y, count)
  query = paste( "(RNA-Seq[Title/Abstract]|RNAseq[Title/Abstract]|RNA seq[Title/Abstract]|RNA sequencing[Title/Abstract]|RNAsequencing[Title/Abstract]|RNA-sequencing[Title/Abstract]) AND ", y, "[PDAT]", sep = "")
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  temp = c(y, count)
  divided  = c(y, temp[2]/tmp[2])
  artCount = rbind(artCount, divided)
  print(artCount)
}

#Create dataframe for the graph:
colnames(artCount) <- c("Year",  "Value")
df <- as.data.frame(artCount)
df[,2] <- as.numeric(as.character(df[,2]))
df[,1] <- as.numeric(as.character(df[,1]))


pdf("All Articles VS. RNA-seq Articles- Divided.pdf", width = 15, height = 10)

theme_set(theme_classic())
ggplot(df, aes(x=df[,1], y=df[,2])) + labs(x = "Years", y = "Published RNA-seq Articles divided by All Published Articles", fill = NULL)+  #, title = "All Published Articles VS. All RNA-seq Published Articles")+
  theme(axis.title = element_text(face = "italic", size = 20),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6), 
        axis.text =   element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm"))+
  #scale_fill_manual(values=m_colors) +
  geom_point(size=6, alpha=0.5)+ geom_path()

dev.off()
