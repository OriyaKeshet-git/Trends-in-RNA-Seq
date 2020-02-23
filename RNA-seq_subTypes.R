# This script search for RNA-seq sub methods through the years, and plot the results

library(rentrez)
library(RISmed)
library(pubmed.mineR)
library(ggplot2)
library(forecast)
library(forcats)
library(dplyr)

types = c('"targeted RNA-seq"[Title/Abstract] | "targeted RNASeq"[Title/Abstract] | "targeted RNA sequencing"[Title/Abstract] | "targeted RNA Seq"[Title/Abstract]'
          ,'"miRNA-seq"[Title/Abstract] | "miRNASeq"[Title/Abstract] | "MicroRNA sequencing"[Title/Abstract] | "miR sequencing"[Title/Abstract] | "miRNA sequencing"[Title/Abstract]'
          , '"qtl-seq"[Title/Abstract] | "qtl seq"[Title/Abstract] | "qtl sequencing"[Title/Abstract]'
          ,'"Single-Cell RNA-Seq"[Title/Abstract] | "Single-Cell RNASeq"[Title/Abstract] | "Single-Cell RNA Seq"[Title/Abstract] | "single cell rna sequencing"[Title/Abstract] | "single-cell rna sequencing"[Title/Abstract]')
types_for_plot = c("targeted RNA-seq" ,"miRNA-seq", "qtl-seq" ,"Single-Cell RNA-Seq","DGE-Seq")  

types_count = c()
# for each year search all the sub-types
for(y in 2010:2018)
{
  tmp = c()
  for(j in 1:length(types))
  {
    query = paste("(",types[j],") AND ", y, "[PDAT]", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(types_for_plot[j],y,count)
    types_count = rbind(types_count, tmp)
  }
}

#create data frame
df= as.data.frame(types_count)
df[,3] = as.numeric(as.character(df[,3]))
names(df) <- c("Types", "Year", "Count")

#plot the results
pdf("RNASEQ sub-methods.pdf", width = 10)
theme_set(theme_classic())
ggplot(data=df, aes(x=Year, y=Count, group=Types, colour=Types)) +
  geom_line() +
  theme(axis.text = element_text(family = "Times", size = 20),
        axis.text.x = element_text(family = "Times", vjust=0.6),
        axis.title = element_text(family = "Times", face = "italic", size = 25),
        plot.title = element_text(family = "Times", face = "bold", size = 30),
        legend.text = element_text(family = "Times", size = 25))+
  scale_fill_manual(values=m_colors)+
  geom_point() +
  labs(title = "RNA-Seq Sub-Methods", y = "Number Of Abstracts")
dev.off()

#------------------------------additional uses for rnaseq------------------------------------------

types = c("(Alternative Splicing[Title/Abstract] or Alternate Splicing[Title/Abstract])", 
          "(gene fusion[title/abstract] or genes fusions[title/abstract] or fused genes[title/abstract] or fused gene[title/abstract] OR gene fusions[Title/Abstract] OR genes fusions[Title/Abstract])",
          "(variants calling[Title/Abstract] or SNV[Title/Abstract])")
types_for_plot = c("Alternative Splicing","gene fusion","variants calling")

# for each year search all additional uses for rna seq

rnaseq_str = "(RNA seq[Title/Abstract] OR RNAseq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract])"
types_count = c()
for(y in 2010:2018)
{
  tmp = c()
  for(j in 1:length(types))
  {
    query = paste(types[j]," AND ",rnaseq_str," AND ",y, "[PDAT]", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(types_for_plot[j],y,count)
    types_count = rbind(types_count, tmp)
  }
}

#plot the results
df= as.data.frame(types_count)
df[,3] = as.numeric(as.character(df[,3]))
names(df) <- c("Types", "Year", "Count")

#plot the results
pdf("RNASEQ additional uses.pdf", width = 10)
theme_set(theme_classic())
ggplot(data=df, aes(x=Year, y=Count, group=Types, colour=Types)) +
  geom_line() +
  theme(axis.text = element_text(family = "Times", size = 20),
        axis.text.x = element_text(family = "Times", vjust=0.6),
        axis.title = element_text(family = "Times", face = "italic", size = 25),
        plot.title = element_text(family = "Times", face = "bold", size = 30),
        legend.text = element_text(family = "Times", size = 25))+
  scale_fill_manual(values=m_colors)+
  geom_point() +
  labs(title = "Additional Uses For RNA-Seq", y = "Number Of Abstracts")
dev.off()

