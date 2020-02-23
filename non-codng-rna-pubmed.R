library(ggplot2)
library(dplyr)
library(rentrez)

# -------------search for rnaseq abstract related to non-coding rna molecules--------------------

types = c("rRNA OR ribosomal RNA OR 12S OR 16s OR 25S OR 25S OR 5S OR 5.8S",
          "tRNA OR mt-tRNA", "lincRNA OR lncRNA", 
          "miRNA OR miR OR microRNA", "piRNA","scRNA","snRNA", "sRNA", "snoRNA","miscRNA",
          "circRNA OR Circular RNA", "Circulating")
types_for_plot = c("rRNA","tRNA", "lincRNA", 
          "microRNA", "piRNA",  "scRNA","snRNA", "sRNA", "snoRNA","miscRNA",
          "circRNA", "circulating RNA")


typesCount = c()
types = gsub("OR", "[Title/Abstract] OR", types)
rnaseq_str = "(RNA seq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract]) "

# for each year, iterate through non-coding types and count the number of abstracts
for(y in 2010:2018)
{
  for(i in 1:length(types))
  {
    query = paste(rnaseq_str," AND (", types[i], "[Title/Abstract]) AND ", y, "[PDAT]", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(y,types_for_plot[i],count)
    typesCount = rbind(typesCount, tmp)
    
  }
  print(typesCount)
}

# create data frame
colnames(typesCount) <- c("Year", "Molecule", "Value")
df <- as.data.frame(typesCount)
df[,3] <- as.numeric(as.character(typesCount[,3]))

m_colors = c("#f2a8bb",
             "#ffd1c4",
             "#fff1c7",
             "#e6ffe7",
             "#a7e8c5",
             "#88b5a1",
             "#82d7c5",
             "#77b3d1",
             "#acdeff",
             "#f7faff",
             "#e6dfff",
             "#b8aae0",
             "#c3a1bd")
# plot the results 
pdf("nonCodingRNA.pdf", width = 15, height = 10)


theme_set(theme_classic())

g<-ggplot(df, aes(x=Year, y=Value)) 
g + geom_bar(aes(fill=fct_reorder(Molecule, Value, sum, desc=TRUE)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6, angle=65), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual("RNA Molecule", values=m_colors) +
  labs(y="Number Of Entries" )
 

dev.off()



