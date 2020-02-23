
# This script calculate for non-coding rna molecules, the number of abstracts realted to 
# specific diseases, and create stacked bar plot for the results.

library(assertr)
library(rentrez)
library(ggplot2)
library(dplyr)
library(forcats)


types = c("rRNA OR ribosomal RNA OR 12S OR 16s OR 25S OR 25S OR 5S OR 5.8S",
          "tRNA OR mt-tRNA", "lincRNA OR lncRNA", 
          "miRNA OR miR OR microRNA", "piRNA","scRNA","snRNA", "sRNA", "snoRNA","miscRNA",
          "circRNA OR Circular RNA", "Circulating")
types_for_plot = c("rRNA","tRNA", "lincRNA", 
                   "microRNA", "piRNA",  "scRNA","snRNA", "sRNA", "snoRNA","misc_RNA",
                   "circRNA", "circulating RNA")

types = gsub("OR", "[Title/Abstract] OR", types) 

# read the diseases table
table = read.csv("C:/Users/Racheli/Documents/mini-project-bioinformatics/cancers-by-Body-Location-System.csv", header = TRUE, check.names = F)
table$Skin[8] = "Sezary Syndrome"
names(table)[1] = "Breast"
t = t(table)
t[t == ""] <- NA

#create query from each row
queryVector = col_concat(t, sep = "  ")
query_no_na <- gsub("NA", "", queryVector)
#trim spaces from the end
query_no_space_vector = trimws(query_no_na, which = "both" )
#add or between terms
queryVector <-gsub("  ", "[Title/Abstract] OR ", query_no_space_vector, fixed=TRUE)


rnaseq_str = "(RNA seq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract] OR rnaseq[Title/Abstract]) "

# count the baseline for each disease
general_disesesCount = c()
for(i in 1:length(queryVector))
{
  query = paste(rnaseq_str," AND (", queryVector[i],"[Title/Abstract])", sep = "")
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  general_disesesCount = c(general_disesesCount, count)
}

# count the baseline for each rna type
general_rna_count = c()
for(i in 1:length(types))
{
  query = paste(rnaseq_str," AND (", types[i],"[Title/Abstract])", sep = "")
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  general_rna_count = c(general_rna_count, count)
}

chi_squre_results = c()
typesCount = c()
# for each rna type, iterate all the diseases and send query combined from disease and rna
for(j in 1:length(types))
{
  counts_per_rna = c()
  #iterate diseases list
  for(i in 1:length(queryVector))
  {
    query = paste(rnaseq_str," AND (", queryVector[i],"[Title/Abstract]) AND (", 
                    types[j], "[Title/Abstract]) ", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(names(queryVector)[i], types_for_plot[j], count)
    counts_per_rna = c(counts_per_rna, count)
    typesCount = rbind(typesCount, tmp)
  }
  # for each rna type, calculate chi-squre test for godness of fit
    chi_squre_table = rbind(counts_per_rna, general_disesesCount)
    chi_squre_table = chi_squre_table[,chi_squre_table[2,] != 0]
     result = chisq.test(chi_squre_table, simulate.p.value = T)
     chi_squre_results = c(chi_squre_results, result$p.value)
}


# create data frame from the results
colnames(typesCount) <- c("Cancer", "RNA", "Value")
df <- as.data.frame(typesCount)
df[,3] <- as.numeric(as.character(typesCount[,3]))

# list of colors for plot
m_colors = c("#f2a8bb",
             "#ffd1c4",
             "#e7ae8e",
             "#e1cd96",
             "#fff1c7",
             "#9fb28a",
             "#e6ffe7",
             "#a7e8c5",
             "#88b5a1",
             "#82d7c5",
             "#77b3d1",
             "#acdeff",
             "#e6dfff",
             "#b8aae0",
             "#c3a1bd")

# create stacked bar plot, and save it to pdf file
pdf("nonCoding_vs_Diseases.pdf", width = 15, height = 10)

theme_set(theme_classic())

g<-ggplot(df, aes(x=RNA, y=Value))
g + geom_bar(aes(fill=fct_reorder(Cancer, Value, sum, desc=TRUE)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6, angle=65), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual("Body Location/System", values=m_colors) +
  labs(y="Number Of Entries")

dev.off()

