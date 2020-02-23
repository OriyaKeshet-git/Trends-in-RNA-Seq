library(rentrez)
library(RISmed)
library(pubmed.mineR)
library(ggplot2)
library(forecast)

#---------search the SRA database for RNA data through the years------------

types_count = c() #for RNA entries
all_count = c() # for the relation between RNA entries to all entries  

for(y in 2010:2018)
{
  tmp = c()
  for(j in 1:length(types))
  {
    query = paste("(\"biomol rna\"[Properties]) AND ", y, "[PDAT]", sep = "")
    count = entrez_search(db="sra", term=query, retmax=0)$count
    query2 = paste( y, "[PDAT]", sep = "") # number of publication for specific year
    count2 = entrez_search(db="sra", term=query2, retmax=0)$count
    tmp = c(y,count)
    tmp2 = c(y, count/count2)
    types_count = rbind(types_count, tmp)
    all_count = rbind(all_count, tmp2)
  }
}
# create data frame
df= as.data.frame(types_count)
df[,2] = as.numeric(as.character(df[,2]))
names(df) <- c( "Year", "Count")

df2= as.data.frame(all_count)
df2[,2] = as.numeric(as.character(df2[,2]))
names(df2) <- c( "Year", "Count")

# -----------------------------plot the results------------------------------- 

# for rna entries
pdf("SRA total count.pdf ", width = 13, height = 10)
theme_set(theme_classic())
ggplot(data=df, aes(x=Year, y=Count)) +
  geom_line() +
  geom_point() + 
  theme(axis.title = element_text(size = 25),
        plot.subtitle = element_text(size = 27),
        axis.text =   element_text(size = 20),
        plot.title = element_text(face = "bold", size = 30))+
        scale_fill_manual(values="e6bb9e") +
  labs(title = "SRA", y = "Number Of Entries", x = "Year", subtitle = "RNA Entries Per Year")
dev.off()

# for RNA entries / All entries
pdf("SRA.pdf ", width = 13, height = 10)
theme_set(theme_classic())
ggplot(data=df2, aes(x=Year, y=Count)) +
  geom_line() +
  geom_point() + 
  theme(axis.title = element_text(size = 25),
        plot.subtitle = element_text(size = 27),
        axis.text =   element_text(size = 20),
        plot.title = element_text(face = "bold", size = 30))+
  scale_fill_manual(values="ffe7d0") +
  labs(title = "SRA", y = "Number Of Entries", x = "Year", subtitle = "RNA Entries/All Entries ")
dev.off()
