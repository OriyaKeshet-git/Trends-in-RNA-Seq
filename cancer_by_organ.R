
# Count the number of rnaseq abstracts related to cancer types, and plot the results.
# The cancer types are grouped by body location or system.
#-------cancer types list from :https://www.cancer.gov/types/by-body-location-----------------

library(assertr)
library(rentrez)

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
queryVector <-gsub("  ", " OR ", query_no_space_vector, fixed=TRUE)

rnaseq_str = "(RNA seq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract] OR rnaseq[Title/Abstract]) "

# for each year, iterate all cancer types and count the number of relevent abstracts
typesCount = c()
for(y in 2010:2018)
{
  for(i in 1:length(queryVector))
  {
    query = paste(rnaseq_str," AND (", queryVector[i],") AND ", y, "[PDAT]", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(y,names(queryVector)[i],count)
    typesCount = rbind(typesCount, tmp)
  }
  print(typesCount)
}

# create data frame
colnames(typesCount) <- c("Year", "Type", "Value")
df <- as.data.frame(typesCount)
df[,3] <- as.numeric(as.character(typesCount[,3]))


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

# plot the results
pdf("cancerByOrgan.pdf", width = 13, height = 10)

theme_set(theme_classic())

g<-ggplot(df, aes(x=Year, y=Value)) 
g + geom_bar(aes(fill=fct_reorder(Type, Value, sum, desc=TRUE)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6, angle=65), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual("Body Location/System", values=m_colors) +
  labs(y="Number Of Entries" )


dev.off()




