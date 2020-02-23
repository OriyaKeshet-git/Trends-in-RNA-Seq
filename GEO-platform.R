library(ggplot2)

#------------ Count the number of uses of different array platforms for GEO entries -------------------
# each platform has sevrel accesion numbers, so the table contain column for each type with all accessions
table <- read.csv("C:/Users/Racheli/Documents/mini-project-bioinformatics/platform-types.csv", header = T)

#clean the data
for(i in 1:ncol(table))
{
  table[,i] = paste0(table[,i], " ")
}
table[table == " "] <- NA

names = names(table)


typesCount = c()
for (y in 2010:2018)  
{
  # read the file with the GEO results
  fileName = paste("gds_result", y, ".txt", sep = "")
  file = readLines(fileName)
  #seperate the file to entries, each row is an entry
  k = cumsum(file == '')
  records = by(file, k, paste, collapse = '\n')
  # extract only 'Expression profiling by array' entries
  exprProf = records[grepl("Expression profiling by array", records, fixed = TRUE)]
  # for each column(all accessions numbers related to a pltform)
  for(i in 1:dim(table)[2])
  {
    tmp = c()
    index = sapply(as.vector(table[,i]), function(x) if(!is.na(x)) 
      length(exprProf[grep(x, exprProf, fixed = T)]))
    count = sum(unlist(index)) # sum the number of entries has one of the platform accesion
    tmp = c(y, names[i], count)
    print(tmp)
    typesCount = rbind (typesCount, tmp)
  }
  print(typesCount)
}

# 
a = typesCount
a = a[a[,2] != "Exiqon",]
# create data frame
colnames(a) <- c("Year","Company", "Value")
df <- as.data.frame(table)
df[,4] <- as.numeric(as.character(df[,4]))

m_colors = c("#df4a7a", "#c97b7a", "#de5137","#d08935", "#a78d57","#d2d23e","#67993f",
             "#76d854","#529477","#6387d7","#777ba7","#b159e0","#d6aad9","#bd6cac",
             "#db49ba")

#plot the results
pdf("GEO-platforms.pdf", width = 13, height = 10)

theme_set(theme_classic())

g <- ggplot(df, aes(x = Year, y = Value))
g + geom_bar(aes(fill=fct_reorder(Company, Value, sum, desc=TRUE)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(family = "Times", face = "italic", size = 25),
        plot.title = element_text(family = "Times", face = "bold", size = 30),
        axis.text.x = element_text(angle=0, vjust=0.6), 
        axis.text =   element_text(family = "Times",size = 20),
        legend.text = element_text(family = "Times",size = 22),
        legend.title = element_text(family = "Times",size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual(values=m_colors) +
  labs(title="Expression Profiling By Array Platforms", 
       subtitle="", fill = "Companys", y = "Number Of Entries", x = "Year",
       family = "Times", face = "italic") 

dev.off()

