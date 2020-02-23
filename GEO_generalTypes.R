library(ggplot2)
library(CGPfunctions)
library(forcats)
types = c("Expression profiling","Genome binding/occupancy", "Genome variation profiling",
          "Methylatio profailing", "Non-coding RNA profiling", "Other", "Protein profiling",
          "SNP genotyping", "Third-party reanalysis")


typesCount = c()

for (y in 2010:2018)  
{
  fileName = paste("gds_result", y, ".txt", sep = "")
  file = read.delim(fileName, header = F, sep = "\n")
  
  tmp = c()
  for(i in 1:length(types))
  {
    count = sum(grepl(types[i], file[,1], fixed = T))
    tmp = c(y, types[i], count)
    typesCount = rbind (typesCount, tmp)
  }
  
}

colnames(typesCount) <- c("Year","Type", "Value")
df <- as.data.frame(typesCount)
df[,3] <- as.numeric(as.character(typesCount[,3]))


m_colors = c("#f2a8bb",
             "#ffd1c4",
             "#fff1c7",
             "#e6ffe7",
             "#82d7c5",
             "#77b3d1",
             "#acdeff",
             "#f7faff",
             "#e6dfff",
             "#b8aae0")
             
pdf("GEOTypes_stackedBarChart.pdf", width = 10)

theme_set(theme_classic())

p<-ggplot(df, aes(x=Year, y=Count, group=Exp_Type)) 
p + geom_line(aes(color=Exp_Type), size=2)+
  geom_point(aes(color=Exp_Type), size=4) +
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6, angle=65), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_color_manual("Experiment Type", values=m_colors) +
  labs(y="Number Of Entries" )


dev.off()


#-------------make slope charts for all types -----------------------------

table = read.csv("C:/Users/Racheli/Documents/mini-project-bioinformatics/GEO_types.csv", header = F)
table = table[, -c(1,2,3,4,5,7,8,9)]

expProf = as.vector(table[,1])[1:3]
genBin = as.vector(table[,2])[1:4]
genVar = as.vector(table[,3])[1:4]
meth = as.vector(table[,4])[1:4]
nonCoding = as.vector(table[,5])[1:3]
proProf = as.vector(table[1:2,6])

expProfCount = c()
genBinCount = c()
genVarCount = c()
methCount = c()
nonCodingCount = c()
proProfCount = c()

for (y in 2010:2018)  
{
  fileName = paste("gds_result", y, ".txt", sep = "")
  file = read.delim(fileName, header = F, sep = "\n")
  
  #Expression profiling
  tmp = c()
  for(i in 1:length(expProf))
  {
    
    count = sum(grepl(expProf[i], file[,1], fixed = TRUE))
    if(expProf[i] != "Expression profiling by array" &&  expProf[i] != "Expression profiling by high throughput sequencing")
    {
      tmp = c(y, "Other", count)
    }
    else
    {
      tmp = c(y, expProf[i], count)
      
    }
    expProfCount = rbind (expProfCount, tmp)
  }
  
  #Genome Binding
  tmp = c()
  for(i in 1:length(genBin))
  {
    count = sum(grepl(genBin[i], file[,1], fixed = TRUE))
    tmp = c(y, genBin[i], count)
    genBinCount = rbind (genBinCount, tmp)
  }
  
  #Genome varation
  tmp = c()
  for(i in 1:length(genVar))
  {
    count = sum(grepl(genVar[i], file[,1], fixed = TRUE))
    tmp = c(y, genVar[i], count)
    genVarCount = rbind (genVarCount, tmp)
  }
  
  #Methylation
  tmp = c()
  for(i in 1:length(meth))
  {
    count = sum(grepl(meth[i], file[,1], fixed = TRUE))
    tmp = c(y, meth[i], count)
    methCount = rbind (methCount, tmp)
  }
  
  #non-coding RNA
  tmp = c()
  for(i in 1:length(nonCoding))
  {
    count = sum(grepl(nonCoding[i], file[,1], fixed = TRUE))
    tmp = c(y, nonCoding[i], count)
    nonCodingCount = rbind (nonCodingCount, tmp)
  }
  
  #protein profiling
  tmp = c()
  for(i in 1:length(proProf))
  {
    count = sum(grepl(proProf[i], file[,1], fixed = TRUE))
    tmp = c(y, proProf[i], count)
    proProfCount = rbind (proProfCount, tmp)
  }
  
}


colnames(expProfCount) <- c("Year","Type", "Value")
ep <- as.data.frame(expProfCount)
ep[,3] <- as.numeric(ep[,3])

colnames(genBinCount) <- c("Year","Type", "Value")
gb <- as.data.frame(genBinCount)
gb[,3] <- as.numeric(gb[,3])

colnames(genVarCount) <- c("Year","Type", "Value")
gv <- as.data.frame(genVarCount)
gv[,3] <- as.numeric(gv[,3])

colnames(methCount) <- c("Year","Type", "Value")
mt <- as.data.frame(methCount)
mt[,3] <- as.numeric(mt[,3])

colnames(nonCodingCount) <- c("Year","Type", "Value")
nc <- as.data.frame(nonCodingCount)
nc[,3] <- as.numeric(nc[,3])

colnames(proProfCount) <- c("Year","Type", "Value")
pp <- as.data.frame(proProfCount)
pp[,3] <- as.numeric(pp[,3])


m_colors = 
  c("#f2a8bb",
   "#77b3d1",
   "#b8aae0",
   "#c3a1bd")

pdf("GEOTypes_detailed.pdf", width = 15, height = 10)


theme_set(theme_classic())

g <- ggplot(ep, aes(x = Year, y = Value))
g + geom_bar(aes(fill=fct_reorder(Type, Value, sum, desc=TRUE)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6, angle=65), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm")) +
  scale_fill_manual("Experiment Type",values=m_colors) + 
labs(title="Expression Profiling", 
     subtitle="")
dev.off()


theme_set(theme_classic())

# From on a categorical column variable
g <- ggplot(ep, aes(x = Year, y = Value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Expression Profiling", 
       subtitle="")

# From on a categorical column variable
g <- ggplot(gb, aes(x = Year, y = Value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Genome Binding/Occupancy", 
       subtitle="")

# From on a categorical column variable
g <- ggplot(gv, aes(x = year, y = value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Genome Variation", 
       subtitle="")

# From on a categorical column variable
g <- ggplot(mt, aes(x = Year, y = Value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Methylation", 
       subtitle="")

# From on a categorical column variable
g <- ggplot(nc, aes(x = Year, y = Value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Non-Coding RNA", 
       subtitle="")

# From on a categorical column variable
g <- ggplot(pp, aes(x = Year, y = Value))
g + geom_bar(aes(fill=Type), width = 0.5, stat="identity") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6),
        axis.text.y = element_blank()) +
  scale_fill_brewer(palette="Set3") +
  labs(title="Protein Profiling", 
       subtitle="")

dev.off()




