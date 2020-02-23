# This script compare between cancer precentges in population to cancer in RNASEQ reaserches 
# and plot grouped bar plot for the results
library(rentrez)
library(assertr)
library(ggplot2)
library(ggthemes)


# The function calcultes the number of abstracts related to cancer
all_cancer_abs_num <- function()
{
  table = read.csv("C:/Users/Racheli/Documents/mini-project-bioinformatics/NEW filterd cancers.csv", header = TRUE, check.names = F)
  cancer_names = as.character(table[,1])
  count_all_abs = 0  
  count_rnaseq_abs = 0
  y = 2018
  rnaseq_str = "(RNA seq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract] OR rnaseq[Title/Abstract]) "
  # for each cancer type, calculte the number of realted abstracts,
  # for all articels and for rnaseq articles
  for(i in 1:length(cancer_names))
    {
      query = paste(rnaseq_str, "AND (", cancer_names[i],") AND ", y, "[PDAT]", sep = "")
      count_rnaseq_abs = count_rnaseq_abs + entrez_search(db="pubmed", term=query, retmax=0)$count
      query = paste("(", cancer_names[i],") AND ", y, "[PDAT]", sep = "")
      count_all_abs = count_all_abs + entrez_search(db="pubmed", term=query, retmax=0)$count
    }
  
  return(cbind(count_all_abs, count_rnaseq_abs))
}

# read the table with the data about cancer in popoulation
table = read.csv("C:/Users/Racheli/Documents/mini-project-bioinformatics/cancersPopulation_2018.csv", header = F)
cancer_types = gsub(" OR ", "[Title/Abstract] | ",table[,1])

# The numbers are from the table(psted it because of technical problems)
general_population_count = c(2093876,2088849,1800977,1276106,1033701,841080, 572034,569847,
  567233,549393,509590,458918,437033,403262,382069,354864,296851,295414,287723,219420,
  177422,159985,129079,92887,80608,79990,71105,52799,48541,44235,41799,34475,30443,17600)
general_population_prec = table[,3]

# calculte the number of abstracts related to all cancer types
counts = all_cancer_abs_num()
all_cancer_abs_num = counts[,1]
all_cancer_rnaseq_abs_num = counts[,2]

# for each cancer type in the table from population, calculte the number of related abstracts
cancer_abs_count = c()
for(j in 1:length(cancer_types))
{
  query = paste("(",cancer_types[j], "[Title/Abstract])", "AND 2018 [PDAT]",sep = " ") 
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  cancer_abs_count = c(cancer_abs_count, count)
}

# for each cancer type in the table from population, calculte the number of related rnaseq abstracts

rnaseq_str = "(RNA seq[Title/Abstract] OR RNA-seq[Title/Abstract] OR rna sequencing[Title/Abstract] OR rnaseq[Title/Abstract]) "
cancer_rnaseq_abs_count = c()
for(i in 1:length(cancer_types))
{
  query = paste("(",cancer_types[i], "[Title/Abstract])", "AND 2018 [PDAT] AND ",rnaseq_str ,sep = " ")
  count = entrez_search(db="pubmed", term=query, retmax=0)$count
  cancer_rnaseq_abs_count = c(cancer_rnaseq_abs_count, count)
}

#chi squre for population sickness in 2018
chi_squre_table = rbind(cancer_rnaseq_abs_count, cancer_abs_count)
chi_squre_table = chi_squre_table[,chi_squre_table[2,] != 0]
result = chisq.test(chi_squre_table, simulate.p.value = T)

#chi squre test for abstracts number
chi_squre_table = rbind(cancer_rnaseq_abs_count, general_population_count)
chi_squre_table = chi_squre_table[,chi_squre_table[2,] != 0]
result = chisq.test(chi_squre_table, simulate.p.value = T)

#chi squre between number of abstracts to popoulation sickness
chi_squre_table = rbind(cancer_abs_count, general_population_count)
chi_squre_table = chi_squre_table[,chi_squre_table[2,] != 0]
result = chisq.test(chi_squre_table, simulate.p.value = T)

#---------------------------------------------------------------------------------------

# calculte the precentges of cancer related abstracts
cancer_abs_prec = round((cancer_abs_count / all_cancer_abs_num * 100), digits = 3)
cancer_abs_rnaseq_prec = round((cancer_rnaseq_abs_count / all_cancer_rnaseq_abs_num * 100), digits = 3)

# read cancer types for plot
cancer_types_for_plot = read.csv("~/mini-project-bioinformatics/cancer_types_for_plot.csv", header = T)
cancer_types_for_plot = as.character(cancer_types_for_plot[,1])

# create data for plot 1, for cancers more reserched from their precentges in population
data_for_plot1 = c()
# iterate all cancer types
for(i in 1:length(cancer_types))
{
  if(cancer_abs_rnaseq_prec[i] > general_population_prec[i])
  {
    tmp = c(cancer_types_for_plot[i], 'RNASEQ', cancer_abs_rnaseq_prec[i])
    tmp1 = c(cancer_types_for_plot[i], 'Population', general_population_prec[i])
    data_for_plot1 = rbind(data_for_plot1, tmp)
    data_for_plot1 = rbind(data_for_plot1, tmp1)
  }
  
}
# create data frame
m_df = as.data.frame(data_for_plot1)
m_df[,3] <- as.numeric(as.character(m_df[,3]))
names(m_df) = c("cancer", "where", "prec")

# create data for plot 2, for cancers with higher precentges in population compared to reaserch
data_for_plot2 = c()
for(i in 1:length(cancer_types))
{
  if(cancer_abs_rnaseq_prec[i] < general_population_prec[i])
  {
    tmp = c(cancer_types_for_plot[i], 'RNASEQ', cancer_abs_rnaseq_prec[i])
    tmp1 = c(cancer_types_for_plot[i], 'Population', general_population_prec[i])
    data_for_plot2 = rbind(data_for_plot2, tmp)
    data_for_plot2 = rbind(data_for_plot2, tmp1)
  }
  
}
# create data frame
m_df1 = as.data.frame(data_for_plot2)
m_df1[,3] <- as.numeric(as.character(m_df1[,3]))
names(m_df1) = c("cancer", "where", "prec")

# -------------------------plot grouped bar chart for the results---------------------------

pdf("cancerVspopulation.pdf", width = 13, height = 10)
ggplot(m_df, aes(fill=where, y=prec, x=cancer)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.title = element_text(family = "Times", face = "italic", size = 25),
        plot.title = element_text(family = "Times", face = "bold", size = 30),
        axis.text.x = element_text(angle=65, vjust=0.6), 
        plot.subtitle =  element_text(family = "Times", size = 27),
        axis.text =   element_text(family = "Times",size = 20),
        legend.text = element_text(family = "Times",size = 22),
        legend.title = element_text(family = "Times",size = 25),
        legend.key.size = unit(0.95, "cm")) +
  labs(title="", 
       subtitle="Occurrence of cancer in population and researches related to RNA-Seq", fill = "", y = "precentges", x = "cancer",
       family = "Times", face = "italic") 
dev.off()


pdf("cancerVspopulation1.pdf", width = 14, height = 10)
ggplot(m_df1, aes(fill=where, y=prec, x=cancer)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.title = element_text(family = "Times", face = "italic", size = 25),
        plot.title = element_text(family = "Times", face = "bold", size = 30),
        plot.subtitle = element_text(family = "Times", size = 27),
        axis.text.x = element_text(angle=90, vjust=0.6), 
        axis.text =   element_text(family = "Times",size = 20),
        legend.text = element_text(family = "Times",size = 22),
        legend.title = element_text(family = "Times",size = 25),
        legend.key.size = unit(0.95, "cm")) +
  labs(title="", 
       subtitle="Occurrence of cancer in population and researches related to RNA-Seq", fill = "", y = "precentges", x = "cancer",
       family = "Times", face = "italic") 
dev.off()
