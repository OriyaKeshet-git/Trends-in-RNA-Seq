library(ggplot2)
library(dplyr)
library(rentrez)
#----------------create Scatterplot graph for six main diseases categories with Linear Regression---------------------------

types = c("cancer|prostate cancer|breast cancer|lymphoma|leukemic|leukemia|osteosarcoma|breast cancer|prostate cancers|ovarian cancer|colorectal carcinoma|lung cancer|
          B-cell lymphoma|Hodgkin lymphoma|melanoma|sarcoma|lymphoma entities|carcinogenesis|Follicular lymphoma|carcinoma|hepatocarcinogenesis|glioblastoma|glioblastoma tumors|
          esophageal cancer|sarcomas|liver cancer|carcinomas|clear cell renal cell carcinomas|squamous cell carcinomas|gastric cancer|Chronic myeloid leukemia|colon cancer|
          lymphomas|gliomas|Buccal mucosal cancer|neuroblastoma tumour|cutaneous melanoma metastases|myeloid leukemia|glioma|Horn cancer|inflammatory breast cancer|colorectal adenocarcinoma|
          acute myeloid leukemia|glioblastomas|squamous cell carcinoma|gastric carcinoma|lung adenocarcinoma|hepatocellular carcinoma|non-small-cell lung cancer (NSCLC) tumors|colon adenocarcinoma|
          adenocarcinoma|carcinogenic|pancreatic adenocarcinoma cancer|Hepatocellular carcinoma (HCC) tumors|clear-cell renal cell carcinoma|colon carcinoma|non-small-cell lung cancer|
          melanoma tumors|tuberculosis|renal cell carcinoma|Burkitt lymphoma|pancreatic ductal adenocarcinoma|clear cell renal cell carcinoma|melanoma eDMRs|non-small cell lung cancer|
          Multiple Cancers|lung squamous cell carcinoma|malignant osteosarcoma|osteosarcomas|Esophageal Squamous Cell Carcinoma|cancers|colorectal cancer|pancreatic carcinogenesis|
          Glioblastoma Atlas|gloiblastoma|neuroblastoma|head and neck cancers|neck squamous cell carcinoma|melanomas|oral squamous cell carcinomas|BRCA1 deficient|tumor|tumors|primary tumors|
          breast tumor|ovarian tumors|tumour|aggressive tumours|Metastasis-competent circulating tumour|
          SCC|NSCLC|CRC|HCC|sarcomatoid RCC|RCC|GBM",
          
          "NMD|pleiotropic developmental defects|polyploidy|chromosomal aneuploidies|Paleopolyploidy|allopolyploidy|hepatic RDD|rosai dorfman disease|rd|PK|hyperdiploid/hypotriploid|
          autosomal dominant hereditary cataracts|osteogenesis|Huntington's disease|GS|genetic disorders|genetic disease|genetic diseases|genetic interaction|genetic abnormalities|inherited diseases|
          NPC|genetic",
          
          "neuro|brain|neurodegenerative diseases|neurodegenerative brain diseases|neurodegenerative disease|Parkinson's disease|Alzheimer's disease|Alzheimer's|AD|AD temporal lobe|LB|epilepsy|
          neurological disorders|neuropsychiatric diseases|neuropsychiatric disorders|psychiatric disorders|necrotic lesions of cerebral infarction|cerebral infarction|neurodegeneration|
          neuro-developmental disorders|neurological diseases|autism, schizophrenia, bipolar disorder|autism|schizophrenia and bipolar disorder|Autism Spectrum Disorders|schizophrenia|TICs", 
          
          "bacteria|bacterial|Salmonella-infected|TSSs|columnaris disease|bacterial infection|bacteremia|bacterial disease infections|bacteria infection|parasitic nematodes and rhizobial bacteria|Tuberculosis Annotation Jamboree|
          TSS|aureus infection",
          
          "blood|atherosclerosis|MDS|Myelodysplastic|myelodysplastic syndrome|cotton SE|SE",
          
          "outoimmune|immune|immono|rheumatoid arthritis|amyotrophic lateral sclerosis|multiple sclerosis|psoriasis|autoimmune disease|autoimmune disease systemic lupus erythematosus|systemic lupus erythematosus|LPS"
          
          
)

#Initialize count:
diseaseCount = c()

#Find all articles with RNA related term in their Title/Abstract:
for(y in 2010:2018)
{
  
  for(i in 1:length(types))
  {
    query = paste("(RNA-Seq[Title/Abstract]|RNAseq[Title/Abstract]|RNA seq[Title/Abstract]|RNA sequencing[Title/Abstract]|RNAsequencing[Title/Abstract]|RNA-sequencing[Title/Abstract]) AND (", types[i], ") AND ", y, "[PDAT]", sep = "")
    count = entrez_search(db="pubmed", term=query, retmax=0)$count
    tmp = c(y,types[i],count)
    diseaseCount = rbind(diseaseCount, tmp)
    
  }
  print(diseaseCount)
}

#Create tags for each category:
types = c(
  "Cancer", 
  "Genetic Diseases", 
  "Brain diseases", 
  "Bacterial Diseases", 
  "Hematologic Diseases",
  "Autoimmune Diseases")

#create a dataframe for the plot:
colnames(diseaseCount) <- c("Year", "Diseases", "Value")
diseaseCount[,2] <- types
df <- as.data.frame(diseaseCount)
df[,3] <- as.numeric(as.character(df[,3]))
df[,1] <- as.numeric(as.character(df[,1]))
a = c(0.8898, 0.9853, 0.9625, 0.9513, 0.9779, 0.9358)



pdf("Diseases Scatterplot.pdf", width = 15, height = 10)

theme_set(theme_classic())
ggplot(df, aes(x=df[,1], y=df[,3], color=Diseases)) + labs(x = "Years", y = "Number Of Articles", fill = NULL, color = "Diseases" )+ #, title = "Linear Regression per Diseases")+
geom_point(size=6, alpha=0.5)+ geom_smooth(method=lm,  se=FALSE)+
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.key.size = unit(0.9, "cm"))
dev.off()
