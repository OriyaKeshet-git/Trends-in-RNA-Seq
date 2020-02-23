library(ggplot2)
library(dplyr)
library(rentrez)
library(ggsci)
#----------------create stacked bar graph for 18 diseases categories---------------------------

types = c("cancer OR prostate cancer OR breast cancer OR lymphoma OR leukemic OR leukemia OR osteosarcoma OR breast cancer OR prostate cancers OR ovarian cancer OR colorectal carcinoma OR lung cancer OR 
                  B-cell lymphoma OR Hodgkin lymphoma OR melanoma OR sarcoma OR lymphoma entities OR carcinogenesis OR Follicular lymphoma OR carcinoma OR hepatocarcinogenesis OR glioblastoma OR glioblastoma tumors OR 
          esophageal cancer OR sarcomas OR liver cancer OR carcinomas OR clear cell renal cell carcinomas OR squamous cell carcinomas OR gastric cancer OR Chronic myeloid leukemia OR colon cancer OR 
          lymphomas OR gliomas OR Buccal mucosal cancer OR neuroblastoma tumour OR cutaneous melanoma metastases OR myeloid leukemiaORglioma OR Horn cancer OR inflammatory breast cancer OR colorectal adenocarcinoma OR 
          acute myeloid leukemia OR glioblastomas OR squamous cell carcinoma OR gastric carcinoma OR lung adenocarcinoma OR hepatocellular carcinoma OR non-small-cell lung cancer (NSCLC) tumors OR colon adenocarcinoma OR 
          adenocarcinoma OR carcinogenic OR pancreatic adenocarcinoma cancer OR Hepatocellular carcinoma (HCC) tumors OR clear-cell renal cell carcinoma OR colon carcinoma OR non-small-cell lung cancer OR 
          melanoma tumors OR tuberculosis OR renal cell carcinoma OR Burkitt lymphoma OR pancreatic ductal adenocarcinoma OR clear cell renal cell carcinoma OR melanoma eDMRs OR non-small cell lung cancer OR 
          Multiple Cancers OR lung squamous cell carcinoma OR malignant osteosarcoma OR osteosarcomas OR Esophageal Squamous Cell Carcinoma OR cancers OR colorectal cancer OR pancreatic carcinogenesis OR 
          Glioblastoma AtlasORgloiblastoma OR neuroblastoma OR head and neck cancers OR neck squamous cell carcinoma OR melanomas OR oral squamous cell carcinomas OR BRCA1 deficient OR tumor OR tumors OR primary tumors OR 
          breast tumor OR ovarian tumors OR tumour OR aggressive tumours OR Metastasis-competent circulating tumour OR 
          SCC OR NSCLC OR CRC OR HCC OR sarcomatoid RCC OR RCC OR GBM [Title/Abstract]",
          
          "NMD OR pleiotropic developmental defects OR polyploidy OR chromosomal aneuploidies OR Paleopolyploidy OR allopolyploidy OR hepatic RDD OR rosai dorfman disease OR rd OR PK OR hyperdiploid/hypotriploid OR 
          autosomal dominant hereditary cataracts OR osteogenesis OR Huntington's disease OR GS OR genetic disorders OR genetic disease OR genetic diseases OR genetic interaction OR genetic abnormalities OR inherited diseases OR 
          NPC OR genetic [Title/Abstract]",
          
          "neuro OR brain OR neurodegenerative diseases OR neurodegenerative brain diseases OR neurodegenerative disease OR Parkinson's disease OR Alzheimer's disease OR Alzheimer's OR AD OR AD temporal lobe OR LB OR epilepsy OR 
          neurological disorders OR neuropsychiatric diseases OR neuropsychiatric disorders OR psychiatric disorders OR necrotic lesions of cerebral infarction OR cerebral infarction OR neurodegeneration OR 
          neuro-developmental disorders OR neurological diseases OR autism, schizophrenia, bipolar disorder OR autism OR schizophrenia and bipolar disorder OR Autism Spectrum Disorders OR schizophrenia OR TICs [Title/Abstract]", 
          
          "bacteria OR bacterial OR Salmonella-infected OR TSSs OR columnaris disease OR bacterial infection OR bacteremia OR bacterial disease infections OR bacteria infection OR parasitic nematodes and rhizobial bacteria OR Tuberculosis Annotation Jamboree OR 
          TSS OR aureus infection [Title/Abstract]",
          
          "viral OR malaria OR malaria parasite OR HCV OR HCV infection OR metabolic impact of HCV infection OR EBV OR Dengue viruses OR Dengue Fever Vector OR viral infections OR viral infection OR virus infection [Title/Abstract]",
          
          "liver OR liver disease OR liver cirrhosis OR liver toxicity [Title/Abstract]",
          
          "blood OR atherosclerosis OR MDS OR Myelodysplastic OR myelodysplastic syndrome OR cotton SE OR SE [Title/Abstract]",
          
          "heart OR cardiac OR cardio OR cardiac hypertrophy OR cardiac septation defects OR myocardial ischemia OR cardiovascular diseases OR myocardial infarction OR cardiac dysfunction OR Congenital heart disease OR PHD [Title/Abstract]", 
          
          "mitochondrial OR mitochondrial dysfunction OR mitochondria dysfunction [Title/Abstract]",
          
          "lung OR raspiratory OR URD OR lung disease OR lung infection OR lung diseases OR aeruginosa causes chronic lung infection OR respiratory syndrome OR pneumococcal OR chronic obstructive pulmonary disease OR pulmonary diseases OR 
          asthma [Title/Abstract]",
          
          "gastro OR gastric OR gastrointestinal disease OR gastrointestinal OR inflammatory bowel disease OR metabolic diseases OR metabolic syndrome OR hepatic metabolism and liver disease OR metabolic disorder OR metabolic disease OR colonic polyp lesion [Title/Abstract]",
          
          "eye OR retiba OR retinal disease OR diabetic retinopathy and retinal inflammation OR retinal leukostasis [Title/Abstract]", 
          
          "Zn deficiency OR N deficiency OR diabetes OR diabetes-induced memory deficits [Title/Abstract]", 
          
          "parasite OR parasitic OR parasitic diseases [Title/Abstract]", 
          
          "depression OR depressive disorder [Title/Abstract]", 
          
          "outoimmune OR immune OR immono OR rheumatoid arthritis OR amyotrophic lateral sclerosis OR multiple sclerosis OR psoriasis OR autoimmune disease OR autoimmune disease systemic lupus erythematosus OR systemic lupus erythematosus OR LPS [Title/Abstract]",
          
          "kidney OR chronic kidney disease[Title/Abstract]", 
          
          "fungal meningitis OR meningitis OR fungal infection OR fungal [Title/Abstract]"
)
myCount = c('16', '33', '73', '174', '293', '488', '679', '954', '1369' )
#Initialize count:
disease_Count = c()
j=1
#Find all RNA related articles for each category per year:
for(y in 2010:2018)
{
  for(i in 1:length(types))
  {
    if (i==1) {
        myT <- gsub('OR', "[Title/Abstract] OR ", types[i])
        query = paste("(RNA-Seq[Title/Abstract]|RNAseq[Title/Abstract]|RNA seq[Title/Abstract]|RNA sequencing[Title/Abstract]|RNAsequencing[Title/Abstract]|RNA-sequencing[Title/Abstract]) AND (", myT, ") AND ", y, "[PDAT]", sep = "")
        print(query)
        count = myCount[j]
        tmp = c(y,myT,count)
        disease_Count = rbind(disease_Count, tmp)
      i= i+1
      j=j+1
    }
    else {
      myT <- gsub('OR', "[Title/Abstract] OR ", types[i])
      query = paste("(RNA-Seq[Title/Abstract]|RNAseq[Title/Abstract]|RNA seq[Title/Abstract]|RNA sequencing[Title/Abstract]|RNAsequencing[Title/Abstract]|RNA-sequencing[Title/Abstract]) AND (", myT, ") AND ", y, "[PDAT]", sep = "")
      count = entrez_search(db="pubmed", term=query, retmax=0)$count
      tmp = c(y,myT,count)
      disease_Count = rbind(disease_Count, tmp)
    }
  }
  print(disease_Count)
}

#disease_Count = rbind(disease_Count, tmp)

#Create tags for the categories:
types = c(
  "Cancer", 
  "Genetic Diseases", 
  "Brain diseases", 
  "Bacterial Diseases", 
  "Viral Diseases", 
  "Liver diseases",
  "Hematologic Diseases",
  "Cardiac Diseases", 
  "Mitocondrial Diseases",
  "Lungs Diseases", 
  "Gastric Diseases", 
  "Eye Diseases", 
  "Deficiencies", 
  "Parasitic Diseases", 
  "Mental Diseases", 
  "Autoimmune Diseases", 
  "Kidney Diseases", 
  "Fungal Diseases")

#Create dataframe for the graph:
colnames(disease_Count) <- c("Year", "Diseases", "Value")
disease_Count[,2] <- types
df <- as.data.frame(disease_Count)
df[,3] <- as.numeric(as.character(df[,3]))
df[,1] <- as.numeric(as.character(df[,1]))

m_colors = c("#f2a8bb",
             "#ffd1c4",
             "#e7ae8e",
             "#c5a68e",
             "#e1cd96",
             "#fff1c7",
             "#9fb28a",
             "#e6ffe7",
             "#a7e8c5",
             "#88b5a1",
             "#82d7c5",
             "#64c8d3",
             "#77b3d1",
             "#acdeff",
             "#f7faff",
             "#e6dfff",
             "#b8aae0",
             "#c3a1bd")

pdf("Diseases Stacked Graph.pdf", width = 15, height = 10)

theme_set(theme_classic())
g <- ggplot(df, aes(x = Year, y = Value, fill = Diseases))
g + geom_bar(alpha = .9 , size=0.01, colour="grey", aes(fill=forcats::fct_reorder(Diseases, Value)), width = 0.5, stat="identity") + 
  theme(axis.title = element_text(face = "italic", size = 27),
        plot.title = element_text( face = "bold", size = 30),
        axis.text.x = element_text( vjust=0.6), 
        axis.text =   element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))+
  scale_fill_manual(values=m_colors) +
  labs(
       subtitle="",
       y = "Number Of Articles")

dev.off()
