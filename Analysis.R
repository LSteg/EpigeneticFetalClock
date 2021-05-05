#Code to paper "Novel epigenetic clock for fetal brain development predicts prenatal developmental age for cellular stem cell models and derived neurons."
author: "Leonard C Steg"
date: "05.05.2021"

setwd("/mnt/data1/LeoSteg/FetalClock")

library(tidyverse)
library(wateRmelon)  
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(xlsx)
library(readxl)
library(ggpubr)
library(rstatix)
library(reshape)
library(caTools)

########################################################################
##### Setup of Training & Testing datasets                         #####
########################################################################

##### Input: 3 QCed and Normalized datasets (1. Helen Spiers data (our lab, previously published but QCed once again with more strigent procedure)/ 2. data from Jaffe Group (Baltimore) / 3. Novel fetal brain samples from our lab)

### Big Brain Dataset - Helen Spiers fetal samples

load("/mnt/data1/BrainData/NormalisedData.Rdata")

pheno.spiers <- full.pheno[full.pheno$Age < 0,]
range(pheno.spiers$Age)
# -0.7041096 -0.2630137

table(pheno.spiers$Sex)
#F  M 
#68 86 

rownames(pheno.spiers) <- pheno.spiers$Sentrix_Full
betas.spiers <- betas.dasen[,rownames(pheno.spiers)]
dim(betas.spiers)
#419671    154

pheno.spiers[,"dpc"] <- pheno.spiers$Age*365+280

### Baltimore / Jaffe Dataset

load("/mnt/data1/EmmaW/QC_datasets/Leiber_Jaffe_GSE74193/QC/GSE74193_Normalised.rdat")
pheno.jaffe <- pheno[pheno$age..in.years..ch1 < 0,]
range(pheno.jaffe$age..in.years..ch1)
# -0.498630 -0.383561

table(pheno.jaffe$sex..clinical.gender..ch1)
#F  M 
#35 29  

pheno.jaffe[,"dpc"] <- pheno.jaffe$age..in.years..ch1*365+280

betas.jaffe <- betas[,rownames(pheno.jaffe)]
dim(betas.jaffe)
#446570     64

### Novel Fetal Brain Samples Dataset

load("/mnt/data1/Array_Projects/Bit_plate/Nicks_Samples/Nick_Samples_Normalised.rdat")
Nick_Brays_Sample_sheet <- read_csv("Nick_Brays_Sample_sheet.csv")

Nick_Brays_Sample_sheet <- Nick_Brays_Sample_sheet[!is.na(Nick_Brays_Sample_sheet$Age..PCW.),]
Nick_Brays_Sample_sheet <- Nick_Brays_Sample_sheet[Nick_Brays_Sample_sheet$Sample_ID%in%SampleSheet$Sample_ID,]
Nick_Brays_Sample_sheet <- Nick_Brays_Sample_sheet[-26,]

red_nick <- Nick_Brays_Sample_sheet[c("Sample_ID", "Age..PCW.", "Sex")]
colnames(red_nick) <- c("Sample_ID", "PCW", "Sex1")

pheno.bray <- left_join(SampleSheet, red_nick, by = "Sample_ID")
pheno.bray$Age..PCW. <- pheno.bray$PCW
pheno.bray$PCW <- NULL
pheno.bray$Sex <- pheno.bray$Sex1
pheno.bray$Sex1 <- NULL
rownames(pheno.bray) <- pheno.bray$Basename

betas.bray <- betas[,rownames(pheno.bray)]
dim(betas.bray)
#793459     40

#Adding Cell proportion estimates
setwd("/mnt/data1/reference_files/CETs")
library(cets)
load("cetsBrain.rda")
load("cetsDilution.rda")

modelIdx <- list(neuron = pdBrain$celltype == "N", glia = pdBrain$celltype == "G")
refProfile <- getReference(brain, modelIdx)
head(refProfile)
pheno.bray$prop <- estProportion(betas.bray, profile = refProfile)
setwd("/mnt/data1/LeoSteg/FetalClock")

pheno.bray$Age..PCW. <- as.numeric(as.character(pheno.bray$Age..PCW.))
pheno.bray$dpc <- pheno.bray$Age..PCW.*7
pheno.bray$Age <- (pheno.bray$dpc-280)/365

range(pheno.bray$Age)
#-0.5369863 -0.4027397
table(pheno.bray$Sex)
#F  M 
#14 26 

### Finding overlapping CpGs of all 3 datasets and EPIC manifest

epicManifest<-read.csv("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7)
x <- intersect(intersect(rownames(betas.spiers), intersect(rownames(betas.jaffe), epicManifest$IlmnID)),rownames(betas.bray))
length(x)
# 385069

betas.spiers <- betas.spiers[x,]
dim(betas.spiers)
# 385069    154
betas.jaffe <- betas.jaffe[x,]
dim(betas.jaffe)
# 385069     64
betas.bray <- betas.bray[x,]
dim(betas.bray)
# 385069     40

### Creating common pheno files with all needed data. Including different epigenetic clocks (MTC, GAC, CPC)

colnamespheno <- c("SampleID", "Sentrix_Full", "Sex", "Age", "dpc", "BrainRegion", "prop")

pheno.spiers <- pheno.spiers[,c(1,14,6,5,25,2,24)]
colnames(pheno.spiers) <- colnamespheno
rownames(pheno.spiers) <- pheno.spiers[,2]

pheno.jaffe <- pheno.jaffe[,c(59,75,74,57,83,9,63)]
colnames(pheno.jaffe) <- colnamespheno
rownames(pheno.jaffe) <- pheno.jaffe[,2]

pheno.bray <- pheno.bray[,c(1,25,9,10,39, 14,38)]
colnames(pheno.bray) <- colnamespheno
rownames(pheno.bray) <- pheno.bray[,2]

pheno.spiers$Cohort <- rep("BigBrain", nrow(pheno.spiers))
pheno.jaffe$Cohort <- rep("Baltimore", nrow(pheno.jaffe))
pheno.bray$Cohort <- rep("bray", nrow(pheno.bray))

# MTC
pheno.spiers$Horvath <- agep(betas.spiers)
pheno.jaffe$Horvath <- agep(betas.jaffe)
pheno.bray$Horvath <- agep(betas.bray)

pheno.spiers$Horvathdpc <- pheno.spiers$Horvath*365+280
pheno.jaffe$Horvathdpc <- pheno.jaffe$Horvath*365+280
pheno.bray$Horvathdpc <- pheno.bray$Horvath*365+280

# GAC
datClock <- read.csv("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Gest_age/cgprobesGApredictor.csv",as.is=T)
GACoef <- datClock$CoefficientTraining
names(GACoef) <- datClock$CpGmarker

pheno.spiers$GA <- agep(betas.spiers, coeff = GACoef, method = "hannum")+GACoef[1]
pheno.spiers$GA <- (pheno.spiers$GA-40)/52
pheno.spiers$GAdpc <- pheno.spiers$GA*365+280

pheno.jaffe$GA <- agep(betas.jaffe, coeff = GACoef, method = "hannum")+GACoef[1]
pheno.jaffe$GA <- (pheno.jaffe$GA-40)/52
pheno.jaffe$GAdpc <- pheno.jaffe$GA*365+280

pheno.bray$GA <- agep(betas.bray, coeff = GACoef, method = "hannum")+GACoef[1]
pheno.bray$GA <- (pheno.bray$GA-40)/52
pheno.bray$GAdpc <- pheno.bray$GA*365+280

# CPC
placCpGs<-read.csv("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Plac_age/Placental_coeff.csv", header = T, stringsAsFactors = F)
placCoeff <- placCpGs[, 3]
names(placCoeff)<-placCpGs[,1]

pheno.spiers$Placdpc<-(agep(betas.spiers, coef=placCoeff, method="hannum") + 13.062)*7
pheno.spiers$PlacAge <- (pheno.spiers$Placdpc-280)/365

pheno.jaffe$Placdpc<-(agep(betas.jaffe, coef=placCoeff, method="hannum") + 13.062)*7
pheno.jaffe$PlacAge <- (pheno.jaffe$Placdpc-280)/365

pheno.bray$Placdpc<-(agep(betas.bray, coef=placCoeff, method="hannum") + 13.062)*7
pheno.bray$PlacAge <- (pheno.bray$Placdpc-280)/365

### Split datasets into 75% / 25% portions

set.seed(202)
sample.train <- sample.split(pheno.spiers$Sentrix_Full, SplitRatio = .75)
pheno.spiers.75 <- subset(pheno.spiers, sample.train == TRUE)
pheno.spiers.25 <- subset(pheno.spiers, sample.train == FALSE)
betas.spiers.75 <- betas.spiers[,rownames(pheno.spiers.75)]
betas.spiers.25 <- betas.spiers[,rownames(pheno.spiers.25)]

set.seed(203)
sample.test <- sample.split(pheno.jaffe$Sentrix_Full, SplitRatio = .75)
pheno.jaffe.75 <- subset(pheno.jaffe, sample.test == TRUE)
pheno.jaffe.25 <- subset(pheno.jaffe, sample.test == FALSE)
betas.jaffe.75 <- betas.jaffe[,rownames(pheno.jaffe.75)]
betas.jaffe.25 <- betas.jaffe[,rownames(pheno.jaffe.25)]

set.seed(204)
sample.bray <- sample.split(pheno.bray$Sentrix_Full, SplitRatio = .75)
pheno.bray.75 <- subset(pheno.bray, sample.bray == TRUE)
pheno.bray.25 <- subset(pheno.bray, sample.bray == FALSE)
betas.bray.75 <- betas.bray[,rownames(pheno.bray.75)]
betas.bray.25 <- betas.bray[,rownames(pheno.bray.25)]

pheno.train <- rbind(pheno.spiers.75, pheno.jaffe.75, pheno.bray.75)
betas.train <- cbind(betas.spiers.75, betas.jaffe.75, betas.bray.75)
identical(rownames(pheno.train), colnames(betas.train))
#TRUE

pheno.test <- rbind(pheno.spiers.25, pheno.jaffe.25, pheno.bray.25)
betas.test <- cbind(betas.spiers.25, betas.jaffe.25, betas.bray.25)
identical(rownames(pheno.test), colnames(betas.test))
#TRUE

save(pheno.train, betas.train, file = "data.train.Rdata")
save(pheno.test, betas.test, file = "data.test.Rdata")

########################################################################
##### Setup of validation dataset                                  #####
########################################################################

##### New fetal brain data from the lab 

load("/mnt/data1/Array_Projects/Schizophrenia_fetal/FetalBrainData/FetalBrainQC/March2021/Normalised_FetalBrain_v4.rdat")
pheno.new <- SampleSheet
betas.new <- betas

### Cleanup
#Remove cerebellum samples
pheno.new <- pheno.new[pheno.new$Tissue_Type == "Prefrontal Cortex",]
#Create new column with age in pcw as numeric
pheno.new$pcw <- as.numeric(str_replace(pheno.new$Age, "pcw", ""))
#Remove samples from Institute "Oxford", as 3 of 4 samples have ages after pregnancy (up to pcw 456) 
pheno.new <- pheno.new[pheno.new$Institute != "Oxford",]

#Remove all other samples with PCW older than 40
pheno.new <- pheno.new[pheno.new$pcw <= 40,]

pheno.new$dpc <- pheno.new$pcw *7

betas.new <- betas.new[,rownames(pheno.new)]
dim(betas.new)
#807898     96

identical(rownames(pheno.new), colnames(betas.new))
# TRUE

### Run epigenetic clocks

##### MTC
pheno.new$Horvath <- agep(betas.new, method = "horvath")*365+280

##### GAC
GestCpGs<- read.csv ("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Gest_age/cgprobesGApredictor.csv",stringsAsFactor=F, header=T)
GestCoeff<-GestCpGs[,2]
names(GestCoeff)<-GestCpGs[,1]
pheno.new$GestPred<-((agep(betas.new, coef=GestCoeff, method="hannum")+41.726)*7)

##### CPC
placCpGs<-read.csv("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Plac_age/Placental_coeff.csv", header = T, stringsAsFactors = F)
placCoeff <- placCpGs[, 3]
names(placCoeff)<-placCpGs[,1]
pheno.new$PlacPred<-((agep(betas.new, coef=placCoeff, method="hannum") + 13.062)*7)

##### FBC
source("/mnt/data1/LeoSteg/FetalClock/Scripts/FetalClockFunction.R")
dir <- "/mnt/data1/LeoSteg/FetalClock/"
pheno.new$FetalPred <- FetalClock(betas.new,  dir = dir)
# No missing probes

### Load in other validation datasets

load("/mnt/data1/LeoSteg/FetalClock/data.colunga.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/data.chatterton.Rdata")

# Merge and set up dataframe

pheno.val <- as.data.frame(matrix(NA, nrow = sum(nrow(pheno.new), nrow(pheno.chatterton), nrow(pheno.colunga)), ncol = 9))

colnames(pheno.val) <- c("Basename/Accession", "Dataset", "Sex", "Tissue", "dpc", "Horvath","GestPred", "PlacPred", "FetalPred")

pheno.val$`Basename/Accession` <- c(as.character(pheno.new$Basename), pheno.chatterton$Geo_Accession, pheno.colunga$Geo_Accession)
pheno.val$Dataset <- c(rep("Exeter Fetal 2", nrow(pheno.new)), rep("chatterton", nrow(pheno.chatterton)), rep("colunga", nrow(pheno.colunga)))
pheno.val$Sex <- c(as.character(pheno.new$Sex), pheno.chatterton$Sex, pheno.colunga$Sex)
pheno.val$Tissue <- c(as.character(pheno.new$Tissue_Type), pheno.chatterton$Tissue, pheno.colunga$Tissue)
pheno.val$dpc <- c(pheno.new$dpc, pheno.chatterton$dpc, pheno.colunga$dpc)
pheno.val$Horvath <- c(pheno.new$Horvath, pheno.chatterton$Horvath, pheno.colunga$Horvath)
pheno.val$GestPred <- c(pheno.new$GestPred, pheno.chatterton$GestPred, pheno.colunga$GestPred)
pheno.val$PlacPred <- c(pheno.new$PlacPred, pheno.chatterton$PlacPred, pheno.colunga$PlacPred)
pheno.val$FetalPred <- c(pheno.new$FetalPred, pheno.chatterton$FetalPred, pheno.colunga$FetalPred)

rownames(pheno.val) <- pheno.val$`Basename/Accession`

save(pheno.val,file = "data.val.Rdata")

########################################################################
##### Exploration of datasets (Histogram = Figure S1)              #####
########################################################################

load("/mnt/data1/LeoSteg/FetalClock/data.test.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/data.train.Rdata")

table(pheno.train$Sex)
#F   M 
#83 110 
table(pheno.test$Sex)
#F  M 
#34 31 
table(pheno.val$Sex)
#F  M
#61 58

median(pheno.train$dpc)
#[1]  99
median(pheno.test$dpc)
#[1]  99
median(pheno.val$dpc)
#[1]  112

range(pheno.train$dpc)
#[1]  37 184
range(pheno.test$dpc)
#[1]  23 153
range(pheno.val$dpc)
#[1]  42 280

histtrain <- ggplot(data = pheno.train, aes(x = dpc))+
  geom_histogram(bins = 20, color = "black", fill = "grey") +
  theme_cowplot(16) + 
  ylim(0,25) +
  annotate("text",  x = min(pheno.train$dpc)-10, y = 25, hjust = 0,vjust = 1 , label = paste("n = ", nrow(pheno.train)), 
           fontface = "bold", size = 4) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

histtest <- ggplot(data = pheno.test, aes(x = dpc))+
  geom_histogram(bins = 20, color = "black", fill = "grey") +
  theme_cowplot(16) + 
  ylim(0,25) +
  annotate("text",  x = min(pheno.test$dpc)-10, y = 25, hjust = 0,vjust = 1 , label = paste("n = ", nrow(pheno.test)), 
           fontface = "bold", size = 4) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

histval <- ggplot(data = pheno.val, aes(x = dpc))+
  geom_histogram(bins = 20, color = "black", fill = "grey") +
  theme_cowplot(16) + 
  ylim(0,25) +
  annotate("text",  x = min(pheno.val$dpc)-10, y = 25, hjust = 0,vjust = 1 , label = paste("n = ", nrow(pheno.val)), 
           fontface = "bold", size = 4) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

plot <- plot_grid(histtrain, histtest, histval, align = "vh" , labels = c("A Training Data", "B Testing Data", "C Validation Data"), label_x = -0.15, nrow = 1)

x.grob1 <- textGrob("Chronological Age", 
                    gp=gpar(fontface="bold", col="black", fontsize=15))
x.grob2 <- textGrob("in days post-conception", 
                    gp=gpar(fontface = "italic", col="black", fontsize=13))
y.grob1 <- textGrob("Number of samples", 
                    gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)

#add to plot
pdf("figures/histogram.all.pdf")
grid.arrange(arrangeGrob(y.grob1,plot, nrow = 1, widths = c(0.04, 0.92,0.04)), nrow = 4, x.grob1, x.grob2, heights = c(0.93,0.04,0.03, 0.65))
dev.off() 

########################################################################
##### Training of the new fetal Clock                              #####
########################################################################

library(glmnet)

### Read in data
load("/mnt/data1/LeoSteg/FetalClock/data.test.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/data.train.Rdata")

##### Training the Clock without transformation
### Run Elastic Net Regression
# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
alpha<-0.5
glmnet.Training.CV = cv.glmnet(t(betas.train), pheno.train$dpc,  nfolds=10,alpha=alpha,family="gaussian") ## CV = cross-validation
# The definition of the lambda parameter:
lambda.glmnet.Training = glmnet.Training.CV$lambda.min
lambda.glmnet.Training
# 3.272082
# Fit the elastic net predictor to the training data
glmnet.Training = glmnet(t(betas.train),pheno.train$dpc, family="gaussian", alpha=0.5, nlambda=100)
# Arrive at an estimate of of DNAmAge
DNAmdpcBrainTraining=predict(glmnet.Training,t(betas.train),type="response",s=lambda.glmnet.Training)
DNAmdpcBrainTesting=predict(glmnet.Training,t(betas.test),type="response",s=lambda.glmnet.Training)
# Save files
write.csv(DNAmdpcBrainTraining,"training_dpc_prediction.csv", row.names=F)
write.csv(DNAmdpcBrainTesting,"testing_dpc_prediction.csv", row.names=F)

### Extract Coefs
tmp_coeffs <- coef(glmnet.Training.CV, s = "lambda.min")
myCoeff<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
write.table(myCoeff,"FetalClock_DPC_coefficients.txt", row.names=F, col.names=T, quote=F)

#### Renaming file of Coefficients to fetalclock_coefficients.txt plus setting up a script with the function for sourcing into other scripts
#Function name = fetalclock(betas, ageinyears = FALSE, dir)

### Create reference file
ref <- apply(betas.train, 1, mean)
save(ref, file = "fetalclock_reference.rdat")

source("/mnt/data1/LeoSteg/FetalClock/Scripts/FetalClockFunction.R")
dir <- "/mnt/data1/LeoSteg/FetalClock/"

########################################################################
##### Applying fetal brain clock to testing data (Figure 1)        #####
########################################################################

compareStats<-function(data){
  colnames(data)<-c("ID","age.raw","Horvath", "GA","PlacClock", "FetalClock")
  data<-data[complete.cases(data$age.raw),]
  data<-data[complete.cases(data$Horvath),]
  corv<-c(round(cor(data[,2],data[,3]),2),round(cor(data[,2],data[,4]),2),round(cor(data[,2],data[,5]),2),round(cor(data[,2],data[,6]),2))
  rmse<-c(round(sqrt(mean((data[,2]-data[,3])^2)),2),round(sqrt(mean((data[,2]-data[,4])^2)),2),round(sqrt(mean((data[,2]-data[,5])^2)),2),round(sqrt(mean((data[,2]-data[,6])^2)),2))
  stats<-matrix(ncol=4,nrow=2)
  colnames(stats)<-c("Horvath's Clock","Gestational Clock","PlacClock", "FetalClock")
  rownames(stats)<-c("Correlation","RMSE")
  stats[1,]<-corv
  stats[2,]<-rmse
  print(stats)
}

pheno.test$FetalPred <- fetalclock(betas.test, dir = dir)

testing <- pheno.test[,c(1,5,10,12,13,15)]
colnames(testing) <- c("sample", "Age", "Horvath","Gestational", "Placental", "Fetal")

summarytest<-compareStats(testing)
write.xlsx(summarytest, file = "summarytest.xlsx")

ggFetaltest <-ggplot(data = pheno.test,aes(x=dpc,y=FetalPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.test$dpc), y = max(pheno.test$FetalPred), hjust = 0,vjust = 1 , label = paste("r = ", summarytest[1,4], "\nRMSE = ", summarytest[2,4]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

ggHorvathtest <-ggplot(data = pheno.test,aes(x=dpc,y=Horvathdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.test$dpc), y = max(pheno.test$Horvathdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytest[1,1], "\nRMSE = ", summarytest[2,1]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggGAtest <-ggplot(data = pheno.test,aes(x=dpc,y=GAdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.test$dpc), y = max(pheno.test$GAdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytest[1,2], "\nRMSE = ", summarytest[2,2]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggPlactest <-ggplot(data = pheno.test,aes(x=dpc,y=Placdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.test$dpc), y = max(pheno.test$Placdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytest[1,3], "\nRMSE = ", summarytest[2,3]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
                
plot <- plot_grid(ggFetaltest,ggHorvathtest,ggGAtest,ggPlactest, align = "vh" , labels = c("A FBC", "B MTC", "C GAC", "D CPC"))
    
x.grob1 <- textGrob("Chronological Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
x.grob2 <- textGrob("in days post-conception", 
                   gp=gpar(fontface = "italic", col="black", fontsize=13))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
                     
#add to plot
pdf("figures/plotgrid.testdata.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)), nrow = 3, x.grob1, x.grob2, heights = c(0.93,0.04,0.03))
dev.off()

########################################################################
##### Applying fetal brain clock to training data (Figure S2)      #####
########################################################################

pheno.train$FetalPred <- fetalclock(betas.train, dir = dir)

training <- pheno.train[,c(1,5,10,12,13,15)]
colnames(training) <- c("sample", "Age", "Horvath","Gestational", "Placental", "Fetal")

summarytrain<-compareStats(training)
write.xlsx(summarytrain, file = "summarytrain.xlsx")

ggFetaltrain <-ggplot(data = pheno.train,aes(x=dpc,y=FetalPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.train$dpc)-10, y = max(pheno.train$FetalPred), hjust = 0,vjust = 1 , label = paste("r = ", summarytrain[1,4], "\nRMSE = ", summarytrain[2,4]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

ggHorvathtrain <-ggplot(data = pheno.train,aes(x=dpc,y=Horvathdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.train$dpc)-10, y = max(pheno.train$Horvathdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytrain[1,1], "\nRMSE = ", summarytrain[2,1]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggGAtrain <-ggplot(data = pheno.train,aes(x=dpc,y=GAdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.train$dpc)-10, y = max(pheno.train$GAdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytrain[1,2], "\nRMSE = ", summarytrain[2,2]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggPlactrain <-ggplot(data = pheno.train,aes(x=dpc,y=Placdpc) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.train$dpc)-10, y = max(pheno.train$Placdpc), hjust = 0,vjust = 1 , label = paste("r = ", summarytrain[1,3], "\nRMSE = ", summarytrain[2,3]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

                
plot <- plot_grid(ggFetaltrain,ggHorvathtrain,ggGAtrain,ggPlactrain, align = "vh" , labels = c("A FBC", "B MTC", "C GAC", "D CPC"))
    
x.grob1 <- textGrob("Chronological Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
x.grob2 <- textGrob("in days post-conception", 
                   gp=gpar(fontface = "italic", col="black", fontsize=13))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
                     
#add to plot
pdf("figures/plotgrid.traindata.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)), nrow = 3, x.grob1, x.grob2, heights = c(0.93,0.04,0.03))
dev.off() 

########################################################################
##### Applying fetal brain clock to validation data (Figure 2)     #####
########################################################################


valdf <- pheno.val[,c(1,5,6,7,8,9)]
colnames(valdf) <- c("sample", "Age", "Horvath","Gestational", "Placental", "Fetal")

summaryval<-compareStats(valdf)
write.xlsx(summaryval, file = "Revision/Validation/summaryvalidation.xlsx")

ggFetalval <-ggplot(data = pheno.val,aes(x=dpc,y=FetalPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.val$dpc)-10, y = max(pheno.val$FetalPred), hjust = 0,vjust = 1 , label = paste("r = ", summaryval[1,4], "\nRMSE = ", summaryval[2,4]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )

plot <- plot_grid(ggFetalval, labels = c("Validation Dataset"))

x.grob1 <- textGrob("Chronological Age", 
                    gp=gpar(fontface="bold", col="black", fontsize=15))
x.grob2 <- textGrob("in days post-conception", 
                    gp=gpar(fontface = "italic", col="black", fontsize=13))
y.grob1 <- textGrob("Predicted Age", 
                    gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                    gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)

#add to plot
pdf("figures/plotgrid.valdata.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)), nrow = 3, x.grob1, x.grob2, heights = c(0.93,0.04,0.03))
dev.off() 


########################################################################
##### Applying fetal brain clock to adult data (Figure S3)         #####
########################################################################

##### BDR

load("/mnt/data1/BDR/QC/FINAL/BDR_DNAm_FINAL.rdat")
dim(betas)
#800916   1221

pheno.bdr <- pheno
betas.bdr <- betas

##### MTC
pheno.bdr$Horvath <- agep(betas.bdr, method = "horvath")

##### GAC
GestCpGs<- read.csv ("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Gest_age/cgprobesGApredictor.csv",stringsAsFactor=F, header=T)
GestCoeff<-GestCpGs[,2]
names(GestCoeff)<-GestCpGs[,1]
pheno.bdr$GestPred<-(((agep(betas.bdr, coef=GestCoeff, method="hannum")+41.726)*7)-280)/365

##### CPC
placCpGs<-read.csv("/mnt/data1/Jenny/EPIC_QC/iPSCsThroughDiff/wo2samples/Age_clocks/Plac_age/Placental_coeff.csv", header = T, stringsAsFactors = F)
placCoeff <- placCpGs[, 3]
names(placCoeff)<-placCpGs[,1]
pheno.bdr$PlacPred<-(((agep(betas.bdr, coef=placCoeff, method="hannum") + 13.062)*7)-280)/365

##### FBC
pheno.bdr$FetalPred <- fetalclock(betas.bdr, ageinyears = TRUE, dir = dir)

pheno.bdr <- pheno.bdr[!is.na(pheno.bdr$Age),]

bdrdf <- pheno.bdr[,c(3,10,15,16,17,18)]
colnames(bdrdf) <- c("sample", "Age", "Horvath","Gestational", "Placental", "Fetal")

summarybdr<-compareStats(bdrdf)
write.xlsx(summarybdr, file = "Adult/summarybdr.xlsx")

ggFetalbdr <-ggplot(data = pheno.bdr,aes(x=Age,y=FetalPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.bdr$Age)-10, y = max(pheno.bdr$FetalPred), hjust = 0,vjust = 1 , label = paste("r = ", summarybdr[1,4], "\nRMSE = ", summarybdr[2,4]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )


ggHorvathbdr <-ggplot(data = pheno.bdr,aes(x=Age,y=Horvath) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.bdr$Age)-10, y = max(pheno.bdr$Horvath), hjust = 0,vjust = 1 , label = paste("r = ", summarybdr[1,1], "\nRMSE = ", summarybdr[2,1]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggGAbdr <-ggplot(data = pheno.bdr,aes(x=Age,y=GestPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.bdr$Age)-10, y = max(pheno.bdr$GestPred), hjust = 0,vjust = 1 , label = paste("r = ", summarybdr[1,2], "\nRMSE = ", summarybdr[2,2]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )
        
ggPlacbdr <-ggplot(data = pheno.bdr,aes(x=Age,y=PlacPred) )+
  geom_abline(intercept=0,slope=1) +
  geom_point()+
  theme_cowplot(16) + 
  annotate("text",  x = min(pheno.bdr$Age)-10, y = max(pheno.bdr$PlacPred), hjust = 0,vjust = 1 , label = paste("r = ", summarybdr[1,3], "\nRMSE = ", summarybdr[2,3]), 
           fontface = "bold", size = 3.3) +
  theme(axis.line = element_line(colour = "black"),
        axis.title=element_blank(),
        plot.title=element_blank(),
        axis.text.x=element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.margin=unit(c(2,0,0.1,0.4), "lines") )


plot <- plot_grid(ggFetalbdr,ggHorvathbdr,ggGAbdr,ggPlacbdr, align = "vh" , labels = c("A FBC", "B MTC", "C GAC", "D CPC"))
    
x.grob1 <- textGrob("Chronological Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
x.grob2 <- textGrob("in years", 
                   gp=gpar(fontface = "italic", col="black", fontsize=13))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in years", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
                     
#add to plot
pdf("Adult/plotgrid.bdr.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)), nrow = 3, x.grob1, x.grob2, heights = c(0.93,0.04,0.03))
dev.off() 

range(pheno.bdr$FetalPred)
# -0.4425534 -0.2945308
median(pheno.bdr$FetalPred)
# -0.3680768
mean(pheno.bdr$FetalPred)
# -0.3683025
range(pheno.bdr$PlacPred)
# -0.20923277 -0.01666141
median(pheno.bdr$PlacPred)
# -0.1231471
mean(pheno.bdr$PlacPred)
# -0.1240357
range(pheno.bdr$GestPred)
# -0.30906345 -0.03306837
median(pheno.bdr$GestPred)
# -0.1897878
mean(pheno.bdr$GestPred)
# -0.190318

########################################################################
##### Testing for Sex effect                                       #####
########################################################################


### Testing dataset
summary(lm(FetalPred~dpc+Sex+dpc*Sex, data = pheno.test))
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 15.00372    9.25956   1.620 0.110316    
#dpc          0.83647    0.08325  10.048 1.49e-14 ***
#SexM        41.48709   11.62507   3.569 0.000706 ***
#dpc:SexM    -0.40261    0.10828  -3.718 0.000438 ***

pdf("figures/TestSexInteractionTest.pdf")
par(oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,2,2) + 0.1)
plot(pheno.test$FetalPred~pheno.test$dpc, col = c(as.factor(pheno.test$Sex)), pch = 16)
abline(lm(FetalPred~dpc, data = pheno.test[pheno.test$Sex == "F",]))
abline(lm(FetalPred~dpc, data = pheno.test[pheno.test$Sex == "M",]), col = "red")
abline(lm(FetalPred~dpc, data = pheno.test), col = "blue")
legend(x = "bottomright",legend = c("Female", "Male","Sum"), pch = 16, cex = 0.8, pt.cex = 1.5,col = c("black","red","blue"))
title(ylab = "Predicted Age",
      outer = TRUE, line = 1.5, cex = 1.5)
title(xlab = "Chronological Age",
      outer = TRUE, line = 1, cex = 1.5)
title(main = "Effect of Sex in Test data")
dev.off()   

### Taking out the two extreme values driving the model

pheno.test.red <- pheno.test[pheno.test$dpc>35 & pheno.test$dpc < 150,]

summary(lm(FetalPred~dpc+Sex+dpc*Sex, data = pheno.test.red))
#           Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 15.00372    8.18763   1.832   0.0719 .  
#dpc          0.83647    0.07361  11.363   <2e-16 ***
#SexM        20.24844   11.38790   1.778   0.0805 .  
#dpc:SexM    -0.18910    0.10754  -1.758   0.0839 .  

pdf("figures/TestSexInteractiontestWoOutliers.pdf")
par(oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,2,2) + 0.1)
plot(pheno.test.red$FetalPred~pheno.test.red$dpc, col = c(as.factor(pheno.test.red$Sex)), pch = 16)
abline(lm(FetalPred~dpc, data = pheno.test.red[pheno.test.red$Sex == "F",]))
abline(lm(FetalPred~dpc, data = pheno.test.red[pheno.test.red$Sex == "M",]), col = "red")
abline(lm(FetalPred~dpc, data = pheno.test.red), col = "blue")
legend(x = "bottomright",legend = c("Female", "Male","Sum"), pch = 16, cex = 0.8, pt.cex = 1.5,col = c("black","red","blue"))
title(ylab = "Predicted Age",
      outer = TRUE, line = 1.5, cex = 1.5)
title(xlab = "Chronological Age",
      outer = TRUE, line = 1, cex = 1.5)
title(main = "Effect of Sex in Reduced test data")
dev.off()  

### Validation dataset

pheno.val.sex <- pheno.val[complete.cases(pheno.val),]
pheno.val.sex <- pheno.val.sex[pheno.val.sex$Sex == "M" | pheno.val.sex$Sex == "F",]

summary(lm(FetalPred~dpc+Sex+dpc*Sex, data = pheno.val.sex))
#Coefficients:
#             Estimate    Std. Error    t   value Pr(>|t|)    
#(Intercept)  73.88269    2.85924  25.840  < 2e-16 ***
#dpc           0.35090    0.02293  15.305  < 2e-16 ***
#SexM        -14.25647    4.96774  -2.870  0.00489 ** 
#dpc:SexM      0.14036    0.04389   3.198  0.00179 ** 

# Significant interaction between sex and age on the predicted age with the FBC. Seems to be driven by outliers.


pdf("figures/ValSexInteractionTestFullData.pdf")
par(oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,2,2) + 0.1)
plot(pheno.val.sex$FetalPred~pheno.val.sex$dpc, col = c(as.factor(pheno.val.sex$Sex))+1, pch = 16)
abline(lm(FetalPred~dpc, data = pheno.val.sex[pheno.val.sex$Sex == "F",]), col = "red")
abline(lm(FetalPred~dpc, data = pheno.val.sex[pheno.val.sex$Sex == "M",]), col = "green")
abline(lm(FetalPred~dpc, data = pheno.val.sex), col = "black")
legend(x = "bottomright",legend = c("Female", "Male","Sum"), pch = 16, cex = 0.8, pt.cex = 1.5,col = c("red","green","black"))
title(ylab = "Predicted Age",
      outer = TRUE, line = 1.5, cex = 1.5)
title(xlab = "Chronological Age",
      outer = TRUE, line = 1, cex = 1.5)
title(main = "Effect of Sex in Reduced Validation data")
dev.off() 

#Removing 3 samples with dpc > 200
summary(lm(FetalPred~dpc+Sex+dpc*Sex, data = pheno.val.sex[pheno.val.sex$dpc <= 185,]))
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 60.368179   3.129229  19.292   <2e-16 ***
#dpc          0.488594   0.028180  17.338   <2e-16 ***
#SexM        -0.741960   4.640434  -0.160    0.873    
#dpc:SexM     0.002668   0.042313   0.063    0.950   


pdf("figures/ValSexInteractionTestNoOutliers.pdf")
par(oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,2,2) + 0.1)
plot(pheno.val.sex$FetalPred~pheno.val.sex$dpc, col = c(as.factor(pheno.val.sex$Sex))+1, pch = 16)
abline(lm(FetalPred~dpc, data = pheno.val.sex[pheno.val.sex$Sex == "F",]), col = "red")
abline(lm(FetalPred~dpc, data = pheno.val.sex[pheno.val.sex$Sex == "M",]), col = "green")
abline(lm(FetalPred~dpc, data = pheno.val.sex), col = "black")
legend(x = "bottomright",legend = c("Female", "Male","Sum"), pch = 16, cex = 0.8, pt.cex = 1.5,col = c("red","green","black"))
title(ylab = "Predicted Age",
      outer = TRUE, line = 1.5, cex = 1.5)
title(xlab = "Chronological Age",
      outer = TRUE, line = 1, cex = 1.5)
title(main = "Effect of Sex in Reduced Validation data")
dev.off() 


########################################################################
##### Setup ipsc - neuron datasets into one dataframe                   #####
########################################################################

##### Input: 5 datasets, all QCed and Normalized, with 4 epigenetic clocks already applied (3 datasets from our lab, supplemented with 2 publically available datasets)

load("/mnt/data1/LeoSteg/FetalClock/IPSCs/data.imm.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/IPSCs/data.price.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/IPSCs/Nazor_Data/data.nazor.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/IPSCs/FS_Data/data.fs.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/IPSCs/Sultanov_Data/data.sultanov.Rdata")

ipsc_summary <- as.data.frame(matrix(NA, nrow = 98, ncol = 7))

colnames(ipsc_summary) <- c("Basename / Sample", "Cell_State", "Cohort", "Horvath", "GestPred", "PlacPred", "FetalPred")

ipsc_summary[1:14,1] <- as.character(pheno.imm[,8])
ipsc_summary[1:14,2] <- as.character(pheno.imm[,10])
ipsc_summary[1:14,4] <- pheno.imm$Horvath
ipsc_summary[1:14,5] <- pheno.imm$GestPred
ipsc_summary[1:14,6] <- pheno.imm$PlacPred
ipsc_summary[1:14,7] <- pheno.imm$FetalPred
ipsc_summary[1:14,3]<- "Imm"

pheno.price <- pheno.price[pheno.price$Cell.type != "Keratinocyte",]
ipsc_summary[15:32,1] <- as.character(pheno.price$Basename)
ipsc_summary[15:32,2] <- as.character(pheno.price$Cell.type)
ipsc_summary[15:32,4] <- pheno.price$Horvath
ipsc_summary[15:32,5] <- pheno.price$GestPred
ipsc_summary[15:32,6] <- pheno.price$PlacPred
ipsc_summary[15:32,7] <- pheno.price$FetalPred
ipsc_summary[15:32,3] <- "Price"

pheno.nazor <- pheno.nazor[pheno.nazor$Cell_State == "iPSC" | pheno.nazor$Cell_State == "NPC",]
ipsc_summary[33:65,1] <- as.character(pheno.nazor$Geo_Accession)
ipsc_summary[33:65,2] <- pheno.nazor$Cell_State
ipsc_summary[33:65,4] <- pheno.nazor$Horvath
ipsc_summary[33:65,5] <- pheno.nazor$GestPred
ipsc_summary[33:65,6] <- pheno.nazor$PlacPred
ipsc_summary[33:65,7] <- pheno.nazor$FetalPred
ipsc_summary[33:65,3] <- "Nazor"

ipsc_summary[66:77,1] <- as.character(pheno.sultanov$Sample)
ipsc_summary[66:77,2] <- pheno.sultanov$Cell_state
ipsc_summary[66:77,4] <- pheno.sultanov$Horvath
ipsc_summary[66:77,5] <- pheno.sultanov$GestPred
ipsc_summary[66:77,6] <- pheno.sultanov$PlacPred
ipsc_summary[66:77,7] <- pheno.sultanov$FetalPred
ipsc_summary[66:77,3] <- "Sultanov"

ipsc_summary[78:98,1] <- as.character(pheno.fs$Sample)
ipsc_summary[78:98,2] <- pheno.fs$Cell_State
ipsc_summary[78:98,4] <- pheno.fs$Horvath
ipsc_summary[78:98,5] <- pheno.fs$GestPred
ipsc_summary[78:98,6] <- pheno.fs$PlacPred
ipsc_summary[78:98,7] <- pheno.fs$FetalPred
ipsc_summary[78:98,3] <- "F. - S."


ipsc_summary[ipsc_summary$Cell_State == "IPS" |ipsc_summary$Cell_State == "IPSC" ,2] <- "iPSC"
ipsc_summary[ipsc_summary$Cell_State == "N-D37" |ipsc_summary$Cell_State == "N-D58" ,2] <- "Neuron"

ipsc_summary$Tag <- paste(ipsc_summary$Cell_State,"_",ipsc_summary$Cohort, sep = "")

ipsc_summary$Tag <- factor(ipsc_summary$Tag,levels = c("iPSC_Imm", "iPSC_Price", "iPSC_F. - S.", "iPSC_Sultanov", "NPC_Imm", "NPC_Nazor", "Neuron_Imm","Neuron_Price", "Neuron_F. - S.", "Neuron_Sultanov"))
ipsc_summary$Horvath <- as.numeric(ipsc_summary$Horvath)
ipsc_summary$FetalPred <- as.numeric(ipsc_summary$FetalPred)

ipsc_summary$TagSimple <- "Neuron"
ipsc_summary[ipsc_summary$Cell_State == "iPSC","TagSimple"] <- "iPSC"
ipsc_summary[ipsc_summary$Cell_State == "NPC","TagSimple"] <- "NPC"
table(ipsc_summary$TagSimple)
#iPSC Neuron    NPC 
#30     48      4 
ipsc_summary$TagSimple <- factor(ipsc_summary$TagSimple, levels = c("iPSC", "NPC", "Neuron"))

save(ipsc_summary, file = "IPSCs/ipsc_summary.Rdata")

########################################################################
##### Plotting per dataset (Figure 3A)                             #####
########################################################################

cells.imm <- ipsc_summary[ipsc_summary$Cohort == "Imm",]
cells.price <- ipsc_summary[ipsc_summary$Cohort == "Price",]
cells.price$TagSimple <- factor(cells.price$TagSimple, levels = c("iPSC", "Neuron"))
cells.nazor <- ipsc_summary[ipsc_summary$Cohort == "Nazor",]
cells.nazor$TagSimple <- factor(cells.nazor$TagSimple, levels = c( "Neuron"))
cells.fs <- ipsc_summary[ipsc_summary$Cohort == "F. - S.",]
cells.fs$TagSimple <- factor(cells.fs$TagSimple, levels = c("iPSC", "Neuron"))
cells.sultanov<- ipsc_summary[ipsc_summary$Cohort == "Sultanov",]
cells.sultanov$TagSimple <- factor(cells.sultanov$TagSimple, levels = c("iPSC", "Neuron"))

stat.test.imm <-aov(FetalPred~TagSimple,data = cells.imm) %>% tukey_hsd() 
ypos.imm <- add_y_position(data = cells.imm, formula = FetalPred ~ TagSimple, test = tukey_hsd(aov(FetalPred~TagSimple,data = cells.imm)))

ggImm <- ggplot(data = cells.imm,aes(x=TagSimple,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, aes(fill = TagSimple)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_brewer(palette = "Set2") +
  stat_pvalue_manual(stat.test.imm, label = "p.adj", y.position = ypos.imm$y.position+c(-8,-4,0), size = 3)
  
    
stat.test.price <-t_test(FetalPred~TagSimple,data = cells.price)
ypos.price <- add_y_position(data = cells.price, formula = FetalPred ~ TagSimple, test = t_test(FetalPred~TagSimple,data = cells.price))

ggPrice <- ggplot(data = cells.price,aes(x=TagSimple,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(1,3)], aes(fill = TagSimple)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  stat_pvalue_manual(stat.test.price, label = "p", y.position = ypos.price$y.position+3, size = 3)


ggNazor <- ggplot(data = cells.nazor,aes(x=TagSimple,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(3)], aes(fill = TagSimple), width = 0.4) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
  

stat.test.fs <-t_test(FetalPred~TagSimple,data = cells.fs)
ypos.fs <- add_y_position(data = cells.fs, formula = FetalPred ~ TagSimple, test = t_test(FetalPred~TagSimple,data = cells.fs))

ggFS <- ggplot(data = cells.fs,aes(x=TagSimple,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(1,3)], aes(fill = TagSimple)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  stat_pvalue_manual(stat.test.fs, label = "p", y.position = ypos.fs$y.position+3, size = 3)
  

stat.test.sultanov <-t_test(FetalPred~TagSimple,data = cells.sultanov)
ypos.sultanov <- add_y_position(data = cells.sultanov, formula = FetalPred ~ TagSimple, test = t_test(FetalPred~TagSimple,data = cells.sultanov))

ggSultanov <- ggplot(data = cells.sultanov,aes(x=TagSimple,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(1,3)], aes(fill = TagSimple)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  stat_pvalue_manual(stat.test.sultanov, label = "p", y.position = ypos.sultanov$y.position+3, size = 3)


plot <- plot_grid(ggImm, ggPrice, ggNazor, ggFS, ggSultanov, align = "vh" , labels = c("A Imm", "B Price", "C Nazor", "D F.- S.", "E Sultanov"))

x.grob <- textGrob("Cell Stage", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
#add to plot
pdf("iPSCs/AllDataFetal.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, bottom = x.grob, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)))
dev.off() 

########################################################################
##### Plotting merged cellular data (Figure3B) and lme model       #####
########################################################################

library(lme4)
library(lmerTest)

ipsc_summary$D1 <- 0
ipsc_summary[ipsc_summary$Cell_State == "NPC", "D1"] <- 1 
ipsc_summary$D2 <- 0
ipsc_summary[ipsc_summary$Cell_State == "Neuron", "D2"] <- 1 


coef(summary(lmerTest::lmer(ipsc_summary$FetalPred ~ ipsc_summary$D1 + ipsc_summary$D2+(1|ipsc_summary$Cohort))))
#                   Estimate Std. Error        df   t value     Pr(>|t|)
#(Intercept)     72.810219   2.555235  4.635174 28.494525 2.194738e-06
#ipsc_summary$D1  4.332468   2.654613 93.559110  1.632052 1.060300e-01
#ipsc_summary$D2 13.346693   1.756687 94.757342  7.597650 2.126568e-11

multilevel_stats <- as.data.frame(matrix(NA, ncol = 6, nrow = 1))
rownames(multilevel_stats) <- c("FetalPred")
colnames(multilevel_stats) <- c("NPC Estimate", "NPC SE", "NPC p value", "Neuron Estimate", "Neuron SE", "Neuron p value")

multilevel_stats[1,c(1:3)] <- coef(summary(lmerTest::lmer(ipsc_summary$FetalPred ~ ipsc_summary$D1 + ipsc_summary$D2+(1|ipsc_summary$Cohort))))[2,c(1,2,5)]
multilevel_stats[1,c(4:6)] <- coef(summary(lmerTest::lmer(ipsc_summary$FetalPred ~ ipsc_summary$D1 + ipsc_summary$D2+(1|ipsc_summary$Cohort))))[3,c(1,2,5)]

ipsc_summary$Cohort <- factor(ipsc_summary$Cohort, levels = c("Imm", "Price", "Nazor", "F. - S.", "Sultanov"))

multilevel_stats <- read_excel("multilevel_stats.xlsx")

stat.test.lme <- as.data.frame(stat.test.imm[1:2,])
stat.test.lme$p.adj <- signif(as.numeric(c(multilevel_stats[1,3],multilevel_stats[1,6])),3)
  
ggiPSCSumm <-ggplot(data = ipsc_summary,aes(x=TagSimple,y=FetalPred, color = Cohort) )+
  geom_boxplot(outlier.shape = NA, color = "black", fill = "#dfdfdf" ) +
  geom_jitter(size=2, alpha=1, width = 0.15) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_brewer(palette = "Set1") + 
  stat_pvalue_manual(stat.test.lme, label = "p.adj", y.position = c(100,105), size = 3)

    
x.grob <- textGrob("Cell Stage", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
                     
legend <- legendGrob(labels = levels(ipsc_summary$Cohort), vgap = 0.1,pch = 16, gp = gpar(fontsize = 11,col = brewer.pal(5, "Set1"))  )                 
#add to plot
pdf("IPSCs/iPSCSummFBCColored.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,ggiPSCSumm,legend, bottom = x.grob, nrow = 1, widths = c(0.04, 0.03, 0.81,0.12)))
dev.off()

########################################################################
##### Setup ipsc-MN dataset and plotting  (Figure 3C)              #####
########################################################################

load("mnt/data1/LeoSteg/FetalClock/IPSCs/ipsc_mn_summary.Rdata")

stat.test.mn <-t_test(FetalPred~Cell_State,data = ipsc_mn_summary)
ypos.mn <- add_y_position(data = ipsc_mn_summary, formula = FetalPred ~ Cell_State, test = t_test(FetalPred~Cell_State,data = ipsc_mn_summary))

ggMN <- ggplot(data = ipsc_mn_summary,aes(x=Cell_State,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(1,3)], aes(fill = Cell_State)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  stat_pvalue_manual(stat.test.mn, label = "p", y.position = ypos.mn$y.position+3, size = 3)

plot <- plot_grid(ggMN, labels = c("iPSC - Motor neuron"))

x.grob <- textGrob("Cell Stage", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
#add to plot
pdf("IPSCs/ipscMN.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, bottom = x.grob, nrow = 1, widths = c(0.04, 0.03, 0.86,0.07)))
dev.off() 

########################################################################
##### Setup ESC - neuron dataset and plotting  (Figure 3D)         #####
########################################################################


load("/mnt/data1/LeoSteg/FetalClock/IPSCs/Kim_Data/data.kim.Rdata")
load("/mnt/data1/LeoSteg/FetalClock/IPSCs/Nazor_Data/data.nazor.Rdata")

esc_summary <- as.data.frame(matrix(NA, nrow = 27, ncol = 7))

colnames(esc_summary) <- c("Basename / Sample", "Cell_State", "Cohort", "Horvath", "GestPred", "PlacPred", "FetalPred")

pheno.nazor <- pheno.nazor[pheno.nazor$Cell_State == "ESC" | pheno.nazor$Cell_State == "ESC NPC",]
esc_summary[1:21,1] <- as.character(pheno.nazor$Geo_Accession)
esc_summary[1:21,2] <- pheno.nazor$Cell_State
esc_summary[1:21,4] <- pheno.nazor$Horvath
esc_summary[1:21,5] <- pheno.nazor$GestPred
esc_summary[1:21,6] <- pheno.nazor$PlacPred
esc_summary[1:21,7] <- pheno.nazor$FetalPred
esc_summary[1:21,3] <- "Nazor"

esc_summary[22:27,1] <- as.character(pheno.kim$Sample)
esc_summary[22:27,2] <- as.character(pheno.kim$Cell_State)
esc_summary[22:27,4] <- pheno.kim$Horvath
esc_summary[22:27,5] <- pheno.kim$GestPred
esc_summary[22:27,6] <- pheno.kim$PlacPred
esc_summary[22:27,7] <- pheno.kim$FetalPred
esc_summary[22:27,3] <- "Kim"

esc_summary[esc_summary$Cell_State == "hES" ,2] <- "ESC"
esc_summary[esc_summary$Cell_State == "ESC NPC",2] <- "NPC"

table(esc_summary$Cell_State)
#ESC Neuron    NPC 
#21      2      4  

esc_summary$Horvath <- as.numeric(esc_summary$Horvath)
esc_summary$FetalPred <- as.numeric(esc_summary$FetalPred)
esc_summary$Cell_State <- factor(esc_summary$Cell_State, levels = c("ESC", "NPC", "Neuron"))

save(esc_summary, file = "IPSCs/ESC/esc_summary.Rdata")

cells.nazor <- esc_summary[esc_summary$Cohort == "Nazor",]
cells.nazor$Cell_State <- factor(cells.nazor$Cell_State, levels = c("ESC", "NPC"))
cells.kim <- esc_summary[esc_summary$Cohort == "Kim",]

stat.test.kim <-aov(FetalPred~Cell_State,data = cells.kim) %>% tukey_hsd() 
ypos.kim <- add_y_position(data = cells.kim, formula = FetalPred ~ Cell_State, test = tukey_hsd(aov(FetalPred~Cell_State,data = cells.kim)))

ggKim <- ggplot(data = cells.kim,aes(x=Cell_State,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, aes(fill = Cell_State)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_fill_brewer(palette = "Set2") +
  stat_pvalue_manual(stat.test.kim, label = "p.adj", y.position = ypos.kim$y.position +c(0,-2,-4), size = 3)
  
stat.test.nazor <-t_test(FetalPred~Cell_State,data = cells.nazor)
ypos.nazor <- add_y_position(data = cells.nazor, formula = FetalPred ~ Cell_State, test = t_test(FetalPred~Cell_State,data = cells.nazor))

ggNazor <- ggplot(data = cells.nazor,aes(x=Cell_State,y=FetalPred) )+
  geom_boxplot(outlier.shape = NA, fill =brewer.pal(3, "Set2")[c(1,2)], aes(fill = Cell_State)) +
  geom_jitter(color="black", size=1, alpha=1, width = 0.1) +
  theme_cowplot(16) +
  theme(legend.position = "none",panel.border = element_blank(),
        axis.text=element_text(size=12),plot.margin=unit(c(2,0,0.1,0.4), "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  stat_pvalue_manual(stat.test.nazor, label = "p", y.position = ypos.nazor$y.position+3, size = 3)
  
plot <- plot_grid(ggNazor, ggKim, align = "vh" , labels = c("A Nazor", "B Kim"))

x.grob <- textGrob("Cell Stage", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))
y.grob1 <- textGrob("Predicted Age", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot = 90)
y.grob2 <- textGrob("in days post-conception", 
                     gp=gpar(fontface = "italic", col="black", fontsize=13), rot = 90)
#add to plot
pdf("IPSCs/ESCNeuron.pdf")
grid.arrange(arrangeGrob(y.grob1, y.grob2,plot, nrow = 1,widths = c(0.04, 0.03, 0.86,0.07)), nrow = 3, x.grob, heights = c(0.46,0.04,0.5))
dev.off() 


