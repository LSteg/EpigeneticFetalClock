#####################################################################################################################################################
####### FETAL BRAIN CLOCK                                                                                                                         ###
### Function to predict fetal age of brain samples and for age predictions in cellular stem cell models and their derived neuronal/cortical cells ###
### It was build using elastic net regression as described by Horvath (2013)                                                                      ###
### Training data consisted of 193 fetal brain samples with DNA methylation quantified using 450K or EPIC methylation array data                  ###
### Output of clock in days post-conception (default, ageinyears = FALSE) or in years after birth (ageinyears = TRUE)                             ### 
###                                                                                                                                               ###
### Author: Leonard Steg                                                                                                                          ###
#####################################################################################################################################################

FetalClock <- function(betas, # betas matrix; rownames = cpgs, colnames = sample IDs
                       ageinyears = FALSE,
                       dir = "/mnt/data1/reference_files/FetalBrainClock/") {  # directory with coefficients (fetalclock_coefficients_wo_intercept.txt) and reference values (fetalclock_reference.rdat)
  
  ### Read in coefficients
  coef<-read.table(paste0(dir,"fetalclock_coefficients_wo_intercept.txt",sep=""),stringsAsFactor=F,header=T)
  
  ### Find overlap between betas and coefficients
  
  overlap<-coef[which(coef$probe %in% rownames(betas)),]
  if (nrow(overlap) < nrow(coef) ){
    print("Some probes of the coefficients are missing in the betas. We will need to impute values here - The final predictions will be less accurate")
  } else {
    print("All the probes overlap between your data and the clock probes - No need for imputing missing values")
  }
  
  ### IF probes are missing: 
  ### Add reference betas for missing values -  imputation method adapted from:  https://github.com/qzhang314/DNAm-based-age-predictor   
  
  if (length(overlap) < nrow(coef)) {
    ## Transform betas to be cpg in col
    betas<-t(betas)
    
    ###########  Read in ref data and match
    load(paste0(dir,"fetalclock_reference.rdat",sep=""))
    ref<-ref[which(names(ref) %in% coef$probe) , drop=F]
    
    betas<-betas[,colnames(betas)%in%names(ref)]
    if(ncol(betas)<length(ref)){
      missprobe<-setdiff(names(ref),colnames(betas))
      refmiss<-ref[missprobe]
      refmiss<-matrix(refmiss,ncol=length(missprobe),nrow=nrow(betas),byrow=T)
      refmiss<-as.data.frame(refmiss)
      colnames(refmiss)<-missprobe
      rownames(refmiss)<-rownames(betas)
      betas<-cbind(betas,refmiss)
    }
    
    betas<-betas[,names(ref)]     ########### match betas
    
    ### Replace missing probes with reference values
    
    ## Impute function 
    imputeNA<-function(betas){
      betas[is.na(betas)]<-mean(betas,na.rm=T)
      return(betas)
    }
    
    ## Apply function 
    betasNona<-apply(betas,2,function(x) imputeNA(x))  
    
    
    ## Tranform betas - CpG in row
    betas<-t(betasNona)
  
    
    ### Age prediction
    
    coef<-coef[match(rownames(betas), coef$probe),]
    pred<-coef$coef%*%betas-58.2413413292108
    
    
    ## ELSE (No impuation for missing values needed):
    
  } else {
    
    ### Direct Age prediction
    
    coef<-coef[match(rownames(betas), coef$probe),]
    pred<-coef$coef%*%betas-58.2413413292108
    
    
  }
  
  ### IF output in years transform from dpc into years after birth
  
  if (ageinyears == T) {
    pred = (pred-280)/365
  } 
  t(pred)
}
                         
                       
                       