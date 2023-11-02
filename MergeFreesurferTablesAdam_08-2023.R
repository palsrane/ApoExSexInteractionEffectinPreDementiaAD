# Changes made by June 2023
# Changed the way duplicates were removed using RID and VISCODE, to RID and EXAMDATE
# Updated to latest data files for ADNI3
# Changed output file from allnewfreesurfer_4 to allnewfreesurfer_5

MergeSubjectTime <- function(df1, df2, mergecol, timecol1, timecol2) {
  
  
  m_data <- lapply(intersect(df1[[mergecol]], df2[[mergecol]]),function(id) {
    d1   <- subset(df1,df1[[mergecol]] == id)
    d2   <- subset(df2,df2[[mergecol]] == id)
    
    d1$indices <- sapply(d1[[timecol1]],function(d) which.min(abs(d2[[timecol2]] - d)))
    d2$indices <- 1:nrow(d2)
    
    base::merge(d1, d2, by = c(mergecol, 'indices'), all.x=FALSE, all.y=FALSE)
  })
  mergeddata <- do.call(rbind, m_data)
  mergeddata$indices <- NULL
  return(mergeddata)
}


ApoeIndADNI <- function(df) {
  indvar <- c()
  drop.rows <- c()
  na.drop <- which(is.na(df$APGEN1), arr.ind = TRUE)
  if(length(na.drop) != 0) {
    df <- df[-na.drop,]
  }
  for( i in 1:nrow(df)) {
    a <- df["APGEN1"][i,]
    b <- df["APGEN2"][i,]
    geneval <- c(a, b)
    #if(a == 2 | b == 2) {
    # drop.rows <- append(drop.rows, i)
    #}
    if(sum(geneval == 4) == 0) { 
      indvar <- append(indvar, 0)
    } else if(sum(geneval == 4) == 1){
      indvar <- append(indvar, 1)
    } else {
      indvar <- append(indvar, 2)
    }
  }
  df["APOE4Ind"] <- indvar
  return(df)
}

ApoeInd2ADNI <- function(df) {
  indvar <- c()
  drop.rows <- c()
  na.drop <- which(is.na(df$APGEN1), arr.ind = TRUE)
  if(length(na.drop) != 0) {
    df <- df[-na.drop,]
  }
  for( i in 1:nrow(df)) {
    a <- df["APGEN1"][i,]
    b <- df["APGEN2"][i,]
    geneval <- c(a, b)
    if(sum(geneval == 2) == 0) { 
      indvar <- append(indvar, 0)
    } else if(sum(geneval == 2) == 1){
      indvar <- append(indvar, 1)
    } else {
      indvar <- append(indvar, 2)
    }
  }
  df["APOE2Ind"] <- indvar
  return(df)
}

FindIndices <- function(data, trainvar, trainlevel) {
  #
  # Constructs vector of rownames or indices based on subsetting parameters
  # Helpful for GroupModel if interested in training models based
  #  on subset of data
  #
  # Args:
  #  data: full dataframe
  #  trainvar: string of variable name (must be factor variable) 
  #  trainlevel: string of factor level of interest
  #
  # Returns:
  #  Vector of rownames or indices
  #
  index      <- which(data[[trainvar]] == trainlevel, 
                      arr.ind = TRUE)
  index      <- as.vector(index)
  return(index)
}



library(ADNIMERGE)

#load freesurfer tables
# ADNI1 Data  - UCSFFSX_11_02_15 #Check
# ADNI1 3T    - UCSFFSX51_ADNI1_3T_02_01_16 #Check
# ADNI2/GO    - UCSFFSX51_ADNI_11_08_19 #Check
# ADNI3       -UCSFFSX_12_13_21 #Check

#6/22/2023 Changed the 'duplicate removal'. Earlier it was done with RID and VISCODE, but decided to change it to RID and EXAMDATE

#load adni1 1.5T
adni1_1.5         <- read.csv("C:/Projects/APOESex/FreesurferData/UCSFFSX_11_02_15.csv")#("ADNI1_1.5T.csv")
adni1_1.5$COLPROT <- "ADNI1"

### Check QC before dropping
adni1_1.5 = adni1_1.5[adni1_1.5$OVERALLQC=="Pass"| adni1_1.5$OVERALLQC=="Hippocampus Only" |((adni1_1.5$OVERALLQC=="" | adni1_1.5$OVERALLQC=="Partial") & adni1_1.5$TEMPQC=="Pass"),]
adni1_1.5 = adni1_1.5[!duplicated(adni1_1.5[,c("RID", "EXAMDATE")]),] #"VISCODE"

#load adni1 3T

adni1_3           <- read.csv("C:/Projects/APOESex/FreesurferData/UCSFFSX51_ADNI1_3T_02_01_16.csv")#("ADNI1_3T.csv")
adni1_3$COLPROT   <- "ADNI1_3T"

### Check QC before dropping
adni1_3 = adni1_3[adni1_3$OVERALLQC=="Pass"| adni1_3$OVERALLQC=="Hippocampus Only"| (adni1_3$OVERALLQC=="Partial" & adni1_3$TEMPQC=="Pass"),]

#drop duplicates
adni1_3 <- adni1_3[!duplicated(adni1_3[,c("RID", "EXAMDATE")]),] #"VISCODE"

#load adni2
# adni2_go_old          <- read.csv("C:/Projects/ADNI/FreesurferData/UCSFFSX51_11_08_19.csv") 
adni2_go       <- read.csv("C:/Projects/APOESex/FreesurferData/UCSFFSX51FINAL_11_08_19.csv") #same as UCSFFSL51_03_01_22.csv

adni2_go$VISCODE  <- adni2_go$VISCODE2

### Check QC before dropping
adni2_go = adni2_go[adni2_go$OVERALLQC=="Pass"| adni2_go$OVERALLQC=="Hippocampus Only"| ((adni2_go$OVERALLQC=="Partial"|adni2_go$OVERALLQC=="") & adni2_go$TEMPQC=="Pass"),]

#drop nv site
adni2_go <- adni2_go[adni2_go$COLPROT!="nv",]
#drop duplicates
adni2_go <- adni2_go[!duplicated(adni2_go[,c("RID", "EXAMDATE")]),] #"VISCODE"

#load adni3
# adni3             <- read.csv("C:/Projects/ADNI/FreesurferData/UCSFFSX6_12_13_21.csv")
# adni3             <- read.csv("C:/Projects/ADNI/FreesurferData/UCSFFSX6_01_21_22.csv")
# adni3             <- read.csv("C:/Projects/APOESex/FreesurferData/UCSFFSX6_06_07_22_ADNI3.csv")
adni3             <- read.csv("C:/Projects/APOESex/FreesurferData/UCSFFSX6_08_17_22_22Jun2023.csv")
adni3$VISCODE     <- adni3$VISCODE2
### Checking QC before. Dropping OVERALLQC="", as those scans are not verified by Marc as ok.
adni3 = adni3[adni3$OVERALLQC=="Pass" | adni3$OVERALLQC=="Hippocampus Only" | (adni3$OVERALLQC=="Partial"& adni3$TEMPQC=="Pass"),]

#drop duplicates
adni3 <- adni3[!duplicated(adni3[,c("RID", "EXAMDATE")]),] #"VISCODE"


#merge freesurfer tables
colnames.keeps <- Reduce(intersect,list(colnames(adni1_1.5),colnames(adni1_3),colnames(adni2_go),colnames(adni3)))
adni1_1.5 <- adni1_1.5[ ,colnames.keeps]
adni1_3   <- adni1_3[   ,colnames.keeps]
adni2_go  <- adni2_go[  ,colnames.keeps]
adni3     <- adni3[     ,colnames.keeps]


fulladni <- rbind(adni1_1.5,
                      adni1_3,
                      adni2_go,
                      adni3)

rm(adni1_1.5,adni1_3,adni2_go,adni3)

#merge with adnimerge table by subject and closest EXAMDATE
fulladni$EXAMDATE <- as.POSIXct(fulladni$EXAMDATE)

adnimerge.subset  <- ADNIMERGE::adnimerge[,c("RID", "VISCODE","AGE", "M", "EXAMDATE")]
adnimerge.hardmerge <- ADNIMERGE::adnimerge[,c("RID", "PTGENDER","PTEDUCAT")]
adnimerge.hardmerge <- adnimerge.hardmerge[!duplicated(adnimerge.hardmerge$RID),]
adnimerge.subset$EXAMDATE <- as.POSIXct(adnimerge.subset$EXAMDATE)

fulladni                  <- merge(fulladni, adnimerge.hardmerge, by="RID")
fulladni                  <- MergeSubjectTime(fulladni, adnimerge.subset, "RID", "EXAMDATE", "EXAMDATE")


fulladni           <- fulladni[order(fulladni$RID, fulladni$M, decreasing = FALSE), ]
rownames(fulladni) <- 1:nrow(fulladni)

alldiagnosis<-ADNIMERGE::dxsum
alldiagnosis$USERDATE[is.na(alldiagnosis$USERDATE) & !is.na(alldiagnosis$EXAMDATE)]=alldiagnosis$EXAMDATE[is.na(alldiagnosis$USERDATE) & !is.na(alldiagnosis$EXAMDATE)]
alldiagnosis$EXAMDATE[is.na(alldiagnosis$EXAMDATE) & !is.na(alldiagnosis$USERDATE)]=alldiagnosis$USERDATE[is.na(alldiagnosis$EXAMDATE) & !is.na(alldiagnosis$USERDATE)]

alldiagnosis<-alldiagnosis[!is.na(alldiagnosis$USERDATE) & !is.na(alldiagnosis$EXAMDATE),]
alldiagnosis$dt=as.Date(alldiagnosis$USERDATE)-as.Date(alldiagnosis$EXAMDATE)
alldiagnosis<-alldiagnosis[!is.na(alldiagnosis$dt),]

alldiagnosis$USERDATE[alldiagnosis$dt>0]<-alldiagnosis$EXAMDATE[alldiagnosis$dt>0] #keep the smaller of the EXAMDATE and USERDATE as USEDATE which will be used to align the diagnosis info to other parameters

alldiagnosis<- alldiagnosis[,c("RID", "USERDATE", "DIAGNOSIS")]
alldiagnosis$USERDATE <- as.POSIXct(alldiagnosis$USERDATE)
fulladni <- MergeSubjectTime(fulladni, alldiagnosis, "RID", "EXAMDATE.x", "USERDATE")

#drop duplicates
fulladni <- fulladni[order(fulladni$RID,fulladni$EXAMDATE.x),]
fulladni <- fulladni[!duplicated(fulladni[,c("RID", "EXAMDATE.x", "COLPROT")]),] #"VISCODE.x" replaced with "EXAMDATE.x"

apoe.adni1   <- apoeres[,c("RID", "APGEN1", "APGEN2")]
apoe.adni2go <- apoego2[,c("RID", "APGEN1", "APGEN2")]
apoe.adni3   <- apoe3[,c("RID", "APGEN1", "APGEN2")]
adni.apoe    <- rbind(apoe.adni1, apoe.adni2go, apoe.adni3)
adni.apoe    <- adni.apoe[!duplicated(adni.apoe$RID),]
fulladni     <- merge(fulladni, adni.apoe, by="RID")


#0 for no e4, 1 for 1 e4 allele, 2 for 2 e4 alleles
fulladni     <- ApoeIndADNI(fulladni)
fulladni     <- ApoeInd2ADNI(fulladni)


# Correct age variable since it is originally baseline age in adnimerge
fulladni$AGE <- fulladni$AGE + (fulladni$M / 12)
fulladni     <- fulladni[order(fulladni$RID, fulladni$AGE, decreasing = FALSE), ] 

#fix classes
fulladni$STUDY    <- fulladni$COLPROT
fulladni$STUDY    <- factor(fulladni$STUDY)
fulladni$AGE      <- as.numeric(fulladni$AGE)
fulladni$PTGENDER <- factor(fulladni$PTGENDER)
fulladni <- fulladni[-which(fulladni$STUDY==""),]
fulladni$STUDY    <- factor(fulladni$STUDY)
fulladni$APOE4Ind  <- factor(fulladni$APOE4Ind)
fulladni$APOE2Ind  <- factor(fulladni$APOE2Ind)

#keep just important features and only complete rows

# fulladni <- fulladni[fulladni$OVERALLQC!="Fail",] #Fail cases have already been removed


features <- c( "ST29SV", "ST88SV", "ST24CV", "ST32CV", "ST40CV", 
               "ST26CV", "ST83CV", "ST91CV", "ST99CV", "ST85CV",#hippocampus and meta roi 
               "ST54CV", "ST113CV","ST14CV", "ST73CV") #Left rostral ACC, Right rostral ACC, Left Caudal ACC,Right Caudal ACC
#covs     <- c("AGE", "PTGENDER","PTEDUCAT", "APOE4Ind", "APOE2Ind", "STUDY")
covs     <- c("AGE", "PTGENDER","PTEDUCAT", "APGEN1", "APGEN2", "STUDY")

features.icv <- c("ST10CV")
covs.icv     <- c("PTGENDER", "STUDY")
fulladni <- fulladni[order(fulladni$RID, fulladni$AGE, decreasing = FALSE),]
completerows <- complete.cases(fulladni[,c(features, covs, "DIAGNOSIS")])
fulladni <- fulladni[completerows,]

fulladni$Baseline <- !duplicated(fulladni[,c("RID", "STUDY")])
fulladni.baseline <- subset(fulladni, Baseline==TRUE)
train.rows <- FindIndices(fulladni.baseline, "DIAGNOSIS", "CN") #FindIndices function is not present

fulladni.features.baseline <- fulladni.baseline[,features]
fulladnicovs <- fulladni.baseline[,covs]
fulladni.icv.features <-as.data.frame(fulladni.baseline[,features.icv])
colnames(fulladni.icv.features) <- "ICV"
fulladni.icv.covs     <- fulladni.baseline[,c("PTGENDER", "STUDY")]
long.features <- fulladni[,features]
long.cov      <- fulladni[,covs]
long.icv  <- as.data.frame(fulladni[,features.icv])
colnames(long.icv) <- "ST10CV"
long.covs.icv <- fulladni[,covs.icv]

fulladni.features.baseline <- fulladni.features.baseline[train.rows,]
fulladnicovs <- fulladnicovs[train.rows,]
fulladni.icv.features <- as.data.frame(fulladni.icv.features[train.rows,])
colnames(fulladni.icv.features) <- "ICV"
fulladni.icv.covs <- fulladni.icv.covs[train.rows,]
all.data.list <- list("baseline.im" = fulladni.features.baseline,
                      "baseline.covs" = fulladnicovs,
                      "baseline.icv.feats"=fulladni.icv.features,
                      "baseline.icv.cov" = fulladni.icv.covs,
                      "long.features" = long.features,
                      "long.covs" = long.cov,
                      "long.icv" =long.icv,
                      "long.cov.icv" =long.covs.icv,
                      "train.indicies" = train.rows,
                      "fulladni" = fulladni)
saveRDS(all.data.list, "C:/Projects/APOESex/FreesurferData/allnewfreesurfer_5.rds")



