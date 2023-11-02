#Changes made in July 2023
#Changed ApoE to APOE
#Changed output file from freesurferharmnewadni3_trial_2.csv to freesurferharmnewadni3_trial_3.csv

feature.correction <- function(training.data,  data, formula, cr.feat1, feat) {
  model <- lm(formula = as.formula(formula), data = training.data)
  mean.val1 <- mean(training.data[[cr.feat1]])
  mean.val1 <- rep(mean.val1, nrow(data))
  coef.correction1 <- model$coefficients[[cr.feat1]]
  new.feat <- data[[feat]]
  new.feat <- (new.feat - (coef.correction1*data[[cr.feat1]]))
  new.feat <- new.feat  + (coef.correction1*mean.val1)
  return(new.feat)
}

freesurfer.data.list <- readRDS("C:/Projects/APOESex/FreesurferData/allnewfreesurfer_5.rds")
baseline.ims.fs <- freesurfer.data.list$baseline.im
baseline.ims.fs$ST10CV <- NULL
baseline.covs.fs <- freesurfer.data.list$baseline.covs
baseline.ims.icv.fs <- freesurfer.data.list$baseline.icv.feats
baseline.ims.icv.covs.fs <- freesurfer.data.list$baseline.icv.cov
training.indicies.fs <- freesurfer.data.list$train.indicies
long.ims.fs <- freesurfer.data.list$long.features
long.ims.fs$ST10CV <- NULL
long.covs.fs  <- freesurfer.data.list$long.covs
long.ims.icv.fs <- freesurfer.data.list$long.icv
long.covs.icv.fs <- freesurfer.data.list$long.cov.icv
full.fs          <- freesurfer.data.list$fulladni
#harmonize baseline icv
icv_harm.fs <- ComGamHarm(feature.data   = baseline.ims.icv.fs,
                       covar.data        = baseline.ims.icv.covs.fs,
                       eb=FALSE)

baseline.icv.harm <- as.numeric(t(icv_harm.fs$harm.results))
baseline.covs.fs$ICV_harmonized <- baseline.icv.harm

long.icv.harm <- t(ApplyHarm(long.ims.icv.fs, long.covs.icv.fs, 
                                           icv_harm.fs$stan.dict, icv_harm.fs$shift.scale.params, 
                                           icv_harm.fs$models.list))

# harmonize baseline features
features_harm.fs <- ComGamHarm(feature.data = baseline.ims.fs,
                          covar.data        = baseline.covs.fs,
                          smooth.terms      = "AGE",
                          k                 = 9,
                          eb                = TRUE)


long.covs.fs$ICV_harmonized <- as.numeric(long.icv.harm)
long.features.harm <- as.data.frame(t(ApplyHarm(long.ims.fs, long.covs.fs, 
                                                features_harm.fs$stan.dict, features_harm.fs$shift.scale.params, 
                                                features_harm.fs$models.list)))

colnames(long.features.harm) <- paste(colnames(long.features.harm), "_harmonized", sep="")
long.features.harm$ICV_harmonized <- as.numeric(long.icv.harm)
full.fs <- cbind(long.features.harm, full.fs)
full.fs$ICV <- full.fs$ST10CV

# # grouping a subject's scans by the ADNI arm and taking the median of harmonized ICV 
# # This was not used and it is ok to just use the median of all the records of harmonized_ICV per subject as during 
# # harmonization of ICV, the effect of different sites/ADNI arms is already nutralized
# full.fs$RIDCOLPROT=interaction(full.fs$RID,full.fs$COLPROT)
# uniqueRIDCOLPROT=unique(full.fs$RIDCOLPROT)
# full.fs$ICV_harmonized_median=c(0)
# for (RIDCOLPROT in uniqueRIDCOLPROT) {
#   full.fs$ICV_harmonized_median[full.fs$RIDCOLPROT==RIDCOLPROT]=median(full.fs$ICV_harmonized[full.fs$RIDCOLPROT==RIDCOLPROT])
# }
# # full.fs=subset(full.fs, select = -c(RIDCOLPROT))


full.fs.split <- split(full.fs, full.fs$RID)
full.fs.combine <- list()
for(i in 1:length(full.fs.split)) {
  subj <- full.fs.split[[i]]
  subj$ICV_harmonized_median <- median(subj$ICV_harmonized)
  full.fs.combine[[i]] <- subj
}
full.fs <- do.call(rbind, full.fs.combine)

### ICV_adjust values
icv.adj.frame.init <- data.frame(matrix(nrow = nrow(full.fs)))
training.data.fc <- subset(full.fs, Baseline==TRUE & DIAGNOSIS=="CN")
correction.features <- colnames(training.data.fc)[grep("_harmonized", colnames(training.data.fc))]
correction.features <- correction.features[-which(correction.features=="ICV_harmonized")]
correction.features <- correction.features[-which(correction.features=="ICV_harmonized_median")]

for(i in 1:length(correction.features)) {
  corrected.name <- correction.features[i]
  print(correction.features[i])
  formula.fs <- paste(corrected.name, "~ AGE + PTGENDER + APOE4Ind + APOE2Ind +PTEDUCAT + ICV_harmonized_median", sep="")
  corrected.name.icv.adj <- paste(corrected.name, "_icv_adjusted", sep="")
  corrected.feat <- feature.correction(training.data.fc, full.fs, formula.fs, "ICV_harmonized_median", corrected.name)
  icv.adj.frame.init[corrected.name.icv.adj] <- corrected.feat
}
icv.adj.frame.init[,1] <- NULL #What is this doing?

full.fs <- cbind(icv.adj.frame.init, full.fs)

full.fs$LeftMeta_harmonized_icv_adj <- rowSums2(as.matrix(full.fs[,c("ST24CV_harmonized_icv_adjusted", "ST32CV_harmonized_icv_adjusted", 
                                                           "ST40CV_harmonized_icv_adjusted", "ST26CV_harmonized_icv_adjusted")]))

full.fs$RightMeta_harmonized_icv_adj <- rowSums2(as.matrix(full.fs[,c("ST83CV_harmonized_icv_adjusted", "ST91CV_harmonized_icv_adjusted", 
                                                                     "ST99CV_harmonized_icv_adjusted", "ST85CV_harmonized_icv_adjusted")]))

full.fs$RightAnteriorCingulate_harmonized_icv_adj <- rowSums2(as.matrix(full.fs[,c("ST113CV_harmonized_icv_adjusted", "ST73CV_harmonized_icv_adjusted")]))

full.fs$LeftAnteriorCingulate_harmonized_icv_adj <- rowSums2(as.matrix(full.fs[,c("ST54CV_harmonized_icv_adjusted","ST14CV_harmonized_icv_adjusted")]))

full.fs$LeftMeta_harmonized <- rowSums2(as.matrix(full.fs[,c("ST24CV_harmonized", "ST32CV_harmonized", 
                                                                     "ST40CV_harmonized", "ST26CV_harmonized")]))

full.fs$RightMeta_harmonized <- rowSums2(as.matrix(full.fs[,c("ST83CV_harmonized", "ST91CV_harmonized", 
                                                                      "ST99CV_harmonized", "ST85CV_harmonized")]))

full.fs$RightAnteriorCingulate_harmonized <- rowSums2(as.matrix(full.fs[,c("ST113CV_harmonized", "ST73CV_harmonized")]))

full.fs$LeftAnteriorCingulate_harmonized <- rowSums2(as.matrix(full.fs[,c("ST54CV_harmonized","ST14CV_harmonized")]))


full.fs$LeftMeta <- rowSums2(as.matrix(full.fs[,c("ST24CV", "ST32CV", 
                                                             "ST40CV", "ST26CV")]))

full.fs$RightMeta <- rowSums2(as.matrix(full.fs[,c("ST83CV", "ST91CV", 
                                                              "ST99CV", "ST85CV")]))

full.fs$RightAnteriorCingulate <- rowSums2(as.matrix(full.fs[,c("ST113CV", "ST73CV")]))

full.fs$LeftAnteriorCingulate <- rowSums2(as.matrix(full.fs[,c("ST54CV","ST14CV")]))


full.fs$EXAMDATE <- full.fs$EXAMDATE.x
full.fs$VISCODE <- full.fs$VISCODE.x
full.fs$EXAMDATE.x <- full.fs$EXAMDATE.y <- full.fs$VISCODE.x <- full.fs$VISCODE.y <- full.fs$X <- NULL


write.csv(full.fs, "freesurferharmnewadni3_trial_3.csv")

