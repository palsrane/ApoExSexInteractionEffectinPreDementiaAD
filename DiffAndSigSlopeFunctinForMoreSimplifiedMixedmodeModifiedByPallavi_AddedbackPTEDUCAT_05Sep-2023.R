#Changes made in July 2023
# All ApoE (varibles/mentions) were changed to APOE to accomodate changes made in main code.
# Changed background color of all plots to white from grey


#Changes made since last version
# Some how loading of mgcv library before the others is resulting in changes in 
# the syntax for predict in main function and predict.gam function used in Calculate derivative fucntion
# I had to change 'exclude_terms =' from predict function that calculated the 
# main lpmatrix to just 'exclude ='
# And also had to include the 'exclude=' portion in the predict.gam function that was 
# being used inside the CalculateDerivatives function
library(patchwork)
library(stringr)
#Slope change / comparing model predictions
`%notin%` <- Negate(`%in%`)

# New color scheme July 2023

#Calculate difference between two predicitons from the lpmatrix prediction
CalculateDifference <- function(m, xp, r1, r2,c1,c2,age) {
  #linear  predictor matrix
  
  #take difference in rows of data relating to that smooth spline
  X <- xp[r1, ] - xp[r2,]
  # #zero out all columns except for smooth spline of interest
  # cols=c(c1,c2)
  # X[,colnames(X)[colnames(X) %notin% cols]] <- 0
  
  # Following is recommended as per Koychev paper and referenced blog https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
  # But not for us as our parametric coefficients are modeled by group and we can't ignore them.
  # X[,cols[!grepl('s\\(', cols)]] <- 0 #All parametric coefficients 0
  # print(cols[!grepl('s\\(', cols)])
  dif <- X %*% coef(m)
  
  #get row wise standard error and CI's
  se <- sqrt(rowSums((X %*% vcov(m, unconditional = TRUE)) * X))
  crit <- qt(.05/2, df.residual(m), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  ret.data <- data.frame("diff"  = dif,
                         "se"    = se,
                         "upper" = upr,
                         "lower" = lwr)
  
  #find points that are significantly different from 0 (CI does not contain zero)
  sig.init <- c()
  sig.age <-c()
  for(i in 1:nrow(ret.data)) {
    if(0 >= ret.data["lower"][i,] & 0 <= ret.data["upper"][i,]) {
      sig.init <- append(sig.init, 0)
      sig.age <- append (sig.age,NA)
    } else {
      sig.init <- append(sig.init, 1)
      sig.age <- append (sig.age,age[i])
    }
  }
  ret.data$Significant <- factor(sig.init)
  ret.data$Significance <- sig.init
  ret.data$Sig.Age <-sig.age
  ret.data$Age <- age
  return(ret.data)
}

#Plot the difference
PlotDifference <- function(data, ylab, title) {
  gg <- ggplot(data, aes(x=Age, y=diff)) + geom_point(aes(colour=Significant), size=4) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3)
  gg <- gg+ylab(ylab) + labs(title = title) + scale_color_manual(breaks = c("0", "1"), values=c("#E65D54", "#005555"))

  if (length(data$Sig.Age[!is.na(data$Sig.Age)])>0){
    bounds<-data.frame(c(data$Sig.Age[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=1,fill=NA))],data$Sig.Age[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=-1,fill=NA))]))
    names(bounds)="sigage"
    bounds$sig=c(data$diff[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=1,fill=NA))],data$diff[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=-1,fill=NA))])
    gg=gg+ geom_vline(xintercept = bounds$sigage, linetype="dashed",size=3)+
      geom_text(size=20,data=bounds,aes(x=sigage-0.75,y=min(data$lower)+ abs(mean(data$upper)),label=round(sigage,digits=2),angle=90))
  }
  gg<- gg+theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60),title = element_text(size=65))+
    theme(
      panel.background = element_rect(fill = "white", colour = "#6D9EC1",
                                      size = 2, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour ="black")
    )
    
  return(gg)
}

CalculateDerivative <- function(m, eps, pdat, rows, cols, age) {
  
  #Finite differencing with respect to AGE
  pdat_eps <- pdat
  pdat_eps$Agediff <- pdat_eps$Agediff + eps
  
  X0 <- predict.gam(m,exclude = c("s(RID)", "s(RID,Agediff)"), newdata = pdat, type="lpmatrix")
  X1 <- predict.gam(m, exclude = c("s(RID)", "s(RID,Agediff)"),newdata = pdat_eps, type="lpmatrix")
  
  #Only keep spline of interest
  X0 <- X0[rows,]
  X1 <- X1[rows,]
  
  #Finite differencing
  Xp <- (X1 - X0) / eps
  # #zero out all columns except for smooth spline of interest
  # Xp[,colnames(Xp)[colnames(Xp) %notin% cols]] <- 0
  
  # Following is recommended as per Koychev paper and referenced blog https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
  # But my results seem wrong if I do
  # Xp[,cols[!grepl('s\\(', cols)]] <- 0 #All parametric coefficients 0
  # print(cols[!grepl('s\\(', cols)])
  der <- Xp %*% coef(m)
  
  #Standard error and CI
  se <- rowSums(Xp %*% m$Vp * Xp)^.5
  crit <- qt(.975, df.residual(m))
  lower <- der - (se*crit)
  upper <- der + (se*crit)
  ret.data <- data.frame("Derivative" = der,
                         "lower" = lower,
                         "upper" = upper)
  sig.init <- c()
  sig.age <- c()
  for(i in 1:nrow(ret.data)) {
    if(0 >= ret.data["lower"][i,] & 0 <= ret.data["upper"][i,]) {
      sig.init <- append(sig.init, 0)
      sig.age <- append(sig.age,NA)
    } else {
      sig.init <- append(sig.init, 1)
      sig.age <- append(sig.age,age[i])
    }
  }
  
  ret.data$Significant <- factor(sig.init)
  ret.data$Significance <-sig.init
  ret.data$Sig.Age <- sig.age
  ret.data$Age <- age
  return(ret.data)
}

PlotDerivative <- function(data, ylab, title) {
  gg <- ggplot(data, aes(x=Age, y=Derivative)) + geom_point(aes(colour=Significant), size=4) + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3)
  gg <- gg+ylab(ylab) + labs(title = title) + scale_color_manual(breaks = c("0", "1"), values=c("#E65D54", "#005555"))
  if (length(data$Sig.Age[!is.na(data$Sig.Age)])>0){
    bounds<-data.frame(c(data$Sig.Age[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=1,fill=NA))],data$Sig.Age[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=-1,fill=NA))]))
    names(bounds)="sigage"
    bounds$sig=c(data$Derivative[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=1,fill=NA))],data$Derivative[!is.na(data$Sig.Age) & is.na(data.table::shift(data$Sig.Age,n=-1,fill=NA))])
    gg=gg+ geom_vline(xintercept = bounds$sigage, linetype="dashed",size=3)+
      geom_text(size=20,data=bounds,aes(x=sigage-0.75,y=min(data$lower)+ abs(mean(data$upper)),label=round(sigage,digits=2),angle=90))
    
  }
  gg<- gg+
    theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60),title = element_text(size=65))+
    theme(
      panel.background = element_rect(fill = "white", colour = "#6D9EC1",
                                      size = 2, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour ="black")
    )
  return(gg)
}


CompareSplineSigSlope<-function(MOI,ModelType,tmpdd,m1,m2,Covarlist,refSex4m1)
{
  tmpdd1=tmpdd #from MixedModelWMH&ICVeffect_Long&CrossSectionalData.R
  tmpdd1$RID <- factor(tmpdd1$RID)
  
  MOIFilename=str_replace_all(MOI," ","")

  #creating a new sample  new dataset with same values that were used in our predictions
  pdat2 <- expand.grid(Age          = seq(min(tmpdd1$Age), max(tmpdd1$Age), length = 100),
                       APOE         = unique(tmpdd1$APOE),
                       Sex          = unique(tmpdd1$Sex),#comment this if using APOE only effects model
                       RID          = tmpdd1$RID[1]) #Necessary to be in dataset for prediction but coefficient get set to 0
  pdat2$Agediff = pdat2$Age-min(pdat2$Age)
  if (grepl("CDSOB",Covarlist)){pdat2$CDSOB=c(0)}
  if (grepl("WML",Covarlist)){pdat2$WML=c(0)}
  if (grepl("WMH",Covarlist)){pdat2$WMH=c(0)}
  if (grepl("wmlICV",Covarlist)){pdat2$wmlICV=c(mean(tmpdd1$wmlICV,na.rm=T))}
  if (grepl("Centiloid",Covarlist)){pdat2$Centiloid=20}
  if (grepl("PTEDUCAT",Covarlist)){pdat2$PTEDUCAT=c(16)}
  if (grepl("Edudiff",Covarlist)){pdat2$Edudiff=c(0)}
  
  # pdat2$int_sex_apoe=paste0(pdat2$Sex,'.',pdat2$APOE) Do not have this dummy group interaction term in the model any more.
  
  # APOE = unique(tmpdd1$APOE) creats x3 records, 1 for each APOE group 
  # & Sex = unique(tmpdd1$Sex) creates x2, 1 for each gender during the expand.grid  
  # the $int_sex_apoe is just a combo of the two which needs to present when we are
  # using model where that tmpdd1 and m2 itself had that parameter
  
  # #Modified to use my own model from previous code
  # predictions <-data.frame(predict.gam(m2, exclude = c("s(RID)", "s(RID,Age)"), newdata=pdat2))
  # #newdata is probably what I was using as values 
  # names(predictions)='m2_predict'
  # predictions$Sex=pdat2$Sex
  # predictions$Age=pdat2$Age
  # predictions$APOE=pdat2$APOE
  # plot.apoe <- ggplot(predictions, aes(x=Age, y=m2_predict, linetype=Sex)) + geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE) + facet_wrap(~APOE) + ylab(MOI) + labs(title=paste(MOI,"By APOE"))
  # plot.sex <- ggplot(predictions, aes(x=Age, y=m2_predict, colour=APOE)) + geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE) + facet_wrap(~Sex)+ ylab(MOI) + labs(title=paste(MOI,"By Sex"))
  
  print("Here at m2")
  #Using my own model by adding type = 'lpmatrix' to the predictions
  m2_predict_lp<-predict(m2,exclude = c('s(RID)','s(RID,Agediff)'), # don't need c('s(RID)','s(RID,Age)') as s(RID) term is not there to begin with
                         newdata = pdat2,type = 'lpmatrix')
  print("Here at m2")
  age.prediction.frame <- seq(min(tmpdd1$Age), max(tmpdd1$Age), length = 100) 
  print(colnames(m2_predict_lp))
  #Actually we do not really need this as columns not belonging to a particular group are zero in rows belonging to that group anyways
  #This is incorrect for the model without the mock interaction group term int_sex_apoe
  #So code where columns were zeroed out in CalculateDerivative and calculateDifference functions are commented out. 
  f.e2.cols <- colnames(m2_predict_lp)[grepl('F.E3/E2', colnames(m2_predict_lp))]
  f.e3.cols <- colnames(m2_predict_lp)[grepl('F.E3/E3', colnames(m2_predict_lp))]
  f.e4.cols <- colnames(m2_predict_lp)[grepl('F.E3/E4', colnames(m2_predict_lp))]
  
  m.e2.cols <- colnames(m2_predict_lp)[grepl('M.E3/E2', colnames(m2_predict_lp))]
  m.e3.cols <- colnames(m2_predict_lp)[grepl('M.E3/E3', colnames(m2_predict_lp))]
  m.e4.cols <- colnames(m2_predict_lp)[grepl('M.E3/E4', colnames(m2_predict_lp))]
  
  #m2_perdict ipmatrix is created in the same order as tmpdd
  f.e2.rows <- which(pdat2$APOE=="E3/E2" & pdat2$Sex=="F")
  f.e3.rows <- which(pdat2$APOE=="E3/E3" & pdat2$Sex=="F")
  f.e4.rows <- which(pdat2$APOE=="E3/E4" & pdat2$Sex=="F")
  
  m.e2.rows <- which(pdat2$APOE=="E3/E2" & pdat2$Sex=="M")
  m.e3.rows <- which(pdat2$APOE=="E3/E3" & pdat2$Sex=="M")
  m.e4.rows <- which(pdat2$APOE=="E3/E4" & pdat2$Sex=="M")
  
  #The way to print from lpmatrix of preditions----
  #Single group
  XFE4=m2_predict_lp[f.e4.rows,] #rows corresponding to 
  XFE4=data.frame(XFE4%*%coef(m2))
  names(XFE4)="predictions"
  XFE4$Age= seq(min(tmpdd1$Age), max(tmpdd1$Age), length = 100)
  print(ggplot(XFE4, aes(x=Age, y=predictions)) + geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE) )
  rm(XFE4)
  
  # Couple of groups
  X = m2_predict_lp[c(f.e4.rows,f.e2.rows),]
  X=data.frame(X%*%coef(m2))
  names(X)="predictions"
  X$Age= pdat2[c(f.e4.rows,f.e2.rows),c("Age")]
  X$Sex=pdat2[c(f.e4.rows,f.e2.rows),c("Sex")]
  X$APOE=pdat2[c(f.e4.rows,f.e2.rows),c("APOE")]
  print(ggplot(X, aes(x=Age,y=predictions, color = APOE, linetype = Sex)) +
          geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE)+ #se=true or false does not matter because in the pdat2 matrix there is only a single subject per group at any given Age point.
          labs(title=paste(MOI,"prediction within ",dxtext,"for Female E4 and E2"),
               x ="Age (years)", y = MOI))
  rm(X)
  
  #All sub-groups
  X=data.frame(drop (m2_predict_lp %*% coef(m2)))
  names(X)="predictions"
  X$se <- sqrt(rowSums((m2_predict_lp %*% vcov(m2, unconditional = TRUE)) * m2_predict_lp))
  X$crit <- qt(.05/2, df.residual(m2), lower.tail = FALSE)
  X$upr=X$predictions+(X$crit * X$se)
  X$lwr=X$predictions-(X$crit * X$se)
  X$Age=pdat2$Age
  X$Sex=pdat2$Sex
  X$APOE=pdat2$APOE
  print(ggplot(X, aes(x=Age,y=predictions, color = APOE, linetype = Sex)) +
          # geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + #This requires fit and se.fit columns to be present
          geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE)+
          # geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1) + #alpha is transperance level of teh ribbon
          labs(title=paste(MOI,"prediction within ",dxtext),
               x ="Age (years)", y = MOI)+ facet_wrap(vars(Sex)))
  rm(X)
  
  #Starting the process of calculating differences
  diff_fe3.fe2 <- CalculateDifference(m2, m2_predict_lp, f.e3.rows, f.e2.rows, f.e3.cols, f.e2.cols, age.prediction.frame)
  diff_fe3.fe4 <- CalculateDifference(m2, m2_predict_lp, f.e3.rows, f.e4.rows, f.e3.cols, f.e4.cols, age.prediction.frame)
  diff_fe2.fe4 <- CalculateDifference(m2, m2_predict_lp, f.e2.rows, f.e4.rows, f.e2.cols, f.e4.cols, age.prediction.frame)
  
  diff_me3.me2 <- CalculateDifference(m2, m2_predict_lp, m.e3.rows, m.e2.rows, m.e3.cols, m.e2.cols, age.prediction.frame)
  diff_me3.me4 <- CalculateDifference(m2, m2_predict_lp, m.e3.rows, m.e4.rows, m.e3.cols, m.e4.cols, age.prediction.frame)
  diff_me2.me4 <- CalculateDifference(m2, m2_predict_lp, m.e2.rows, m.e4.rows, m.e2.cols, m.e4.cols, age.prediction.frame)
  
  diff_me2.fe2 <- CalculateDifference(m2, m2_predict_lp, m.e2.rows, f.e2.rows, m.e2.cols, f.e2.cols, age.prediction.frame)
  diff_me3.fe3 <- CalculateDifference(m2, m2_predict_lp, m.e3.rows, f.e3.rows, m.e3.cols, f.e3.cols, age.prediction.frame)
  diff_me4.fe4 <- CalculateDifference(m2, m2_predict_lp, m.e4.rows, f.e4.rows, m.e4.cols, f.e4.cols, age.prediction.frame)
  
  
  plot_diff.f.e3.e2 <- PlotDifference(diff_fe3.fe2, paste(MOI,"Difference"), "Female E3 - Female E2")
  # print(diff_fe2.fe3$Sig.Age)
  plot_diff.f.e3.e4 <- PlotDifference(diff_fe3.fe4, paste(MOI,"Difference"), "Female E3 - Female E4")
  # print(diff_fe3.fe4$Sig.Age)
  plot_diff.f.e2.e4 <- PlotDifference(diff_fe2.fe4, paste(MOI,"Difference"), "Female E2 - Female E4")
  # print(diff_fe2.fe4$Sig.Age)
  jpeg(paste0(MOIFilename,ModelType,"ModelBetweenFemaleDiff.jpg"),width=6000,height = 1500)
  print(plot_diff.f.e3.e2 + plot_diff.f.e3.e4 + plot_diff.f.e2.e4)
  dev.off()
  
  
  plot_diff.m.e3.e2 <- PlotDifference(diff_me3.me2, paste(MOI,"Difference"), "Male E3 - Male E2")
  # print(diff_me2.me3$Sig.Age)
  plot_diff.m.e3.e4 <- PlotDifference(diff_me3.me4, paste(MOI,"Difference"), "Male E3 - Male E4")
  # print(diff_me3.me4$Sig.Age)
  plot_diff.m.e2.e4 <- PlotDifference(diff_me2.me4, paste(MOI,"Difference"), "Male E2 - Male E4")
  # print(diff_me2.me4$Sig.Age)
  jpeg(paste0(MOIFilename,ModelType,"ModelBetweenMaleDiff.jpg"),width=6000,height = 1500)
  print(plot_diff.m.e3.e2 + plot_diff.m.e3.e4 + plot_diff.m.e2.e4)
  dev.off()
  
  plot_diff.me2.fe2 <- PlotDifference(diff_me2.fe2, paste(MOI,"Difference"), "Male E2 - Female E2")
  # print(diff_me2.fe2$Sig.Age)
  plot_diff.me3.fe3 <- PlotDifference(diff_me3.fe3, paste(MOI,"Difference"), "Male E3 - Female E3")
  # print(diff_me3.fe3$Sig.Age)
  plot_diff.me4.fe4 <- PlotDifference(diff_me4.fe4, paste(MOI,"Difference"), "Male E4 - Female E4")
  # print(diff_me4.fe4$Sig.Age)
  jpeg(paste0(MOIFilename,ModelType,"ModelWithinAPOEBetweenSexDiff.jpg"),width=6000,height = 1500)
  print(plot_diff.me3.fe3 + plot_diff.me2.fe2 + plot_diff.me4.fe4)
  dev.off()
  
  
  ##### Calculating derivative of Age smooths across groups
  
  ### 
  deriv_fe2 <- CalculateDerivative(m2, 0.001, pdat2, f.e2.rows, f.e2.cols, age.prediction.frame)
  deriv_fe3 <- CalculateDerivative(m2, 0.001, pdat2, f.e3.rows, f.e3.cols, age.prediction.frame)
  deriv_fe4 <- CalculateDerivative(m2, 0.001, pdat2, f.e4.rows, f.e4.cols, age.prediction.frame)
  
  deriv_me2 <- CalculateDerivative(m2, 0.001, pdat2, m.e2.rows, m.e2.cols, age.prediction.frame)
  deriv_me3 <- CalculateDerivative(m2, 0.001, pdat2, m.e3.rows, m.e3.cols, age.prediction.frame)
  deriv_me4 <- CalculateDerivative(m2, 0.001, pdat2, m.e4.rows, m.e4.cols, age.prediction.frame)
  
  plot_deriv_fe2 <- PlotDerivative(deriv_fe2, paste(MOI,"Derivative"), "Female E2")
  plot_deriv_fe3 <- PlotDerivative(deriv_fe3, paste(MOI,"Derivative"), "Female E3")
  plot_deriv_fe4 <- PlotDerivative(deriv_fe4, paste(MOI,"Derivative"), "Female E4")
  
  jpeg(paste0(MOIFilename,ModelType,"ModelFemaleDerivatives.jpg"),width=6000,height = 1500)
  print(plot_deriv_fe3 + plot_deriv_fe2 + plot_deriv_fe4)
  dev.off()
  
  plot_deriv_me2 <- PlotDerivative(deriv_me2, paste(MOI,"Derivative"), "Male E2")
  plot_deriv_me3 <- PlotDerivative(deriv_me3, paste(MOI,"Derivative"), "Male E3")
  plot_deriv_me4 <- PlotDerivative(deriv_me4, paste(MOI,"Derivative"), "Male E4")
  jpeg(paste0(MOIFilename,ModelType,"ModelMaleDerivatives.jpg"),width=6000,height = 1500)
  print(plot_deriv_me3 + plot_deriv_me2 + plot_deriv_me4)
  dev.off()
  
  # m1----
  #creating a new sample  new dataset with same values that were used in our predictions
  pdat1 <- expand.grid(Age          = seq(min(tmpdd1$Age), max(tmpdd1$Age), length = 100),
                       APOE         = unique(tmpdd1$APOE),
                       Sex          = refSex4m1,
                       RID          = tmpdd1$RID[1]) #Necessary to be in dataset for prediction but coefficient get set to 0
  pdat1$Agediff = pdat1$Age-min(pdat1$Age)
  if (grepl("CDSOB",Covarlist)){pdat1$CDSOB=c(0)}
  if (grepl("WML",Covarlist)){pdat1$WML=c(0)}
  if (grepl("WMH",Covarlist)){pdat1$WMH=c(0)}
  if (grepl("wmlICV",Covarlist)){pdat1$wmlICV=c(mean(tmpdd1$wmlICV,na.rm=T))}
  if (grepl("Centiloid",Covarlist)){pdat1$Centiloid=20}
  if (grepl("PTEDUCAT",Covarlist)){pdat1$PTEDUCAT=c(16)}
  if (grepl("Edudiff",Covarlist)){pdat1$Edudiff=c(0)}
  
  #Using my own model by adding type = 'lpmatrix' to the predictions
  print("Here at m1")
  m1_predict_lp<-predict(m1,exclude = c('s(RID)','s(RID,Agediff)'), 
                         newdata = pdat1, type = 'lpmatrix')
  age.prediction.frame <- seq(min(tmpdd1$Age), max(tmpdd1$Age), length = 100) 
  
  print("Here at m1")
  #columns belonging to each APOE group
  e2.cols <- colnames(m1_predict_lp)[grepl('E3/E2', colnames(m1_predict_lp))]
  e3.cols <- colnames(m1_predict_lp)[grepl('E3/E3', colnames(m1_predict_lp))]
  e4.cols <- colnames(m1_predict_lp)[grepl('E3/E4', colnames(m1_predict_lp))]
  # e2.cols <- colnames(m1_predict_lp)[grepl('E3/E2', colnames(m1_predict_lp)) & grepl('Age', colnames(m1_predict_lp))]
  # e3.cols <- colnames(m1_predict_lp)[grepl('E3/E3', colnames(m1_predict_lp)) & grepl('Age', colnames(m1_predict_lp))]
  # e4.cols <- colnames(m1_predict_lp)[grepl('E3/E4', colnames(m1_predict_lp)) & grepl('Age', colnames(m1_predict_lp))]
  # Somehow, not having the covariate terms columns doesn't make a difference
  
  #m1_perdict ipmatrix is created in the same order as tmpdd
  e2.rows <- which(pdat1$APOE=="E3/E2")
  e3.rows <- which(pdat1$APOE=="E3/E3")
  e4.rows <- which(pdat1$APOE=="E3/E4")
  
  #The way to print from lpmatrix of preditions----
  #Single group
  XE4=m1_predict_lp[e4.rows,] #rows corresponding to 
  XE4=data.frame(XE4%*%coef(m1))
  names(XE4)="predictions"
  XE4$Age= age.prediction.frame
  print(ggplot(XE4, aes(x=Age, y=predictions)) + geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE) )
  rm(XE4)
  
  # Couple of groups
  X = m1_predict_lp[c(e4.rows,e2.rows),]
  X=data.frame(X%*%coef(m1))
  names(X)="predictions"
  X$Age= pdat1[c(e4.rows,e2.rows),c("Age")]
  X$APOE=pdat1[c(e4.rows,e2.rows),c("APOE")]
  print(ggplot(X, aes(x=Age,y=predictions, color = APOE)) +
          geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE)+ #se=true or false does not matter because in the pdat1 matrix there is only a single subject per group at any given Age point.
          labs(title=paste(MOI,"prediction within ",dxtext,"for E4 and E2"),
               x ="Age (years)", y = MOI))
  rm(X)
  
  #All sub-groups
  X=data.frame(drop (m1_predict_lp %*% coef(m1)))
  names(X)="predictions"
  X$se <- sqrt(rowSums((m1_predict_lp %*% vcov(m1, unconditional = TRUE)) * m1_predict_lp))
  X$crit <- qt(.05/2, df.residual(m1), lower.tail = FALSE)
  X$upr=X$predictions+(X$crit * X$se)
  X$lwr=X$predictions-(X$crit * X$se)
  X$Age=pdat1$Age
  X$APOE=pdat1$APOE
  print(ggplot(X, aes(x=Age,y=predictions, color = APOE)) +
          # geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + #This requires fit and se.fit columns to be present
          geom_smooth(method = "gam", formula = y ~ s(x), se=FALSE)+
          # geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1) + #alpha is transperance level of teh ribbon
          labs(title=paste(MOI,"prediction within ",dxtext),
               x ="Age (years)", y = MOI))
  rm(X)
  
  #Starting the process of calculating differences
  diff_e3.e2 <- CalculateDifference(m1, m1_predict_lp, e3.rows, e2.rows, e3.cols, e2.cols, age.prediction.frame)
  diff_e3.e4 <- CalculateDifference(m1, m1_predict_lp, e3.rows, e4.rows, e3.cols, e4.cols, age.prediction.frame)
  diff_e2.e4 <- CalculateDifference(m1, m1_predict_lp, e2.rows, e4.rows, e2.cols, e4.cols, age.prediction.frame)
  
  print("plot_diff.e2.e3")
  plot_diff.e3.e2 <- PlotDifference(diff_e3.e2, paste(MOI,"Difference"), "E3 - E2")
  # print(diff_fe2.fe3$Sig.Age)
  plot_diff.e3.e4 <- PlotDifference(diff_e3.e4, paste(MOI,"Difference"), "E3 - E4")
  # print(diff_fe3.fe4$Sig.Age)
  plot_diff.e2.e4 <- PlotDifference(diff_e2.e4, paste(MOI,"Difference"), "E2 - E4")
  # print(diff_fe2.fe4$Sig.Age)
  
  jpeg(paste0(MOIFilename,ModelType,"ModelBetweenAPOEDiff.jpg"),width=6000,height = 1500)
  print(plot_diff.e3.e2 + plot_diff.e3.e4 + plot_diff.e2.e4)
  dev.off()
  
  
  ##### Calculating derivative of Age smooths across groups ### 
  deriv_e2 <- CalculateDerivative(m1, 0.001, pdat1, e2.rows, e2.cols, age.prediction.frame)
  deriv_e3 <- CalculateDerivative(m1, 0.001, pdat1, e3.rows, e3.cols, age.prediction.frame)
  deriv_e4 <- CalculateDerivative(m1, 0.001, pdat1, e4.rows, e4.cols, age.prediction.frame)
  
  plot_deriv_e2 <- PlotDerivative(deriv_e2, paste(MOI,"Derivative"), "E2")
  plot_deriv_e3 <- PlotDerivative(deriv_e3, paste(MOI,"Derivative"), "E3")
  plot_deriv_e4 <- PlotDerivative(deriv_e4, paste(MOI,"Derivative"), "E4")
  
  jpeg(paste0(MOIFilename,ModelType,"ModelAPOEDerivatives.jpg"),width=6000,height = 1500)
  print(plot_deriv_e3 + plot_deriv_e2 + plot_deriv_e4)
  dev.off()
  
  rm(tmpdd,tmpdd1,m1,m2,Covarlist)
  
  allenv<-mget(ls())
  
  saveRDS(allenv,paste0(MOIFilename,ModelType,"ModelallDerivSigSlopeObjects.rds"))
  return(list("deriv_e2"=deriv_e2,"deriv_e3"=deriv_e3,"deriv_e4"=deriv_e4,"deriv_fe2"=deriv_fe2,"deriv_me2"=deriv_me2,"deriv_fe3"=deriv_fe3,"deriv_me3"=deriv_me3,"deriv_fe4"=deriv_fe4,"deriv_me4"=deriv_me4))
}
