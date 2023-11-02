#Changes made Sep 2023
# Added back PTEDUCAT (i.e. education level) to the model because the differences between the APOE groups and between the sex groups were significant.
# Adding NwithRepeatedMeasures, AvgTimeBetweenVisitsMean, and TotalFollowupDuration to individual dataset demographics 

#Changes made Aug 2023
# CDSOBAgeBin for ref for other agebin distributions (Changed how different elements of the plot are displayed and also the plot background)
# Added function for plotting data scatter and a function for plotting predictions
# Decided to combine the whole and the sex-wise predictions in one variable in the plotting function so that a uniform y axis is applied.
# Created a new table for participant records, and also a table for whole demographics
# Added a new function to create demographics table for each dataset
# Added a new function to plot data by age bins. It plots two plots: 1 for APOE distribution and 1 for Sex distribution for each dataset
# Added an addition section at the end to generate box plots that look similar to the ones above but for APOE and sex distribution in the overall dataset

#Changes made in July 2023
# -Changed the way plots are produced by
#   changing colors to black/gray for E3, Light/dark cyan for E2, and pink/red for E4
#   and also adding the significant change duration to the main plot itself with the 
#   darker of the colors representing the significant durations
# -Changed ApoE to APOE to correctly represent the allele and not the protein
# -Changed to latest data files
# -Changed to use to teh earliest of the EXAMDATE and USERDATE as the actual EXAMDATE
#   in cases where both were present. 
# -Added additional APOE info that was made available for ADNI3 data
library(mgcv)
library(plyr)
library(nlme)
library(lme4)
library(ggplot2)
library(plsdepot)
library(stringr)
library(tidymv)
library(olsrr)
library(MatchIt)
library(tidyr)
library(patchwork)
library(ADNIMERGE)
library(grid)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(VennDiagram)
library(venn)
library(ggpolypath)
# install.packages("remotes")
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
library(magick)
library(ggplotify)

rm(list = ls())
`%notin%` <- Negate(`%in%`)
setwd("C:/Projects/APOESex")

# cdr.sob=0
dx = c('CN','MCI') #Not including AD as very almost no apoe2 carriers with AD
dxtext=dx[1]
if(length(dx)>1){
  for (i in 2:length(dx)){dxtext= paste(dxtext, dx[i], sep=" & ")} #added the '_' as separator for file names
}

m1zscores<-matrix(nrow=0,ncol=5)
colnames(m1zscores)<-c("Age","APOE","z","Significance","Measure")

m2zscores<-matrix(nrow=0,ncol=6)
colnames(m2zscores)<-c("Age","APOE","Sex","z","Significance","Measure")

demog_info<-as.data.frame(matrix(nrow=0,ncol=25))
colnames(demog_info)=c("Measure","APOESexgroup","NParticipants","NTotRecords","NCUBaselineRecords","NMCIBaselineRecords","AgeBaselineMean","AgeBaselineStd","AgeMin","AgeMax","EduMean","EduStd","PercentWhite","Nwhite","PercentAsian","NAsian","PercentBlack","NBlack","PercentHispanic","NHispanic","AgeMean","AgeStd","NPartWRepMeasures","NCURecords","NMCIRecords")

Participant_Records=data.frame(RID=0)

#Age limits for plots
x_min=55 #minimum Age plotted
x_max=92 #95 #maximum Age plotted

#Colors for plots
cols <- c("0E3/E2"="lightcyan","1E3/E2"="darkcyan","0E3/E3"="darkgray","1E3/E3"="black","0E3/E4"="pink","1E3/E4"="red") 
# cols_datascatter <- c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red")
# cols_datascatter <- c("E3/E2"="cyan","E3/E3"="gray","E3/E4"="pink") 

plot_theme =theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1", size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour ="black"))
  

# # Old code If I want to test Sex effects----
# m1 <- bam(formula = reg ~ s(Age, bs='cs', by= Sex) + 
#             s(RID,Age, bs="re")+ #Accounting for random intercept and random slope by individual subjects' repeated measures
#             Sex+
#             CDSOB:Sex + PTEDUCAT:Sex, 
#           data = tmpdd, method = "REML")
# m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Age)'), 
#                           values=list(PTEDUCAT=c(16), CDSOB=cdr.sob, 
#                                                                 +                                     AmyloidPos=0, 
#                                                                 +                                     WML=0,Centiloid=20,RID=NULL))
# m1_predict$Sex=factor(m1_predict$Sex)
# m1_predict$Sex = relevel(m1_predict$Sex,ref='M')
# ggp1=ggplot(m1_predict, aes(Age, fit, color = Sex)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("PACC prediction within ",dxtext,'with Sex effects'), 
#        x ="Age (years)", y = 'PACC')
# ggp1 #----

plotdatascatter <-function(tmpdd,plt_labs){
  tmpdd$APOE_Colors=paste0('1',tmpdd$APOE)
  ggp=ggplot(tmpdd, aes(x = Age, y = reg,
                    group = RID, #Switching to Cross sectional analysis
                    colour = APOE_Colors)) +
    geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) + #,size=1.1
    scale_color_manual(values=cols) + scale_x_continuous(limits = c(x_min, x_max)) +
    geom_point()+ #aes(size=3)
    plt_labs+
    theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 2, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="black"),
      strip.background = element_blank(),
      strip.placement = "outside")  
  print(paste("Agerange",min(tmpdd$Age),max(tmpdd$Age)))
  return(ggp)
}

plot_m1andm2predicts <- function(m1_predict,m2_predict){
  tmp=m1_predict
  tmp$Sex="Total"
  tmp=rbind(tmp,m2_predict)
  tmp$Sex<-relevel(factor(tmp$Sex),ref = 'Total')
  tmp$APOE_lineCol=paste0('0',tmp$APOE)
  SignificanceBounds<-as.data.frame(matrix(nrow=0,ncol=3))
  names(SignificanceBounds)=c("APOE","Sex","Age")
  
  #Marking one continuously significant-change section as one group to draw a single line through it which does not connect to the next significant section
  tmp$grp=0
  grpcount=0
  for (tmpSex in levels(tmp$Sex)) {
    for (tmpAPOE in levels(tmp$APOE)) {
      idx<- c(0,diff(tmp$Significance[tmp$APOE==tmpAPOE & tmp$Sex==tmpSex])) #this creates an array of same length as tmp for givensex and APOE and marks whereever the significance chges with a non-zero number
      i2 <- c(1,which(idx != 0), length(tmp$Significance[tmp$APOE==tmpAPOE & tmp$Sex==tmpSex])+1)
      tmp$grp[tmp$APOE==tmpAPOE & tmp$Sex==tmpSex]=rep((grpcount+1):(grpcount+length(diff(i2))), diff(i2))
      grpcount=max(tmp$grp[tmp$APOE==tmpAPOE & tmp$Sex==tmpSex])
      print(paste(tmpSex,tmpAPOE,tmp$Age[tmp$APOE==tmpAPOE & tmp$Sex==tmpSex][idx!=0]))#Ages where the plots change significance levels
      for (count in unique(tmp$grp[tmp$Significance==1 & tmp$APOE==tmpAPOE & tmp$Sex==tmpSex])){
        SignificanceBounds[nrow(SignificanceBounds)+1,]<-c(tmpAPOE,tmpSex,min(tmp$Age[tmp$grp==count]))
        SignificanceBounds[nrow(SignificanceBounds)+1,]<-c(tmpAPOE,tmpSex,max(tmp$Age[tmp$grp==count]))
      }
      
    }
    
  }
  SignificanceBounds$Age=as.numeric(SignificanceBounds$Age)
  SignificanceBounds$VLineColor=paste0("1",SignificanceBounds$APOE)
  SignificanceBounds$Sex<-relevel(factor(SignificanceBounds$Sex),ref = 'Total')
  SignificanceBounds$APOE<-relevel(factor(SignificanceBounds$APOE),ref = 'E3/E3')
  
  ggp=ggplot(tmp, aes(x = Age, y = fit, colour = APOE_lineCol)) + #, size=Significance , linetype=Sex
    geom_line(aes(size=4))+ geom_smooth_ci(size = 4, ci_alpha=0.1) + scale_color_manual(values=cols) +  
    # geom_point(data=tmp[tmp$Significance!=0,],aes(colour = SigAPOE,size=4))+
    geom_line(data=tmp[tmp$Significance!=0,],aes(colour = SigAPOE,size=4, group=grp))+
    plt_labs + scale_x_continuous(limits = c(x_min, x_max))+ 
    facet_wrap(vars(Sex),labeller = labeller(Sex=c(Total="Full Dataset", M="Males",F="Females"))) +
    theme(legend.position="none",axis.text = element_text(size = 60), 
          axis.title = element_text(size = 60),strip.text = element_text(size =60),
          panel.background = element_rect(fill = "white", colour = "#6D9EC1", size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="white"), #changed from "black"
          strip.background = element_blank(),
          plot.margin = unit(c(1,0.25,1,3),"inches"))+
    labs(x=NULL,y=NULL) #making x and Y axis labels null for the sake of creating a uniform display
  
  YaxisBottom=layer_scales(ggp)$y$range$range[1]
  ggp=ggp+ geom_vline(data=SignificanceBounds,aes(xintercept = Age,colour =VLineColor,size=2), linetype="dashed")+
    geom_text_repel(size=15,data=SignificanceBounds,aes(x=Age-0.75,y=YaxisBottom,label=round(Age,digits=2),angle=90,colour =VLineColor),nudge_x = 0.5,nudge_y = 0) #Changes geom_text to geom_text_repel to avoid overlaps
  
  
  
  return(ggp)
  
}

Save_demog_distribution_withintmpdd<-function(tmpdd,measure){
  tmpdd1=tmpdd[,c("RID","EXAMDATE","Age","APOE","Sex","DXSUM","NumOfScans","PTETHCAT","PTRACCAT","PTEDUCAT","AvgYrsBtnVisits","TotalFollowupYrs")]
  tmpdd1=tmpdd1[order(tmpdd1$RID,tmpdd1$EXAMDATE),]
  tmpbaseline=tmpdd1[!duplicated(tmpdd1$RID),]
  tmpbaseline <- transform(tmpbaseline, WhiteOrNot= ifelse(PTRACCAT==5, 1, 0))
  tmpbaseline <- transform(tmpbaseline, BlackOrNot= ifelse(PTRACCAT==4, 1, 0))
  tmpbaseline <- transform(tmpbaseline, AsianOrNot= ifelse(PTRACCAT==2, 1, 0))
  tmpbaseline <- transform(tmpbaseline, HispanicOrNot= ifelse(PTETHCAT==1, 1, 0))
  
  ret_demog_info=as.data.frame(matrix(nrow=1,ncol = 12))
  ret_demog_info[1,]=c(paste0(measure," dataset"),paste0("Total N=",length(tmpbaseline$RID),", N longitudinal records=",length(tmpbaseline$RID[tmpbaseline$NumOfScans>1])," (",round(length(tmpbaseline$RID[tmpbaseline$NumOfScans>1])*100/length(tmpbaseline$RID),1),"%), Visits/repeated participants =",min(tmpbaseline$NumOfScans[tmpbaseline$NumOfScans>1]),"-",max(tmpbaseline$NumOfScans[tmpbaseline$NumOfScans>1]),", Avg yrs betn visits=",round(mean(tmpbaseline$AvgYrsBtnVisits[tmpbaseline$NumOfScans>1]),2),", Total followup duration=",min(tmpbaseline$TotalFollowupYrs[tmpbaseline$NumOfScans>1]),"-",max(tmpbaseline$TotalFollowupYrs[tmpbaseline$NumOfScans>1])," yrs"),"","","","","","","","","","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("APOE groups",paste0("E3/E3 (N=",length(unique(tmpdd1$RID[tmpdd1$APOE=="E3/E3"])),")"),"","",paste0("E3/E2 (N=",length(unique(tmpdd1$RID[tmpdd1$APOE=="E3/E2"])),")"),"","",paste0("E3/E4 (N=",as.character(length(unique(tmpdd1$RID[tmpdd1$APOE=="E3/E4"]))),")"),"","","E3vsE2 \n M:F","E3vsE4 \n M:F")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("","Female","Male","p","Female","Male","p","Female","Male","p","p","p")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("N # Participants",
                                            length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),"",
                                            length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),"",
                                            length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),"",
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE!="E3/E4",c("Sex","APOE")]))$p,digits=2),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE!="E3/E2",c("Sex","APOE")]))$p,digits=2))
  ret_demog_info[nrow(ret_demog_info)+1,]=c("Total # Records",
                                            length(tmpdd1$RID[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="F"]),length(tmpdd1$RID[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="M"]),"",
                                            length(tmpdd1$RID[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="F"]),length(tmpdd1$RID[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="M"]),"",
                                            length(tmpdd1$RID[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="F"]),length(tmpdd1$RID[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="M"]),"","","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("Age (years; baseline);\n mean(SD) [Range]",
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="F"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="F"]),2),"]"),
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="M"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E3"&tmpdd1$Sex=="M"]),2),"]"),
                                            format(t.test(Age ~ Sex, data = tmpbaseline[tmpbaseline$APOE=="E3/E3",])$p.value,digits=2),
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="F"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="F"]),2),"]"),
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="M"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E2"&tmpdd1$Sex=="M"]),2),"]"),
                                            format(t.test(Age ~ Sex, data = tmpbaseline[tmpbaseline$APOE=="E3/E2",])$p.value,digits=2),
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="F"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="F"]),2),"]"),
                                            paste0(round(mean(tmpbaseline$Age[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"(",round(stats::sd(tmpbaseline$Age[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),")\n[",round(min(tmpdd1$Age[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="M"]),2)," - ",round(max(tmpdd1$Age[tmpdd1$APOE=="E3/E4"&tmpdd1$Sex=="M"]),2),"]"),
                                            format(t.test(Age ~ Sex, data = tmpbaseline[tmpbaseline$APOE=="E3/E4",])$p.value,digits=2),
                                            "","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("Education (years);\n mean(SD) [Range]",
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),1),"]"),
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),1),"]"),
                                            format(t.test(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"],tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"])$p.value,digits=2),
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),1),"]"),
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),1),"]"),
                                            format(t.test(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"],tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"])$p.value,digits=2),
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),1),"]"),
                                            paste0(round(mean(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),1),"(",round(stats::sd(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),1),")\n[",round(min(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),1)," - ",round(max(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),1),"]"),
                                            format(t.test(tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"],tmpbaseline$PTEDUCAT[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"])$p.value,digits=2),
                                            "","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("%White (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E3",c("Sex","WhiteOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E2",c("Sex","WhiteOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==5]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==5]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E4",c("Sex","WhiteOrNot")]))$p,scientific = FALSE,digits=2),
                                            "","") 
  #PTRACCAT: 1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander;
  #          4=Black or African American; 5=White; 6=More than one race; 7=Unknown
  #PTETHCAT: 1=Hispanic or Latino; 2=Not Hispanic or Latino; 3=Unknown
  ret_demog_info[nrow(ret_demog_info)+1,]=c("% Black/African American (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E3",c("Sex","BlackOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E2",c("Sex","BlackOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==4]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==4]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E4",c("Sex","BlackOrNot")]))$p,scientific = FALSE,digits=2),
                                            "","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("% Asian (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E3",c("Sex","AsianOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E2",c("Sex","AsianOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTRACCAT==2]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTRACCAT==2]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E4",c("Sex","AsianOrNot")]))$p,scientific = FALSE,digits=2),
                                            "","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("% Hispanic or Latinx (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E3",c("Sex","HispanicOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E2",c("Sex","HispanicOrNot")]))$p,scientific = FALSE,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$PTETHCAT==1]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$PTETHCAT==1]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E4",c("Sex","HispanicOrNot")]))$p,scientific = FALSE,digits=2),
                                            "","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("% CU at baseline (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"]),")"),
                                            "",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"]),")"),
                                            "",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="CN"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="CN"]),")"),
                                            "","","")
  ret_demog_info[nrow(ret_demog_info)+1,]=c("% MCI at baseline (N)",
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E3"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E3",c("Sex","DXSUM")]))$p,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E2"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E2",c("Sex","DXSUM")]))$p,digits=2),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="F"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            paste0(round(length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"])*100/length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"]),2),"% (",length(tmpbaseline$RID[tmpbaseline$APOE=="E3/E4"&tmpbaseline$Sex=="M"&tmpbaseline$DXSUM=="MCI"]),")"),
                                            format(fisher.test(table(tmpbaseline[tmpbaseline$APOE=="E3/E4",c("Sex","DXSUM")]))$p,digits=2),
                                            "","")
  write.csv(ret_demog_info,paste0(measure,"_demographics.csv"))
  rm(tmpdd1,tmpbaseline,ret_demog_info)
  
  
}

demog_distribution<-function(tmpdd,measure){
  # c("Measure","APOESexgroup","NTotRecords","NParticipants","NPartWRepMeasures",
  # "NCURecords","NMCIRecords","AgeMin","AgeMax","AgeMean","AgeStd","AgeBaselineMean","AgeBaselineStd","EduMean","EduStd",
  # "PercentWhite","PercentBlack","PercentHispanic")
  tmpdd1=tmpdd[,c("RID","EXAMDATE","Age","APOE","Sex","DXSUM","NumOfScans","PTETHCAT","PTRACCAT","PTEDUCAT")]
  ret_demog_info=as.data.frame(demog_info)
  ret_demog_info=ret_demog_info[0,]

  for (AP in c("Total","E3/E3","E3/E2","E3/E4")) {
    for (Sx in c("Total","F","M")) {
      if(AP=="Total") {
        if(Sx=="Total") demog_data=tmpdd1 else demog_data=tmpdd1[tmpdd1$Sex==Sx,]
      }
      if(AP!="Total") {
        if(Sx=="Total") demog_data=tmpdd[tmpdd1$APOE==AP,] else demog_data=tmpdd1[tmpdd1$APOE==AP & tmpdd1$Sex==Sx,]
      }
      tmp_demog_info=as.data.frame(demog_info)
      tmp_demog_info=tmp_demog_info[0,]
      tmp_demog_info=tmp_demog_info[1,]
      demog_data=demog_data[order(demog_data$RID,demog_data$EXAMDATE),]
      tmp_demog_info$Measure=measure
      tmp_demog_info$APOESexgroup=paste0(AP,Sx)
      tmp_demog_info$NTotRecords=length(demog_data$RID)
      tmp_demog_info$NParticipants=length(unique(demog_data$RID))
      tmp_demog_info$NPartWRepMeasures=length(unique(demog_data$RID[demog_data$NumOfScans>1]))
      tmp_demog_info$NCURecords=length(demog_data$RID[demog_data$DXSUM=="CN"])
      tmp_demog_info$NMCIRecords=length(demog_data$RID[demog_data$DXSUM=="MCI"])
      tmp_demog_info$NCUBaselineRecords=length(demog_data$RID[demog_data$DXSUM=="CN" & !duplicated(demog_data$RID)])
      tmp_demog_info$NMCIBaselineRecords=length(demog_data$RID[demog_data$DXSUM=="MCI"& !duplicated(demog_data$RID)])
      tmp_demog_info$AgeMin=min(demog_data$Age)
      tmp_demog_info$AgeMax=max(demog_data$Age)
      tmp_demog_info$AgeMean=mean(demog_data$Age)
      tmp_demog_info$AgeStd=stats::sd(demog_data$Age)
      tmp_demog_info$AgeBaselineMean=mean(demog_data$Age[!duplicated(demog_data$RID)])
      tmp_demog_info$AgeBaselineStd=stats::sd(demog_data$Age[!duplicated(demog_data$RID)])
      tmp_demog_info$EduMean=mean(demog_data$PTEDUCAT[!duplicated(demog_data$RID)])
      tmp_demog_info$EduStd=stats::sd(demog_data$PTEDUCAT[!duplicated(demog_data$RID)])
      # if(AP=="Total" & Sx=="Total") {
      #PTRACCAT: 1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander;
      #          4=Black or African American; 5=White; 6=More than one race; 7=Unknown
      #PTETHCAT: 1=Hispanic or Latino; 2=Not Hispanic or Latino; 3=Unknown
      tmp_demog_info$NWhite =length(unique(demog_data$RID[demog_data$PTRACCAT==5]))
      tmp_demog_info$PercentWhite =(length(unique(demog_data$RID[demog_data$PTRACCAT==5]))*100)/length(unique(demog_data$RID))
      tmp_demog_info$NAsian =length(unique(demog_data$RID[demog_data$PTRACCAT==2]))
      tmp_demog_info$PercentAsian =(length(unique(demog_data$RID[demog_data$PTRACCAT==2]))*100)/length(unique(demog_data$RID))
      tmp_demog_info$NBlack =length(unique(demog_data$RID[demog_data$PTRACCAT==4]))
      tmp_demog_info$PercentBlack =(length(unique(demog_data$RID[demog_data$PTRACCAT==4]))*100)/length(unique(demog_data$RID))
      tmp_demog_info$NHispanic =length(unique(demog_data$RID[demog_data$PTETHCAT==1]))
      tmp_demog_info$PercentHispanic =(length(unique(demog_data$RID[demog_data$PTETHCAT==1]))*100)/length(unique(demog_data$RID))
      # }
      ret_demog_info=rbind(ret_demog_info,tmp_demog_info)
      
    }
    
  }
  
  return(ret_demog_info)
    
}

PrintN <- function(dset){
  print(paste0("N Total ", dim(dset)[1]))
  print(paste0("F E3 ",dim(dset[dset$Sex=='F' & dset$APOE=='E3/E3',])[1]))
  print(paste0("F E2 ",dim(dset[dset$Sex=='F' & dset$APOE=='E3/E2',])[1]))
  print(paste0("F E4 ",dim(dset[dset$Sex=='F' & dset$APOE=='E3/E4',])[1]))
  print(paste0("M E3 ", dim(dset[dset$Sex=='M' & dset$APOE=='E3/E3',])[1]))
  print(paste0("M E2 ", dim(dset[dset$Sex=='M' & dset$APOE=='E3/E2',])[1]))
  print(paste0("M E4 ", dim(dset[dset$Sex=='M' & dset$APOE=='E3/E4',])[1]))}

PrintAmyloidPosN <- function(dset){
  print("Ignore amyloid positivity as different datasets are not integrated")
  dset=dset[order(dset$RID,dset$EXAMDATE),]
  dset1=dset[!duplicated(dset$RID),]
  
  print(paste0("Group TotalN %AmyloidPos AgeMean AgeStdDev EduMean EduStdDev CDSOBMean CDSOBStdDev CN:MCI %MCI %MCIbyAPOE"))
  print(paste0("Total ",length(dset$RID),' ',round(length(dset$RID[dset$AmyloidPos==1])*100/length(dset$RID),digits=2),' ',round(mean(dset$Age,na.rm=TRUE),digits=2),
               ' ',round(stats::sd(dset$Age,na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT,na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT,na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB,na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB,na.rm=TRUE),digits=2),' ',length(dset$RID[dset$DXSUM=='CN']),':',length(dset$RID[dset$DXSUM=='MCI']),' ',round(length(dset$RID[dset$DXSUM=='MCI'])*100/length(dset$RID),digits=2)))
  print(paste0("F_E3 ",length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3']),' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E3']),digits=2),' ',round(length(dset$RID[dset$APOE=='E3/E3' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$APOE=='E3/E3']),digits=2)))
  print(paste0("M_E3 ",length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3']),' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E3']),digits=2)))
  print(paste0("F_E2 ",length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2']),' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E2']),digits=2),' ',round(length(dset$RID[dset$APOE=='E3/E2' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$APOE=='E3/E2']),digits=2)))
  print(paste0("M_E2 ",length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2']),' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E2']),digits=2)))
  print(paste0("F_E4 ",length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4']),' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='F' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='F' & dset$APOE=='E3/E4']),digits=2),' ',round(length(dset$RID[dset$APOE=='E3/E4' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$APOE=='E3/E4']),digits=2)))
  print(paste0("M_E4 ",length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4']),' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4' & dset$AmyloidPos==1])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4']),digits=2),
               ' ',round(mean(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$Age[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(mean(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$PTEDUCAT[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset$CDSOB[dset$Sex=='M' & dset$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4' & dset$DXSUM=='CN']),':',length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4' & dset$DXSUM=='MCI']),
               ' ',round(length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4' & dset$DXSUM=='MCI'])*100/length(dset$RID[dset$Sex=='M' & dset$APOE=='E3/E4']),digits=2)))
  
  
  cat("\n Baseline only \n")
  print(paste0("Group FirstRecordN FirstRecord%AmyloidPos FirstRecordAgeMean FirstRecordAgeStdDev EducationMean EducationStdDev CDSOBMean CDSOBStdDev CN:MCI %MCI %MCIbyAPOE"))
  print(paste0("Total ",length(dset1$RID),' ',round(length(dset1$RID[dset1$AmyloidPos==1])*100/length(dset1$RID),digits=2),' ',round(mean(dset1$Age,na.rm=TRUE),digits=2),
               ' ',round(stats::sd(dset1$Age,na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT,na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT,na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB,na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB,na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$DXSUM=='MCI']),' ',round(length(dset1$RID[dset1$DXSUM=='MCI'])*100/length(dset1$RID),digits=2)))
  print(paste0("F_E3 ",length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3']),' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E3']),digits=2),' ',round(length(dset1$RID[dset1$APOE=='E3/E3' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$APOE=='E3/E3']),digits=2)))
  print(paste0("M_E3 ",length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3']),' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E3'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E3']),digits=2)))
  print(paste0("F_E2 ",length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2']),' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E2']),digits=2),' ',round(length(dset1$RID[dset1$APOE=='E3/E2' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$APOE=='E3/E2']),digits=2)))
  print(paste0("M_E2 ",length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2']),' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E2'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E2']),digits=2)))
  print(paste0("F_E4 ",length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4']),' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='F' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='F' & dset1$APOE=='E3/E4']),digits=2),' ',round(length(dset1$RID[dset1$APOE=='E3/E4' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$APOE=='E3/E4']),digits=2)))
  print(paste0("M_E4 ",length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4']),' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4' & dset1$AmyloidPos==1])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4']),digits=2),
               ' ',round(mean(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$Age[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(mean(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$PTEDUCAT[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),
               ' ',round(mean(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',round(stats::sd(dset1$CDSOB[dset1$Sex=='M' & dset1$APOE=='E3/E4'],na.rm=TRUE),digits=2),' ',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4' & dset1$DXSUM=='CN']),':',length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4' & dset1$DXSUM=='MCI']),
               ' ',round(length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4' & dset1$DXSUM=='MCI'])*100/length(dset1$RID[dset1$Sex=='M' & dset1$APOE=='E3/E4']),digits=2)))
  
  
  }

#Plotting Scatter of the Dataset for use with Full ADNI dataset
PlotScatter <- function(dset, Xtitle, Ytitle){
  if('M'%in% levels(factor(dset$Sex))){
    
    print(ggplot(dset[dset$Sex=='M',],aes(Age,reg,color=APOE))+geom_point(size=3, alpha =0.8)+
            labs(title=paste(Ytitle, " APOE effect \n within ",dxtext," for Males"),x =Xtitle, y = Ytitle) +
            xlim(min(dset$Age),max(dset$Age))+ylim(min(dset$reg),max(dset$reg)))
  }
  if('F'%in% levels(factor(dset$Sex))){
    print(ggplot(dset[dset$Sex=='F',],aes(Age,reg,color=APOE))+geom_point(size=3, shape=17, alpha =0.8)+
            labs(title=paste(Ytitle, " APOE effect \n within ",dxtext, 'for Females'),x =Xtitle, y = Ytitle) +
            xlim(min(dset$Age),max(dset$Age))+ylim(min(dset$reg),max(dset$reg)))
  }
  # #Scatter For the whole dataset
  # ggplot(dset,aes(Age,reg,color=APOE,shape=DXSUM,size=Sex))+
  #   geom_point()+scale_shape_manual(values=c(17,3))+scale_size_manual(values = c(2,4))+
  #   labs(title=paste(Ytitle, " APOE effect \n within CN and MCI for Females"),x =Xtitle,
  #        y = Ytitle)+ facet_wrap(vars(Sex))
}

plotAPOEnSexinAgeBins<-function(tmpdd,measure){
  APOE_cols=c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red") 
  dset=tmpdd
  dset=dset[,c("RID","APOE","Age","Sex")]
  dset$AgeBin=NA
  dset$AgeBin[dset$Age<65]="<65"
  dset$AgeBin[dset$Age>=65 & dset$Age<75]="65-75"
  dset$AgeBin[dset$Age>=75 & dset$Age<85]="75-85"
  dset$AgeBin[dset$Age>=85]=">=85"
  dset=dset[!is.na(dset$AgeBin),]
  dset=dset[order(dset$RID,dset$Age),]
  dset1=dset[!duplicated(dset$RID),]
  
  dmp=data.frame(expand.grid(AgeBin=c("<65","65-75","75-85",">=85"),APOE=c('E3/E3','E3/E2','E3/E4'),Count=0,CountType=c("N","NumRecords"))) #6065,
  dmp$AgeBin=factor(dmp$AgeBin,levels = c("<65","65-75","75-85",">=85"))

  for (agebin in unique(dmp$AgeBin)){
    for (apoe in unique(dmp$APOE)){
      dmp$Count[dmp$AgeBin==agebin & dmp$APOE==apoe & dmp$CountType=="NumRecords"]=length(dset$RID[dset$AgeBin==agebin & dset$APOE==apoe])
      dmp$Count[dmp$AgeBin==agebin & dmp$APOE==apoe & dmp$CountType=="N"]=length(dset1$RID[dset1$AgeBin==agebin & dset1$APOE==apoe])
      
    }
  }
  dmp$Percent[dmp$CountType=="NumRecords"]=dmp$Count[dmp$CountType=="NumRecords"]*100/length(dset$RID)
  dmp$Percent[dmp$CountType=="N"]=dmp$Count[dmp$CountType=="N"]*100/length(dset1$RID)
  
  
  ggp1= ggplot(data=dmp)+
    geom_bar_pattern(aes(x=CountType,y=Count,fill=APOE,pattern=CountType),stat = "identity",width = 1)+
    labs(y = '', x = 'Age (Years)')+
    scale_fill_manual(values = APOE_cols)+ 
    scale_pattern_manual(values=c("none","stripe"))+
    theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
          panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
          panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
          strip.background = element_blank(),strip.placement = "outside",strip.text = element_text(size=60),
          axis.text.x = element_blank(),axis.text.y = element_text(size = 60),axis.title = element_text(size=60))+
    facet_grid(~AgeBin,switch = "both") 

  #Sex in Age bins
  
  dmp=data.frame(expand.grid(AgeBin=c("<65","65-75","75-85",">=85"),Sex=c(levels(factor(dset$Sex))),Count=0,CountType=c("N","NumRecords"))) #6065,
  dmp$AgeBin=factor(dmp$AgeBin,levels = c("<65","65-75","75-85",">=85"))
  # dmp$APOE_cols=paste0(dmp$APOE,dmp$CountType)
  
  for (agebin in unique(dmp$AgeBin)){
    for (sx in unique(dmp$Sex)){
      dmp$Count[dmp$AgeBin==agebin & dmp$Sex==sx & dmp$CountType=="NumRecords"]=length(dset$RID[dset$AgeBin==agebin & dset$Sex==sx])
      dmp$Count[dmp$AgeBin==agebin & dmp$Sex==sx & dmp$CountType=="N"]=length(dset1$RID[dset1$AgeBin==agebin & dset1$Sex==sx])
      
    }
  }
  dmp$Percent[dmp$CountType=="NumRecords"]=dmp$Count[dmp$CountType=="NumRecords"]*100/length(dset$RID)
  dmp$Percent[dmp$CountType=="N"]=dmp$Count[dmp$CountType=="N"]*100/length(dset1$RID)
  
  
  ggp2=ggplot(data=dmp)+
    geom_bar_pattern(aes(x=CountType,y=Count,fill=Sex,pattern=CountType),stat = "identity",width = 1)+
    labs(y = '', x = 'Age (Years)')+
    scale_fill_manual(values = c("M"="lightblue","F"="lightpink"))+
    scale_pattern_manual(values=c("none","stripe"))+
    theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
          panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
          panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
          strip.background = element_blank(),strip.placement = "outside",strip.text = element_text(size=60),
          axis.text = element_blank(),axis.title = element_text(size=60))+
    facet_grid(~AgeBin,switch = "both") 
  

  jpeg(paste0(measure,"_Sex&APOEinAgeBins.jpg"),width=3000,height = 1500)
  print(ggp1+ggp2)
  dev.off()

  frink <- image_read(paste0(measure,"_Sex&APOEinAgeBins.jpg"))
  frink<-image_annotate(frink, paste("Measure:",measure,"    \n# Participants:",length(dset1$RID),"   \n# Total Records:",length(dset$RID)), color = "darkblue", size = 40,location = "+225+100",boxcolor = "white")
  print(frink)
  paste(image_write(frink,paste0(measure,"_Sex&APOEinAgeBins.jpg")))
  rm(dset,dset1,dmp,ggp1,ggp2)
}

compareTimelines<-function(refE,tstE,FileName){ 
  # This function from the other file eventually
}

#Getting Data from different files from Loni ADNI
getAPOE <-function(){
  apoe=read.csv('APOERES_1.csv')
  apoe1=read.csv('APOE_ADNI3.csv') #Send by Loni rep, not yet on ADNI Loni website
  
  apoe=apoe[, c('RID','APGEN1','APGEN2')]
  apoe$APOE4=NA
  apoe$APOE4[apoe$APGEN1==4 | apoe$APGEN2==4]=1
  apoe$APOE4[apoe$APGEN1==4 & apoe$APGEN2==4]=2
  apoe$APOE4[apoe$APGEN1!=4 & apoe$APGEN2!=4]=0
  apoe$APOE2=NA
  apoe$APOE2[apoe$APGEN1==2 | apoe$APGEN2==2]=1
  apoe$APOE2[apoe$APGEN1==2 & apoe$APGEN2==2]=2
  apoe$APOE2[apoe$APGEN1!=2 & apoe$APGEN2!=2]=0
  apoe$APOE = 'NA'
  apoe$APOE[apoe$APOE2==0 & apoe$APOE4==0] = 'E3/E3'
  apoe$APOE[apoe$APOE2==1 & apoe$APOE4==0] = 'E3/E2'
  apoe$APOE[apoe$APOE2==0 & apoe$APOE4==1] = 'E3/E4'
  # apoe$APOE[apoe$APOE2==2 & apoe$APOE4==0] = 'E2/E2'
  # apoe$APOE[apoe$APOE2==0 & apoe$APOE4==1] = 'E3/E4'
  # 
  
  apoe=apoe[apoe$APOE!='NA',]
  apoe$APOE = factor(apoe$APOE)
  apoe$APOE <- relevel(apoe$APOE, ref='E3/E3')
  apoe=apoe[!duplicated(apoe[c("RID")]),c('RID','APOE')]
  
  apoe1=apoe1[apoe1$RID %notin% apoe$RID,]
  apoe1$APOE=paste0(apoe1$APOE.Value.1,"/",apoe1$APOE.Value.2)
  apoe1$APOE[apoe1$APOE=='E2/E3']='E3/E2'
  apoe1$APOE[apoe1$APOE=='E4/E3']='E3/E4'
  apoe1=apoe1[!duplicated(apoe1[c("RID")]) & (apoe1$APOE=='E3/E3' | apoe1$APOE=='E3/E2' | apoe1$APOE=='E3/E4'),c('RID','APOE')]
  apoe=rbind(apoe,apoe1)
  return(apoe)
}

getCDR <- function(){
  cdr=read.csv('CDR_11Jun2023.csv') 
  cdr=cdr[(cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))| (cdr$USERDATE!=''& !is.na(cdr$USERDATE)),]
  cdr$EXAMDATE[(cdr$EXAMDATE=='' | is.na(cdr$EXAMDATE))& (cdr$USERDATE!='' & !is.na(cdr$USERDATE))]=cdr$USERDATE[(cdr$EXAMDATE=='' | is.na(cdr$EXAMDATE))& (cdr$USERDATE!='' & !is.na(cdr$USERDATE))]
  cdr$USERDATE[(cdr$USERDATE=='' | is.na(cdr$USERDATE))& (cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))]=cdr$EXAMDATE[(cdr$USERDATE=='' | is.na(cdr$USERDATE))& (cdr$EXAMDATE!='' & !is.na(cdr$EXAMDATE))]
  
  cdr$dt=as.Date(cdr$EXAMDATE,format="%Y-%m-%d")-as.Date(cdr$USERDATE,format="%Y-%m-%d")
  cdr$EXAMDATE[cdr$dt>0]=cdr$USERDATE[cdr$dt>0]

    colnames(cdr)[colnames(cdr)=='USERDATE']='CDR.DATE'
  cdr=cdr[cdr$CDGLOBAL!=-1 & !is.na(cdr$CDGLOBAL),]
  cdr$CDSOB = cdr$CDMEMORY + cdr$CDORIENT + cdr$CDJUDGE + cdr$CDCOMMUN + cdr$CDHOME + cdr$CDCARE
  cdr=cdr[,c('RID','CDR.DATE','CDSOB')]
  cdr <- cdr[order(cdr$RID, cdr$CDR.DATE),]
  return(cdr)
}

getCentiloidData <- function(){
#Old Data   # av45=read.csv('UCBERKELEYAV45.csv') ----
#   # av45$AV45.DATE=av45$EXAMDATE
#   # av45=av45[,c('RID','AV45.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_WHOLECEREBNORM","SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF",
#   #              "SUMMARYSUVR_COMPOSITE_REFNORM","SUMMARYSUVR_COMPOSITE_REFNORM_0.79CUTOFF")]
#   av45=read.csv('UCBERKELEYAV45_04_26_22.csv')
#   # #av45=read.csv('UCBERKELEYAV45_11_16_21.csv')
#   av45$AV45.DATE=av45$EXAMDATE
#   av45=av45[,c('RID','AV45.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_WHOLECEREBNORM","SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF")]
#   #Removed ,"SUMMARYSUVR_COMPOSITE_REFNORM","SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF" as it was never used.
#   
#   # fbb=read.csv('UCBERKELEYFBB.csv')
#   # fbb$FBB.DATE=fbb$EXAMDATE
#   # fbb=fbb[,c('RID','FBB.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_WHOLECEREBNORM",
#   #            "SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF","SUMMARYSUVR_CEREBGREYNORM",
#   #            "SUMMARYSUVR_CEREBGREYNORM_1.33CUTOFF")]
#   fbb=read.csv('UCBERKELEYFBB_04_26_22.csv')
#   # #fbb=read.csv('UCBERKELEYFBB_11_16_21.csv')
#   fbb$FBB.DATE=fbb$EXAMDATE
#   fbb=fbb[,c('RID','FBB.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_WHOLECEREBNORM",
#              "SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF")]
#   
#   #Removed ,"SUMMARYSUVR_CEREBGREYNORM","SUMMARYSUVR_CEREBGREYNORM_1.33CUTOFF" as it was never used
#   d=merge(av45,fbb,by=c('RID','VISCODE','VISCODE2'),all=T)
#   d$VISCODE[d$VISCODE2!=""]=d$VISCODE2[d$VISCODE2!=""]  #Copying VISCODE2 to VISCODE as long as VISCODE2 was not empty. VISCODE2 seems to have more levels
#   # rm(av45,fbb)
#   d$Centiloid=196.9*d$SUMMARYSUVR_WHOLECEREBNORM.x-196.03
#   d$Centiloid[(is.na(d$Centiloid))]=159.08*d$SUMMARYSUVR_WHOLECEREBNORM.y[is.na(d$Centiloid)]-151.65
#   d$Cent.DATE=d$AV45.DATE
#   d$Cent.DATE[is.na(d$Cent.DATE)]=d$FBB.DATE[is.na(d$Cent.DATE)]
#   d$AmyloidPos=d$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF
#   d$AmyloidPos[is.na(d$AmyloidPos)]=d$SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF[is.na(d$AmyloidPos)]
#   d = d[!is.na(d$Centiloid),c('RID','Cent.DATE','Centiloid','AmyloidPos')]
#   d <- d[order(d$RID, d$Cent.DATE),]
#   d <- d[!duplicated(d[c("RID","Cent.DATE")]),] #Keep all unique centiloid scan visits
#   return(d)
   
   # The new AV45 and FBB files do not have the SUMMARYSUVR_WHOLECEREBNORM, ----
   # and instead have the SUMMARY_SUVR column. Assuming they are the same.
   # ROIs tested in Rubinski paper medial-orbitofrontal, precuneus, posterior cingulate, 
   # inferior parietal, inferior temporal and parahippocampal gyrus
   av45=read.csv("UCBERKELEYAV45_8mm_02_17_23_11Jun2023.csv")
   colnames(av45)[colnames(av45)=='EXAMDATE']='AV45.DATE'
   av45=av45[,c('RID','AV45.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF",
                "SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
                "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]
   
   av45[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
           "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
           "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]=196.9* av45[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
                                                                                   "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                                                                                   "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]-196.03
   fbb=read.csv("UCBERKELEYFBB_8mm_02_17_23_11Jun2023.csv")
   fbb$FBB.DATE=fbb$EXAMDATE
   fbb=fbb[,c('RID','FBB.DATE','VISCODE','VISCODE2',"SUMMARYSUVR_COMPOSITE_REFNORM_0.74CUTOFF",
              "SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
              "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
              "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]
   fbb[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
          "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
          "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]=159.08*fbb[,c("SUMMARY_SUVR", "CTX_MEDIALORBITOFRONTAL_SUVR","CTX_PRECUNEUS_SUVR",
                                                                                 "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                                                                                 "CTX_INFERIORTEMPORAL_SUVR","CTX_PARAHIPPOCAMPAL_SUVR")]-151.65
   
   Centiloid=merge(av45,fbb,by=c('RID','VISCODE','VISCODE2'),all=T)
   Centiloid$VISCODE2[Centiloid$VISCODE2==""]=Centiloid$VISCODE[Centiloid$VISCODE2==""]  #IF VISCODE2 is empty, copying to VISCODE to it. VISCODE2 seems to have more levels, hence maintaining it.
   Centiloid$Centiloid=Centiloid$SUMMARY_SUVR.x #SUMMARY_SUVR.y added below
   # Centiloid$MOF_Centiloid=Centiloid$CTX_MEDIALORBITOFRONTAL_SUVR.x
   # Centiloid$MOF_Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$CTX_MEDIALORBITOFRONTAL_SUVR.y[is.na(Centiloid$Centiloid)]
   # Centiloid$PRECU_Centiloid=Centiloid$CTX_PRECUNEUS_SUVR.x
   # Centiloid$PRECU_Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$CTX_PRECUNEUS_SUVR.y[is.na(Centiloid$Centiloid)]
   # Centiloid$PCC_Centiloid=Centiloid$CTX_POSTERIORCINGULATE_SUVR.x
   # Centiloid$PCC_Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$CTX_POSTERIORCINGULATE_SUVR.y[is.na(Centiloid$Centiloid)]
   # Centiloid$INFP_Centiloid=Centiloid$CTX_INFERIORPARIETAL_SUVR.x
   # Centiloid$INFP_Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$CTX_INFERIORPARIETAL_SUVR.y[is.na(Centiloid$Centiloid)]
   # Centiloid$PHC_Centiloid=Centiloid$CTX_PARAHIPPOCAMPAL_SUVR.x
   # Centiloid$PHC_Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$CTX_PARAHIPPOCAMPAL_SUVR.y[is.na(Centiloid$Centiloid)]
   Centiloid$Cent.DATE=Centiloid$AV45.DATE
   Centiloid$Cent.DATE[is.na(Centiloid$Centiloid)]=Centiloid$FBB.DATE[is.na(Centiloid$Centiloid)]
   Centiloid$AmyloidPos=Centiloid$SUMMARYSUVR_COMPOSITE_REFNORM_0.78CUTOFF
   Centiloid$AmyloidPos[is.na(Centiloid$Centiloid)]=Centiloid$SUMMARYSUVR_COMPOSITE_REFNORM_0.74CUTOFF[is.na(Centiloid$Centiloid)]
   Centiloid$Centiloid[is.na(Centiloid$Centiloid)]= Centiloid$SUMMARY_SUVR.y[is.na(Centiloid$Centiloid)]
   
   Centiloid = Centiloid[!is.na(Centiloid$Centiloid),c('RID','Cent.DATE','Centiloid','AmyloidPos')]
   Centiloid <- Centiloid[order(Centiloid$RID, Centiloid$Cent.DATE),]
   Centiloid <- Centiloid[!duplicated(Centiloid[c("RID","Cent.DATE")]),] #Keep all unique centiloid scan visits
   rm(av45,fbb)
   return(Centiloid)
   }

getWML <- function(){
  # #wml_old=read.csv('ADNI_UCD_WMH_09_01_20_1.csv') #UCDevis
  # #wml=read.csv('ADNI_UCD_WMH_12_13_21.csv')
  wml=read.csv('ADNI_UCD_WMH_05_02_22.csv')
  wml=wml[order(wml$RID,wml$EXAMDATE),]
  wml$VISCODE=wml$VISCODE2
  wml=wml[wml$STATUS=="",] #Status with "NO3DT1" and "Bad Flair" were removed
  colnames(wml)[colnames(wml)=='TOTAL_BRAIN']='wmlICV'
  colnames(wml)[colnames(wml)=='TOTAL_WMH']='WMH'
  # wml$WML=wml$WMH/wml$wmlICV*100
  wml=wml[,c('RID','EXAMDATE','WMH','wmlICV','COLPROT')]
  wml=wml[!is.na(wml$WMH) & !is.na(wml$wmlICV) & !is.na(wml$EXAMDATE),]
  wml=wml[order(wml$RID, wml$EXAMDATE),]
  wml=wml[!duplicated(wml[c("RID","EXAMDATE")]),]
  colnames(wml)[colnames(wml)=='EXAMDATE']='wml.DATE'
  wml$RIDCOLPROT=interaction(wml$RID,wml$COLPROT)
  uniqueRIDCOLPROT=unique(wml$RIDCOLPROT)
  wml$medianwmlICVbyCOLPROT=wml$wmlICV #just a place holder
  for (RIDCOLPROT in uniqueRIDCOLPROT) {
    wml$medianwmlICVbyCOLPROT[wml$RIDCOLPROT==RIDCOLPROT]=median(wml$wmlICV[wml$RIDCOLPROT==RIDCOLPROT])
  }
  wml$WML=wml$WMH/wml$medianwmlICVbyCOLPROT *100 #This will be the covariate in the rest of the measures
  # wml=subset(wml, select = -c(COLPROT,RIDCOLPROT))
  return(wml)
}

getDemog <- function(){
  dem=read.csv('PTDEMOG_3.csv') #as of 5/30/2023
  dem=dem[!is.na(dem$PTGENDER) & !is.na(dem$PTEDUCAT) & dem$PTGENDER!=-4,]
  colnames(dem)[colnames(dem)=='PTGENDER']='Sex'
  dem$Sex=factor(dem$Sex)
  levels(dem$Sex)=c('M','F') #initially 1 and 2
  dem$refDOB=as.Date(paste0(as.character(dem$PTDOBYY),'-',as.character(dem$PTDOBMM),'-01'),format="%Y-%m-%d") #added 01 of birthday month as reference date to calculate  at different tests
  colnames(dem)[colnames(dem)=='USERDATE']='dem.DATE'
  dem=dem[,c('RID','PTEDUCAT','Sex','dem.DATE','refDOB',"PTRACCAT","PTETHCAT")]
  #PTRACCAT: 1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander;
  #          4=Black or African American; 5=White; 6=More than one race; 7=Unknown
  #PTETHCAT: 1=Hispanic or Latino; 2=Not Hispanic or Latino; 3=Unknown
  
  return(dem)
}

getDXSUM <- function(){
  dxsm=read.csv('DXSUM_PDXCONV_ADNIALL_14Jun2023.csv')
  # dxsm$dt=as.Date(dxsm$EXAMDATE,format="%Y-%m-%d")-as.Date(dxsm$USERDATE,format="%Y-%m-%d")
  # dxsm$EXAMDATE[dxsm$EXAMDATE!='' & dxsm$USERDATE!='' & dxsm$dt>0]=dxsm$USERDATE[dxsm$EXAMDATE!='' & dxsm$USERDATE!='' & dxsm$dt>0]
  # dxsm$EXAMDATE[dxsm$EXAMDATE=='' & dxsm$USERDATE!='']=dxsm$USERDATE[dxsm$EXAMDATE=='' & dxsm$USERDATE!='']
  dxsm=dxsm[(dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))| (dxsm$USERDATE!=''& !is.na(dxsm$USERDATE)),]
  dxsm$EXAMDATE[(dxsm$EXAMDATE=='' | is.na(dxsm$EXAMDATE)) & (dxsm$USERDATE!='' & !is.na(dxsm$USERDATE))]=dxsm$USERDATE[(dxsm$EXAMDATE=='' | is.na(dxsm$EXAMDATE)) & (dxsm$USERDATE!='' & !is.na(dxsm$USERDATE))]
  dxsm$USERDATE[(dxsm$USERDATE=='' | is.na(dxsm$USERDATE))& (dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))]=dxsm$EXAMDATE[(dxsm$USERDATE=='' | is.na(dxsm$USERDATE))& (dxsm$EXAMDATE!='' & !is.na(dxsm$EXAMDATE))]
  dxsm$dt=as.Date(dxsm$EXAMDATE,format="%Y-%m-%d")-as.Date(dxsm$USERDATE,format="%Y-%m-%d")
  dxsm$EXAMDATE[dxsm$dt>0]=dxsm$USERDATE[dxsm$dt>0]
  
  colnames(dxsm)[colnames(dxsm)=='EXAMDATE']='DX.DATE'
  dxsm=dxsm[!is.na(dxsm$DXCHANGE) | !is.na(dxsm$DIAGNOSIS),]
  dxsm$DXSUM=NA
  dxsm$DXSUM[dxsm$DIAGNOSIS==1 | dxsm$DXCHANGE==1 | dxsm$DXCHANGE==7 | dxsm$DXCHANGE==9]='CN'
  dxsm$DXSUM[dxsm$DIAGNOSIS==2 | dxsm$DXCHANGE==2 | dxsm$DXCHANGE==4 | dxsm$DXCHANGE==8]='MCI'
  dxsm$DXSUM[dxsm$DIAGNOSIS==3 | dxsm$DXCHANGE==3 | dxsm$DXCHANGE==5 | dxsm$DXCHANGE==6]='AD'
  dxsm$DXSUM<-factor(dxsm$DXSUM)
  dxsm=dxsm[!is.na(dxsm$DX.DATE) & !is.na(dxsm$DXSUM),c('RID','DX.DATE','DXSUM')]
  return(dxsm)
}

# PlotDataBins <- function(dset){ #----
#   
#   Calculator <- function(tmpd,sexname)
#   {
#     tmpd$AgeBin=NA
#     #tmpd$AgeBin[tmpd$Age>=60 & tmpd$Age<65]=6065
#     tmpd$AgeBin[tmpd$Age>=65 & tmpd$Age<70]=6570
#     tmpd$AgeBin[tmpd$Age>=70 & tmpd$Age<75]=7075
#     tmpd$AgeBin[tmpd$Age>=75 & tmpd$Age<80]=7580
#     tmpd$AgeBin[tmpd$Age>=80 & tmpd$Age<=85]=8085
#     tmpd=tmpd[!is.na(tmpd$AgeBin),]
#     
#     dmp=data.frame(expand.grid(AgeBin=c(6570,7075,7580,8085),DXSUM=c(levels(factor(tmpd$DXSUM))),APOE=c('E3/E3','E3/E2','E3/E4'),Count=0,Percent=0,DXAP=NA,TotPop=0)) #6065,
#     for (AgeCount in levels(factor(dmp$AgeBin))) 
#     {
#       for (APOECount in levels(factor(dmp$APOE))) 
#       {
#         for (DX in levels(factor(dmp$DXSUM))) 
#         {
#           dmp$Count[dmp$DXSUM==DX & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=length(tmpd$RID[tmpd$DXSUM==DX & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])
#           
#         }
#         dmp$TotPop[dmp$APOE == APOECount & dmp$AgeBin==as.numeric(AgeCount)]= sum(dmp$Count[dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount])
#         dmp$Percent[dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount] = dmp$Count[dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]*100/dmp$TotPop[dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]
#         
#         # dmp$Count[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=length(tmpd$RID[tmpd$DXSUM=='CN' & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])
#         # dmp$Count[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=length(tmpd$RID[tmpd$DXSUM=='MCI' & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])
#         # dmp$Percent[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=dmp$Count[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]*100/(dmp$Count[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]+dmp$Count[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount])
#         # dmp$Percent[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=dmp$Count[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]*100/(dmp$Count[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]+dmp$Count[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount])
#         # dmp$DXAP[dmp$DXSUM=='CN' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=paste0('CN_',APOECount)
#         # dmp$DXAP[dmp$DXSUM=='MCI' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]=paste0('MCI_',APOECount)
#       }
#     }
#     dmp$Percent[is.na(dmp$Percent)]=0
#     dmp$DXAP=paste0(dmp$DXSUM,'_',dmp$APOE)
#     dmp$AgeBin=factor(dmp$AgeBin)
#     dmp$DXSUM=factor(dmp$DXSUM)
#     
#     dmp$APOE=factor(dmp$APOE) 
#     dmp <- with(dmp, dmp[order(AgeBin, DXSUM, APOE),])
#     cols <- c("CN_E3/E3"="#FF3300","MCI_E3/E3"="#CC3300","AD_E3/E3"="#990000","CN_E3/E2"="#66CC00","MCI_E3/E2"="#009933","AD_E3/E2"="#006600","CN_E3/E4"="#99CCFF","MCI_E3/E4"="#3399FF","AD_E3/E4"="#003366") 
#     print(ggplot(data=dmp, aes(x=APOE, y=Count, fill=DXAP, width=0.5)) + geom_bar(stat="identity") + scale_fill_manual(values=cols) +
#             labs(title=paste(" N within AgeBins for",sexname)) + facet_grid(~AgeBin))
#     print(ggplot(data=dmp, aes(x=APOE, y=Percent, fill=DXAP, width=0.5)) + geom_bar(stat="identity") + scale_fill_manual(values=cols) + 
#             labs(title=paste(" N% within AgeBins for",sexname)) + facet_grid(~AgeBin))
#     rm(dmp,tmpd)
#   }
#   
#   if('M'%in% levels(factor(dset$Sex))){
#     Calculator(dset[dset$Sex=='M',c("RID","APOE","Age","DXSUM")],'Males')
#     
#   }
#   if('F'%in% levels(factor(dset$Sex))){
#     Calculator(dset[dset$Sex=='F',c("RID","APOE","Age","DXSUM")],'Females')
#     
#   }
# }
# 
# PlotDataBinByYear <- function(tmpdd,datasetname){
#   Calculate <-function(tmpd,datasetname)
#   {
#     tmpd$Agebin=floor(tmpd$Age)
#     AgeBin=unique(tmpd$Agebin)
#     dmp<-data.frame(expand.grid(AgeBin=c(AgeBin),APOE=c('E3/E3','E3/E2','E3/E4'),DXSUM=c(levels(factor(tmpd$DXSUM))), AgebinMcount=0,AgebinFcount=0,AgebinTotcount=0))
#     # print(dmp)
#     for (AgeCount in AgeBin) 
#     {
#       for (APOECount in levels(factor(dmp$APOE))) 
#       {
#         for (DX in levels(factor(dmp$DXSUM)))
#         {
#           dmp$AgebinMcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$DXSUM==DX]=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$DXSUM==DX & tmpd$Sex=='M'])
#           dmp$AgebinFcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$DXSUM==DX]=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$DXSUM==DX & tmpd$Sex=='F'])
#         }
#       }
#     }
#     dmp$AgebinTotcount=dmp$AgebinTotcount+dmp$AgebinFcount
#     dmp$AgeBin=factor(dmp$AgeBin)
#     dmp$APOE=factor(dmp$APOE)
#     dmp$APOE=relevel(dmp$APOE,ref = "E3/E3")
#     dmp$DXSUM=factor(dmp$DXSUM)
#     dmp$DXSUM=relevel(dmp$DXSUM,ref = "CN")
#     dmp <- with(dmp, dmp[order(AgeBin, APOE,DXSUM),])
#     # print(dmp)
#     
#     dmp$DXAP=paste0(dmp$DXSUM,'_',dmp$APOE)
#     cols <- c("CN_E3/E3"="#FF3300","MCI_E3/E3"="#990000","CN_E3/E2"="#66CC00","MCI_E3/E2"="#006600","CN_E3/E4"="#3399FF","MCI_E3/E4"="#003366") #"AD_E3/E3"="#CC3300","AD_E3/E2"="#009933","AD_E3/E4"="#99CCFF"
#     APOEcols <- c("E3/E3"="#FF3300","E3/E2"="#66CC00","E3/E4"="#3399FF") #"AD_E3/E3"="#CC3300","AD_E3/E2"="#009933","AD_E3/E4"="#99CCFF"
#     DXcols <- c("CN"="#336600","MCI"="#660033")
#     #Plotting total distribution
#     ggp1=ggplot(data=dmp, aes(x=AgeBin, y=AgebinTotcount, fill=DXAP, width=0.5)) + 
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) + scale_color_manual(values=DXcols) +
#       geom_smooth(aes(group=DXAP, color=DXSUM), se=FALSE)+
#       labs(title=paste(" N in AgeBins in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp2=ggplot(data=dmp, aes(x=AgeBin, y=AgebinMcount, fill=DXAP, width=0.5)) + 
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) + scale_color_manual(values=DXcols) +
#       geom_smooth(aes(group=DXAP, color=DXSUM), se=FALSE)+
#       labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp3=ggplot(data=dmp, aes(x=AgeBin, y=AgebinFcount, fill=DXAP, width=0.5)) + 
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) + scale_color_manual(values=DXcols) +
#       geom_smooth(aes(group=DXAP, color=DXSUM), se=FALSE)+
#       labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset")) + facet_grid(~APOE)
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0(datasetname,"DatasetNDistributionBarPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
# 
#     # ggp1=ggplot(data=dmp[dmp$Sex=='Total',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#     #   geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#     #   labs(title=paste(" N in AgeBins in",datasetname,"dataset"))
#     # ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#     #   geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#     #   labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset"))
#     # ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#     #   geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#     #   labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset"))
#     # print(ggp1+ ggp2 + ggp3)
#     # jpeg(paste0(datasetname,"DatasetNDistributionPlots.jpg"),width=2000,height = 500)
#     # print(ggp1+ ggp2 + ggp3)
#     # dev.off()
#   }
#   Calculate(tmpdd,datasetname)
# }
# 
# PlotDataBinByYearByDiagnosis <- function(tmpdd,datasetname){
#   Calculate <-function(tmpd,datasetname)
#   {
#     tmpd$Agebin=floor(tmpd$Age)
#     AgeBin=unique(tmpd$Agebin)
#     dmp<-data.frame(expand.grid(AgeBin=c(AgeBin),APOE=c('E3/E3','E3/E2','E3/E4'),Sex=c('BothSexes','M','F'), AgebinCNcount=0,AgebinMCIcount=0,Agebincount=0)) #6065,,Sex=c(levels(factor(tmpd$Sex)))
#     for (AgeCount in AgeBin) 
#     {
#       for (APOECount in levels(factor(dmp$APOE))) 
#       {
#         dmp$AgebinCNcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$Sex=='BothSexes']=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$DXSUM=='CN'])
#         dmp$AgebinMCIcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$Sex=='BothSexes']=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$DXSUM=='MCI'])
#         
#         for (SexCount in levels(factor(dmp$Sex)))
#         {
#           if (SexCount!="BothSexes")
#           {
#             dmp$AgebinCNcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$Sex==SexCount]=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$Sex==SexCount & tmpd$DXSUM=='CN'])
#             dmp$AgebinMCIcount[dmp$AgeBin==AgeCount & dmp$APOE==APOECount & dmp$Sex==SexCount]=length(tmpd$RID[tmpd$Agebin==AgeCount & tmpd$APOE==APOECount & tmpd$Sex==SexCount & tmpd$DXSUM=='MCI'])
#           }
#         }
#       }
#     }
#     dmp$Agebincount=dmp$AgebinCNcount+dmp$AgebinMCIcount
#     print(dmp)
#     dmp$AgeBin=factor(dmp$AgeBin)
#     dmp$Sex=factor(dmp$Sex)
#     dmp$APOE=factor(dmp$APOE)
#     dmp$APOE=relevel(dmp$APOE,ref = "E3/E3")
#     dmp$Sex=relevel(dmp$Sex,ref = "BothSexes")
#     dmp <- with(dmp, dmp[order(AgeBin, APOE,Sex),])
#     cols <- c("E3/E3"="#CC3300","E3/E2"="#009933","E3/E4"="#3399FF")
#     
#     #Plotting total distribution
#     ggp1=ggplot(data=dmp[dmp$Sex=='BothSexes',], aes(x=AgeBin, y=Agebincount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste(" N in AgeBins in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(x=AgeBin, y=Agebincount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(x=AgeBin, y=Agebincount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset")) + facet_grid(~APOE)
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0(datasetname,"DatasetNDistributionBarPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     ggp1=ggplot(data=dmp[dmp$Sex=='Total',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in",datasetname,"dataset"))
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset"))
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset"))
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0(datasetname,"DatasetNDistributionPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     
#     
#     #Plotting CN distribution
#     ggp1=ggplot(data=dmp[dmp$Sex=='BothSexes',], aes(x=AgeBin, y=AgebinCNcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("CN N in AgeBins in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(x=AgeBin, y=AgebinCNcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("CN N in AgeBins in Males in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(x=AgeBin, y=AgebinCNcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("CN N in AgeBins in Females in",datasetname,"dataset")) + facet_grid(~APOE)
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0("CNin",datasetname,"DatasetNDistributionBarPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     ggp1=ggplot(data=dmp[dmp$Sex=='Total',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in",datasetname,"dataset"))
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset"))
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset"))
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0("CNin",datasetname,"DatasetNDistributionPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     
#     #Plotting MCI distribution
#     ggp1=ggplot(data=dmp[dmp$Sex=='BothSexes',], aes(x=AgeBin, y=AgebinMCIcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("MCI N in AgeBins in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(x=AgeBin, y=AgebinMCIcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("MCI N in AgeBins in Males in",datasetname,"dataset")) + facet_grid(~APOE)
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(x=AgeBin, y=AgebinMCIcount, fill=APOE, width=0.5)) + #fill=ASbin
#       geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#       labs(title=paste("MCI N in AgeBins in Females in",datasetname,"dataset")) + facet_grid(~APOE)
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0("MCIin",datasetname,"DatasetNDistributionBarPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     ggp1=ggplot(data=dmp[dmp$Sex=='Total',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in",datasetname,"dataset"))
#     ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Males in",datasetname,"dataset"))
#     ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(AgeBin, Agebincount, color=APOE, group=APOE)) +
#       geom_smooth() + scale_color_manual(values=cols) + scale_fill_manual(values=cols)+
#       labs(title=paste(" N in AgeBins in Females in",datasetname,"dataset"))
#     print(ggp1+ ggp2 + ggp3)
#     jpeg(paste0("MCIin",datasetname,"DatasetNDistributionPlots.jpg"),width=2000,height = 500)
#     print(ggp1+ ggp2 + ggp3)
#     dev.off()
#     rm(dmp,tmpd,ggp1,ggp2,ggp3)
#   }
#   Calculate(tmpdd,datasetname)
# }
# 
# PlotCDSOBBins <- function(tmpd,datasetname){
#   tmpd=tmpd[,c("RID","Age","APOE","Sex","CDSOB")]
#   tmpd$AgeBin=NA
#   #tmpd$AgeBin[tmpd$Age>=60 & tmpd$Age<65]=6065
#   tmpd$AgeBin[tmpd$Age>=65 & tmpd$Age<70]=6570
#   tmpd$AgeBin[tmpd$Age>=70 & tmpd$Age<75]=7075
#   tmpd$AgeBin[tmpd$Age>=75 & tmpd$Age<80]=7580
#   tmpd$AgeBin[tmpd$Age>=80 & tmpd$Age<=85]=8085
#   tmpd=tmpd[!is.na(tmpd$AgeBin),]
#   
#   dmp=data.frame(expand.grid(AgeBin=c(6570,7075,7580,8085),APOE=c('E3/E3','E3/E2','E3/E4'),Sex=c('Total','M','F'),CDSOBMean=0,CDSOBStd=0)) #6065,,Sex=c(levels(factor(tmpd$Sex)))
#   for (AgeCount in levels(factor(dmp$AgeBin))) 
#   {
#     for (APOECount in levels(factor(dmp$APOE))) 
#     {
#       dmp$CDSOBMean[dmp$Sex=='Total' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount] = mean(tmpd$CDSOB[tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])
#       dmp$CDSOBStd[dmp$Sex=='Total' & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount] = stats::sd(tmpd$CDSOB[tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])/sqrt(length(tmpd$CDSOB[tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount]))
#       for (Sex in levels(factor(tmpd$Sex))) 
#       {
#         if (Sex!="Total")
#         {
#           dmp$CDSOBMean[dmp$Sex==Sex & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]= mean(tmpd$CDSOB[tmpd$Sex==Sex & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])
#           dmp$CDSOBStd[dmp$Sex==Sex & dmp$AgeBin==as.numeric(AgeCount) & dmp$APOE == APOECount]= stats::sd(tmpd$CDSOB[tmpd$Sex==Sex & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount])/sqrt(length(tmpd$CDSOB[tmpd$Sex==Sex & tmpd$AgeBin==as.numeric(AgeCount) & tmpd$APOE == APOECount]))
#           
#         }
#       }
#     }
#   }
#   
#   dmp$ASbin=paste0(dmp$Sex,'_',dmp$APOE)
#   dmp$AgeBin=factor(dmp$AgeBin)
#   dmp$Sex=factor(dmp$Sex)
#   dmp$APOE=factor(dmp$APOE)
#   dmp$APOE=relevel(dmp$APOE,ref = "E3/E3")
#   dmp$Sex=relevel(dmp$Sex,ref = "Total")
#   dmp <- with(dmp, dmp[order(AgeBin, APOE,Sex),])
#   dmp$APOE_Colors=paste0('1',dmp$APOE)
#   
#   # print(dmp)
#   # cols <- c("M_E3/E3"="#FF3300","F_E3/E3"="#CC3300","M_E3/E2"="#66CC00","F_E3/E2"="#009933","M_E3/E4"="#99CCFF","F_E3/E4"="#3399FF","Total_E3/E3"="#990000","Total_E3/E2"="#006600","Total_E3/E4"="#003366") 
#   # cols <- c("E3/E3"="#CC3300","E3/E2"="#009933","E3/E4"="#3399FF") 
#   ggp1=ggplot(data=dmp[dmp$Sex=='Total',], aes(x=AgeBin, y=CDSOBMean, fill=APOE_Colors, width=0.5)) + #fill=ASbin
#     geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#     geom_errorbar(aes(x=AgeBin,ymin=CDSOBMean-CDSOBStd,ymax=CDSOBMean+CDSOBStd)) +
#     labs(title=paste(" CDSOB in AgeBins in",datasetname,"dataset")) + facet_grid(~APOE) +
#     scale_y_continuous(limits = c(0, max(dmp$CDSOBMean)+max(dmp$CDSOBStd)))+
#     scale_x_discrete(breaks=c("6065","6570","7075","7580","8085"),labels=c("60-65","65-70","70-75","75-80","80-85"))+
#     theme(title = element_text(size = 20),
#           legend.position="none",axis.text = element_text(size = 20), axis.text.x = element_text(size =20,angle = 45),
#           axis.title = element_text(size = 20),strip.text.y = element_text(size =20),strip.text.x = element_text(size =20),
#           panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 2, linetype = "solid"),
#           panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="gray"))
#   
#   ggp2=ggplot(data=dmp[dmp$Sex=='M',], aes(x=AgeBin, y=CDSOBMean, fill=APOE_Colors, width=0.5)) + #fill=ASbin
#     geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#     geom_errorbar(aes(x=AgeBin,ymin=CDSOBMean-CDSOBStd,ymax=CDSOBMean+CDSOBStd)) +
#     labs(title=paste(" CDSOB in AgeBins in Males in",datasetname,"dataset")) + facet_grid(~APOE)+
#     scale_y_continuous(limits = c(0, max(dmp$CDSOBMean)+max(dmp$CDSOBStd)))+
#     scale_x_discrete(breaks=c("6065","6570","7075","7580","8085"),labels=c("60-65","65-70","70-75","75-80","80-85"))+
#     theme(title = element_text(size = 20),
#           legend.position="none",axis.text = element_text(size = 20), axis.text.x = element_text(size =20,angle = 45),
#           axis.title = element_text(size = 20),strip.text.y = element_text(size =20),strip.text.x = element_text(size =20),
#           panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 2, linetype = "solid"),
#           panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="gray"))
#   
#   ggp3=ggplot(data=dmp[dmp$Sex=='F',], aes(x=AgeBin, y=CDSOBMean, fill=APOE_Colors, width=0.5)) + #fill=ASbin
#     geom_bar(stat="identity",) + scale_fill_manual(values=cols) +
#     geom_errorbar(aes(x=AgeBin,ymin=CDSOBMean-CDSOBStd,ymax=CDSOBMean+CDSOBStd)) +
#     labs(title=paste(" CDSOB in AgeBins in Females in",datasetname,"dataset")) + facet_grid(~APOE)+
#     scale_y_continuous(limits = c(0, max(dmp$CDSOBMean)+max(dmp$CDSOBStd)))+
#     scale_x_discrete(breaks=c("6065","6570","7075","7580","8085"),labels=c("60-65","65-70","70-75","75-80","80-85"))+
#     theme(title = element_text(size = 20),
#           legend.position="none",axis.text = element_text(size = 20), axis.text.x = element_text(size =20,angle = 45),
#           axis.title = element_text(size = 20),strip.text.y = element_text(size =20),strip.text.x = element_text(size =20),
#           panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 2, linetype = "solid"),
#           panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="gray"))
#   
#   print(ggp1+ ggp2 + ggp3)
#   jpeg(paste0(datasetname,"DatasetCDSOBDistributionPlots.jpg"),width=2000,height = 500)
#   print(ggp1+ ggp2 + ggp3)
#   dev.off()
#   rm(dmp,tmpd)
# }
#----

# WMH ----
# # Changed this to WMH with ICV (called wmlICV here) 
# # For analysis WMH will be used, 
# # and for other measures WML=WMH/median_wmlICV_by_subject will be used as covariate
# wml=getWML()
# colnames(wml)[colnames(wml)=='wml.DATE']='EXAMDATE'
# 
# apoe=getAPOE()
# wml=merge(wml,apoe,by=c('RID')) # removed ,all.x=T as we need APOE for each record
# rm(apoe)
# wml=wml[!is.na(wml$APOE),]
#
# Adni=adnimerge
# colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
# Adni$Sex=factor(Adni$Sex)
# levels(Adni$Sex)=c('F','M')
# Adni$Sex = relevel(Adni$Sex,ref = 'M')
# colnames(Adni)[colnames(Adni)=='DX']='DXSUM'
# colnames(Adni)[colnames(Adni)=='AGE']='Age'
# colnames(Adni)[colnames(Adni)=='EXAMDATE']='Adnimerge.DATE'
# colnames(Adni)[colnames(Adni)=='EXAMDATE.bl']='Adnimerge.DATE.bl'
# levels(Adni$DXSUM)[levels(Adni$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# Adni=Adni[,c("RID","Adnimerge.DATE","Age","Years.bl","Sex","DXSUM","PTEDUCAT","CDSOB")]
# 
# wml=merge(wml,Adni,by=c('RID'),all.x=T)
# wml$dt=abs(as.Date(wml$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(wml$EXAMDATE,format="%Y-%m-%d"))
# wml=wml[!is.na(wml$CDSOB) & !is.na(wml$PTEDUCAT) & !is.na(wml$Age) & !is.na(wml$DXSUM) & !is.na(wml$Sex) ,]
# wml <- wml[wml$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# wml <- wml[order(wml$RID,wml$EXAMDATE,wml$dt),]
# wml <- wml[!duplicated(wml[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# wml$Age=wml$Age+wml$Years.bl+(as.Date(wml$EXAMDATE,format="%Y-%m-%d")-as.Date(wml$Adnimerge.DATE,format="%Y-%m-%d"))/365
# wml=subset(wml, select = -c(Adnimerge.DATE,Years.bl))
# wml$Sex=factor(wml$Sex)
# wml$DXSUM=factor(wml$DXSUM)
# rm(Adni)




#Alternate way of adding demographic info
# dem=getDemog()
# dem=dem[dem$RID %in% unique(wml$RID),]
# wml=merge(wml,dem,by=c('RID'),all.x=T)
# wml$dt=abs(as.Date(wml$dem.DATE,format="%Y-%m-%d")-as.Date(wml$EXAMDATE,format="%Y-%m-%d")) 
# wml <- wml[order(wml$RID,wml$EXAMDATE,wml$dt),]
# wml <- wml[!duplicated(wml[c("RID","EXAMDATE")]),] #Keep the demographic info collected closest in date to the WML EXAMDATE
# wml$Age=as.numeric(as.Date(wml$EXAMDATE,format="%Y-%m-%d")-wml$refDOB)/365 #Age at the time of wml scan
# wml=wml[!is.na(wml$Age),]
# rm(dem)
# wml=subset(wml, select = -c(refDOB,dem.DATE,dt) )
# 
# cdr=getCDR()
# cdr=cdr[cdr$RID %in% unique(wml$RID),]
# wml=merge(wml,cdr,by=c('RID'),all.x=T)
# wml$dt=abs(as.Date(wml$CDR.DATE,format="%Y-%m-%d")-as.Date(wml$EXAMDATE,format="%Y-%m-%d")) 
# wml<-wml[wml$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# wml <- wml[order(wml$RID,wml$EXAMDATE,wml$dt),]
# wml <- wml[!duplicated(wml[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the WML EXAMDATE
# wml=wml[!is.na(wml$CDSOB),]
# wml=subset(wml, select = -c(CDR.DATE))
# rm(cdr)
# 
# dxsm=getDXSUM()
# wml=merge(wml,dxsm,by=c('RID'),all.x=T)
# wml$dt=abs(as.Date(wml$DX.DATE,format="%Y-%m-%d")-as.Date(wml$EXAMDATE,format="%Y-%m-%d"))
# wml<-wml[wml$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# wml <- wml[order(wml$RID,wml$EXAMDATE,wml$dt),]
# wml <- wml[!duplicated(wml[c("RID","EXAMDATE")]),] #Keep the DXSM exam closest in date to the WML EXAMDATE
# wml=subset(wml, select = -c(DX.DATE) )
# wml=wml[!is.na(wml$DXSUM),]
# rm(dxsm)
# 

# # Adding Centiloid/AmyloidPos info
# CentiloidData=getCentiloidData()
# #CentiloidData <- CentiloidData[!duplicated(CentiloidData[c("RID")]),] #Keep the first Centiloid scan
# wml=merge(wml,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# wml$dt=abs(as.Date(wml$EXAMDATE,format="%Y-%m-%d")-as.Date(wml$Cent.DATE,format="%Y-%m-%d")) 
# wml<-wml[wml$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# wml <- wml[order(wml$RID, wml$EXAMDATE, wml$dt),]
# wml <- wml[!duplicated(wml[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the wml scan
# rm(CentiloidData)
# wml=subset(wml, select = -c(Cent.DATE) )
# # PlotDataBins(wml)

# wml$Age = as.numeric(wml$Age)
# WholeDatasetwithAD=wml
# write.csv(wml,"WMHdataNoCovar.csv")
# rm(wml)

WholeDatasetwithAD=read.csv("WMHdataNoCovar1.csv")
WholeDatasetwithAD$reg=WholeDatasetwithAD$WMH
WholeDatasetwithAD$Sex<-relevel(factor(WholeDatasetwithAD$Sex),ref = 'M')
WholeDatasetwithAD$APOE<-relevel(factor(WholeDatasetwithAD$APOE),ref = 'E3/E3')
WholeDatasetwithAD$int_sex_apoe <- interaction(WholeDatasetwithAD$Sex, WholeDatasetwithAD$APOE) #interaction variable that allows teasing out details of individual groups
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
tmpdd$DXSUM=factor(tmpdd$DXSUM)

#tmpdd = tmpdd[tmpdd$Age>=63 & tmpdd$Age<=87,]#This one is not used currently 6/21/2023
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]

tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)
# tmpdd$Edudiff=tmpdd$PTEDUCAT-16 #Since we have been predicting for 16 years of education

#This one is not used currently 6/21
#tmpdd = tmpdd[tmpdd$RID %notin% tmpdd$RID[tmpdd$reg<100],] # Removing 5 outlier datapoints, all from 1 subject
# tmpdd$RID=factor(tmpdd$RID) This is causing error

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'WMH',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   # group = RID, #Switching to Cross sectional analysis
#                   colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) + #, size=1
#   geom_point()+ #aes(size=1)
#   plt_labs + plot_theme
ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"WMH"))

Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],WMH=1),by=c("RID"),all = TRUE)

jpeg("WMHdataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()

plotAPOEnSexinAgeBins(tmpdd,"WMH")
Save_demog_distribution_withintmpdd(tmpdd,"WMH")

# Complex Model without CDSOB
#Assess Centiloid effect in secondary model, and effect ICV does not change by group, so Centiloid and wmlICV:APOE terms were removed
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE) +
            s(RID,Agediff, bs="re")  + s(RID,bs="re")  + #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex + wmlICV, data = tmpdd,method = "REML") #+ Edudiff + Edudiff:APOE 
# m1$model$PTEDUCAT=m1$model$Edudiff+16
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"WMHPTEDUCATCovarModelWOGender.Rds")

m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex)) +
            s(RID, Agediff, bs = "re") + s(RID,bs="re")  + 
            PTEDUCAT + APOE + Sex + APOE:Sex + wmlICV,
          data = tmpdd,
          method = "REML")
# m2$model$PTEDUCAT=m2$model$Edudiff+16
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"WMHPTEDUCATCovarModelWithGender.Rds")



# #predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL,Sex='F',wmlICV=mean(tmpdd$wmlICV,na.rm=T))) #Edudiff=0,Centiloid=20,
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("WMH prediction within ",dxtext,'without Sex effects'), 
#        x ="Age (years)", y = 'WMH')
# ggp1
# 
#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL,wmlICV=mean(tmpdd$wmlICV,na.rm=T)))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("WMH prediction within ",dxtext),
#     x ="Age (years)", y = '')#'WMH')
# ggp2 + facet_wrap(vars(Sex))
# 
# 
# jpeg("WMHPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()




# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)


#Printing progression with significance info 
m1=readRDS("WMHPTEDUCATCovarModelWOGender.Rds")
m2=readRDS("WMHPTEDUCATCovarModelWithGender.Rds")
SigOrNot=CompareSplineSigSlope("WMH","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT,wmlICV",'F')

m1_predict$z=0
m1_predict$z= (m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$Measure="WMH"
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)

ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= (m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$Measure="WMH"
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))

m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])




jpeg("WMHPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("WMH")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# Printing legends separately ----

tmpm2predict=m2_predict
tmpm2predict$SigAPOE<-factor(tmpm2predict$SigAPOE,levels = c("0E3/E3","1E3/E3","0E3/E4","1E3/E4","0E3/E2","1E3/E2"))
ggptmp=ggplot(tmpm2predict, aes(Age, fit, color = SigAPOE)) + #, linetype = Sex
  geom_smooth_ci(size = 4, ci_z=1, ci_alpha=0.001) + scale_color_manual(values=cols, labels=NULL,name=NULL) +
  labs(x ="Age (years)", y = '') #+  theme(legend.text = element_text(size = 60),legend.title = element_text(size = 60),legend.key.size = unit(4,'cm'),legend.direction = "horizontal")

t=get_legend(ggptmp+ theme(legend.text = element_text(size = 60),legend.title = element_text(size = 60),
                            legend.key.size = unit(6,'cm'),legend.direction = "horizontal"))
grid.newpage()
grid.draw(t)
jpeg("Legends.jpg",width=3000)
grid.draw(t) 
dev.off()

rm(ggptmp,t,tmpm2predict)


# Centiloid ----

# d=getCentiloidData()
# colnames(d)[colnames(d)=='Cent.DATE']<-'EXAMDATE'
# #apoe info
# apoe=getAPOE()
# d=merge(d,apoe,by=c('RID'),all.x=T) #Only  keeping RIDs that have the APOE data
# d=d[!is.na(d$APOE),]
# rm(apoe)

# Adni=adnimerge
# colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
# Adni$Sex=factor(Adni$Sex)
# levels(Adni$Sex)=c('F','M')
# Adni$Sex = relevel(Adni$Sex,ref = 'M')
# colnames(Adni)[colnames(Adni)=='DX']='DXSUM'
# colnames(Adni)[colnames(Adni)=='AGE']='Age'
# colnames(Adni)[colnames(Adni)=='EXAMDATE']='Adnimerge.DATE'
# colnames(Adni)[colnames(Adni)=='EXAMDATE.bl']='Adnimerge.DATE.bl'
# levels(Adni$DXSUM)[levels(Adni$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# Adni=Adni[,c("RID","Adnimerge.DATE","Age","Years.bl","Sex","DXSUM","PTEDUCAT","CDSOB")]
# 
# d=merge(d,Adni,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d"))
# d=d[!is.na(d$CDSOB) & !is.na(d$PTEDUCAT) & !is.na(d$Age) & !is.na(d$DXSUM) & !is.na(d$Sex) ,]
# d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# d$Age=d$Age+d$Years.bl+(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Adnimerge.DATE,format="%Y-%m-%d"))/365
# d=subset(d, select = -c(Adnimerge.DATE,Years.bl,dt))
# d$Sex=factor(d$Sex)
# d$DXSUM=factor(d$DXSUM)
# rm(Adni)

##Alternative way of adding demographic info. Decided to use ADNIMERGE on 6/15/2023
# dem=getDemog()
# dem=dem[dem$RID %in% unique(d$RID),]
# dem=dem[!duplicated(dem[,c("RID","refDOB")]),]
# d=merge(d,dem,by=c('RID'),all.x=T)
# rm(dem)
# d$Age=as.numeric(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$refDOB,format="%Y-%m-%d"))/365 #Age at Centiloid scan
# d=d[!is.na(d$Age),]
# d=subset(d, select = -c(refDOB,dem.DATE) )
# 
# cdr=getCDR()
# cdr=cdr[cdr$RID %in% unique(d$RID),]
# d=merge(d,cdr,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$CDR.DATE,format="%Y-%m-%d"))
# d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# d=d[!is.na(d$CDSOB),]
# d=subset(d, select = -c(CDR.DATE,dt) )
# rm(cdr)
# 
# dxsm=getDXSUM()
# dxsm=dxsm[dxsm$RID %in% unique(d$RID),]
# d=merge(d,dxsm,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$DX.DATE,format="%Y-%m-%d"))
# d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the DXSUM closest in date to the Centiloid scan visid=d[!is.na(d$DXSUM),]
# d=subset(d, select = -c(DX.DATE,dt) )
# rm(dxsm)

# # # Not Adding WML info currently
# wml=getWML()
# wml=wml[wml$RID %in% unique(d$RID),]
# d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format=ations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# d=d[!is.na(d$CDSOB),]
# d=subset(d, select = -c(CDR.DATE,dt) )
# rm(cdr)
#


# dd = d[d$APOE!='NA' & !is.na(d$Age) & !is.na(d$APOE) & !is.na(d$Centiloid),]
# # PlotDataBins(d[d$APOE!='NA' & !is.na(d$Age) & !is.na(d$APOE) & !is.na(d$reg),])
# dd <- dd[order(dd$RID,dd$Age),]
# dd <- dd[!duplicated(dd[c("RID","Age")]),]
# write.csv(dd,"CentiloiddataNoCovar2.csv")
# rm(dd)

WholeDatasetwithAD=read.csv("CentiloiddataNoCovar2.csv")#This was done with new ADNIMERGE installed on 8/2/23 and gives 3 extra records compared to previous version
WholeDatasetwithAD$reg =WholeDatasetwithAD$Centiloid
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM%in%dx,]

tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# tmpdd = tmpdd[tmpdd$Age >=63 & tmpdd$Age <=87,] # F E2 had very limited data outside of thei age range.
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
tmpdd=tmpdd[tmpdd$RID %notin% tmpdd$RID[tmpdd$reg>300],] #Removing 1 RID records that has an outlier record
# tmpdd$RID<-factor(tmpdd$RID) Somehow this is messing up the modelling
tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$DXSUM <- factor(tmpdd$DXSUM)
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Age = as.numeric(tmpdd$Age)
# tmpdd$Edudiff=tmpdd$PTEDUCAT-16 #Since we have been predicting for 16 years of education
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}
plt_labs <- labs(y = 'Amyloid (Centiloid units)',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   group = RID, #Switching to Cross sectional analysis
#                   colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) + #,size=1.1
#   geom_point()+ #aes(size=3)
#   scale_color_manual(values=cols_datascatter) + scale_x_continuous(limits = c(x_min, x_max)) +
#   plt_labs
# 
ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"Centiloid"))
Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],AmyloidPET=1),by=c("RID"),all = TRUE)

jpeg("CentiloiddataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()

plotAPOEnSexinAgeBins(tmpdd,"AmyloidPET")
Save_demog_distribution_withintmpdd(tmpdd,"AmyloidPET")


#Chosen ComplexModel without CDSOB
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE) +
            s(RID,Agediff, bs="re") + s(RID,bs="re")  + 
            PTEDUCAT + APOE + Sex,
          data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"CentiloidPTEDUCATCovarModelWOGender.Rds")

m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex)) +
            s(RID, Agediff, bs = "re") + 
            PTEDUCAT + APOE + Sex + APOE:Sex,
                 data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"CentiloidPTEDUCATCovarModelWithGender.Rds")



# #predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,Sex='F',RID=NULL))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age

# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) +
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("Centiloid prediction within",dxtext,'without Sex effects'),
#        x ="Age (years)", y = 'Centiloid')
# ggp1


#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) +
#   theme(legend.position="none",axis.text = element_text(size = 60)) +
#   #scale_linetype_manual(values=c("dotdash", "longdash")) +
#   labs(#title=paste("Centiloid prediction within",dxtext),
#       x ="Age (years)", y = '')#Centiloid')
#
# ggp2 + facet_wrap(vars(Sex))


m1=readRDS("CentiloidPTEDUCATCovarModelWOGender.Rds")
m2=readRDS("CentiloidPTEDUCATCovarModelWithGender.Rds")
SigOrNot=CompareSplineSigSlope("Centiloid","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT","F")

m1_predict$z=0
m1_predict$z= (m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']))/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$APOE_lineCol=paste0('0',m1_predict$APOE)
m1_predict$Measure="Centiloid"

ggplot(m1_predict, aes(x = Age, y = fit, colour = APOE_lineCol)) + 
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 2, ci_alpha=0.1) +  
  plt_labs + scale_x_continuous(limits = c(x_min, x_max)) + geom_point(data=m1_predict[m1_predict$Significance!=0,],aes(colour = SigAPOE,size=1.2))+
  theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1", size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour ="black"))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= (m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$APOE_lineCol=paste0('0',m2_predict$APOE)
m2_predict$Measure="Centiloid"



ggplot(m2_predict, aes(x = Age, y = fit, colour = APOE_lineCol)) + 
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 2, ci_alpha=0.1) +  
  plt_labs + scale_x_continuous(limits = c(x_min, x_max)) + geom_point(data=m2_predict[m2_predict$Significance!=0,],aes(colour = SigAPOE,size=1.2))+
  theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1", size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour ="black"))+facet_wrap(vars(Sex))


# ggp2 = ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance
#   geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) + 
#   plt_labs + scale_x_continuous(limits = c(x_min, x_max))


jpeg("CentiloidPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()


m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])

sink("ModelSummaryoutput.txt",append = TRUE)
print("Amyloid (Centiloid scale)")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

# Volumetric data ----

### d_old=read.csv('harmonized_imaging_for_pallavi.csv')

# #d=read.csv('freesurferharmnewadni3_trial_1.csv')
# #d=read.csv('freesurferharmnewadni3_trial_2.csv')#Includes ADNI3 data from 6/7/2022 
# 
# d=read.csv('freesurferharmnewadni3_trial_3.csv') #Includes ADNI3 data from 6/23/2023
# # d=tidyr::separate(d,EXAMDATE.y,c("EXAMDATE"),sep=" ",remove = TRUE, extra = "drop")
# # #the EXAMDATE.y had an additional time component that was dropped.
# # names(d)[names(d)=="EXAMDATE.x"]="EXAMDATE"
# # d=subset(d, select = -c(EXAMDATE.y))
# d=d[order(d$RID,d$EXAMDATE,d$STUDY),]
# d=d[!duplicated(d[c("RID","EXAMDATE")]),]
# names(d)[names(d)=="DIAGNOSIS"]="DXSUM"
# d=d[,c(sort(names(d)[grepl("RID|EXAMDATE|AGE|Baseline|harmonized|ICV|DX|PTEDUCAT|PTGENDER|APGEN|OVERALLQC", names(d))]))]
# d$DXSUM[d$DXSUM=="Dementia"]="AD"
# names(d)[names(d)=="PTGENDER"]="Sex"
# d$Sex=relevel(factor(substr(d$Sex,1,1)),ref="M")#converting Male to M and Female to F, factoring it and setting M as the reference
# names(d)[names(d)=="AGE"]="Age"

# d$APOE4=NA
# d$APOE4[d$APGEN1==4 | d$APGEN2==4]=1
# d$APOE4[d$APGEN1==4 & d$APGEN2==4]=2
# d$APOE4[d$APGEN1!=4 & d$APGEN2!=4]=0
# d$APOE2=NA
# d$APOE2[d$APGEN1==2 | d$APGEN2==2]=1
# d$APOE2[d$APGEN1==2 & d$APGEN2==2]=2
# d$APOE2[d$APGEN1!=2 & d$APGEN2!=2]=0
# 
# d$APOE = 'NA'
# d$APOE[d$APOE2==0 & d$APOE4==0] = 'E3/E3'
# d$APOE[d$APOE2==1 & d$APOE4==0] = 'E3/E2'
# d$APOE[d$APOE2==0 & d$APOE4==1] = 'E3/E4'
# d=d[d$APOE!='NA' & !is.na(d$APOE),]
# d$APOE = factor(d$APOE)
# d$APOE <- relevel(d$APOE, ref='E3/E3')
# d=subset(d, select = -c(APGEN1,APGEN2,APOE2,APOE4))
# rid=unique(d$RID)
# 
# # disc=read.csv('fressurfer_data_dict.csv')
# # for(i in 1:length(names(d)[grepl("harm", names(d))])){print(paste(names(d)[grepl("harm", names(d))][i],disc$TEXT[disc$FLDNAME == unlist(strsplit(names(d)[grepl("harm", names(d))][i], split='_', fixed=TRUE))[1]]))}
# # rm(disc)
# # #So far I have these volumes (harmonized between sites) and then meta left and right as the last two
# # # [1] "ST103CV_harm Volume (Cortical Parcellation) of RightParahippocampal"
# # # [1] "ST44CV_harm Volume (Cortical Parcellation) of LeftParahippocampal"
# # # [1] "ST29SV_harm Volume (WM Parcellation) of LeftHippocampus"
# # # [1] "ST88SV_harm Volume (WM Parcellation) of RightHippocampus"
# # # [1] "ST24CV_harm Volume (Cortical Parcellation) of LeftEntorhinal"
# # # [1] "ST32CV_harm Volume (Cortical Parcellation) of LeftInferiorTemporal"
# # # [1] "ST40CV_harm Volume (Cortical Parcellation) of LeftMiddleTemporal"
# # # [1] "ST26CV_harm Volume (Cortical Parcellation) of LeftFusiform"
# # # [1] "ST83CV_harm Volume (Cortical Parcellation) of RightEntorhinal"
# # # [1] "ST91CV_harm Volume (Cortical Parcellation) of RightInferiorTemporal"
# # # [1] "ST99CV_harm Volume (Cortical Parcellation) of RightMiddleTemporal"
# # # [1] "ST85CV_harm Volume (Cortical Parcellation) of RightFusiform"
# # # [1] "left.meta_harm "
# # # [1] "right.meta_harm "
# 
# Adni=adnimerge
# colnames(Adni)[colnames(Adni)=='EXAMDATE']='Adnimerge.DATE'
# colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# Adni=Adni[,c("RID","Adnimerge.DATE","CDSOB")]
# 
# d=merge(d,Adni,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d"))
# d = d[!is.na(d$CDSOB) & d$CDSOB != 'NA' & !is.na(d$dt),]
# d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# d=subset(d, select = -c(Adnimerge.DATE))
# d$Sex=factor(d$Sex)
# d$DXSUM=factor(d$DXSUM)
# rm(Adni)


# cdr=getCDR() Decided to go with ADNIMERGE during 8/3/2023
# cdr=cdr[cdr$RID %in% unique(d$RID),]
# d=merge(d,cdr,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$CDR.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d")) 
# d=d[d$dt<=180,]
# d = d[!is.na(d$CDSOB) & d$CDSOB != 'NA' & !is.na(d$dt),]
# d <- d[order(d$RID, d$EXAMDATE, d$dt),]
# d <- d[!duplicated(d[c("RID", "EXAMDATE")]),] #Keep the CDR exam closest in date to the EXAMDATE
# rm(cdr)
# d=subset(d, select = -c(CDR.DATE))
# 
# Not adding Amyloid and WML info to the dataset 8/3/2023
# # Adding AmyloidPos and Centiloid info
# CentiloidData=getCentiloidData()
# CentiloidData=CentiloidData[CentiloidData$RID %in% unique(d$RID),]
# d=merge(d,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d")) 
# d=d[d$dt<=180,]
# d=d[!is.na(d$dt),]
# d <- d[order(d$RID, d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the FDG scan
# rm(CentiloidData)
# d=subset(d, select = -c(Cent.DATE))
# Adding WML info
# wml=getWML()
# wml=wml[wml$RID %in% unique(d$RID),]
# d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$wml.DATE,format="%Y-%m-%d")) 
# d=d[d$dt<=180,]
# d=d[!is.na(d$dt),]
# d <- d[order(d$RID, d$EXAMDATE ,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the WML scan closest in date to the d scan
# rm(wml)
# d=subset(d, select = -c(wml.DATE) )

##### Hippocampus volume from Freesurfer Data ########
# reg=c('Volume (WM Parcellation) of LeftHippocampus','Volume (WM Parcellation) of RightHippocampus')
# # d$reg = rowMeans(d[,paste(as.character(c("ST29SV_harmonized","ST88SV_harmonized")))])/d$ICV_harmonized_median*100
# d$reg = rowMeans(d[,paste(as.character(c("ST29SV_harmonized_icv_adjusted","ST88SV_harmonized_icv_adjusted")))])
# 
# dd = d[d$APOE!='NA' & !is.na(d$Age) #& !is.na(d$Centiloid)#& !is.na(d$APOE) & !is.na(d$AmyloidPos)
#        & !is.na(d$reg) & d$OVERALLQC!="Partial",] #we are only keeping data with OVERALLQC= "Pass" or "Hippocampus Only",
# #OVERALLQC ="" and "FAIL" were removed during harmonization process.
# PlotDataBins(d[!is.na(d$reg) & d$APOE !='NA' & !is.na(d$APOE) & !is.na(d$Age),])
# write.csv(dd,"HippocampalVoldataNoCovar1.csv") #Switched on 8/3/2023

WholeDatasetwithAD=read.csv("HippocampalVoldataNoCovar1.csv")
WholeDatasetwithAD$APOE <- factor(WholeDatasetwithAD$APOE)
WholeDatasetwithAD$APOE=relevel(WholeDatasetwithAD$APOE,ref = 'E3/E3')
WholeDatasetwithAD$DXSUM<-factor(WholeDatasetwithAD$DXSUM)
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
tmpdd <- tmpdd[order(tmpdd$RID, tmpdd$EXAMDATE),]
tmpdd <- tmpdd[!duplicated(tmpdd[c("RID","EXAMDATE")]),]

tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups

#Limiting data by mins and max of individual groups in the dataset
# tmpdd=tmpdd[tmpdd$Age>=62 & tmpdd$Age<=87,] 8/3/2023
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID) this creates errors
tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)
# tmpdd$Edudiff=tmpdd$PTEDUCAT-16 #Since we have been predicting for 16 years of education

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'ICV adjusted Hippocampal Vol (mm^3)',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   group = RID, colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID, size=1.1)) +
#   geom_point(aes(size=3))+ 
#   plt_labs

ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"HippocampalVol"))

if("VolumetricData" %in% names(Participant_Records)){
  Participant_Records$VolumetricData[Participant_Records$RID%in%tmpdd$RID]=1
  if(length(unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]))>0){
    Participant_Records=merge(Participant_Records,data.frame(RID=unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]),VolumetricData=1),by=c("RID","VolumetricData"),all = TRUE)
  }
}else{
  Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],VolumetricData=1),by=c("RID"),all = TRUE)
}

Save_demog_distribution_withintmpdd(tmpdd,"HippocampalVol")

jpeg("HippocampalVoldataScatterPTEDUCATCovars.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()


# plotAPOEnSexinAgeBins(tmpdd,"") These are not saved here because the MetaROI Vol data set in inclusive of these
# Save_demog_distribution_withintmpdd(tmpdd,"")

#Simple Model
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE,k=4) + 
            s(RID,Agediff, bs="re")+ s(RID,bs="re")  + #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex + s(RID,bs="re"), data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"HippocampalVolPTEDUCATCovarModelWOGender_ICV_adjusted.Rds")
summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,Sex='F',RID=NULL))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("Hippocampal Vol prediction within ",dxtext,'without Sex effects'), 
#        x ="Age (years)", y = 'Hippocampal Volume (mm^3)')
# ggp1

# Simple Model
m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex),k=4) + 
            s(RID, Agediff, bs = "re") + 
            PTEDUCAT + APOE + Sex + APOE:Sex, 
          data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"HippocampalVolPTEDUCATCovarModelWithGender_ICV_adjusted.Rds")
# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) +
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("Hippocampal Vol prediction within ",dxtext),
#     x ="Age (years)", y = '')#Harmonized ICV adjusted Hippocampal Volume  (mm^3)')
# ggp2 + facet_wrap(vars(Sex))

# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)


m1=readRDS("HippocampalVolPTEDUCATCovarModelWOGender_ICV_adjusted.Rds")
m2=readRDS("HippocampalVolPTEDUCATCovarModelWithGender_ICV_adjusted.Rds")
SigOrNot=CompareSplineSigSlope("Hippocampal Vol","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT",'F')


m1_predict$z=0
m1_predict$z= -(m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="HippoVol"

ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= -(m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="HippoVol"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))

m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])


jpeg("HippocampalVolPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1 + geom_line(aes(size=1))+
#   theme(legend.position="none",axis.text = element_text(size = 60), 
#         axis.title = element_text(size = 60),strip.text = element_text(size =60)) +
#   (ggp2+facet_wrap(vars(Sex)))+
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("Hippocampal Volume")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

##### MetaROI volume from Freesurfer Data ########

# # MetaROI thickness analysis based on previous code Individual ROIs
# # 'Thickness.ctxlhentorhinal'ST24CV_harm,'Thickness.ctxrhentorhinal'ST83CV_harm,
# # 'Thickness.ctxlhinferiortemporal'ST32CV_harm,'Thickness.ctxrhinferiortemporal'ST91CV_harm,
# # 'Thickness.ctxlhmiddletemporal'ST40CV_harm,'Thickness.ctxrhmiddletemporal'ST99CV_harm,
# # 'Thickness.ctxlhfusiform'ST26CV_harm,'Thickness.ctxrhfusiform'ST85CV_harm
# 
# reg=c('LeftMeta_harmonized_icv_adj','RightMeta_harmonized_icv_adj')
# # "ST24CV_harmonized_icv_adjusted", "ST32CV_harmonized_icv_adjusted",
# # "ST40CV_harmonized_icv_adjusted", "ST26CV_harmonized_icv_adjusted"
# # "ST83CV_harmonized_icv_adjusted", "ST91CV_harmonized_icv_adjusted",
# # "ST99CV_harmonized_icv_adjusted", "ST85CV_harmonized_icv_adjusted"
# d$reg = rowMeans(d[,paste(as.character(reg))])
# dd = d[d$APOE!='NA' & !is.na(d$Age) #& !is.na(d$Centiloid)#& !is.na(d$APOE) & !is.na(d$AmyloidPos)
#        & !is.na(d$reg) & d$OVERALLQC!="Hippocampus Only",] #we are only keeping data with OVERALLQC= "Pass", or "Partial" where TMPQC has passed
# #PlotDataBins(d[!is.na(d$reg) & d$APOE !='NA' & !is.na(d$APOE) & !is.na(d$Age),])
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM=factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$EXAMDATE),]
# dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),]
# 
# write.csv(dd,"MetaROIVoldataNoCovar1.csv") #Written on 8/3/2023

WholeDatasetwithAD=read.csv("MetaROIVoldataNoCovar1.csv")
WholeDatasetwithAD$APOE <- factor(WholeDatasetwithAD$APOE)
WholeDatasetwithAD$APOE=relevel(WholeDatasetwithAD$APOE,ref = 'E3/E3')
WholeDatasetwithAD$DXSUM<-factor(WholeDatasetwithAD$DXSUM)
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
tmpdd <- tmpdd[order(tmpdd$RID, tmpdd$EXAMDATE),]
tmpdd <- tmpdd[!duplicated(tmpdd[c("RID","EXAMDATE")]),]


tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups

#Limiting data by mins and max of individual groups in the dataset
# tmpdd=tmpdd[tmpdd$Age>=62 & tmpdd$Age <=87,] 8/3/2023
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID)This causeserror in model
tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'ICV adjusted MetaROI Vol',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   group = RID, colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID,size=1.1)) +
#   geom_point(aes(size=3))+ 
#   plt_labs

ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"MetaROIVol"))

if("VolumetricData" %in% names(Participant_Records)){
  Participant_Records$VolumetricData[Participant_Records$RID%in%unique(tmpdd$RID)]=1
  if(length(unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]))>0){
    Participant_Records=merge(Participant_Records,data.frame(RID=unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]),VolumetricData=1),by=c("RID","VolumetricData"),all = TRUE)
  }
}else{
  Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],VolumetricData=1),by=c("RID"),all = TRUE)
}

jpeg("MetaROIVoldataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()


plotAPOEnSexinAgeBins(tmpdd,"VolumetricMRI")
Save_demog_distribution_withintmpdd(tmpdd,"MetaROIVolData")



#Simple Model
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE,k=5) + #8/3/2023 added k=5 as otherwise the fit is too wobbly while trying to fit the data.
            s(RID,Agediff, bs="re") + s(RID,bs="re")  +  #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex, 
          data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"MetaROIVolPTEDUCATCovarModelWOGender_ICVadjusted.Rds")
# summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL, Sex='F'))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age

# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("MetaROI Vol prediction within ",dxtext,'without Sex effects with CDSOB'), 
#        x ="Age (years)", y = 'MetaROI Volume (mm^3)')
# ggp1

# Simple Model
m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex),k=5) + 
            s(RID, Agediff, bs = "re") + 
            PTEDUCAT + APOE +Sex+ APOE:Sex, 
          data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"MetaROIVolPTEDUCATCovarModelWithGender_ICVadjusted.Rds")
# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
                        values=list(PTEDUCAT=16,RID=NULL))
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)#Getting back from age difference to actual age

# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("MetaROI Vol predictions within ",dxtext,'with CDSOB'),
#     x ="Age (years)", y = '')#Harmonized ICV adjusted MetaROI Volume (mm^3)')
# ggp2 + facet_wrap(vars(Sex))




m1=readRDS("MetaROIVolPTEDUCATCovarModelWOGender_ICVadjusted.Rds")
m2=readRDS("MetaROIVolPTEDUCATCovarModelWithGender_ICVadjusted.Rds")
SigOrNot=CompareSplineSigSlope("MetaROI Vol","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT","F")


m1_predict$z=0
m1_predict$z= -(m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="MetaROIVol"

ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= -(m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="MetaROIVol"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))

m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])


jpeg("MetaROIVolPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1 + geom_line(aes(size=1))+
#   theme(legend.position="none",axis.text = element_text(size = 60), 
#         axis.title = element_text(size = 60),strip.text = element_text(size =60)) +
#   (ggp2+facet_wrap(vars(Sex)))+
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("MetaROI volume")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

# ##### Anterior Cingulate Cortex Volume ----
# # reg=c('LeftAnteriorCingulate_harmonized_icv_adj','RightAnteriorCingulate_harmonized_icv_adj')
# # #"ST113CV", "ST73CV","ST54CV","ST14CV"
# # 
# # d$reg = rowMeans(d[,paste(as.character(reg))])
# # dd = d[d$APOE!='NA' & !is.na(d$Age) #& !is.na(d$Centiloid)#& !is.na(d$APOE) & !is.na(d$AmyloidPos) 
# #         & !is.na(d$reg) & d$OVERALLQC=="Pass",] #we are only keeping data with OVERALLQC= "Pass"
# # #PlotDataBins(d[!is.na(d$reg) & d$APOE !='NA' & !is.na(d$APOE) & !is.na(d$Age),])
# # dd$APOE <- factor(dd$APOE)
# # dd$DXSUM=factor(dd$DXSUM)
# # dd <- dd[order(dd$RID, dd$EXAMDATE),]
# # dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),]
# # write.csv(dd,"ACCVoldataNoCovar.csv")
# 
# # WholeDatasetwithAD=read.csv("ACCVoldataNoCovar.csv")# ACCVol Not used----
# # WholeDatasetwithAD$APOE <- factor(WholeDatasetwithAD$APOE)
# # WholeDatasetwithAD$APOE=relevel(WholeDatasetwithAD$APOE,ref = 'E3/E3')
# # WholeDatasetwithAD$DXSUM<-factor(WholeDatasetwithAD$DXSUM)
# # tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
# # tmpdd <- tmpdd[order(tmpdd$RID, tmpdd$EXAMDATE),]
# # tmpdd <- tmpdd[!duplicated(tmpdd[c("RID","EXAMDATE")]),]
# # 
# # tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# # 
# # #Limiting data by mins and max of individual groups in the dataset
# # tmpdd=tmpdd[tmpdd$Age <=87 & tmpdd$Age>=62,]
# # tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
# #                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# # # tmpdd$RID<-factor(tmpdd$RID) causes error in model
# # 
# # tmpdd$APOE <- factor(tmpdd$APOE)
# # tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# # tmpdd$Sex=factor(tmpdd$Sex)
# # tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# # tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# # tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# # tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)
# # 
# # tmpdd$NumOfScans=c(1)
# # tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
# # for (id in unique(tmpdd$RID))
# # {
# #   tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
# # }
# # plt_labs <- labs(y = 'ICV adjusted AnteriorCingulate Vol',
# #                  x = 'Age in Years',
# #                  colour = 'APOE')
# # ggp=ggplot(tmpdd, aes(x = Age, y = reg,
# #                   group = RID, colour = APOE)) +
# #   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
# #   geom_point()+ 
# #   plt_labs
# # 
# # ggp + facet_wrap(vars(Sex,APOE))
# # # ggp + facet_wrap(vars(Sex)) 
# # # ggp + facet_wrap(vars(APOE))
# # 
# # jpeg("AnteriorCingulateVoldataScatter.jpg",width=700,height = 500)
# # ggp+facet_wrap(vars(Sex,APOE))
# # dev.off()
# # 
# # PlotCDSOBBins(tmpdd,"ACCVol")
# # 
# # PrintAmyloidPosN(tmpdd)
# # min(tmpdd$Age)
# # max(tmpdd$Age)
# # # PrintN(tmpdd)
# # # PlotScatter(tmpdd,'Age (years)','AnteriorCingulate Vol')
# # # PlotDataBins(tmpdd)
# # 
# # #Simple Model
# # m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE) + 
# #             s(RID,Agediff, bs="re")+ s(RID,bs="re")  + #Accounting for random intercept and random slope by individual subjects' repeated measures
# #             APOE + Sex, 
# #           data = tmpdd, method = "REML")
# # m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
# # saveRDS(m1,"AnteriorCingulateVolNoCovarModelWOGender_ICVadjusted.Rds")
# # summary(m1)
# # 
# # #predict_gam
# # m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
# #                         values=list(Sex='F',RID=NULL))
# # m1_predict$APOE=factor(m1_predict$APOE)
# # m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# # m1_predict$Age=m1_predict$Agediff++min(tmpdd$Age)
# # ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
# #   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
# #   labs(title=paste("Anterior Cingulate Vol prediction within ",dxtext,'without Sex effects with CDSOB'), 
# #        x ="Age (years)", y = 'Harmonized ICV adjusted Anterior Cingulate Volume (mm^3)')
# # ggp1
# # 
# # #Modelling interactions without separating the groups
# # # Simple Model
# # m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex)) + 
# #             s(RID, Agediff, bs = "re") + 
# #             APOE + Sex + APOE:Sex, 
# #           data = tmpdd, method = "REML")
# # m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
# # saveRDS(m2,"AnteriorCingulateVolNoCovarModelWithGender_ICVadjusted.Rds")
# # 
# # summary(m2)
# # 
# # #predict_gam
# # m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
# #                         values=list(RID=NULL))
# # m2_predict$APOE=factor(m2_predict$APOE)
# # m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# # m2_predict$Sex=factor(m2_predict$Sex)
# # m2_predict$Sex = relevel(m2_predict$Sex,'M')
# # m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# # ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
# #   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + labs(title=paste("AnteriorCingulate Vol predictions within ",dxtext,'with CDSOB'), 
# #                                                         x ="Age (years)", y = 'Harmonized ICV adjusted AnteriorCingulate Volume (mm^3)')
# # ggp2 + facet_wrap(vars(Sex))
# # 
# # jpeg("AnteriorCingulateVolNoCovarModelPredictionPlots.jpg",width=2000,height = 500)
# # ggp1+ (ggp2+facet_wrap(vars(Sex)))
# # dev.off()
# # m1=readRDS("AnteriorCingulateVolNoCovarModelWOGender_ICVadjusted.Rds")
# # m2=readRDS("AnteriorCingulateVolNoCovarModelWithGender_ICVadjusted.Rds")
# # SigOrNot=CompareSplineSigSlope("Anterior Cingulate Vol","NoCovar",tmpdd,m1,m2,"","F")
# # 
# # m1_predict$z=0
# # m1_predict$z= -(m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
# # m1_predict$Significance=0
# # m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
# # m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
# # m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
# # ggplot(m1_predict, aes(x = Age, y = z, colour = APOE, size=Significance)) +
# #   geom_line()+ 
# #   plt_labs 
# # m1_predict$Measure="AnteriorCingulateVol"
# # m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])
# # 
# # m2_predict$z=0
# # m2_predict$z= -(m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
# # m2_predict$Significance=0
# # m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
# # m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
# # m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
# # m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
# # m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
# # m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
# # m2_predict$Measure="AnteriorCingulateVol"
# # ggplot(m2_predict, aes(x = Age, y = z, colour = APOE, size=Significance)) +
# #   geom_line()+ 
# #   plt_labs + facet_wrap(vars(Sex))
# # m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])
# # 
# # 
# # AIC(m1)
# # BIC(m1)
# # AIC(m2)
# # BIC(m2)
# # summary(m1)
# # summary(m2)
# 
##### Individual ROI from MetaROI----
# MetaROI thickness analysis based on previous code Individual ROIs
# 'Thickness.ctxlhentorhinal'ST24CV_harm,'Thickness.ctxrhentorhinal'ST83CV_harm,
# 'Thickness.ctxlhinferiortemporal'ST32CV_harm,'Thickness.ctxrhinferiortemporal'ST91CV_harm,
# 'Thickness.ctxlhmiddletemporal'ST40CV_harm,'Thickness.ctxrhmiddletemporal'ST99CV_harm,
# 'Thickness.ctxlhfusiform'ST26CV_harm,'Thickness.ctxrhfusiform'ST85CV_harm


# FDG for Meta-ROI metabolism ----

# Old data file
# # pet=read.csv('UCBERKELEYFDG_05_28_20.csv')
# # colnames(pet)[colnames(pet)=='EXAMDATE']='FDG.DATE'
# # # Now changing different ROI rows to columns
# # d=pet[,c('RID','VISCODE','VISCODE2','FDG.DATE')]
# # d=unique(d)
# # tempreg=pet[pet$ROINAME=='CingulumPost',c('RID','FDG.DATE','MEAN')]
# # colnames(tempreg) <- c('RID','FDG.DATE','FDG.BilateralCingulumPosterior')
# # d=merge(d,tempreg,by=c('RID','FDG.DATE'))
# # tempreg=pet[pet$ROINAME=='Temporal' & pet$ROILAT=='Left',c('RID','FDG.DATE','MEAN')]
# # colnames(tempreg) <- c('RID','FDG.DATE','FDG.LeftTemporal')
# # d=merge(d,tempreg,by=c('RID','FDG.DATE'))
# # tempreg=pet[pet$ROINAME=='Temporal' & pet$ROILAT=='Right',c('RID','FDG.DATE','MEAN')]
# # colnames(tempreg) <- c('RID','FDG.DATE','FDG.RightTemporal')
# # d=merge(d,tempreg,by=c('RID','FDG.DATE'))
# # tempreg=pet[pet$ROINAME=='Angular' & pet$ROILAT=='Left',c('RID','FDG.DATE','MEAN')]
# # colnames(tempreg) <- c('RID','FDG.DATE','FDG.LeftAngular')
# # d=merge(d,tempreg,by=c('RID','FDG.DATE'))
# # tempreg=pet[pet$ROINAME=='Angular' & pet$ROILAT=='Right',c('RID','FDG.DATE','MEAN')]
# # colnames(tempreg) <- c('RID','FDG.DATE','FDG.RightAngular')
# # d=merge(d,tempreg,by=c('RID','FDG.DATE'))
# # rm(pet,tempreg)
# # d <- d[!is.na(d$FDG.BilateralCingulumPosterior) | !is.na(d$FDG.LeftTemporal) | !is.na(d$FDG.RightTemporal)  | !is.na(d$FDG.LeftAngular) | !is.na(d$FDG.RightAngular),]
# # d <- d[order(d$RID, d$FDG.DATE),]
# # # d <- d[!duplicated(d[c('RID','FDG.DATE')]),]
#  

# pet=read.csv("UCBERKELEYFDG_03_02_22.csv")
# pet=read.csv("UCBERKELEYFDG_8mm_02_17_23_11Jun2023.csv")
# pet1=pet[pet$ROINAME =="MetaROI",c("RID","EXAMDATE","MEAN","TOTVOX")]#was metaroi, changed ROI names to current on 8/3/20203 as per latest file
# pet2=pet[pet$ROINAME =="Top50PonsVermis",c("RID","EXAMDATE","MEAN","TOTVOX")]#pons-vermis
# colnames(pet1)[colnames(pet1)=='MEAN']='FDGmetaROI'
# colnames(pet1)[colnames(pet1)=='TOTVOX']='metaROITOTVOX'
# colnames(pet2)[colnames(pet2)=='MEAN']='PonsVermis'
# colnames(pet2)[colnames(pet2)=='TOTVOX']='PonsVermisTOTVOX'
# # pet2$PonsVermisIntensPerVox=pet2$PonsVermis/pet2$PonsVermisTOTVOX
# # pet2$AdjFDGmetaROI=median(pet2$PonsVermisIntensPerVox)-pet2$PonsVermisIntensPerVox
# d <- merge(pet1,pet2,by=c("RID","EXAMDATE"),all.x=T)
# d$AdjFDGmetaROI= d$FDGmetaROI/d$PonsVermis#(d$AdjFDGmetaROI+(d$FDGmetaROI/d$metaROITOTVOX))*d$metaROITOTVOX
# rm(pet, pet1,pet2)
# #names(d)[names(d)=="EXAMDATE"]="FDG.DATE"
# 
# apoe=getAPOE()
# apoe=apoe[apoe$RID %in% unique(d$RID),]
# d=merge(d,apoe,by=c('RID'),all.x=T) #Only  keeping RIDs that have the APOE data
# rm(apoe)
# d=d[!is.na(d$APOE),]
# 
# Adni=adnimerge
# colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
# Adni$Sex=factor(Adni$Sex)
# levels(Adni$Sex)=c('F','M')
# Adni$Sex = relevel(Adni$Sex,ref = 'M')
# colnames(Adni)[colnames(Adni)=='DX']='DXSUM'
# colnames(Adni)[colnames(Adni)=='AGE']='Age'
# colnames(Adni)[colnames(Adni)=='EXAMDATE']='Adnimerge.DATE'
# colnames(Adni)[colnames(Adni)=='EXAMDATE.bl']='Adnimerge.DATE.bl'
# levels(Adni$DXSUM)[levels(Adni$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# Adni=Adni[,c("RID","Adnimerge.DATE","Age","Years.bl","Sex","DXSUM","PTEDUCAT","CDSOB")]
# 
# d=merge(d,Adni,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d"))
# d=d[!is.na(d$CDSOB) & !is.na(d$PTEDUCAT) & !is.na(d$Age) & !is.na(d$DXSUM) & !is.na(d$Sex) ,]
# d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# d$Age=d$Age+d$Years.bl+(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Adnimerge.DATE,format="%Y-%m-%d"))/365
# d=subset(d, select = -c(Adnimerge.DATE,Years.bl))
# d$Sex=factor(d$Sex)
# d$DXSUM=factor(d$DXSUM)
# rm(Adni)

#Alternate way of getting demographic data ----
#Decided on ADNIMERGE instead
# dem=getDemog()
# dem=dem[dem$RID %in% unique(d$RID),]
# d=merge(d,dem,by=c('RID'),all.x=T)
# rm(dem)
# d$Age=as.numeric(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$refDOB,format="%Y-%m-%d"))/365 #Age at FDG scan
# d=d[!is.na(d$Age),]
# d=subset(d, select = -c(refDOB) )
# 
# cdr=getCDR()
# cdr=cdr[cdr$RID %in% unique(d$RID),]
# d=merge(d,cdr,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$CDR.DATE,format="%Y-%m-%d")) #Days since FDG scan
# d=d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d <- d[order(d$RID, d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the FDG scan
# d=d[!is.na(d$CDSOB),]
# d=subset(d, select = -c(CDR.DATE) )
# rm(cdr)
# 
# dxsm=getDXSUM()
# dxsm=dxsm[dxsm$RID %in% unique(d$RID),]
# d=merge(d,dxsm,by=c('RID'),all.x=T)
# d=d[!is.na(d$DXSUM),]
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$DX.DATE,format="%Y-%m-%d")) 
# d<-d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d<-d[!is.na(d$dt),]
# d <- d[order(d$RID, d$EXAMDATE, d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the DXSUM closest in date to the FDG scan
# d=subset(d, select = -c(DX.DATE,dt) )
# rm(dxsm)
# 

# not adding WML and Amyloid info 8/3/2023 ----
# # Adding AmyloidPos info
# CentiloidData=getCentiloidData()
# CentiloidData=CentiloidData[CentiloidData$RID %in% unique(d$RID),]
# d=merge(d,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d")) 
# d=d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d<-d[!is.na(d$dt),]
# d <- d[order(d$RID,d$EXAMDATE, d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the FDG scan
# rm(CentiloidData)
# d=subset(d, select = -c(Cent.DATE))
# 
# wml=getWML()
# wml=wml[wml$RID %in% unique(d$RID),]
# d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$wml.DATE,format="%Y-%m-%d")) 
# d=d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# d=d[!is.na(d$dt),]
# d <- d[order(d$RID,d$EXAMDATE, d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the WML scan closest in date to the FDG scan
# rm(wml)
# d=subset(d, select = -c(wml.DATE))
# # ----

# # d$FDGmetaROI <- (2*d$FDG.BilateralCingulumPosterior 
# #                  + d$FDG.LeftTemporal + d$FDG.RightTemporal 
# #                  + d$FDG.LeftAngular + d$FDG.RightAngular)/6
# # reg = 'FDGmetaROI'
# 
# reg = 'AdjFDGmetaROI'#MetaROI values adjusted with same values that would adjust the Pons_vermis FDG to the median FDG/Voxel values
# d$reg = d[,paste(as.character(reg))]
# dd = d[d$APOE!='NA' & !is.na(d$Age) & !is.na(d$APOE) & !is.na(d$reg),]# & !is.na(d$Centiloid),]

# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# dd <- dd[order(dd$RID,dd$EXAMDATE),]
# dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),]
# dd$Age = as.numeric(dd$Age)
# write.csv(dd,"MetaROIFDGdataNoCovar1.csv") #New file saved on 8/3/2023
# rm(dd)

WholeDatasetwithAD = read.csv("MetaROIFDGdataNoCovar1.csv") #New file saved on 8/3/2023
WholeDatasetwithAD$APOE <- factor(WholeDatasetwithAD$APOE)
WholeDatasetwithAD$APOE=relevel(WholeDatasetwithAD$APOE,ref = 'E3/E3')
WholeDatasetwithAD$DXSUM<-factor(WholeDatasetwithAD$DXSUM)
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
tmpdd <- tmpdd[order(tmpdd$RID, tmpdd$EXAMDATE),]
tmpdd <- tmpdd[!duplicated(tmpdd[c("RID","EXAMDATE")]),]

tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# tmpdd=tmpdd[tmpdd$Age>=62 & tmpdd$Age<=88,]
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]

tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'MetaROI FDG SUVR',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   # group = RID, #Switching to Cross sectional analysis
#                   colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID, size=1.1)) +
#   geom_point(aes(size=3))+ 
#   plt_labs

ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"FDGMetaROI"))
Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],FDGPET=1),by=c("RID"),all = TRUE)


jpeg("MetaROIFDGSUVRdataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()

plotAPOEnSexinAgeBins(tmpdd,"FDGPET")
Save_demog_distribution_withintmpdd(tmpdd,"FDGPET")


# PlotCDSOBBins(tmpdd,"MetaROIFDG")
# 
# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)

# PrintN(tmpdd)
# PlotScatter(tmpdd,'Age (years)','WML')
# PlotDataBins(tmpdd)


#Simple Model
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE) + 
            s(RID,Agediff, bs="re")+ s(RID,bs="re")  + #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex, data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelWOGender.Rds")

# summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL,Sex='F'))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("Adj MetaROI FDG Uptake prediction within ",dxtext,'without Sex effects'), 
#        x ="Age (years)", y = 'MetaROI FDG SUVR' )# adjusted with PonsVermis')
# ggp1

#Simple Model
m2 <- gam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex)) + 
            s(RID, Agediff, bs = "re") + 
            PTEDUCAT + APOE + Sex + APOE:Sex,
          data = tmpdd,
          method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelWithGender.Rds")
# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
                        values=list(PTEDUCAT=16,RID=NULL))
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("Adj MetaROI FDG Uptake prediction within ",dxtext),
#     x ="Age (years)", y = '')#MetaROI FDG SUVR')# adjusted with PonsVermis')
# ggp2 + facet_wrap(vars(Sex))



m1=readRDS("MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelWOGender.Rds")
m2=readRDS("MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelWithGender.Rds")
SigOrNot=CompareSplineSigSlope("MetaROI FDG SUVR","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT","F")#,PTEDUCAT,WML,Centiloid")
# "MetaROI FDG Uptake Adj With PonsVermis"

m1_predict$z=0
m1_predict$z= -(m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="FDG"


ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))


m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= -(m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="FDG"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))

m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])


jpeg("MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1 + geom_line(aes(size=1))+
#   theme(legend.position="none",axis.text = element_text(size = 60), 
#         axis.title = element_text(size = 60),strip.text = element_text(size =60)) +
#   (ggp2+facet_wrap(vars(Sex)))+
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("MetaROI FDG SUVR")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

# reg='FDG.BilateralCingulumPosterior'
# d$reg = d[,paste(as.character(reg))]
# dd = d[d$APOE!='NA' & !is.na(d$Age) & !is.na(d$APOE) & d$DXSUM%in%dx & !is.na(d$reg) & !is.na(d$Centiloid),]
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <-factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$Age),]
# dd <- dd[!duplicated(dd$RID),]
# ceiling(min(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),max,na.rm=T)))
# floor(max(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),min,na.rm=T)))
# PlotScatter(dd,'Age (Years)','Bilateral posterior cingulum metabolism') #This function will need TAU.Age changed to Age
# dd = dd[dd$Age >= floor(max(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),min,na.rm=T))) &
#           dd$Age <= ceiling(min(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),max,na.rm=T))),]
# PrintN(dd)
# b1 <- mgcv::gam(reg ~ s(Age, bs='cs') 
#                 + s(Age, bs='cs', by = interaction(APOE,Sex)) + s(Centiloid, bs='cs') 
#                 + APOE + Sex + PTEDUCAT + CDSOB + Centiloid, 
#                 data = dd, method="REML")#removed +ICV
# b1p <- predict_gam(b1, values=list(Centiloid=20, PTEDUCAT=16, CDSOB=cdr.sob))
# ggp=ggplot(b1p, aes(Age, fit, color = APOE,linetype = Sex)) +
#   geom_smooth_ci( ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("Bilateral posterior cingulum metabolism APOE x Sex effect \n within ",dxtext), 
#        x ="Age (years)", y = 'Bilateral posterior cingulum metabolism')
# ggp + facet_wrap(vars(Sex))
# summary(b1)
#                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           summary(b1)
# rm(d,dd,b1,b1p,ggp)


# Tau Braak Region 1 ----
### d_old=read.csv('ADNImerge-TAU-wCogn-ROIs.csv') 
# ##d=read.csv('UCBERKELEYAV1451_PVC_01_12_22.csv')
# #d=read.csv('UCBERKELEYAV1451_PVC_04_29_22.csv')

# d=read.csv('UCBERKELEYAV1451_PVC_8mm_02_17_23.csv')

# Adni=adnimerge
# colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
# Adni$Sex=factor(Adni$Sex)
# levels(Adni$Sex)=c('F','M')
# Adni$Sex = relevel(Adni$Sex,ref = 'M')
# colnames(Adni)[colnames(Adni)=='DX']='DXSUM'
# colnames(Adni)[colnames(Adni)=='AGE']='Age'
# colnames(Adni)[colnames(Adni)=='EXAMDATE']='Adnimerge.DATE'
# colnames(Adni)[colnames(Adni)=='EXAMDATE.bl']='Adnimerge.DATE.bl'
# levels(Adni$DXSUM)[levels(Adni$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# Adni_demog=Adni[,c("RID","Adnimerge.DATE","Age","Years.bl","Sex","PTEDUCAT")]
# Adni_demog=Adni_demog[order(Adni_demog$RID,Adni_demog$Adnimerge.DATE),]
# Adni_demog=Adni_demog[!duplicated(Adni_demog$RID),]
# Adni_DXSUM_CDSOB=Adni[,c("RID","Adnimerge.DATE","DXSUM","CDSOB")]
# 
# d=merge(d,Adni_demog,by=c('RID')) #,all.x=T
# d$Age=d$Age+d$Years.bl+(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Adnimerge.DATE,format="%Y-%m-%d"))/365
# d=subset(d, select = -c(Adnimerge.DATE,Years.bl))
# 
# 
# #Do not need CDSOB info. Will add DXSUM independently.
# #Adnimerge might not be the best way to add DXSUM and CDSOB info as there are a lot of CDSOB and DXSUM data missing from ADNIMERGE, 
# # and the Adnimerge dates do not necessarily represent the actual EXAMDATEs for those measures
# # d=merge(d,Adni_DXSUM_CDSOB,by=c('RID'),all.x=T)
# # d$dt=abs(as.Date(d$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d"))
# # d=d[!is.na(d$CDSOB) & !is.na(d$PTEDUCAT) & !is.na(d$Age) & !is.na(d$DXSUM) & !is.na(d$Sex) ,]
# # d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# # d$Sex=factor(d$Sex)
# # d$DXSUM=factor(d$DXSUM)
# # rm(Adni)
# # 
# # d$NumOfScans=c(1)
# # for (id in unique(d$RID))
# # {
# #   d$NumOfScans[d$RID==id]=length(d$RID[d$RID==id])
# #   
# # }
# # d=d[!(d$dt>180 & d$NumOfScans==1),]
# # # # Following was performed to see which records with dt>180 to keep
# # d1=d[d$RID %in% d$RID[d$dt>180],]
# # d1$dt=as.Date(d1$Adnimerge.DATE,format="%Y-%m-%d")-as.Date(d1$EXAMDATE,format="%Y-%m-%d")
# # #Based on manual inspection of d1, following were to be removed from d
# # #because the DXSUM and CDSOB could not be assumed/accpted based on the preior or subsequect records of that subject
# # # RID   EXAMDATE
# # # 416 2019-09-24
# # # 1378 2019-04-30
# # # 4198 2016-09-20
# # # 4429 2020-01-09 
# # # 6062 2021-01-14
# # # 6178 2020-09-08 
# # # 6206 2019-11-12
# # # 6251 2021-06-22
# # # 6251 2022-07-13
# # # 6297 2022-04-19
# # # 6347 2021-05-26 Thought this was AD from the record before, the CDSOB score can't be taken from 244 days before
# # # 6385 2021-05-25 
# # # 6493 2022-04-05  
# # # 6629 2022-01-25 
# # # 6837 2021-01-21
# # 
# # 
# # # 4643 2018-08-21 was left because the dxsum remained CN throughout the years, though at this timepoint the CDSOB has gone up to 0.5
# # 
# # d=d[!(d$RID==4198 & d$EXAMDATE=="2016-09-20")&!(d$RID==4301 & d$EXAMDATE=="2019-03-12")&
# #       !(d$RID==6229 & d$EXAMDATE=="2022-04-26")&!(d$RID==6297 & d$EXAMDATE=="2022-04-19")&
# #       !(d$RID==6810 & d$EXAMDATE=="2021-12-01"),]
# 
# 
# ##Adding DXSUM 
# dxsm=getDXSUM()
# dxsm=dxsm[dxsm$RID %in% unique(d$RID),]
# d=merge(d,dxsm,by=c('RID'),all.x=T)
# d$dt=abs(as.Date(d$DX.DATE,format="%Y-%m-%d")-as.Date(d$EXAMDATE,format="%Y-%m-%d"))
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the DXSUM closest in date to the Tau PET scan 
# d=d[!is.na(d$DXSUM),]
# # #based on manual inspection of d1 similar to above, the RIDs with multiple records following 
# # #were eliminated based on whether the diagnosis was stable from the 
# # #time-point before and after the record with dt>180
# # #RID 4198 EXAMDATE 2016-09-20 as the sub diagnosis changed from CN on 2016-02-02 to MCI on 2017-11-09
# # #RID 4301 EXAMDATE 2019-03-12 as last record for DXSUM as MCI was on 2018-03-19
# # #RID 6229 EXAMDATE 2022-04-26 as last record for DXSUM as MCI was on 2021-03-31
# # #RID 6297 EXAMDATE 2022-04-19 as no record for DXSUM was noted after CN record on 2021-04-14 after being MCI on all previous instances
# # #RID 6810 EXAMDATE 2021-12-01 has record for DXSUM on 2021-01-18 was noted as AD and record on 2022-07-11 reports MCI
# # 
# # # Records such as following was included back in
# # #RID 6428 EXAMDATE 2019-06-26 has record as MCI on 2018-07-11 as well as on 2020-06-16
# # # along with Records where dt was >180 but a later DXSUM record reports CN, or an earlier DXSUM reported AD
# # 
# d=d[!(d$RID==722 & d$EXAMDATE=="2015-10-20")& !(d$RID==896  & d$EXAMDATE=="2018-11-14")&
#     !(d$RID==2389 & d$EXAMDATE=="2018-01-16")&!(d$RID==4198 & d$EXAMDATE=="2016-09-20")&
#       !(d$RID==2389 & d$EXAMDATE=="2018-01-16")&!(d$RID==4198 & d$EXAMDATE=="2016-09-20")&
#       !(d$RID==4301 & d$EXAMDATE=="2019-03-12")&!(d$RID==4349 & d$EXAMDATE=="2022-04-26")&
#       !(d$RID==4382 & d$EXAMDATE=="2020-09-01")&!(d$RID==4429 & d$EXAMDATE=="2020-01-09")&
#       !(d$RID==4952 & d$EXAMDATE=="2019-06-05")&!(d$RID==4974 & d$EXAMDATE=="2017-10-12")&
#       !(d$RID==5079 & d$EXAMDATE=="2017-11-15")&!(d$RID==5127 & d$EXAMDATE=="2016-04-20")&
#       !(d$RID==5269 & d$EXAMDATE=="2016-08-16")&!(d$RID==6062 & d$EXAMDATE=="2021-01-14")&
#       !(d$RID==6084 & d$EXAMDATE=="2018-06-12")&!(d$RID==6213 & d$EXAMDATE=="2019-01-30")&
#       !(d$RID==6229 & d$EXAMDATE=="2022-04-26")&!(d$RID==6251 & d$EXAMDATE=="2022-07-13")&
#       !(d$RID==6297 & d$EXAMDATE=="2022-04-19")&!(d$RID==6229 & d$EXAMDATE=="2022-04-26")&
#       !(d$RID==6412 & d$EXAMDATE=="2019-07-24")&!(d$RID==6493 & d$EXAMDATE=="2022-04-05")&
#       !(d$RID==6563 & d$EXAMDATE=="2019-12-19")&!(d$RID==6629 & d$EXAMDATE=="2022-01-25")&
#       !(d$RID==6810 & d$EXAMDATE=="2021-12-01")&!(d$RID==6880 & d$EXAMDATE=="2021-05-12")&
#       !(d$RID==6883 & d$EXAMDATE=="2021-04-08")&!(d$RID==6942 & d$EXAMDATE=="2022-03-22")&
#       !(d$RID==7018 & d$EXAMDATE=="2022-08-17"),]
# 
# d=subset(d, select = -c(DX.DATE,dt) )
# 
# rm(dxsm)
# 
# #Using Adnimerge for demographics and not adding CDR info
# # dem=getDemog()
# # dem=dem[dem$RID %in% unique(d$RID),]
# # dem=dem[!duplicated(dem[,c("RID","refDOB")]),]
# # d=merge(d,dem,by=c('RID'),all.x=T)
# # rm(dem)
# # d$Age=as.numeric(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$refDOB,format="%Y-%m-%d"))/365 #Age at Centiloid scan
# # d=d[!is.na(d$Age),]
# # d=subset(d, select = -c(refDOB,dem.DATE) )
# # 
# # cdr=getCDR() NOT GOING TO ADD CDR AS NOT USED IN MODEL
# # cdr=cdr[cdr$RID %in% unique(d$RID),]
# # d=merge(d,cdr,by=c('RID'),all.x=T)
# # d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$CDR.DATE,format="%Y-%m-%d"))
# # d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the CDR exam closest in date to the Centiloid scan visit
# # d=d[!is.na(d$CDSOB),]
# # #d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# #Instead the EXAMDATES and corresponding CDR.DATEs and CDSOB values were examined. Cases where the CDR values remained unchanged from before Tau EXAMDATE to after EXAMDATE, were given a go.
# #If the CDSOB value changed during that period, as in the following records ,those records were eliminated
# #1378 2019-04-30
# #2200 2016-03-15
# #2234 2019-01-16
# #2389 2018-01-16
# #4198 2017-10-05
# #4216 2020-01-27 
# #... and more
# 
# # d=subset(d, select = -c(CDR.DATE,dt) )
# # rm(cdr)
# # 
# 
# 
# apoe=getAPOE()
# apoe=apoe[apoe$RID %in% unique(d$RID),]
# d=merge(d,apoe,by=c('RID'),all.x=T) #Only  keeping RIDs that have the APOE data
# rm(apoe)
# d=d[!is.na(d$APOE),]
# 
# 
# 
# #not adding amyloid positivity info
# # # d$AmyloidPos = factor(d$AmyloidPos) #no Amyloid positivity info in this file. Will be added later
# # #d=d[d$QC=="",] #No QC column in new file
# # 
# # 
# # CentiloidData=getCentiloidData()
# # d=merge(d,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# # d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d")) 
# # d <- d[d$dt<=180,]
# # d=d[!is.na(d$dt),]
# # d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the d scan
# # rm(CentiloidData)
# # d=subset(d, select = -c(Cent.DATE,dt) )
# # 
# # wml=getWML()
# # d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# # d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$wml.DATE,format="%Y-%m-%d")) 
# # d <- d[d$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
# # d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the WML scan closest in date to the centiloid scan visit
# # rm(wml)
# # d=subset(d, select = -c(wml.DATE) )
# 
# 
# dd = d[d$APOE != 'NA' & !is.na(d$Age),] #& !is.na(d$AmyloidPos),]
# # PlotDataBins(d[d$APOE != 'NA' & !is.na(d$Age) & !is.na(d$AmyloidPos) & !is.na(d$reg),])
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$EXAMDATE),]
# dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),]
# # PlotScatter(dd,'Age (Years)','Braak I region AV1451 SUVR (partial volume corrected)'
# 
# dd$Age=as.numeric(dd$Age) #This is needed for whatever reason.
# write.csv(dd,"TaudataNoCovar1.csv") #Saved on 8/4/2023


WholeDatasetwithAD=read.csv("TaudataNoCovar1.csv")
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]

#regions affected in Braak I stages of AD
# ##reg=c("PVC.AV1451.ctxlhentorhinal","PVC.AV1451.ctxrhentorhinal")
# ##reg=c("AV1451.ctxlhentorhinal","AV1451.ctxrhentorhinal")
# #reg=c("CTX_LH_ENTORHINAL_SUVR","CTX_RH_ENTORHINAL_SUVR")
# #d$reg = rowMeans(d[,paste(as.character(reg))]) #Average of all the ROIs above

tmpdd$reg=tmpdd$CTX_ENTORHINAL_SUVR/tmpdd$INFERIOR_CEREBGM_SUVR#BRAAK1_SUVR=CTX_ENTORHINAL_SUVR

# tmpdd=tmpdd[tmpdd$INFERIOR_CEREBGM_SUVR>0.9,] #This was to remove one weird record that had really high CEREBGM_SUVR
tmpdd=tmpdd[!is.na(tmpdd$reg),]
tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# tmpdd=tmpdd[tmpdd$reg<5,]#removing two outliers RID=4271 reg=8.86598 APOE=E3/E2, and RID=6606 reg=6.3390 APOE=E3/E4
# Tested and approved following on 8/4/2023 
tmpdd=tmpdd[tmpdd$reg<6,]#removed 1 value that seemed seemed like an outlier with RID 4271, APOE=E3/E2
# tmpdd = tmpdd[tmpdd$Age >= 63 & tmpdd$Age <= 87,]
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID)
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'Tau Braak Region I',
                 x = 'Age in Years',
                 colour = 'APOE')
ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"TauBraak1"))

if("TauPET" %in% names(Participant_Records)){
  Participant_Records$TauPET[Participant_Records$RID%in%tmpdd$RID]=1
  if(length(unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]))>0){
    Participant_Records=merge(Participant_Records,data.frame(RID=unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]),TauPET=1),by=c("RID","TauPET"),all = TRUE)
  }
}else{
  Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],TauPET=1),by=c("RID"),all = TRUE)
}

jpeg("TauBraakIdataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()



# PlotCDSOBBins(tmpdd,"TauBraakI") THis is only saved as Tau in the Braak III dataset which has 2 additional subjects

# PrintAmyloidPosN(tmpdd)
# PrintN(tmpdd)
# PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# PlotDataBins(tmpdd)

m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE,k=4) +
            s(RID,Agediff, bs="re")+ s(RID,bs="re")  +  #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex,data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"TauEntorhinalBraakIPTEDUCATCovarModelWOGenderEffects.Rds")


# summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16, Sex='F',RID=NULL))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("Tau Braak I prediction within ",dxtext,'with APOE effects'), 
#        x ="Age (years)", y = 'Tau (AV1451 SUVR) Braak 1 region')#AV1451 PVC SUVR Braak I regions')
# ggp1

m2 <- gam(formula = reg ~ s(Agediff, bs = "cs", 
                            by = interaction(APOE,Sex),k=4) + 
            s(RID, Agediff, bs = "re") + s(RID,bs="re") + 
            PTEDUCAT + APOE + Sex + APOE:Sex, 
          data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"TauEntorhinalBraakIPTEDUCATCovarModelWithGenderEffects.Rds")

# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID,Age)'),length_out = 100,
                        values=list(PTEDUCAT=16, RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
 
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("Braak I Tau prediction within ",dxtext, 'with Sex x APOE effects'), 
#     x ="Age (years)", y = '')#AV1451 PVC SUVR in Braak I regions')
# ggp2 + facet_wrap(vars(Sex))
# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)


m1=readRDS("TauEntorhinalBraakIPTEDUCATCovarModelWOGenderEffects.Rds")
m2=readRDS("TauEntorhinalBraakIPTEDUCATCovarModelWithGenderEffects.Rds")
SigOrNot=CompareSplineSigSlope("TauBraakI","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT",'F')

m1_predict$z=0
m1_predict$z= (m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="TauBraak1"


ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

# To add Tau Positivity threhold in interested ----
#modelling all the CN irrespctive of their APOE group
m1_CN <- bam(formula = reg ~ s(Agediff, bs='cs',k=4) +
               s(RID,Agediff, bs="re")+ s(RID,bs="re")  +  #Accounting for random intercept and random slope by individual subjects' repeated measures
               PTEDUCAT,data = tmpdd[tmpdd$DXSUM=='CN',], method = "REML")
m1_CN$model$Age=m1_CN$model$Agediff+min(tmpdd[tmpdd$DXSUM=='CN',c('Age')])

m1CN_predict<-predict_gam(m1_CN,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                          values=list(PTEDUCAT=16,RID=NULL))
m1CN_predict$Age=m1CN_predict$Agediff+min(tmpdd[tmpdd$DXSUM=='CN',c("Age")])

TauPositivity_threshold=0.95*max(m1CN_predict$fit)

#and then add the following to the ggplot
+ geom_hline(yintercept = TauPositivity_threshold)



m2_predict$z=0
m2_predict$z= (m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="TauBraak1"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))


m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])

jpeg("TauBraakIPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
# ggp1 + geom_line(aes(size=1))+
#   theme(legend.position="none",axis.text = element_text(size = 60), 
#         axis.title = element_text(size = 60),strip.text = element_text(size =60)) +
#   (ggp2+facet_wrap(vars(Sex)))+
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("Tau Braak I Region")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)


#Braak Region 2 Do not use----
# tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]
# 
# #regions affected in Braak II stages of AD is hippocapus and the selectiveness 
# #of this contrast agent in Hippocampus is not good
# 
# tmpdd$reg=tmpdd$BRAAK2_SUVR/tmpdd$INFERIOR_CEREBGM_SUVR
# tmpdd=tmpdd[tmpdd$INFERIOR_CEREBGM_SUVR>0.9,]
# tmpdd=tmpdd[!is.na(tmpdd$reg),]
# tmpdd$APOE <- factor(tmpdd$APOE)
# tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# tmpdd$Sex=factor(tmpdd$Sex)
# tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# tmpdd = tmpdd[tmpdd$Age >= 63 & tmpdd$Age <= 87,]
# tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
#                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# # tmpdd$RID<-factor(tmpdd$RID)
# tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)
# 
# tmpdd$NumOfScans=c(1)
# 
# for (id in unique(tmpdd$RID))
# {
#   tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
#   
# }
# 
# plt_labs <- labs(y = 'Tau Braak Region 2',
#                  x = 'Age in Years',
#                  colour = 'APOE')
# ggp= ggplot(tmpdd, aes(x = Age, y = reg,
#                        # group = RID, #Switching to Cross sectional analysis
#                        colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
#   geom_point()+ 
#   plt_labs
# 
# ggp+ facet_wrap(vars(Sex,APOE))
# # ggp+facet_wrap(vars(Sex)) 
# # ggp+facet_wrap(vars(APOE))
# jpeg("TauBraak2dataScatter.jpg",width=700,height = 500)
# ggp+ facet_wrap(vars(Sex,APOE))
# dev.off()
# 
# # PlotCDSOBBins(tmpdd,"TauBraak2") 
# 
# PrintAmyloidPosN(tmpdd)
# # PrintN(tmpdd)
# # PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# # PlotDataBins(tmpdd)
# 
# m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE) +
#             s(RID,Agediff, bs="re")+s(RID,bs="re")  +  #Accounting for random intercept and random slope by individual subjects' repeated measures
#             APOE + Sex,data = tmpdd, method = "REML")
# m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
# saveRDS(m1,"TauBraak2NoCovarModelWOGenderEffects.Rds")
# 
# 
# summary(m1)
# 
# #predict_gam
# m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
#                         values=list(Sex='F',RID=NULL))
# m1_predict$APOE=factor(m1_predict$APOE)
# m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("Tau Braak 2 prediction within ",dxtext,'with APOE effects'), 
#        x ="Age (years)", y = 'AV1451 PVC SUVR Braak 2 regions')
# ggp1
# 
# 
# m2 <- gam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex)) + 
#             s(RID, Agediff, bs = "re") + 
#             APOE+Sex+APOE:Sex, 
#           data = tmpdd, method = "REML")
# m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
# saveRDS(m2,"TauBraak2NoCovarModelWithGenderEffects.Rds")
# 
# summary(m2)
# 
# #predict_gam
# m2_predict<-predict_gam(m2,exclude_terms = c('s(RID,Age)'),length_out = 100,
#                         values=list(RID=NULL))
# # m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
# m2_predict$APOE=factor(m2_predict$APOE)
# m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# m2_predict$Sex=factor(m2_predict$Sex)
# m2_predict$Sex = relevel(m2_predict$Sex,'M')
# m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + labs(title=paste("Braak 2 Tau prediction within ",dxtext, 'with Sex x APOE effects'), 
#                                                         x ="Age (years)", y = 'AV1451 PVC SUVR in Braak 2 regions')
# ggp2 + facet_wrap(vars(Sex))
# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)
# 
# jpeg("TauBraak2PredictionPlots.jpg",width=2000,height = 500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()
# 
# m1=readRDS("TauBraak2NoCovarModelWOGenderEffects.Rds.Rds")
# m2=readRDS("TauBraak2NoCovarModelWithGenderEffects.Rds.Rds")
# SigOrNot=CompareSplineSigSlope("TauBraak2","NoCovar",tmpdd,m1,m2,"",'F')
# 
# m1_predict$z=0
# m1_predict$z= (m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
# m1_predict$Significance=0
# m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
# m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
# m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
# ggplot(m1_predict, aes(x = Age, y = z, colour = APOE, size=Significance)) +
#   geom_line()+ 
#   plt_labs 
# m1_predict$Measure="TauBraak2"
# m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])
# 
# m2_predict$z=0
# m2_predict$z= (m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
# m2_predict$Significance=0
# m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
# m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
# m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
# m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
# m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
# m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
# m2_predict$Measure="TauBraak2"
# ggplot(m2_predict, aes(x = Age, y = z, colour = APOE, size=Significance)) +
#   geom_line()+ 
#   plt_labs + facet_wrap(vars(Sex))
# m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])
# 
# 
# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)


#Tau Braak Regions 3 & 4----

tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM %in% dx,]

#regions affected in Braak 3&4 stages of AD
# left and right parahippocampal, fusiform, lingual cortices, amygdala, 
# left and right middle temporal, caudal and rostral anterior cingulate, 
# post cingulate, isthmus cingulate, insula, inferior temporal, and 
# temporal pole cortices

reg=c("CTX_PARAHIPPOCAMPAL_SUVR","CTX_FUSIFORM_SUVR","CTX_LINGUAL_SUVR","AMYGDALA_SUVR",
      "CTX_MIDDLETEMPORAL_SUVR", "CTX_CAUDALANTERIORCINGULATE_SUVR",
      "CTX_ROSTRALANTERIORCINGULATE_SUVR", "CTX_POSTERIORCINGULATE_SUVR",
      "CTX_ISTHMUSCINGULATE_SUVR", "CTX_INSULA_SUVR",
      "CTX_INFERIORTEMPORAL_SUVR", "CTX_TEMPORALPOLE_SUVR")
tmpdd$BRAAK34_SUVR=rowMeans(tmpdd[,paste(as.character(reg))]) #Average of all the ROIs above 
tmpdd$reg= tmpdd$BRAAK34_SUVR/tmpdd$INFERIOR_CEREBGM_SUVR
  
 
# tmpdd=tmpdd[tmpdd$INFERIOR_CEREBGM_SUVR>0.9,]

tmpdd=tmpdd[!is.na(tmpdd$reg),]
tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# tmpdd = tmpdd[tmpdd$Age >= 63 & tmpdd$Age <= 87,]
tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID)
tmpdd$Agediff=tmpdd$Age-min(tmpdd$Age)


tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'Tau Braak Region III & IV',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp= ggplot(tmpdd, aes(x = Age, y = reg,
#                        # group = RID, #Switching to Cross sectional analysis
#                        colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID,size=1.1)) +
#   geom_point(aes(size=3))+ 
#   plt_labs

ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))
# ggp+facet_wrap(vars(Sex))
# ggp+facet_wrap(vars(APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"TauBraak3&4"))

if("TauPET" %in% names(Participant_Records)){
  Participant_Records$TauPET[Participant_Records$RID%in%tmpdd$RID]=1
  if(length(unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID])>0)){#if there are RIDs that are not present in the Participant_Records
    Participant_Records=merge(Participant_Records,data.frame(RID=unique(tmpdd$RID[tmpdd$RID %notin% Participant_Records$RID]),TauPET=1),by=c("RID","TauPET"),all = TRUE)
  }
}else{
  Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],TauPET=1),by=c("RID"),all = TRUE)
}



jpeg("TauBraak3&4dataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()

plotAPOEnSexinAgeBins(tmpdd,"TauPET")
Save_demog_distribution_withintmpdd(tmpdd,"TauPET")


# PlotCDSOBBins(tmpdd,"TauBraak3&4") 

# PrintAmyloidPosN(tmpdd)
# PrintN(tmpdd)
# PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# PlotDataBins(tmpdd)

m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE,k=4) +
            s(RID,Agediff, bs="re")+ s(RID,bs="re") + #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex,data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
saveRDS(m1,"TauBraak3&4PTEDUCATCovarModelWOGenderEffects.Rds")


# summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16, Sex='F',RID=NULL))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60), axis.title.y = element_text(size = 60)) +
#   labs(#title=paste("Tau Braak 3 & 4 prediction within ",dxtext,'with APOE effects'), 
#        x ="Age (years)", y = 'Tau (AV1452 SUVR) Braak 3 & 4 regions')
# ggp1


m2 <- gam(formula = reg ~ s(Agediff, bs = "cs", 
                            by = interaction(APOE,Sex),k=4) + 
            s(RID, Agediff, bs = "re") + s(RID,bs="re") + 
            PTEDUCAT + APOE + Sex + APOE:Sex, 
          data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"TauBraak3&4PTEDUCATCovarModelWithGenderEffects.Rds")

# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=list(PTEDUCAT=16,RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("Braak 3&4 Tau prediction within ",dxtext, 'with Sex x APOE effects'), 
#     x ="Age (years)", y = '')# 'AV1451 PVC SUVR in Braak 3&4 regions')
# ggp2 + facet_wrap(vars(Sex))
# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)


m1=readRDS("TauBraak3&4PTEDUCATCovarModelWOGenderEffects.Rds")
m2=readRDS("TauBraak3&4PTEDUCATCovarModelWithGenderEffects.Rds")
SigOrNot=CompareSplineSigSlope("TauBraakIII&IV","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT",'F')

m1_predict$z=0
m1_predict$z= (m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="TauBraak3&4"

ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= (m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="TauBraak3&4"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))

m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])

jpeg("TauBraakIII&IVPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("Tau Braak III & IV regions")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)



#regions affected in Braak III stages of AD
# # reg=c("PVC.AV1451.LeftAmygdala","PVC.AV1451.RightAmygdala",
# #       "PVC.AV1451.ctxlhparahippocampal","PVC.AV1451.ctxrhparahippocampal",
# #       "PVC.AV1451.ctxlhfusiform","PVC.AV1451.ctxrhfusiform",
# #       "PVC.AV1451.ctxlhlingual","PVC.AV1451.ctxrhlingual")
# 
# reg=c("AMYGDALA_SUVR",
#       "CTX_LH_PARAHIPPOCAMPAL_SUVR","CTX_RH_PARAHIPPOCAMPAL_SUVR",
#       "CTX_LH_FUSIFORM_SUVR","CTX_RH_FUSIFORM_SUVR",
#       "CTX_LH_LINGUAL_SUVR","CTX_RH_LINGUAL_SUVR")
# d$reg = rowMeans(d[,paste(as.character(reg))]) #Average of all the ROIs above
# dd = d[d$APOE != 'NA' & !is.na(d$Age) #& !is.na(d$AmyloidPos) 
#        & d$DXSUM%in%dx & !is.na(d$reg),]
# # PlotDataBins(d[d$APOE != 'NA' & !is.na(d$Age) & !is.na(d$AmyloidPos) & !is.na(d$reg),])
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$EXAMDATE),]
# dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),] 
# 
# #Longitudinal Analysis
# dd$NumOfScans=c(1)
# 
# for (id in unique(dd$RID))
# {
#   dd$NumOfScans[dd$RID==id]=length(dd$RID[dd$RID==id])
# }
# 
# tmpdd=dd
# 
# tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# tmpdd$Age=as.numeric(tmpdd$Age) #This is needed for whatever reason. 
# tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
#                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID)
# write.csv(tmpdd,"TauBraakIIIdata.csv")

# tmpdd=read.csv("TauBraakIIIdata.csv")
# tmpdd$APOE <- factor(tmpdd$APOE)
# tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# tmpdd$Sex=factor(tmpdd$Sex)
# tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# plt_labs <- labs(y = 'Tau Braak Region III',
#                  x = 'Age in Years',
#                  colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   # group = RID, #Switching to Cross sectional analysis
#                   colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
#   geom_point()+ 
#   plt_labs
# 
# ggp+ facet_wrap(vars(Sex,APOE))
# # ggp+facet_wrap(vars(Sex)) 
# # ggp+facet_wrap(vars(APOE))
# jpeg("TauBraakIIIdataScatter.jpg",width=700,height = 500)
# ggp+ facet_wrap(vars(Sex,APOE))
# dev.off()
# 
# PlotCDSOBBins(tmpdd,"Tau")
# 
# PrintAmyloidPosN(tmpdd)
# # PrintN(tmpdd)
# # PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# # PlotDataBins(tmpdd)
# 
# m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) +
#             s(RID,Age, bs="re")+ #Accounting for random intercept and random slope by individual subjects' repeated measures
#             APOE + Sex,
#           data = tmpdd, method = "REML")
# m1$model$Age=m1$model$Agediff+min(tmpdd$Age)
# saveRDS(m1,"TauBraakIIINoCovarModelWOGenderEffects.Rds")
# 
# summary(m1)
# 
# #predict_gam
# m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
#                         values=list(RID=NULL, Sex='F'))
# m1_predict$APOE=factor(m1_predict$APOE)
# m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("Braak III Tau prediction within ",dxtext,'with APOE effects'),
#        x ="Age (years)", y = 'Braak III region AV1451 PVC SUVR')
# ggp1
# 
# m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = interaction(APOE,Sex)) + 
#             s(RID, Age, bs = "re") + APOE +Sex + APOE:Sex, 
#           data = tmpdd, method = "REML")
# m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
# saveRDS(m2,"TauBraakIIINoCovarModelWithGenderEffects.Rds")
# 
# 
# summary(m2)
# 
# #predict_gam
# m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
#                         values=list(RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
# m2_predict$APOE=factor(m2_predict$APOE)
# m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# m2_predict$Sex=factor(m2_predict$Sex)
# m2_predict$Sex = relevel(m2_predict$Sex,'M')
# m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("Braak III Tau prediction within ",
#                    dxtext,'with Sex x APOE effects'),
#        x ="Age (years)", y = 'AV1451 PVC SUVR in Braak III regions')
# ggp2 + facet_wrap(vars(Sex))
# PrintAmyloidPosN(tmpdd)
# 
# jpeg("TauBraakIIIPredictionPlots.jpg",width=2000,height = 500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()
# 
# m1=readRDS("TauBraakIIISimpleModelWOGenderEffects.Rds")
# m2=readRDS("TauBraakIIISimpleModelWithGenderEffects.Rds")
# CompareSplineSigSlope("TauBraakIII","Simple",tmpdd,m1,m2,"CDSOB")
# 
# m1=readRDS("TauBraakIIIComplexModelWOGenderEffects.Rds")
# m2=readRDS("TauBraakIIIComplexModelWithGenderEffects.Rds")
# CompareSplineSigSlope("Tau Braak III","Complex",tmpdd,m1,m2,"CDSOB,PTEDUCAT,WML,Centiloid")
# 
# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

#regions affected in Braak IV stage of AD
# # reg=c("PVC.AV1451.ctxlhmiddletemporal","PVC.AV1451.ctxlhcaudalanteriorcingulate","PVC.AV1451.ctxlhrostralanteriorcingulate",
# #       "PVC.AV1451.ctxlhposteriorcingulate","PVC.AV1451.ctxlhisthmuscingulate",
# #       "PVC.AV1451.ctxlhinferiortemporal","PVC.AV1451.ctxlhinsula","PVC.AV1451.ctxlhtemporalpole",
# #       "PVC.AV1451.ctxrhmiddletemporal","PVC.AV1451.ctxrhcaudalanteriorcingulate","PVC.AV1451.ctxrhrostralanteriorcingulate",
# #       "PVC.AV1451.ctxrhposteriorcingulate","PVC.AV1451.ctxrhisthmuscingulate",
# #       "PVC.AV1451.ctxrhinferiortemporal","PVC.AV1451.ctxrhinsula","PVC.AV1451.ctxrhtemporalpole")
# reg=c("CTX_LH_MIDDLETEMPORAL_SUVR","CTX_LH_ROSTRALANTERIORCINGULATE_SUVR",
#       "CTX_LH_CAUDALANTERIORCINGULATE_SUVR", "CTX_LH_POSTERIORCINGULATE_SUVR",
#       "CTX_LH_ISTHMUSCINGULATE_SUVR", "CTX_LH_INFERIORTEMPORAL_SUVR",
#       "CTX_LH_INSULA_SUVR", "CTX_LH_TEMPORALPOLE_SUVR",
#       "CTX_RH_MIDDLETEMPORAL_SUVR","CTX_RH_ROSTRALANTERIORCINGULATE_SUVR",
#       "CTX_RH_CAUDALANTERIORCINGULATE_SUVR", "CTX_RH_POSTERIORCINGULATE_SUVR",
#       "CTX_RH_ISTHMUSCINGULATE_SUVR", "CTX_RH_INFERIORTEMPORAL_SUVR",
#       "CTX_RH_INSULA_SUVR", "CTX_RH_TEMPORALPOLE_SUVR")
# d$reg = rowMeans(d[,paste(as.character(reg))]) #Average of all the ROIs above
# dd = d[d$APOE != 'NA' & !is.na(d$Age) #& !is.na(d$AmyloidPos) 
#        & d$DXSUM%in%dx & !is.na(d$reg),]
# # PlotDataBins(d[d$APOE != 'NA' & !is.na(d$Age) & !is.na(d$AmyloidPos) & !is.na(d$reg),])
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$EXAMDATE),]
# dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),] 
# 
# #Longitudinal Analysis
# dd$NumOfScans=c(1)
# 
# for (id in unique(dd$RID))
# {
#   dd$NumOfScans[dd$RID==id]=length(dd$RID[dd$RID==id])
# }
# 
# tmpdd=dd
# 
# tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# tmpdd$Age=as.numeric(tmpdd$Age) #This is needed for whatever reason. 
# tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
#                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# tmpdd$RID<-factor(tmpdd$RID)
# write.csv(tmpdd,"TauBraakIVdata.csv")

# tmpdd=read.csv("TauBraakIVdata.csv")
# tmpdd$APOE <- factor(tmpdd$APOE)
# tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# tmpdd$Sex=factor(tmpdd$Sex)
# tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# plt_labs <- labs(y = 'Tau Braak Region IV',
#                  x = 'Age in Years',
#                  colour = 'APOE')
# ggp= ggplot(tmpdd, aes(x = Age, y = reg,
#                   # group = RID, #Switching to Cross sectional analysis
#                   colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
#   geom_point()+ 
#   plt_labs
# 
# ggp+ facet_wrap(vars(Sex,APOE))
# # ggp+facet_wrap(vars(Sex)) 
# # ggp+facet_wrap(vars(APOE))
# jpeg("TauBraakIVdataScatter.jpg",width=700,height = 500)
# ggp+ facet_wrap(vars(Sex,APOE))
# dev.off()
# 
# PrintAmyloidPosN(tmpdd)
# # PrintN(tmpdd)
# # PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# # PlotDataBins(tmpdd)
# 
# # m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) +
# #             # s(RID, bs="re") + 
# #             s(RID,Age, bs="re")+ #Accounting for random intercept and random slope by individual subjects' repeated measures
# #             APOE + CDSOB:APOE + PTEDUCAT:APOE + Centiloid:APOE+ WML:APOE,
# #           data = tmpdd, method = "REML")
# # saveRDS(m1,"TauBaakIVComplexModelWOGenderEffects.Rds")
# 
# m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) +
#             # s(RID, bs="re") + 
#             s(RID,Age, bs="re")+ #Accounting for random intercept and random slope by individual subjects' repeated measures
#             APOE + CDSOB:APOE,
#           data = tmpdd, method = "REML")
# saveRDS(m1,"TauBaakIVSimpleModelWOGenderEffects.Rds")
# 
# summary(m1)
# 
# #predict_gam
# m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Age)'), values=list(PTEDUCAT=c(16), CDSOB=cdr.sob, AmyloidPos=0, WML=0, Centiloid=20,RID=NULL))
# m1_predict$APOE=factor(m1_predict$APOE)
# m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("Braak IV Tau prediction within ",dxtext,'with Centiloid without Sex effects'),
#        x ="Age (years)", y = 'Braak IV region AV1451 PVC SUVR')
# ggp1
# #Modelling interactions without separating the groups
# # m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) + 
# #             # s(RID, bs = "re") + 
# #             s(RID, Age, bs = "re") + 
# #             int_sex_apoe + CDSOB:int_sex_apoe + PTEDUCAT:int_sex_apoe + 
# #             Centiloid:int_sex_apoe + WML:int_sex_apoe, 
# #           data = tmpdd, method = "REML")
# # saveRDS(m2,"TauBaakIVComplexModelWithGenderEffects.Rds")
# 
# m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) + 
#             # s(RID, bs = "re") + 
#             s(RID, Age, bs = "re") + 
#             int_sex_apoe + CDSOB:int_sex_apoe, 
#           data = tmpdd, method = "REML")
# saveRDS(m2,"TauBaakIVSimpleModelWithGenderEffects.Rds")
# 
# summary(m2)
# 
# #predict_gam
# m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Age)'), values=list(PTEDUCAT=c(16), CDSOB=cdr.sob, WML=0, Centiloid=20, RID=NULL))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
# m2_predict$APOE=factor(m2_predict$APOE)
# m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# m2_predict$Sex=factor(m2_predict$Sex)
# m2_predict$Sex = relevel(m2_predict$Sex,'M')
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + labs(title=paste("Braak IV Tau prediction within ",dxtext), 
#                                                         x ="Age (years)", y = 'AV1451 PVC SUVR in Braak IV regions')
# ggp2 + facet_wrap(vars(Sex))
# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)
# 
# jpeg("TauBraakIVPredictionPlots.jpg",width=2000,height = 500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()
# 
# m1=readRDS("TauBaakIVSimpleModelWOGenderEffects")
# m2=readRDS("TauBaakIVSimpleModelWithGenderEffects")
# CompareSplineSigSlope("TauBraakIV","Simple",tmpdd,m1,m2,"CDSOB")
# 
# m1=readRDS("TauBaakIVComplexModelWOGenderEffects")
# m2=readRDS("TauBaakIVComplexModelWithGenderEffects")
# CompareSplineSigSlope("TauBraakIV","Complex",tmpdd,m1,m2,"CDSOB,PTEDUCAT,WML,Centiloid")
# 
# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)
# 
# # NA Tau in regions affected in Braak V stage of AD----
# reg=c("PVC.AV1451.ctxlhsuperiorfrontal","PVC.AV1451.ctxlhlateralorbitofrontal","PVC.AV1451.ctxlhmedialorbitofrontal",
#       "PVC.AV1451.ctxlhfrontalpole","PVC.AV1451.ctxlhcaudalmiddlefrontal","PVC.AV1451.ctxlhrostralmiddlefrontal",
#       "PVC.AV1451.ctxlhparsopercularis","PVC.AV1451.ctxlhparsorbitalis","PVC.AV1451.ctxlhparstriangularis","PVC.AV1451.ctxlhlateraloccipital",
#       "PVC.AV1451.ctxlhsupramarginal","PVC.AV1451.ctxlhinferiorparietal","PVC.AV1451.ctxlhsuperiortemporal","PVC.AV1451.ctxlhsuperiorparietal",
#       "PVC.AV1451.ctxlhprecuneus","PVC.AV1451.ctxlhbankssts","PVC.AV1451.ctxlhtransversetemporal",
#       "PVC.AV1451.ctxrhsuperiorfrontal","PVC.AV1451.ctxrhlateralorbitofrontal","PVC.AV1451.ctxrhmedialorbitofrontal",
#       "PVC.AV1451.ctxrhfrontalpole","PVC.AV1451.ctxrhcaudalmiddlefrontal","PVC.AV1451.ctxrhrostralmiddlefrontal",
#       "PVC.AV1451.ctxrhparsopercularis","PVC.AV1451.ctxrhparsorbitalis","PVC.AV1451.ctxrhparstriangularis","PVC.AV1451.ctxrhlateraloccipital",
#       "PVC.AV1451.ctxrhsupramarginal","PVC.AV1451.ctxrhinferiorparietal","PVC.AV1451.ctxrhsuperiortemporal","PVC.AV1451.ctxrhsuperiorparietal",
#       "PVC.AV1451.ctxrhprecuneus","PVC.AV1451.ctxrhbankssts","PVC.AV1451.ctxrhtransversetemporal")
# d$reg = rowMeans(d[,paste(as.character(reg))]) #Average of all the ROIs above
# dd = d[d$APOE != 'NA' & !is.na(d$Age) & !is.na(d$AmyloidPos) & d$DXSUM%in%dx & !is.na(d$reg),]
# PlotDataBins(d[d$APOE != 'NA' & !is.na(d$Age) & !is.na(d$AmyloidPos) & !is.na(d$reg),])
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# dd <- dd[order(dd$RID, dd$Age),]
# dd <- dd[!duplicated(dd$RID),] 
# PlotScatter(dd,'Age (Years)','Braak V region AV1451 SUVR (partial volume corrected)')
# # ceiling(min(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),max,na.rm=T)))
# # floor(max(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),min,na.rm=T)))
# # dd = dd[dd$Age >= floor(max(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),min,na.rm=T))) &
# #           dd$Age <= ceiling(min(tapply(dd[,paste(as.character("Age"))],interaction(dd$APOE,dd$Sex),max,na.rm=T))),]
# dd=dd[dd$reg<2.8,] #removing outliers with reg>2.8, both M, RID=4262 6632
# dd=dd[dd$Age<=85 &dd$Age>=65,]
# PrintN(dd)
# 


# # # CDR-SOB clinician-rated staging method for cognitive dysfunction and functional ability ----
# # #  0 = no dementia, .5 = questionable dementia, 1 = mild dementia, 2 = moderate dementia, and 3 = severe dementia
# # library(ADNIMERGE)
# # d=adnimerge
# # d$Age=d$AGE+d$Years.bl #Age @EXAMDATE = baseline AGE + Years since baseline
# # d=d[!is.na(d$CDRSB) & !is.na(d$PTGENDER) & !is.na(d$PTEDUCAT),]
# # d=subset(d,select=c(RID,Age,MMSE,DX,PTGENDER,PTEDUCAT,EXAMDATE,CDRSB)) #
# # colnames(d)[colnames(d)=='PTGENDER']='Sex'
# # d$Sex=factor(d$Sex)
# # levels(d$Sex)=c('F','M')
# # d$Sex = relevel(d$Sex,ref = 'M')
# # colnames(d)[colnames(d)=='DX']='DXSUM'
# # levels(d$DXSUM)[levels(d$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# # 
# # colnames(d)[colnames(d)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# # 
# # apoe=getAPOE()
# # d=merge(d,apoe,by=c('RID'),all.x=T) #Only  keeping RIDs that have the APOE data
# # rm(apoe)
# # d=d[d$APOE!='NA' & !is.na(d$APOE),]
# # 
# # wml=getWML()
# # wml=wml[wml$RID %in% unique(d$RID), ]
# # d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# # d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$wml.DATE,format="%Y-%m-%d"))
# # d <- d[d$dt<=180,]
# # d=d[!is.na(d$dt),]
# # d <- d[order(d$RID,d$EXAMDATE, d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the WML scan closest in date to the d scan
# # rm(wml)
# # d=subset(d, select = -c(wml.DATE,dt) )
# # 
# # # Adding AmyloidPos info
# # CentiloidData=getCentiloidData()
# # CentiloidData=CentiloidData[CentiloidData$RID %in% unique(d$RID), ]
# # d=merge(d,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# # d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d"))
# # d <- d[d$dt<=180,]
# # d=d[!is.na(d$dt),]
# # d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# # d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the d scan
# # rm(CentiloidData)
# # d=subset(d, select = -c(Cent.DATE,dt) )
# # 
# # reg='CDSOB'
# # d$reg = d[,c(paste(as.character(reg)))]
# # dd = d[d$APOE!='NA' & !is.na(d$Age) & d$DXSUM%in%dx & !is.na(d$reg),]
# #        # & !is.na(d$wmlICV)& !is.na(d$WMH),]
# # # PlotDataBins(d[!is.na(d$reg) & d$APOE !='NA' & !is.na(d$APOE) & !is.na(d$Age),])
# # dd$APOE <- factor(dd$APOE)
# # dd$DXSUM <- factor(dd$DXSUM)
# # dd <- dd[order(dd$RID, dd$Age),]
# # 
# # dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),]
# # 
# # dd$NumOfScans=c(1)
# # dd<-dd[order(dd$RID,dd$EXAMDATE),]
# # for (id in unique(dd$RID)) 
# # {
# #   dd$NumOfScans[dd$RID==id]=length(dd$RID[dd$RID==id])
# # }
# # 
# # tmpdd=dd
# # tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# # #Limiting data by mins and max of individual groups in the dataset
# # tmpdd=tmpdd[tmpdd$Age>=65 & tmpdd$Age<=85,]
# # tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
# #                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# # tmpdd$RID<-factor(tmpdd$RID)
# # write.csv(tmpdd,"CDRSOBdata.csv")
# 
# tmpdd=read.csv("CDRSOBdata.csv")
# tmpdd$APOE <- factor(tmpdd$APOE)
# tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# tmpdd$Sex=factor(tmpdd$Sex)
# tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# plt_labs <- labs(y = "CDR-SOB",
#                  x = 'Age in Years',
#                  colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                       group = RID, colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
#   geom_point()+ 
#   plt_labs
# 
# ggp + facet_wrap(vars(Sex,APOE))
# # ggp + facet_wrap(vars(Sex))
# # ggp + facet_wrap(vars(APOE))
# 
# jpeg("CDRSOBdataScatter.jpg",width=700,height = 500)
# ggp+ facet_wrap(vars(Sex,APOE))
# dev.off()
# 
# PrintAmyloidPosN(tmpdd)
# min(tmpdd$Age)
# max(tmpdd$Age)
# # PrintN(tmpdd)
# # PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# # PlotDataBins(tmpdd)
# 
# # m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) + 
# #             # s(RID, bs="re") + 
# #             s(RID,Age, bs="re") + #Accounting for random intercept and random slope by individual subjects' repeated measures
# #             APOE + Centiloid:APOE +WML:APOE + PTEDUCAT:APOE,
# #           data = tmpdd, method = "REML")
# # saveRDS(m1,"CDRSOBComplexModelWOGenderWithWML.Rds")
# 
# m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) + 
#             # s(RID, bs="re") + 
#             s(RID,Age, bs="re") + #Accounting for random intercept and random slope by individual subjects' repeated measures
#             APOE,# + PTEDUCAT:int_sex_apoe,
#           data = tmpdd, method = "REML")
# saveRDS(m1,"CDRSOBSimpleModelWOGenderWithWML.Rds")
# 
# summary(m1)
# 
# #predict_gam
# m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Age)'),
#                         values=list(RID=NULL, Centiloid=20, PTEDUCAT=16, WML=0))
# m1_predict$APOE=factor(m1_predict$APOE)
# m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
#   labs(title=paste("CDR-SOB prediction within ",dxtext), 
#        x ="Age (years)", y = "CDR-SOB score")
# ggp1
# 
# # #Complex Model
# # m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) +
# #             s(RID, Age, bs = "re") +
# #             int_sex_apoe+
# #             Centiloid:int_sex_apoe + WML:int_sex_apoe + PTEDUCAT:int_sex_apoe,
# #           data = tmpdd, method = "REML")
# # saveRDS(m2,"CDRSOBComplexModelWithGenderWithWML.Rds")
# 
# m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) + 
#             s(RID, Age, bs = "re") +
#             int_sex_apoe,# + PTEDUCAT:int_sex_apoe, 
#           data = tmpdd, method = "REML")
# saveRDS(m2,"CDRSOBSimpleModelWithGenderWithWML.Rds")
# 
# summary(m2)
# 
# #predict_gam
# m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Age)'), 
#                         values=list(RID=NULL, Centiloid=20, PTEDUCAT=16, WML=0))
# m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
# m2_predict$APOE=factor(m2_predict$APOE)
# m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# m2_predict$Sex=factor(m2_predict$Sex)
# m2_predict$Sex = relevel(m2_predict$Sex,'M')
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + labs(title=paste("CDR-SOB predictions within ",dxtext), 
#                                                         x ="Age (years)", y = "CDR-SOB")
# ggp2 + facet_wrap(vars(Sex))
# 
# jpeg("CDRSOBSimpleModelPredictionPlots.jpg",width=2000,height = 500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()
# m1=readRDS("CDRSOBSimpleModelWOGenderWithWML.Rds")
# m2=readRDS("CDRSOBSimpleModelWithGenderWithWML.Rds")
# CompareSplineSigSlope("CDRSOB","Complex",tmpdd,m1,m2,"")
# 
# jpeg("CDRSOBComplexModelPredictionPlots.jpg",width=2000,height = 500)
# ggp1+ (ggp2+facet_wrap(vars(Sex)))
# dev.off()
# m1=readRDS("CDRSOBComplexModelWOGenderWithWML.Rds")
# m2=readRDS("CDRSOBComplexModelWithGenderWithWML.Rds")
# CompareSplineSigSlope("CDRSOB","Complex",tmpdd,m1,m2,"Centiloid, WML,PTEDUCAT")
# 
# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

# PACC (Preclinical Alzheimer's Cognitive Composite) & MMSE (Mini Mental State Exam)----
## This needs to run if recreating datasets for MMSE and PACC
# rm(d)
# d=adnimerge
# d = pacc(d,keepComponents = FALSE)
# d$PACC=d$mPACCdigit
# d$PACC[is.na(d$PACC)]=d$mPACCtrailsB[is.na(d$PACC)]
# d = d[!is.na(d$PACC) & !is.na(d$PTGENDER) & !is.na(d$PTGENDER) & !is.na(d$CDRSB),]
# d$Age=d$AGE+d$Years.bl #Age @EXAMDATE = baseline AGE + Years since baseline
# d=subset(d,select=c(RID,Age,PACC,MMSE,DX,PTGENDER,PTEDUCAT,EXAMDATE,CDRSB)) #
# colnames(d)[colnames(d)=='PTGENDER']='Sex'
# d$Sex=factor(d$Sex)
# levels(d$Sex)=c('F','M')
# d$Sex = relevel(d$Sex,ref = 'M')
# colnames(d)[colnames(d)=='DX']='DXSUM'
# levels(d$DXSUM)[levels(d$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
# 
# colnames(d)[colnames(d)=='CDRSB']="CDSOB" #Confirmed below that the Calculated CDSOB is the same as CDRSB
# 
# apoe=getAPOE()
# d=merge(d,apoe,by=c('RID'),all.x=T) #Only  keeping RIDs that have the APOE data
# rm(apoe)
# d=d[d$APOE!='NA' & !is.na(d$APOE),]
# 
# # Adding AmyloidPos info #Not adding 8/2023
# CentiloidData=getCentiloidData()
# CentiloidData=CentiloidData[CentiloidData$RID %in% unique(d$RID), ]
# d=merge(d,CentiloidData,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$Cent.DATE,format="%Y-%m-%d")) 
# d <- d[d$dt<=180,]
# d=d[!is.na(d$dt),]
# d <- d[order(d$RID,d$EXAMDATE,d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the Centiloid scan closest in date to the d scan
# rm(CentiloidData)
# d=subset(d, select = -c(Cent.DATE,dt) )
# 
# wml=getWML()
# wml=wml[wml$RID %in% unique(d$RID), ]
# d=merge(d,wml,by=c('RID')) # No all.x=T because only keeping the common RIDs
# d$dt=abs(as.Date(d$EXAMDATE,format="%Y-%m-%d")-as.Date(d$wml.DATE,format="%Y-%m-%d")) 
# d <- d[d$dt<=180,]
# d=d[!is.na(d$dt),]
# d <- d[order(d$RID,d$EXAMDATE, d$dt),]
# d <- d[!duplicated(d[c("RID","EXAMDATE")]),] #Keep the WML scan closest in date to the d scan
# rm(wml)
# d=subset(d, select = -c(wml.DATE,dt) )

### mPACC ----
# d$reg=d$PACC
# dd = d[!is.na(d$reg) & !is.na(d$CDSOB) & d$APOE !='NA' & !is.na(d$APOE) &
#          !is.na(d$Age),]
# dd$APOE <- factor(dd$APOE)
# dd$DXSUM <- factor(dd$DXSUM)
# write.csv(dd,"PACCdataNoCovar1.csv")

WholeDatasetwithAD=read.csv("PACCdataNoCovar1.csv")
tmpdd=WholeDatasetwithAD[WholeDatasetwithAD$DXSUM%in%dx,]

tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
#Limiting data by mins and max of individual groups in the dataset
# tmpdd=tmpdd[tmpdd$Age<=88 &tmpdd$Age>=62,]#Based on E2 F and E4 F datasets
tmpdd = tmpdd[tmpdd$Age > floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
                tmpdd$Age < ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
#Running above twice because somehow running it only once does not eliminate all of the straggly outliers.

tmpdd$APOE <- factor(tmpdd$APOE)
tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
tmpdd$Sex=factor(tmpdd$Sex)
tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
tmpdd$Agediff=tmpdd$Age - min(tmpdd$Age)

tmpdd <- tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd <- tmpdd[!duplicated(tmpdd[c("RID","EXAMDATE")]),]

tmpdd$NumOfScans=c(1)
tmpdd<-tmpdd[order(tmpdd$RID,tmpdd$EXAMDATE),]
tmpdd$AvgYrsBtnVisits=c(0)
tmpdd$TotalFollowupYrs=c(0)
for (id in unique(tmpdd$RID)){
  tmpdd$NumOfScans[tmpdd$RID==id]=length(tmpdd$RID[tmpdd$RID==id])
  if (length(tmpdd$RID[tmpdd$RID==id])>1){
    tmpdd$AvgYrsBtnVisits[tmpdd$RID==id]=mean(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365)
    tmpdd$TotalFollowupYrs[tmpdd$RID==id]=round(sum(diff(as.numeric(as.Date(tmpdd$EXAMDATE[tmpdd$RID==id])))/365),2)
  }
}

plt_labs <- labs(y = 'mPACC',
                 x = 'Age in Years',
                 colour = 'APOE')
# ggp=ggplot(tmpdd, aes(x = Age, y = reg,
#                   group = RID, colour = APOE)) +
#   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID,size=1.1)) +
#   geom_point(aes(size=3))+ 
#   plt_labs
ggp=plotdatascatter(tmpdd,plt_labs)

ggp +facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))#+ facet_wrap(vars(Sex,APOE))

demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(tmpdd$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
tmpdd=merge(tmpdd,demog,by=c('RID'),all.x = TRUE)
rm(demog)

demog_info=rbind(demog_info,demog_distribution(tmpdd,"mPACC"))
Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],mPACCData=1),by=c("RID"),all = TRUE)

jpeg("mPACCdataScatter.jpg",width=3000,height = 2100)
ggp+ facet_grid(rows=vars(Sex),cols = vars(APOE),labeller = labeller(Sex=c(M="Males",F="Females")))+ #facet_wrap(vars(Sex,APOE))+ 
  geom_point(aes(size=1.1))+ geom_line(aes(size=1))+
  theme(legend.position="none",axis.text = element_text(size = 60), 
        axis.title = element_text(size = 60),strip.text = element_text(size =60)) 
dev.off()

plotAPOEnSexinAgeBins(tmpdd,"mPACC")
Save_demog_distribution_withintmpdd(tmpdd,"mPACC")


# PrintAmyloidPosN(tmpdd)
# max(tmpdd$Age)
# min(tmpdd$Age)
# PrintN(tmpdd)
# PlotScatter(tmpdd,'Age (years)','MetaROI Vol')
# PlotDataBins(tmpdd)

#Simple Model
m1 <- bam(formula = reg ~ s(Agediff, bs='cs', by= APOE, k=3) + 
            s(RID,Agediff, bs="re")+ s(RID,bs="re")  + #Accounting for random intercept and random slope by individual subjects' repeated measures
            PTEDUCAT + APOE + Sex,data = tmpdd, method = "REML")
m1$model$Age=m1$model$Agediff+min(tmpdd$Age)

saveRDS(m1,"mPACCPTEDUCATCovarModelWOGender.Rds")

# summary(m1)

#predict_gam
m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
                        values=list(PTEDUCAT=16, Sex='M',RID=NULL))
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(tmpdd$Age)


ggplot(m1_predict, aes(Age, fit, color = APOE)) +
  geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) +
  # scale_color_manual(values=cols) + 
  labs(x ="Age (years)", y = 'mPACC Score')



#Simple Model
m2 <- bam(formula = reg ~ s(Agediff, bs = "cs", by = interaction(APOE,Sex),k=3) + 
            s(RID, Agediff, bs = "re") + s(RID, bs = "re") +
            PTEDUCAT + APOE + Sex + APOE:Sex,
          data = tmpdd, method = "REML")
m2$model$Age=m2$model$Agediff+min(tmpdd$Age)
saveRDS(m2,"mPACCPTEDUCATCovarModelWithGender.Rds")

# summary(m2)

#predict_gam
m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100, 
                        values=list(PTEDUCAT=16, RID=NULL))
m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(tmpdd$Age)
# ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
#   geom_smooth_ci(size = 2, ci_z=1, ci_alpha=0.1) + 
#   theme(legend.position="none",axis.text = element_text(size = 60)) + 
#   scale_linetype_manual(values=c("dotdash", "longdash")) + 
#   labs(#title=paste("PACC prediction within ",dxtext),
#     x ="Age (years)", y = '')# 'mPACC')
# ggp2 + facet_wrap(vars(Sex))


m1=readRDS("mPACCPTEDUCATCovarModelWOGender.Rds")
m2=readRDS("mPACCPTEDUCATCovarModelWithGender.Rds")
SigOrNot=CompareSplineSigSlope("mPACC","PTEDUCATCovar",tmpdd,m1,m2,"PTEDUCAT","F")

m1_predict$z=0
m1_predict$z= -(m1_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)
m1_predict$Measure="mPACC"


ggplot(m1_predict, aes(x = Age, y = fit, colour=SigAPOE)) + #colour = APOE, size=Significance
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs +scale_x_continuous(limits = c(x_min, x_max))

m1zscores=rbind(m1zscores,m1_predict[,c("Age","APOE","z","Significance","Measure")])

m2_predict$z=0
m2_predict$z= -(m2_predict$fit-mean(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN']) )/stats::sd(tmpdd$reg[tmpdd$APOE=='E3/E3' & tmpdd$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)
m2_predict$Measure="mPACC"

ggplot(m2_predict, aes(x = Age, y = fit, colour = SigAPOE)) + #, size=Significance , linetype=Sex
  geom_line()+ scale_color_manual(values=cols) + geom_smooth_ci(size = 1.2, ci_alpha=0.1) +
  plt_labs + facet_wrap(vars(Sex)) +scale_x_continuous(limits = c(x_min, x_max))


m2zscores=rbind(m2zscores,m2_predict[,c("Age","APOE","Sex","z","Significance","Measure")])

jpeg("mPACCPTEDUCATCovarModelPredictionPlots.jpg",width=6000,height = 1500)
gg=plot_m1andm2predicts(m1_predict,m2_predict) 
gg
dev.off()

sink("ModelSummaryoutput.txt",append = TRUE)
print("mPACC scores")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()

# AIC(m1)
# BIC(m1)
# AIC(m2)
# BIC(m2)
# summary(m1)
# summary(m2)

# # ### MMSE----
# # # d$reg=d$MMSE
# # # dd = d[d$DXSUM%in%dx & !is.na(d$reg) & !is.na(d$Centiloid) & !is.na(d$WMH) &  
# # #          !is.na(d$wmlICV) & !is.na(d$CDSOB) & d$APOE !='NA' & !is.na(d$APOE) & 
# # #          !is.na(d$Age),]
# # # dd$APOE <- factor(dd$APOE)
# # # dd$DXSUM <- factor(dd$DXSUM)
# # # dd <- dd[order(dd$RID,dd$EXAMDATE),]
# # # dd <- dd[!duplicated(dd[c("RID","EXAMDATE")]),] 
# # # 
# # # dd$NumOfScans=c(1)
# # # dd<-dd[order(dd$RID,dd$EXAMDATE),]
# # # for (id in unique(dd$RID)) 
# # # {
# # #   dd$NumOfScans[dd$RID==id]=length(dd$RID[dd$RID==id])
# # # }
# # # 
# # # tmpdd=dd 
# # # tmpdd$int_sex_apoe <- interaction(tmpdd$Sex, tmpdd$APOE) #interaction variable that allows teasing out details of individual groups
# # # #Limiting data by mins and max of individual groups in the dataset
# # # tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
# # #                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# # # tmpdd=tmpdd[tmpdd$Age<=85 & tmpdd$Age>=65,]# This is an additional step based on data because E2 females only had 1 data point after age 85
# # # tmpdd = tmpdd[tmpdd$Age >= floor(max(tapply(tmpdd$Age,tmpdd$int_sex_apoe,min))) &
# # #                 tmpdd$Age <= ceiling(min(tapply(tmpdd$Age,tmpdd$int_sex_apoe,max))),]
# # # #Running above twice because somehow running it only once does not eliminate all of the straggly outliers.
# # # tmpdd=tmpdd[10<tmpdd$reg,]#eliminating 1 record from E3 group (RID 6016) with MMSE values less than -20
# # # 
# # # tmpdd$RID<-factor(tmpdd$RID)
# # # write.csv(tmpdd,"MMSEdata.csv")
# # 
# # 
# # tmpdd=read.csv("MMSEdata.csv")
# # tmpdd$APOE <- factor(tmpdd$APOE)
# # tmpdd$APOE = relevel(tmpdd$APOE,ref = 'E3/E3')
# # tmpdd$Sex=factor(tmpdd$Sex)
# # tmpdd$Sex = relevel(tmpdd$Sex,ref = 'M')
# # tmpdd$int_sex_apoe <- factor(tmpdd$int_sex_apoe)
# # tmpdd$int_sex_apoe = relevel(tmpdd$int_sex_apoe,ref="M.E3/E3")
# # plt_labs <- labs(y = 'MMSE',
# #                  x = 'Age in Years',
# #                  colour = 'APOE')
# # ggp=ggplot(tmpdd, aes(x = Age, y = reg,
# #                   group = RID, colour = APOE)) +
# #   geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
# #   geom_point()+ 
# #   plt_labs
# # 
# # 
# # ggp + facet_wrap(vars(Sex,APOE))
# # #ggp + facet_wrap(vars(Sex)) 
# # # ggp + facet_wrap(vars(APOE))
# # 
# # jpeg("MMSEdataScatter.jpg",width=700,height = 500)
# # ggp+facet_wrap(vars(Sex,APOE))
# # dev.off()
# # 
# # PrintAmyloidPosN(tmpdd)
# # max(tmpdd$Age)
# # min(tmpdd$Age)
# # 
# # 
# # # m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) + 
# # #             # s(RID, bs="re") + 
# # #             s(RID,Age, bs="re") +  #Accounting for random intercept and random slope by individual subjects' repeated measures
# # #             APOE + CDSOB:APOE + Centiloid:APOE + PTEDUCAT:APOE +WML:APOE,
# # #           data = tmpdd, method = "REML")
# # # saveRDS(m1,"MMSEComplexModelWOGender.Rds")
# # 
# # m1 <- bam(formula = reg ~ s(Age, bs='cs', by= APOE) + 
# #             # s(RID, bs="re") + 
# #             s(RID,Age, bs="re") +  #Accounting for random intercept and random slope by individual subjects' repeated measures
# #             APOE + CDSOB:APOE,
# #           data = tmpdd, method = "REML")
# # saveRDS(m1,"MMSESimpleModelWOGender.Rds")
# # summary(m1)
# # 
# # #predict_gam
# # m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Age)'), 
# #                         values=list(RID=NULL, CDSOB=0, Centiloid=20, PTEDUCAT=16, WML=0))
# # m1_predict$APOE=factor(m1_predict$APOE)
# # m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
# # ggp1=ggplot(m1_predict, aes(Age, fit, color = APOE)) +
# #   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + 
# #   labs(title=paste("MMSE prediction within ",dxtext,'with APOE effects'), 
# #        x ="Age (years)", y = 'MMSE')
# # ggp1
# # 
# # #Modelling interactions without separating the groups
# # # m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) + 
# # #             s(RID, Age, bs = "re") + 
# # #             int_sex_apoe + CDSOB:int_sex_apoe+ Centiloid:int_sex_apoe + 
# # #             WML:int_sex_apoe + PTEDUCAT:int_sex_apoe,
# # #           data = tmpdd, method = "REML")
# # # saveRDS(m2,"MMSEComplexModelWithGender.Rds")
# # 
# # m2 <- bam(formula = reg ~ s(Age, bs = "cs", by = int_sex_apoe) + 
# #             s(RID, Age, bs = "re") + 
# #             int_sex_apoe + CDSOB:int_sex_apoe,
# #           data = tmpdd, method = "REML")
# # saveRDS(m2,"MMSESimpleModelWithGender.Rds")
# # 
# # summary(m2)
# # 
# # #predict_gam
# # m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Age)'), 
# #                         values=list(RID=NULL, CDSOB=0, Centiloid=20, PTEDUCAT=16, WML=0))
# # m2_predict<-separate(m2_predict,int_sex_apoe,c('Sex','APOE'),extra='merge')
# # m2_predict$APOE=factor(m2_predict$APOE)
# # m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
# # m2_predict$Sex=factor(m2_predict$Sex)
# # m2_predict$Sex = relevel(m2_predict$Sex,'M')
# # ggp2=ggplot(m2_predict, aes(Age, fit, color = APOE, linetype = Sex)) +
# #   geom_smooth_ci(size = 1, ci_z=1, ci_alpha=0.1) + labs(title=paste("MMSE prediction within ",dxtext), 
# #                                                         x ="Age (years)", y = 'MMSE')
# # ggp2 + facet_wrap(vars(Sex))
# # 
# # jpeg("MMSESimpleModelPredictionPlots.jpg",width=2000,height = 500)
# # ggp1+ (ggp2+facet_wrap(vars(Sex)))
# # dev.off()
# # 
# # m1=readRDS("MMSESimpleModelWOGender.Rds")
# # m2=readRDS("MMSESimpleModelWithGender.Rds")
# # CompareSplineSigSlope("MMSE","Simple",tmpdd,m1,m2,"CDSOB")
# # 
# # jpeg("MMSEComplexModelPredictionPlots.jpg",width=2000,height = 500)
# # ggp1+ (ggp2+facet_wrap(vars(Sex)))
# # dev.off()
# # m1=readRDS("MMSEComplexModelWOGender.Rds")
# # m2=readRDS("MMSEComplexModelWithGender.Rds")
# # CompareSplineSigSlope("MMSE","Complex",tmpdd,m1,m2,"CDSOB,Centiloid,PTEDUCAT,WML")
# # 
# # AIC(m1)
# # BIC(m1)
# # AIC(m2)
# # BIC(m2)
# # summary(m1)
# # summary(m2)
# 
# Saving all the demographic and z score variables----
write.csv(m1zscores,"m1zscores_9-2023.csv")
write.csv(m2zscores,"m2zscores_9-2023.csv")
write.csv(demog_info,"demographic_distribution.csv")
Participant_Records=Participant_Records[Participant_Records$RID!=0,] #removing the initial 1 moc record
write.csv(Participant_Records,"Participant_Records.csv")

#Plotting Participation Venn Diagram----
#Venn function works for higher order venn diagrams
jj=read.csv("Participant_Records.csv") #or jj=Participant_Records
jj[is.na(jj)]=0
jpeg("Participant_RecordsWithVenn.jpg",width=6000,height = 6000)
venn(jj[,c("TauPET","VolumetricData","WMH","AmyloidPET","FDGPET","mPACCData")], zcolor = "style",sncs = 15,ilcs = 12,ggplot = TRUE)
dev.off()
# #Venn.diagram function only works for 5 or less groups.
# venn.diagram(list(WMH=Participant_Records$RID[!is.na(Participant_Records$WMH)],
#                   AmyloidPET=Participant_Records$RID[!is.na(Participant_Records$AmyloidPET)],
#                   TauPET=Participant_Records$RID[!is.na(Participant_Records$TauPET)],
#                   VolData=Participant_Records$RID[!is.na(Participant_Records$VolumetricData)],
#                   # FDGPET=Participant_Records$RID[!is.na(Participant_Records$FDGPET)],
#                   mPACC=Participant_Records$RID[!is.na(Participant_Records$mPACC)]),
#              filename = "ParticipantVennDiagram.png")


#Plotting APOE and Sex distributions in participant records----
jj=read.csv("Participant_Records.csv") #or jj=Participant_Records
apoe=getAPOE()
jj=merge(jj,apoe,by=c("RID"))
Adni=adnimerge
colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
jj=merge(jj,Adni[,c("RID","Sex")],by=c("RID"),all.x = TRUE)
jj=jj[order(jj$RID,jj$APOE,jj$Sex),]
jj=jj[!duplicated(jj[,c("RID","APOE","Sex")]),]
jj=jj[,c("RID","APOE","Sex")]
jj_dist=data.frame(expand.grid(APOE=c('E3/E3','E3/E2','E3/E4'),Sex=c("Male","Female"),Count=0,Percent_within_APOE=0,APOE_total=0))
for(apoe in unique(jj_dist$APOE)){
  for(sx in unique(jj_dist$Sex)){
    print(paste(apoe,sx))
    length(unique(jj$RID[jj$APOE==apoe & jj$Sex==sx]))
    jj_dist$Count[jj_dist$APOE==apoe & jj_dist$Sex==sx]=length(unique(jj$RID[jj$APOE==apoe & jj$Sex==sx]))
  }
  jj_dist$APOE_total[jj_dist$APOE==apoe]=sum(jj_dist$Count[jj_dist$APOE==apoe])
  jj_dist$Percent_within_APOE[jj_dist$APOE==apoe]=jj_dist$Count[jj_dist$APOE==apoe]*100/sum(jj_dist$Count[jj_dist$APOE==apoe])

}
for(sx in unique(jj_dist$Sex)){

  jj_dist$Sex_total[jj_dist$Sex==sx]=sum(jj_dist$Count[jj_dist$Sex==sx])
  jj_dist$Percent_within_Sex[jj_dist$Sex==sx]=jj_dist$Count[jj_dist$Sex==sx]*100/sum(jj_dist$Count[jj_dist$Sex==sx])
}

gg1=ggplot(data=jj_dist)+
  geom_bar_pattern(aes(x=APOE,y=Count,fill=APOE,pattern=Sex),stat = "identity",width = 0.7)+
  labs(y = '# Participants (N)', x = 'APOE')+
  geom_label(aes(label=Count,x=APOE,y=Count,group=Sex),position = position_stack(vjust = 0.5),color="#000000",size=30,label.padding = unit(1, "lines"),show.legend = NA)+
  geom_text(data=jj_dist[!duplicated(jj_dist$APOE),],aes(label=APOE_total,x=APOE,y=APOE_total,vjust=-0.2),color="#000000",size=40)+
  geom_label(aes(label="A",x=0.5,y=0.9*max(APOE_total),group=APOE),size=40)+
  scale_fill_manual(values = c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red"))+
  scale_pattern_manual(values=c('none','wave'))+
  theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
        panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
        panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
        axis.text = element_text(size=70),axis.title = element_text(size=80),
        plot.margin = unit(c(1,0.25,1,1),"inches")) #axis.text = element_blank(),

gg2=ggplot(data=jj_dist)+
  geom_bar_pattern(aes(x=APOE,y=Percent_within_APOE,fill=APOE,pattern=Sex),stat = "identity",width = 0.7)+
  labs(y = 'Percentage ', x = 'APOE')+
  geom_label(aes(label=paste0(round(Percent_within_APOE,digits = 1),"%"),x=APOE,y=Percent_within_APOE,group=Sex),
             position = position_stack(vjust = 0.5),color="#000000",size=30,label.padding = unit(1, "lines"),show.legend = NA)+
  geom_label(aes(label="B",x=0.5,y=90,group=APOE),size=40)+
  scale_fill_manual(values = c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red"))+
  scale_pattern_manual(values=c('none','wave'))+scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
        panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
        panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
        axis.text = element_text(size=70), axis.title = element_text(size=80))

gg3=ggplot(data=jj_dist)+
  geom_bar_pattern(aes(x=Sex,y=Count,fill=APOE,pattern=Sex),stat = "identity",width = 0.5)+
  labs(y = '# Participants (N)', x = 'Sex')+
  geom_label(aes(label=Count,x=Sex,y=Count,group=APOE),position = position_stack(vjust = 0.5),color="#000000",size=30,label.padding = unit(1, "lines"),show.legend = NA)+
  geom_text(data=jj_dist[!duplicated(jj_dist$Sex),],aes(label=Sex_total,x=Sex,y=Sex_total,vjust=-0.2),color="#000000",size=40)+
  geom_label(aes(label="C",x=0.5,y=0.9*max(Sex_total),group=APOE),size=40)+
  scale_fill_manual(values = c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red"))+
  scale_pattern_manual(values=c('none','wave'))+
  theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
        panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
        panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
        axis.text = element_text(size=70),axis.title = element_text(size=80),
        plot.margin = unit(c(1,0.25,1,1),"inches")) #axis.text = element_blank(),

gg4=ggplot(data=jj_dist)+
  geom_bar_pattern(aes(x=Sex,y=Percent_within_Sex,fill=APOE,pattern=Sex),stat = "identity",width = 0.5)+
  labs(y = 'Percentage ', x = 'Sex')+
  geom_label(aes(label=paste0(round(Percent_within_Sex,digits = 1),"%"),x=Sex,y=Percent_within_Sex,group=APOE),
             position = position_stack(vjust = 0.5),color="#000000",size=30,label.padding = unit(1, "lines"),show.legend = NA)+
  geom_label(aes(label="D",x=0.5,y=90,group=APOE),size=40)+
  scale_fill_manual(values = c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red"))+
  scale_pattern_manual(values=c('none','wave'))+scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme(legend.position="none", #legend.text = element_text(size = 60), legend.title = element_text(size=60),
        panel.background = element_rect(fill = "white", colour = "#6D9EC1",size = 0, linetype = "solid"),
        panel.grid.major = element_line(size = 0, linetype = 'solid',colour ="black"),
        axis.text = element_text(size=70), axis.title = element_text(size=80))


jpeg("Participant_RecordsSexAPOEdistribution.jpg",width=6000,height = 6000)
gg1+gg2+gg3+gg4
dev.off()

ggtmp=ggplot(data=jj_dist)+
  geom_bar_pattern(aes(x=APOE,y=Percent_within_APOE,fill=APOE,pattern=Sex),
                   stat = "identity",width = 0.8,
                   pattern_density = 0.5, pattern_spacing = 0.5,
                   pattern_key_scale_factor = 0.1)+
  labs(y = 'Percent ', x = 'APOE')+
  scale_fill_manual(values = c("E3/E2"="darkcyan","E3/E3"="black","E3/E4"="red"))+
  scale_pattern_manual(values=c('none','wave'))+scale_y_continuous(labels = scales::percent_format(scale = 1))+
  theme(legend.position="bottom")


t=get_legend(ggtmp+ theme(legend.text = element_text(size = 60),legend.title = element_text(size = 60),
                           legend.key.size = unit(6,'cm'),legend.direction = "horizontal"))
grid.newpage()
grid.draw(t)
jpeg("Participant_distribution_Legends.jpg",width=3000)
grid.draw(t)
dev.off()

rm(ggptmp,t,tmpm2predict)




# Full Adni dataset with Racial and APOE info added to it----
# This is only to visualize the full dataset. 
# Not using this to add Adnimerge data to the rest of the measures because the 
# Adnimerge date does not exactly match up the DXSUM or CDSOB dates from their independent dataset files.
Adni=adnimerge
colnames(Adni)[colnames(Adni)=='PTGENDER']='Sex'
Adni$Sex=factor(Adni$Sex)
levels(Adni$Sex)=c('F','M')
Adni$Sex = relevel(Adni$Sex,ref = 'M')
colnames(Adni)[colnames(Adni)=='DX']='DXSUM'
levels(Adni$DXSUM)[levels(Adni$DXSUM)=='Dementia']='AD'#Changing Dementia to AD
colnames(Adni)[colnames(Adni)=='AGE']='Age.bl'
Adni$Age=Adni$Age.bl+Adni$Years.bl
colnames(Adni)[colnames(Adni)=='CDRSB']="CDSOB" #Confirmed that the Calculated CDSOB is the same as CDRSB
Adni=Adni[,c("RID","EXAMDATE","Age","Sex","DXSUM","PTEDUCAT","CDSOB")]
Adni=Adni[order(Adni$RID,Adni$EXAMDATE),]
Adni=Adni[!is.na(Adni$Age)&!is.na(Adni$Sex)&!is.na(Adni$DXSUM),]
demog=getDemog() #Adding Ethnicity and Racial info to the dataset
demog=demog[demog$RID%in%unique(Adni$RID),c("RID","PTRACCAT","PTETHCAT")]
demog=demog[(!is.na(demog$PTRACCAT)|!is.na(demog$PTETHCAT))&!duplicated(demog$RID),]
Adni=merge(Adni,demog,by=c('RID'),all.x = TRUE) #PTRACCAT and PTETHCAT will remain NA for those that the information is not available on
rm(demog)

apoe=getAPOE()
apoe=apoe[apoe$RID%in%unique(Adni$RID),]
Adni=merge(Adni,apoe,by=c('RID')) # removed ,all.x=T as we need APOE for each record

Adni=Adni[!is.na(Adni$APOE),]
Full_Participant_Records=read.csv("Participant_Records.csv") #This file actually gets created at the very end of going through all individual datasets and tells us RIDs that participated in at least one of the studies

Adni<-Adni[order(Adni$RID,Adni$EXAMDATE),]
Adni$NumOfScans=1
Adni=Adni[!duplicated(Adni$RID),] #Removing duplicated records, as the repeated measures are for a given measure, and this is only for rest of the demographic info

Adni=Adni[Adni$RID %in% Full_Participant_Records$RID,] 
Adni$AvgYrsBtnVisits=c(0)
Adni$TotalFollowupYrs=c(0)
Save_demog_distribution_withintmpdd(Adni,"Combined Dataset")
demog_info=rbind(demog_info,demog_distribution(Adni,"Combined Dataset"))
write.csv(demog_info,"demographic_distribution.csv")

rm(Adni,apoe,Full_Participant_Records)


# Creating a combined image of all Prediction Plots ----
AmyloidPlot<-image_read('CentiloidPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
Tau1Plot<-image_read('TauBraakIPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
Tau34Plot<-image_read('TauBraakIII&IVPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
HippoVolPlot<-image_read('HippocampalVolPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
MetaROIVolPlot<-image_read('MetaROIVolPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
FDGPlot<-image_read('MetaROIFDGSUVRAdjWithPonsPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
WMHPlot<-image_read('WMHPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
mPACCPlot<-image_read('mPACCPTEDUCATCovarModelPredictionPlots.jpg')%>% image_ggplot()
Legends <- image_read('PredictionPlotLegends.jpg')%>% image_ggplot()


jpeg("AllPredictionPlotsWggarrange1.jpg",width=12000,height = 6490)

#AllPredictionPlotsWggarrange
# ggpubr::ggarrange(as.grob(AmyloidPlot),as.grob(MetaROIVolPlot),
#                   as.grob(Tau1Plot),as.grob(FDGPlot),
#                   as.grob(Tau34Plot),as.grob(WMHPlot),
#                   as.grob(HippoVolPlot),as.grob(mPACCPlot),
#                   as.grob(Legends),
#                   labels=c("A","B","C","D","E","F","G","H","Legends:",""),
#                  ncol = 2,nrow=5)

# AllPredictionPlotsWggarrange1
ggpubr::ggarrange(AmyloidPlot, Tau1Plot,Tau34Plot,Legends,
                  HippoVolPlot,MetaROIVolPlot, FDGPlot,WMHPlot,
                  labels=c("A","B","C","","D","E","F","G"),
                  nrow=2,ncol=4)

# AllPredictionPlotsWgridExtra
# gridExtra::grid.arrange(gridExtra::arrangeGrob(as.grob(AmyloidPlot),as.grob(MetaROIVolPlot),ncol=2),
#                         gridExtra::arrangeGrob(as.grob(Tau1Plot),as.grob(FDGPlot),ncol=2),
#                         gridExtra::arrangeGrob(as.grob(AmyloidPlot),as.grob(MetaROIVolPlot),ncol=2),
#                         gridExtra::arrangeGrob(as.grob(Tau34Plot),as.grob(WMHPlot),ncol=2),
#                         gridExtra::arrangeGrob(as.grob(HippoVolPlot),as.grob(mPACCPlot),ncol=2),
#                         gridExtra::arrangeGrob(as.grob(Legends),ncol=1))
dev.off()


# plotting z- scores ----
m1zscores=read.csv('m1zscores_8-2023.csv')
m1zscores$APOE=factor(m1zscores$APOE)
m1zscores$APOE=relevel(m1zscores$APOE,ref='E3/E3')
m1zscores$Significant=factor(m1zscores$Significance)

m2zscores=read.csv('m2zscores_8-2023.csv')
m2zscores$APOE=factor(m2zscores$APOE)
m2zscores$APOE=relevel(m2zscores$APOE,ref='E3/E3')
m2zscores$Significant=factor(m2zscores$Significance)

# m1zscores=m1zscores[m1zscores$Measure!="TauBraak2",]
# m2zscores=m2zscores[m2zscores$Measure!="TauBraak2",]


zscorecols <- c("AnteriorCingulateVol"="#999999","Centiloid"="#FF00FF","FDG"="#377EB8",
          "HippoVol"="#4DAF4A","MetaROIVol"="#66FF33", "PACC"="#0000CC",
          "TauBraak1"="#E41A1C", "TauBraak2"="#A65628", "TauBraak3&4"="#FF7F00", 
          "WMH"="#F781BF") 

zscorecols1 <- c("Centiloid"="#FF00FF","FDG"="#377EB8",
           "MetaROIVol"="#66FF33", "PACC"="#0000CC",
           "TauBraak1"="#E41A1C", "TauBraak3&4"="#FF7F00")



jpeg("APOE_z-scoreplots.jpg",width=6000,height = 2400)
ggplot(m1zscores, aes(x = Age, y = z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=4)+
  scale_color_manual(breaks= unique(m1zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  labs(y = 'z-scores',x = 'Age (Years)',title = 'Z-scores') +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()

# jpeg("APOE_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# ggplot(m1zscores[m1zscores$z>=0,], aes(x = Age, y = z)) +
#   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
#   scale_color_manual(breaks= unique(m1zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
#   labs(y = 'z-scores',x = 'Age in Years',title = 'Z-scores') +
#   facet_wrap(vars(APOE))
# dev.off()
# 
# 
# # jpeg("APOE_Selected_z-scoreplots.jpg",width=2100,height = 800)
# # ggplot(m1zscores[m1zscores$Measure=='Centiloid' | m1zscores$Measure=='TauBraak1' 
# #                  | m1zscores$Measure=='TauBraak3&4'| m1zscores$Measure=='FDG' | 
# #                    m1zscores$Measure=='MetaROIVol' | m1zscores$Measure=='PACC',],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# # 
# # jpeg("APOE_Selected_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# # ggplot(m1zscores[m1zscores$z>=0 & (m1zscores$Measure=='Centiloid' | m1zscores$Measure=='TauBraak1' 
# #                  | m1zscores$Measure=='TauBraak3&4'| m1zscores$Measure=='FDG' | 
# #                    m1zscores$Measure=='MetaROIVol' | m1zscores$Measure=='PACC'),],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure, group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+ guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# 
# minmeasures<-expand.grid(unique(m1zscores$Measure),unique(m1zscores$APOE))
# names(minmeasures)=c("Measure","APOE")
# m1zscores$adj_z=0
# for (eachAPOE in unique(m1zscores$APOE)) {
#   for (eachMeasure in unique(m1zscores$Measure)) {
#     minmeasures$minVal[minmeasures$Measure==eachMeasure & minmeasures$APOE==eachAPOE]= -min(m1zscores$z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE])
#     minmeasures$Age[minmeasures$Measure==eachMeasure & minmeasures$APOE==eachAPOE]= m1zscores$Age[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE & m1zscores$z==min(m1zscores$z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE])]
#     m1zscores$adj_z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE]=m1zscores$z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE]-min(m1zscores$z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE])
#     # print(paste(eachAPOE,eachMeasure,as.character(min(m1zscores$z[m1zscores$Measure==eachMeasure & m1zscores$APOE==eachAPOE]))))
#   }
# }
# 
# jpeg("APOE_ZeroShifted_z-scores.jpg",width=2100,height = 800)
# ggplot(m1zscores,aes(x = Age, y = adj_z)) +
#   geom_line(aes(color=Measure, group=Measure, alpha=Significance), size=1.1)+ 
#   # geom_point(data=minmeasures,aes(x=Age, y=minVal, color=Measure))+
#   scale_color_manual(breaks= unique(m1zscores$Measure),values = zscorecols)+ guides(alpha="none", size="none") +
#   labs(y = 'Zero shifted z-scores',x = 'Age in Years',title = 'Zero shifted Z-scores') + 
#   
#   facet_wrap(vars(APOE))
# dev.off()
# 
# # APOE X Sex analysis
# 
# # Females
jpeg("Female_z-scoreplots.jpg",width=6000,height = 2400)
ggplot(m2zscores[m2zscores$Sex=='F',], aes(x = Age, y = z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=4)+
  scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+
  guides(alpha="none", size="none") +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  labs(y = 'z-scores',x = 'Age (Years)',title = 'Female Z-scores') +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()
# 
# jpeg("Female_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# ggplot(m2zscores[m2zscores$Sex=='F' & m2zscores$z>=0,], aes(x = Age, y = z)) +
#   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
#   scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
#   labs(y = 'z-scores',x = 'Age in Years',title = 'Female Z-scores') +
#   facet_wrap(vars(APOE))
# dev.off()
# 
# 
# # jpeg("Female_Selected_z-scoreplots.jpg",width=2100,height = 800)
# # ggplot(m2zscores[m2zscores$Sex=='F' & (m2zscores$Measure=='Centiloid' | m2zscores$Measure=='TauBraak1' 
# #                  | m2zscores$Measure=='TauBraak3&4'| m2zscores$Measure=='FDG' | 
# #                    m2zscores$Measure=='MetaROIVol' | m2zscores$Measure=='PACC'),],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Female Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# # 
# # jpeg("Female_Selected_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# # ggplot(m2zscores[m2zscores$Sex=='F' & m2zscores$z>=0 & (m2zscores$Measure=='Centiloid' | m2zscores$Measure=='TauBraak1' 
# #                                    | m2zscores$Measure=='TauBraak3&4'| m2zscores$Measure=='FDG' | 
# #                                      m2zscores$Measure=='MetaROIVol' | m2zscores$Measure=='PACC'),],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure, group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Female Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# 
# #Males
jpeg("Male_z-scoreplots.jpg",width=6000,height = 2400)
ggplot(m2zscores[m2zscores$Sex=='M',], aes(x = Age, y = z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=4)+
  scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+
  guides(alpha="none", size="none") +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  labs(y = 'z-scores',x = 'Age (Years)',title = 'Male Z-scores') +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()
# 
# jpeg("Male_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# ggplot(m2zscores[m2zscores$Sex=='M' & m2zscores$z>=0,], aes(x = Age, y = z)) +
#   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
#   scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+
#   guides(alpha="none", size="none") +
#   labs(y = 'z-scores',x = 'Age in Years',title = 'Male Z-scores') +
#   facet_wrap(vars(APOE))
# dev.off()
# 
# 
# # jpeg("Male_Selected_z-scoreplots.jpg",width=2100,height = 800)
# # ggplot(m2zscores[m2zscores$Sex=='M' & (m2zscores$Measure=='Centiloid' | m2zscores$Measure=='TauBraak1' 
# #                                        | m2zscores$Measure=='TauBraak3&4'| m2zscores$Measure=='FDG' | 
# #                                          m2zscores$Measure=='MetaROIVol' | m2zscores$Measure=='PACC'),],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Male Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# # 
# # jpeg("Male_Selected_z-scoreplots_withZeroMin.jpg",width=2100,height = 800)
# # ggplot(m2zscores[m2zscores$Sex=='M' & m2zscores$z>=0 & (m2zscores$Measure=='Centiloid' | m2zscores$Measure=='TauBraak1' 
# #                                                         | m2zscores$Measure=='TauBraak3&4'| m2zscores$Measure=='FDG' | 
# #                                                           m2zscores$Measure=='MetaROIVol' | m2zscores$Measure=='PACC'),],
# #        aes(x = Age, y = z)) +
# #   geom_line(aes(color=Measure, group=Measure, alpha=Significance), size=1.1)+ 
# #   scale_color_manual(values = zscorecols1)+guides(alpha="none", size="none") +
# #   labs(y = 'z-scores',x = 'Age in Years',title = 'Male Z-scores') +
# #   facet_wrap(vars(APOE))
# # dev.off()
# 
# m2zscores$adj_z=0
# for (eachAPOE in unique(m2zscores$APOE)) {
#   for (eachMeasure in unique(m2zscores$Measure)) {
#     for (eachSex in unique(m2zscores$Sex)) {
#       m2zscores$adj_z[m2zscores$Measure==eachMeasure & m2zscores$APOE==eachAPOE & m2zscores$Sex==eachSex]=m2zscores$z[m2zscores$Measure==eachMeasure & m2zscores$APOE==eachAPOE & m2zscores$Sex==eachSex]-min(m2zscores$z[m2zscores$Measure==eachMeasure & m2zscores$APOE==eachAPOE & m2zscores$Sex==eachSex])
#     }
#   }
# }
# 
# jpeg("Female_ZeroShifted_z-scores.jpg",width=2100,height = 800)
# ggplot(m2zscores[m2zscores$Sex=='F',], aes(x = Age, y = adj_z)) +
#   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
#   scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
#   labs(y = 'Zero shifted z-scores',x = 'Age in Years',title = 'Zero Shifted Female Z-scores') +
#   facet_wrap(vars(APOE))
# dev.off()
# 
# jpeg("Male_ZeroShifted_z-scores.jpg",width=2100,height = 800)
# ggplot(m2zscores[m2zscores$Sex=='M',], aes(x = Age, y = adj_z)) +
#   geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
#   scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
#   labs(y = 'Zero shifted z-scores',x = 'Age in Years',title = 'Zero Shifted Male Z-scores') +
#   facet_wrap(vars(APOE))
# dev.off()
# 


m1zscores$GlobalAdj_z=0
for (eachMeasure in unique(m1zscores$Measure)) {
  m1zscores$GlobalAdj_z[m1zscores$Measure==eachMeasure]=m1zscores$z[m1zscores$Measure==eachMeasure]-min(m1zscores$z[m1zscores$Measure==eachMeasure])
}

jpeg("APOE_GlobalZeroShifted_z-scores.jpg",width=6000,height = 2400)
ggplot(m1zscores,aes(x = Age, y = GlobalAdj_z)) +
  geom_line(aes(color=Measure, group=Measure, alpha=Significance), size=4)+ 
  # geom_point(data=minmeasures,aes(x=Age, y=minVal, color=Measure))+
  scale_color_manual(breaks= unique(m1zscores$Measure),values = zscorecols)+ guides(alpha="none", size="none") +
  labs(y = 'Zero shifted z-scores',x = 'Age (Years)',title = 'Z-scores shifted by min of each measure') + 
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60), strip.text.x = element_text(size =60)) +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()


m2zscores$GlobalAdj_z=0
for (eachMeasure in unique(m2zscores$Measure)) {
  m2zscores$GlobalAdj_z[m2zscores$Measure==eachMeasure]=m2zscores$z[m2zscores$Measure==eachMeasure]-min(m2zscores$z[m2zscores$Measure==eachMeasure])
}
jpeg("Female_GlobalZeroShifted_z-scores.jpg",width=6000,height = 2400)
ggplot(m2zscores[m2zscores$Sex=='F',], aes(x = Age, y = GlobalAdj_z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=4)+ 
  scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
  labs(y = 'Zero shifted z-scores',x = 'Age (Years)',title = 'Female Z-scores shifted by min of each measure') +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()

jpeg("Male_GlobalZeroShifted_z-scores.jpg",width=6000,height = 2400)
ggplot(m2zscores[m2zscores$Sex=='M',], aes(x = Age, y = GlobalAdj_z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=4)+ 
  scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
  labs(y = 'Zero shifted z-scores',x = 'Age (Years)',title = 'Male Z-scores shifted by min of each measure') +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  facet_wrap(vars(APOE))+ scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
dev.off()

# Printing z-score plot legends
ggptmp=ggplot(m2zscores[m2zscores$Sex=='M',], aes(x = Age, y = GlobalAdj_z)) +
  geom_line(aes(color=Measure,group=Measure, alpha=Significance), size=1.1)+ 
  scale_color_manual(breaks= unique(m2zscores$Measure),values = zscorecols)+guides(alpha="none", size="none") +
  labs(y = 'Zero shifted z-scores',x = 'Age in Years',title = 'Male Z-scores shifted by min of each measure') +
  theme(legend.position="none",axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  facet_wrap(vars(APOE))+
  theme(legend.text = element_text(size = 60),legend.title = element_text(size = 60),legend.key.size = unit(2,'cm'),legend.direction = "horizontal")
t=get_legend(ggptmp)
grid.newpage()
grid.draw(t)
jpeg("z-scoreLegends.jpg",width=6000)
grid.draw(t) 
dev.off()
rm(ggptmp,t)

ggptmp=ggplot(tmpdd, aes(x = Age, y = reg,
                  group = RID, colour = APOE)) +
  geom_line(data=tmpdd[tmpdd$NumOfScans>1,],aes(group=RID)) +
  geom_point()+ 
  plt_labs +
  theme(axis.text = element_text(size = 60), axis.title = element_text(size = 60),strip.text.x = element_text(size =60)) +
  theme(legend.text = element_text(size = 100),legend.title = element_text(size = 100),legend.key.size = unit(3,'cm'),legend.direction = "horizontal")+
  guides(colour = guide_legend(override.aes = list(size=5)))
t=cowplot::get_legend(ggptmp)
grid.newpage()
grid.draw(t)
jpeg("DataScatterLegends.jpg",width=6000)
grid.draw(t) 
dev.off()
rm(ggptmp,t)




