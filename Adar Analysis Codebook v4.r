###  This Codebook is for the analysis of the Adar samples that we received.  The sequencing was done at U of M.          ###
###  Sequences were processed in mothur and analyzed using the following code in R by John R. Erb-Downward on January 2, 2018.   ###	
###  Input files were generated in Mothur v.1.39.0.  Followed SOP, SILVA database for align and RDP for classification.  ###
###  This is version 2

library(vegan)
library(gplots)
library(mvabund)
sem<-function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
geo.mean<-function(x, na.rm=TRUE){
  	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}
library(RColorBrewer)

#Here is also a useful helper function for figuring out variance that will work with ggplot
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

   #Geometric Mean functiongeo.mean<-function(x, na.rm=TRUE){
  geo.mean<-function(x, na.rm=TRUE){
  	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
	}
      
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
         geo.mean=geo.mean(xx[[col]],na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)       
        #geo.sd=10^(sd(decostand(xx[[col]],"log"),na.rm=na.rm))
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = paste("Mean.",measurevar,sep="")))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
	datac$geo.mean[which(datac$geo.mean==1)]<-0 #Geometric Mean will take 0 abundance and transform to exactly 1.  This looks for values that are exactly 1 and sets them back to 0.
	
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
	
    return(datac)
}


##########################################################################################################

 data.raw<-read.table("/Users/dunard/Dropbox/Adar_ananlysis/Adar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared", sep="\t", header=T,row.names=2)
data.trim<-data.raw[,-c(1,2)]
data.trim<-as.matrix(data.trim)
dim(data.trim)
# [1]  727 9050

rownames(data.trim)
#   [1] "AE_1_P1"                    "AE_2_P1"                    "AE_3_P1"                    "AE_4_P1"                    "AE_5_P2"                    "AE_6_P2"                   
#   [7] "AE_7_P2"                    "AE_8_P2"                    "Adar_3057_CBAL_1_a_P1"      "Adar_3057_CBAL_1_b_P1"      "Adar_3057_CBAL_1_c_P1"      "Adar_3057_EBC_a_P1"        
#  [13] "Adar_3057_EBC_b_P1"         "Adar_3057_EBC_c_P1"         "Adar_3057_LBAL_a_P1"        "Adar_3057_LBAL_b_P1"        "Adar_3057_LBAL_c_P1"        "Adar_3057_NC_P1"           
#  [19] "Adar_3057_Nasal_a_P1"       "Adar_3057_Nasal_b_P1"       "Adar_3057_Nasal_c_P1"       "Adar_3057_Oral_a_P1"        "Adar_3057_Oral_b_P1"        "Adar_3057_Oral_c_P1"       
#  [25] "Adar_3057_PreWash_a_P1"     "Adar_3057_PreWash_b_P1"     "Adar_3057_PreWash_c_P1"     "Adar_3057_RBAL_a_P1"        "Adar_3057_RBAL_b_P1"        "Adar_3057_RBAL_c_P1"       
#  [31] "Adar_3116_CBAL_1_a_P1"      "Adar_3116_CBAL_1_b_P1"      "Adar_3116_CBAL_1_c_P1"      "Adar_3116_CBAL_2_a_P1"      "Adar_3116_CBAL_2_b_P1"      "Adar_3116_CBAL_2_c_P1"     
#  [37] "Adar_3116_EBC_a_P1"         "Adar_3116_EBC_b_P1"         "Adar_3116_EBC_c_P1"         "Adar_3116_Filter_1_a_P1"    "Adar_3116_Filter_1_b_P1"    "Adar_3116_Filter_1_c_P1"   
#  [43] "Adar_3116_Filter_2_a_P1"    "Adar_3116_Filter_2_b_P1"    "Adar_3116_Filter_2_c_P1"    "Adar_3116_LBAL_a_P1"        "Adar_3116_LBAL_b_P1"        "Adar_3116_LBAL_c_P1"       
#  [49] "Adar_3116_NC_P1"            "Adar_3116_Nasal_a_P1"       "Adar_3116_Nasal_b_P1"       "Adar_3116_Nasal_c_P1"       "Adar_3116_Oral_a_P1"        "Adar_3116_Oral_b_P1"       
#  [55] "Adar_3116_Oral_c_P1"        "Adar_3116_PreWash_a_P1"     "Adar_3116_PreWash_b_P1"     "Adar_3116_PreWash_c_P1"     "Adar_3116_RBAL_a_P1"        "Adar_3116_RBAL_b_P1"       
#  [61] "Adar_3116_RBAL_c_P1"        "Adar_3172_CBAL_1_a_P1"      "Adar_3172_CBAL_1_b_P1"      "Adar_3172_CBAL_1_c_P1"      "Adar_3172_CBAL_2_a_P1"      "Adar_3172_CBAL_2_b_P1"     
#  [67] "Adar_3172_CBAL_2_c_P1"      "Adar_3172_EBC_a_P1"         "Adar_3172_EBC_b_P1"         "Adar_3172_EBC_c_P1"         "Adar_3172_Filter_1_a_P1"    "Adar_3172_Filter_1_b_P1"   
#  [73] "Adar_3172_Filter_1_c_P1"    "Adar_3172_Filter_2_a_P1"    "Adar_3172_Filter_2_b_P1"    "Adar_3172_Filter_2_c_P1"    "Adar_3172_LBAL_a_P1"        "Adar_3172_LBAL_b_P1"       
#  [79] "Adar_3172_LBAL_c_P1"        "Adar_3172_NC_P1"            "Adar_3172_Nasal_a_P1"       "Adar_3172_Nasal_b_P1"       "Adar_3172_Nasal_c_P1"       "Adar_3172_Oral_a_P1"       
#  [85] "Adar_3172_Oral_b_P1"        "Adar_3172_Oral_c_P1"        "Adar_3172_PreWash_a_P1"     "Adar_3172_PreWash_b_P1"     "Adar_3172_PreWash_c_P1"     "Adar_3172_RBAL_a_P1"       
#  [91] "Adar_3172_RBAL_b_P1"        "Adar_3172_RBAL_c_P1"        "Adar_3183_CBAL_1_a_P1"      "Adar_3183_CBAL_1_b_P1"      "Adar_3183_CBAL_1_c_P1"      "Adar_3183_CBAL_2_a_P1"     
#  [97] "Adar_3183_CBAL_2_b_P1"      "Adar_3183_CBAL_2_c_P1"      "Adar_3183_EBC_a_P1"         "Adar_3183_EBC_b_P1"         "Adar_3183_EBC_c_P1"         "Adar_3183_Filter_1_a_P1"   
# [103] "Adar_3183_Filter_1_b_P1"    "Adar_3183_Filter_1_c_P1"    "Adar_3183_Filter_2_a_P1"    "Adar_3183_Filter_2_b_P1"    "Adar_3183_Filter_2_c_P1"    "Adar_3183_LBAL_a_P1"       
# [109] "Adar_3183_LBAL_b_P1"        "Adar_3183_LBAL_c_P1"        "Adar_3183_NC_P1"            "Adar_3183_Nasal_a_P1"       "Adar_3183_Nasal_b_P1"       "Adar_3183_Nasal_c_P1"      
# [115] "Adar_3183_Oral_a_P1"        "Adar_3183_Oral_b_P1"        "Adar_3183_Oral_c_P1"        "Adar_3183_PreWash_a_P1"     "Adar_3183_PreWash_b_P1"     "Adar_3183_PreWash_c_P1"    
# [121] "Adar_3183_RBAL_a_P1"        "Adar_3183_RBAL_b_P1"        "Adar_3183_RBAL_c_P1"        "Adar_3194_CBAL_1_a_P1"      "Adar_3194_CBAL_1_b_P1"      "Adar_3194_CBAL_1_c_P1"     
# [127] "Adar_3194_CBAL_2_a_P1"      "Adar_3194_CBAL_2_b_P1"      "Adar_3194_CBAL_2_c_P1"      "Adar_3194_EBC_a_P1"         "Adar_3194_EBC_b_P1"         "Adar_3194_EBC_c_P1"        
# [133] "Adar_3194_Filter_1_a_P1"    "Adar_3194_Filter_1_b_P1"    "Adar_3194_Filter_1_c_P1"    "Adar_3194_Filter_2_a_P1"    "Adar_3194_Filter_2_b_P1"    "Adar_3194_Filter_2_c_P1"   
# [139] "Adar_3194_LBAL_a_P1"        "Adar_3194_LBAL_b_P1"        "Adar_3194_LBAL_c_P1"        "Adar_3194_NC_P1"            "Adar_3194_Nasal_a_P1"       "Adar_3194_Nasal_b_P1"      
# [145] "Adar_3194_Nasal_c_P1"       "Adar_3194_Oral_a_P1"        "Adar_3194_Oral_b_P1"        "Adar_3194_Oral_c_P1"        "Adar_3194_PreWash_a_P1"     "Adar_3194_PreWash_b_P1"    
# [151] "Adar_3194_PreWash_c_P1"     "Adar_3194_RBAL_a_P1"        "Adar_3194_RBAL_b_P1"        "Adar_3194_RBAL_c_P1"        "Adar_3219_CBAL_1_a_P1"      "Adar_3219_CBAL_1_b_P1"     
# [157] "Adar_3219_CBAL_1_c_P1"      "Adar_3219_CBAL_2_a_P1"      "Adar_3219_CBAL_2_b_P1"      "Adar_3219_CBAL_2_c_P1"      "Adar_3219_EBC_Blank_a_P1"   "Adar_3219_EBC_Blank_b_P1"  
# [163] "Adar_3219_EBC_Blank_c_P1"   "Adar_3219_EBC_a_P1"         "Adar_3219_EBC_b_P1"         "Adar_3219_EBC_c_P1"         "Adar_3219_Feces_Blank_a_P1" "Adar_3219_Feces_Blank_b_P1"
# [169] "Adar_3219_Feces_Blank_c_P1" "Adar_3219_Feces_P2"         "Adar_3219_Filter_1_a_P1"    "Adar_3219_Filter_1_b_P1"    "Adar_3219_Filter_1_c_P1"    "Adar_3219_Filter_2_a_P1"   
# [175] "Adar_3219_Filter_2_b_P1"    "Adar_3219_Filter_2_c_P1"    "Adar_3219_LBAL_a_P1"        "Adar_3219_LBAL_b_P1"        "Adar_3219_LBAL_c_P1"        "Adar_3219_NC_P1"           
# [181] "Adar_3219_Nasal_Blank_a_P1" "Adar_3219_Nasal_Blank_b_P1" "Adar_3219_Nasal_Blank_c_P1" "Adar_3219_Nasal_a_P1"       "Adar_3219_Nasal_b_P1"       "Adar_3219_Nasal_c_P1"      
# [187] "Adar_3219_Oral_a_P1"        "Adar_3219_Oral_b_P1"        "Adar_3219_Oral_c_P1"        "Adar_3219_PreWash_a_P1"     "Adar_3219_PreWash_b_P1"     "Adar_3219_PreWash_c_P1"    
# [193] "Adar_3219_RBAL_a_P1"        "Adar_3219_RBAL_b_P1"        "Adar_3219_RBAL_c_P1"        "Adar_3264_CBAL_1_a_P1"      "Adar_3264_CBAL_1_b_P1"      "Adar_3264_CBAL_1_c_P1"     
# [199] "Adar_3264_CBAL_2_a_P1"      "Adar_3264_CBAL_2_b_P1"      "Adar_3264_CBAL_2_c_P1"      "Adar_3264_EBC_Blank_a_P1"   "Adar_3264_EBC_Blank_b_P1"   "Adar_3264_EBC_Blank_c_P1"  
# [205] "Adar_3264_EBC_a_P1"         "Adar_3264_EBC_b_P1"         "Adar_3264_EBC_c_P1"         "Adar_3264_Feces_Blank_a_P1" "Adar_3264_Feces_Blank_b_P1" "Adar_3264_Feces_Blank_c_P1"
# [211] "Adar_3264_Filter_1_a_P1"    "Adar_3264_Filter_1_b_P1"    "Adar_3264_Filter_1_c_P1"    "Adar_3264_Filter_2_a_P1"    "Adar_3264_Filter_2_b_P1"    "Adar_3264_Filter_2_c_P1"   
# [217] "Adar_3264_LBAL_a_P1"        "Adar_3264_LBAL_b_P1"        "Adar_3264_LBAL_c_P1"        "Adar_3264_NC_P1"            "Adar_3264_Nasal_Blank_a_P1" "Adar_3264_Nasal_Blank_b_P1"
# [223] "Adar_3264_Nasal_Blank_c_P1" "Adar_3264_Nasal_a_P1"       "Adar_3264_Nasal_b_P1"       "Adar_3264_Nasal_c_P1"       "Adar_3264_Oral_a_P1"        "Adar_3264_Oral_b_P1"       
# [229] "Adar_3264_Oral_c_P1"        "Adar_3264_PreWash_a_P1"     "Adar_3264_PreWash_b_P1"     "Adar_3264_PreWash_c_P1"     "Adar_3264_RBAL_a_P1"        "Adar_3264_RBAL_b_P1"       
# [235] "Adar_3264_RBAL_c_P1"        "Adar_3301_CBAL_1_a_P1"      "Adar_3301_CBAL_1_b_P1"      "Adar_3301_CBAL_1_c_P1"      "Adar_3301_CBAL_2_a_P1"      "Adar_3301_CBAL_2_b_P1"     
# [241] "Adar_3301_CBAL_2_c_P1"      "Adar_3301_EBC_Blank_a_P1"   "Adar_3301_EBC_Blank_b_P1"   "Adar_3301_EBC_Blank_c_P1"   "Adar_3301_EBC_a_P1"         "Adar_3301_EBC_b_P1"        
# [247] "Adar_3301_EBC_c_P1"         "Adar_3301_Feces_Blank_a_P1" "Adar_3301_Feces_Blank_b_P1" "Adar_3301_Feces_Blank_c_P1" "Adar_3301_Filter_1_a_P1"    "Adar_3301_Filter_1_b_P1"   
# [253] "Adar_3301_Filter_1_c_P1"    "Adar_3301_Filter_2_a_P1"    "Adar_3301_Filter_2_b_P1"    "Adar_3301_Filter_2_c_P1"    "Adar_3301_LBAL_a_P1"        "Adar_3301_LBAL_b_P1"       
# [259] "Adar_3301_LBAL_c_P1"        "Adar_3301_NC_P1"            "Adar_3301_Nasal_Blank_a_P1" "Adar_3301_Nasal_Blank_b_P1" "Adar_3301_Nasal_Blank_c_P1" "Adar_3301_Nasal_a_P1"      
# [265] "Adar_3301_Nasal_b_P1"       "Adar_3301_Nasal_c_P1"       "Adar_3301_Oral_a_P1"        "Adar_3301_Oral_b_P1"        "Adar_3301_Oral_c_P1"        "Adar_3301_PreWash_a_P1"    
# [271] "Adar_3301_PreWash_b_P1"     "Adar_3301_PreWash_c_P1"     "Adar_3301_RBAL_a_P1"        "Adar_3301_RBAL_b_P1"        "Adar_3301_RBAL_c_P1"        "Adar_3356_CBAL_1_a_P1"     
# [277] "Adar_3356_CBAL_1_b_P1"      "Adar_3356_CBAL_1_c_P1"      "Adar_3356_CBAL_2_a_P1"      "Adar_3356_CBAL_2_b_P1"      "Adar_3356_CBAL_2_c_P1"      "Adar_3356_EBC_Blank_a_P1"  
# [283] "Adar_3356_EBC_Blank_b_P1"   "Adar_3356_EBC_Blank_c_P1"   "Adar_3356_EBC_a_P1"         "Adar_3356_EBC_b_P1"         "Adar_3356_EBC_c_P1"         "Adar_3356_Feces_Blank_a_P1"
# [289] "Adar_3356_Feces_Blank_b_P1" "Adar_3356_Feces_Blank_c_P1" "Adar_3356_Filter_1_a_P1"    "Adar_3356_Filter_1_b_P1"    "Adar_3356_Filter_1_c_P1"    "Adar_3356_Filter_2_a_P1"   
# [295] "Adar_3356_Filter_2_b_P1"    "Adar_3356_Filter_2_c_P1"    "Adar_3356_LBAL_a_P1"        "Adar_3356_LBAL_b_P1"        "Adar_3356_LBAL_c_P1"        "Adar_3356_NC_P1"           
# [301] "Adar_3356_Nasal_Blank_a_P1" "Adar_3356_Nasal_Blank_b_P1" "Adar_3356_Nasal_Blank_c_P1" "Adar_3356_Nasal_a_P1"       "Adar_3356_Nasal_b_P1"       "Adar_3356_Nasal_c_P1"      
# [307] "Adar_3356_Oral_a_P1"        "Adar_3356_Oral_b_P1"        "Adar_3356_Oral_c_P1"        "Adar_3356_PreWash_a_P1"     "Adar_3356_PreWash_b_P1"     "Adar_3356_PreWash_c_P1"    
# [313] "Adar_3356_RBAL_a_P1"        "Adar_3356_RBAL_b_P1"        "Adar_3356_RBAL_c_P1"        "Adar_3378_CBAL_1_a_P1"      "Adar_3378_CBAL_1_b_P1"      "Adar_3378_CBAL_1_c_P1"     
# [319] "Adar_3378_CBAL_2_a_P1"      "Adar_3378_CBAL_2_b_P1"      "Adar_3378_CBAL_2_c_P1"      "Adar_3378_EBC_Blank_a_P1"   "Adar_3378_EBC_Blank_b_P1"   "Adar_3378_EBC_Blank_c_P1"  
# [325] "Adar_3378_EBC_a_P1"         "Adar_3378_EBC_b_P1"         "Adar_3378_EBC_c_P2"         "Adar_3378_Feces_Blank_a_P1" "Adar_3378_Feces_Blank_b_P1" "Adar_3378_Feces_Blank_c_P2"
# [331] "Adar_3378_Feces_P2"         "Adar_3378_Filter_1_a_P1"    "Adar_3378_Filter_1_b_P1"    "Adar_3378_Filter_1_c_P2"    "Adar_3378_Filter_2_a_P1"    "Adar_3378_Filter_2_b_P1"   
# [337] "Adar_3378_Filter_2_c_P2"    "Adar_3378_LBAL_a_P1"        "Adar_3378_LBAL_b_P1"        "Adar_3378_LBAL_c_P1"        "Adar_3378_NC_P1"            "Adar_3378_NC_P2"           
# [343] "Adar_3378_Nasal_Blank_a_P1" "Adar_3378_Nasal_Blank_b_P1" "Adar_3378_Nasal_Blank_c_P2" "Adar_3378_Nasal_a_P1"       "Adar_3378_Nasal_b_P1"       "Adar_3378_Nasal_c_P2"      
# [349] "Adar_3378_Oral_a_P1"        "Adar_3378_Oral_b_P1"        "Adar_3378_Oral_c_P2"        "Adar_3378_PreWash_a_P1"     "Adar_3378_PreWash_b_P1"     "Adar_3378_PreWash_c_P2"    
# [355] "Adar_3378_RBAL_a_P1"        "Adar_3378_RBAL_b_P1"        "Adar_3378_RBAL_c_P1"        "Adar_3389_CBAL_1_a_P2"      "Adar_3389_CBAL_1_b_P2"      "Adar_3389_CBAL_1_c_P2"     
# [361] "Adar_3389_EBC_Blank_a_P2"   "Adar_3389_EBC_Blank_b_P2"   "Adar_3389_EBC_Blank_c_P2"   "Adar_3389_EBC_Control_a_P2" "Adar_3389_EBC_Control_b_P2" "Adar_3389_EBC_Control_c_P2"
# [367] "Adar_3389_EBC_a_P2"         "Adar_3389_EBC_b_P2"         "Adar_3389_EBC_c_P2"         "Adar_3389_Feces_Blank_a_P2" "Adar_3389_Feces_Blank_b_P2" "Adar_3389_Feces_Blank_c_P2"
# [373] "Adar_3389_Feces_P2"         "Adar_3389_Filter_1_a_P2"    "Adar_3389_Filter_1_b_P2"    "Adar_3389_Filter_1_c_P2"    "Adar_3389_Filter_2_a_P2"    "Adar_3389_Filter_2_b_P2"   
# [379] "Adar_3389_Filter_2_c_P2"    "Adar_3389_LBAL_a_P2"        "Adar_3389_LBAL_b_P2"        "Adar_3389_LBAL_c_P2"        "Adar_3389_NC_P2"            "Adar_3389_Nasal_Blank_a_P2"
# [385] "Adar_3389_Nasal_Blank_b_P2" "Adar_3389_Nasal_Blank_c_P2" "Adar_3389_Nasal_a_P2"       "Adar_3389_Nasal_b_P2"       "Adar_3389_Nasal_c_P2"       "Adar_3389_Oral_a_P2"       
# [391] "Adar_3389_Oral_b_P2"        "Adar_3389_Oral_c_P2"        "Adar_3389_PreWash_a_P2"     "Adar_3389_PreWash_b_P2"     "Adar_3389_PreWash_c_P2"     "Adar_3389_RBAL_a_P2"       
# [397] "Adar_3389_RBAL_b_P2"        "Adar_3389_RBAL_c_P2"        "Adar_3482_CBAL_1_a_P2"      "Adar_3482_CBAL_1_b_P2"      "Adar_3482_CBAL_1_c_P2"      "Adar_3482_EBC_a_P2"        
# [403] "Adar_3482_EBC_b_P2"         "Adar_3482_EBC_c_P2"         "Adar_3482_Feces_P2"         "Adar_3482_Filter_1_a_P2"    "Adar_3482_Filter_1_b_P2"    "Adar_3482_Filter_1_c_P2"   
# [409] "Adar_3482_Filter_2_a_P2"    "Adar_3482_Filter_2_b_P2"    "Adar_3482_Filter_2_c_P2"    "Adar_3482_LBAL_a_P2"        "Adar_3482_LBAL_b_P2"        "Adar_3482_LBAL_c_P2"       
# [415] "Adar_3482_NC_P2"            "Adar_3482_Nasal_a_P2"       "Adar_3482_Nasal_b_P2"       "Adar_3482_Nasal_c_P2"       "Adar_3482_Oral_a_P2"        "Adar_3482_Oral_b_P2"       
# [421] "Adar_3482_Oral_c_P2"        "Adar_3482_PreWash_a_P2"     "Adar_3482_PreWash_b_P2"     "Adar_3482_PreWash_c_P2"     "Adar_3482_RBAL_a_P2"        "Adar_3482_RBAL_b_P2"       
# [427] "Adar_3482_RBAL_c_P2"        "Adar_3493_CBAL_1_a_P2"      "Adar_3493_CBAL_1_b_P2"      "Adar_3493_CBAL_1_c_P2"      "Adar_3493_EBC_Blank_a_P2"   "Adar_3493_EBC_Blank_b_P2"  
# [433] "Adar_3493_EBC_Blank_c_P2"   "Adar_3493_EBC_a_P2"         "Adar_3493_EBC_b_P2"         "Adar_3493_EBC_c_P2"         "Adar_3493_Feces_P2"         "Adar_3493_Filter_1_a_P2"   
# [439] "Adar_3493_Filter_1_b_P2"    "Adar_3493_Filter_1_c_P2"    "Adar_3493_Filter_2_a_P2"    "Adar_3493_Filter_2_b_P2"    "Adar_3493_Filter_2_c_P2"    "Adar_3493_LBAL_a_P2"       
# [445] "Adar_3493_LBAL_b_P2"        "Adar_3493_LBAL_c_P2"        "Adar_3493_NC_P2"            "Adar_3493_Nasal_Blank_a_P2" "Adar_3493_Nasal_Blank_b_P2" "Adar_3493_Nasal_Blank_c_P2"
# [451] "Adar_3493_Nasal_a_P2"       "Adar_3493_Nasal_b_P2"       "Adar_3493_Nasal_c_P2"       "Adar_3493_Oral_a_P2"        "Adar_3493_Oral_b_P2"        "Adar_3493_Oral_c_P2"       
# [457] "Adar_3493_PreWash_a_P2"     "Adar_3493_PreWash_b_P2"     "Adar_3493_PreWash_c_P2"     "Adar_3493_RBAL_a_P2"        "Adar_3493_RBAL_b_P2"        "Adar_3493_RBAL_c_P2"       
# [463] "Adar_3518_CBAL_1_a_P2"      "Adar_3518_CBAL_1_b_P2"      "Adar_3518_CBAL_1_c_P2"      "Adar_3518_EBC_a_P2"         "Adar_3518_EBC_b_P2"         "Adar_3518_EBC_c_P2"        
# [469] "Adar_3518_Feces_P2"         "Adar_3518_Filter_1_a_P2"    "Adar_3518_Filter_1_b_P2"    "Adar_3518_Filter_1_c_P2"    "Adar_3518_Filter_2_a_P2"    "Adar_3518_Filter_2_b_P2"   
# [475] "Adar_3518_Filter_2_c_P2"    "Adar_3518_LBAL_a_P2"        "Adar_3518_LBAL_b_P2"        "Adar_3518_LBAL_c_P2"        "Adar_3518_NC_P2"            "Adar_3518_Nasal_a_P2"      
# [481] "Adar_3518_Nasal_b_P2"       "Adar_3518_Nasal_c_P2"       "Adar_3518_Oral_a_P2"        "Adar_3518_Oral_b_P2"        "Adar_3518_Oral_c_P2"        "Adar_3518_PreWash_a_P2"    
# [487] "Adar_3518_PreWash_b_P2"     "Adar_3518_PreWash_c_P2"     "Adar_3518_RBAL_a_P2"        "Adar_3518_RBAL_b_P2"        "Adar_3518_RBAL_c_P2"        "Adar_3552_CBAL_1_a_P2"     
# [493] "Adar_3552_CBAL_1_b_P2"      "Adar_3552_CBAL_1_c_P2"      "Adar_3552_CBAL_2_a_P2"      "Adar_3552_CBAL_2_b_P2"      "Adar_3552_CBAL_2_c_P2"      "Adar_3552_EBC_1_a_P2"      
# [499] "Adar_3552_EBC_1_b_P2"       "Adar_3552_EBC_1_c_P2"       "Adar_3552_EBC_2_a_P2"       "Adar_3552_EBC_2_b_P2"       "Adar_3552_EBC_2_c_P2"       "Adar_3552_Feces_P2"        
# [505] "Adar_3552_Filter_1_a_P2"    "Adar_3552_Filter_1_b_P2"    "Adar_3552_Filter_1_c_P2"    "Adar_3552_Filter_2_a_P2"    "Adar_3552_Filter_2_b_P2"    "Adar_3552_Filter_2_c_P2"   
# [511] "Adar_3552_LBAL_a_P2"        "Adar_3552_LBAL_b_P2"        "Adar_3552_LBAL_c_P2"        "Adar_3552_NC_P2"            "Adar_3552_Nasal_a_P2"       "Adar_3552_Nasal_b_P2"      
# [517] "Adar_3552_Nasal_c_P2"       "Adar_3552_Oral_a_P2"        "Adar_3552_Oral_b_P2"        "Adar_3552_Oral_c_P2"        "Adar_3552_PreWash_a_P2"     "Adar_3552_PreWash_b_P2"    
# [523] "Adar_3552_PreWash_c_P2"     "Adar_3552_RBAL_a_P2"        "Adar_3552_RBAL_b_P2"        "Adar_3552_RBAL_c_P2"        "Adar_3585_CBAL_1_a_P2"      "Adar_3585_CBAL_1_b_P2"     
# [529] "Adar_3585_CBAL_1_c_P2"      "Adar_3585_CBAL_2_a_P2"      "Adar_3585_CBAL_2_b_P2"      "Adar_3585_CBAL_2_c_P2"      "Adar_3585_Feces_P2"         "Adar_3585_LBAL_a_P2"       
# [535] "Adar_3585_LBAL_b_P2"        "Adar_3585_LBAL_c_P2"        "Adar_3585_NC_P2"            "Adar_3585_Nasal_a_P2"       "Adar_3585_Nasal_b_P2"       "Adar_3585_Nasal_c_P2"      
# [541] "Adar_3585_Oral_a_P2"        "Adar_3585_Oral_b_P2"        "Adar_3585_Oral_c_P2"        "Adar_3585_PreWash_a_P2"     "Adar_3585_PreWash_b_P2"     "Adar_3585_PreWash_c_P2"    
# [547] "Adar_3585_RBAL_a_P2"        "Adar_3585_RBAL_b_P2"        "Adar_3585_RBAL_c_P2"        "Adar_3600_CBAL_1_a_P2"      "Adar_3600_CBAL_1_b_P2"      "Adar_3600_CBAL_1_c_P2"     
# [553] "Adar_3600_EBC_a_P2"         "Adar_3600_EBC_b_P2"         "Adar_3600_EBC_c_P2"         "Adar_3600_Feces_P2"         "Adar_3600_Filter_1_a_P2"    "Adar_3600_Filter_1_b_P2"   
# [559] "Adar_3600_Filter_1_c_P2"    "Adar_3600_Filter_2_a_P2"    "Adar_3600_Filter_2_b_P2"    "Adar_3600_Filter_2_c_P2"    "Adar_3600_LBAL_a_P2"        "Adar_3600_LBAL_b_P2"       
# [565] "Adar_3600_LBAL_c_P2"        "Adar_3600_NC_P2"            "Adar_3600_Nasal_a_P2"       "Adar_3600_Nasal_b_P2"       "Adar_3600_Nasal_c_P2"       "Adar_3600_Oral_a_P2"       
# [571] "Adar_3600_Oral_b_P2"        "Adar_3600_Oral_c_P2"        "Adar_3600_PreWash_a_P2"     "Adar_3600_PreWash_b_P2"     "Adar_3600_PreWash_c_P2"     "Adar_3600_RBAL_a_P2"       
# [577] "Adar_3600_RBAL_b_P2"        "Adar_3600_RBAL_c_P2"        "Adar_3611_CBAL_1_a_P2"      "Adar_3611_CBAL_1_b_P2"      "Adar_3611_CBAL_1_c_P2"      "Adar_3611_CBAL_2_a_P2"     
# [583] "Adar_3611_CBAL_2_b_P2"      "Adar_3611_CBAL_2_c_P2"      "Adar_3611_EBC_a_P2"         "Adar_3611_EBC_b_P2"         "Adar_3611_EBC_c_P2"         "Adar_3611_Feces_P2"        
# [589] "Adar_3611_Filter_1_a_P2"    "Adar_3611_Filter_1_b_P2"    "Adar_3611_Filter_1_c_P2"    "Adar_3611_Filter_2_a_P2"    "Adar_3611_Filter_2_b_P2"    "Adar_3611_Filter_2_c_P2"   
# [595] "Adar_3611_LBAL_a_P2"        "Adar_3611_LBAL_b_P2"        "Adar_3611_LBAL_c_P2"        "Adar_3611_NC_P2"            "Adar_3611_Nasal_a_P2"       "Adar_3611_Nasal_b_P2"      
# [601] "Adar_3611_Nasal_c_P2"       "Adar_3611_Oral_a_P2"        "Adar_3611_Oral_b_P2"        "Adar_3611_Oral_c_P2"        "Adar_3611_PreWash_a_P2"     "Adar_3611_PreWash_b_P2"    
# [607] "Adar_3611_PreWash_c_P2"     "Adar_3611_RBAL_a_P2"        "Adar_3611_RBAL_b_P2"        "Adar_3611_RBAL_c_P2"        "Adar_3644_CBAL_1_a_P2"      "Adar_3644_CBAL_1_b_P2"     
# [613] "Adar_3644_CBAL_1_c_P2"      "Adar_3644_CBAL_2_a_P2"      "Adar_3644_CBAL_2_b_P2"      "Adar_3644_CBAL_2_c_P2"      "Adar_3644_EBC_a_P2"         "Adar_3644_EBC_b_P2"        
# [619] "Adar_3644_EBC_c_P2"         "Adar_3644_Feces_P2"         "Adar_3644_Filter_1_a_P2"    "Adar_3644_Filter_1_b_P2"    "Adar_3644_Filter_1_c_P2"    "Adar_3644_Filter_2_a_P2"   
# [625] "Adar_3644_Filter_2_b_P2"    "Adar_3644_Filter_2_c_P2"    "Adar_3644_LBAL_a_P2"        "Adar_3644_LBAL_b_P2"        "Adar_3644_LBAL_c_P2"        "Adar_3644_NC_P2"           
# [631] "Adar_3644_Nasal_a_P2"       "Adar_3644_Nasal_b_P2"       "Adar_3644_Nasal_c_P2"       "Adar_3644_Oral_a_P2"        "Adar_3644_Oral_b_P2"        "Adar_3644_Oral_c_P2"       
# [637] "Adar_3644_PreWash_a_P2"     "Adar_3644_PreWash_b_P2"     "Adar_3644_PreWash_c_P2"     "Adar_3644_RBAL_a_P2"        "Adar_3644_RBAL_b_P2"        "Adar_3644_RBAL_c_P2"       
# [643] "Adar_3769_Feces_P2"         "Adar_3840_CBAL_1_a_P2"      "Adar_3840_CBAL_1_b_P2"      "Adar_3840_CBAL_1_c_P2"      "Adar_3840_Feces_P2"         "Adar_3840_LBAL_a_P2"       
# [649] "Adar_3840_LBAL_b_P2"        "Adar_3840_LBAL_c_P2"        "Adar_3840_NC_P2"            "Adar_3840_Nasal_a_P2"       "Adar_3840_Nasal_b_P2"       "Adar_3840_Nasal_c_P2"      
# [655] "Adar_3840_Oral_a_P2"        "Adar_3840_Oral_b_P2"        "Adar_3840_Oral_c_P2"        "Adar_3840_PreWash_a_P2"     "Adar_3840_PreWash_b_P2"     "Adar_3840_PreWash_c_P2"    
# [661] "Adar_3840_RBAL_a_P2"        "Adar_3840_RBAL_b_P2"        "Adar_3840_RBAL_c_P2"        "Adar_3851_CBAL_1_a_P2"      "Adar_3851_CBAL_1_b_P2"      "Adar_3851_CBAL_1_c_P2"     
# [667] "Adar_3851_CBAL_2_a_P2"      "Adar_3851_CBAL_2_b_P2"      "Adar_3851_CBAL_2_c_P2"      "Adar_3851_LBAL_a_P2"        "Adar_3851_LBAL_b_P2"        "Adar_3851_LBAL_c_P2"       
# [673] "Adar_3851_NC_P2"            "Adar_3851_Nasal_a_P2"       "Adar_3851_Nasal_b_P2"       "Adar_3851_Nasal_c_P2"       "Adar_3851_Oral_a_P2"        "Adar_3851_Oral_b_P2"       
# [679] "Adar_3851_Oral_c_P2"        "Adar_3851_PreWash_a_P2"     "Adar_3851_PreWash_b_P2"     "Adar_3851_PreWash_c_P2"     "Adar_3851_RBAL_a_P2"        "Adar_3851_RBAL_b_P2"       
# [685] "Adar_3851_RBAL_c_P2"        "BLANK1_P1"                  "BLANK2_P1"                  "BLANK3_P1"                  "BLANK4_P1"                  "EMPTY1_P2"                 
# [691] "IsoCtrl_1_a_P1"             "IsoCtrl_1_b_P1"             "IsoCtrl_1_c_P1"             "IsoCtrl_2_a_P1"             "IsoCtrl_2_b_P1"             "IsoCtrl_2_c_P1"            
# [697] "IsoCtrl_3_a_P1"             "IsoCtrl_3_b_P1"             "IsoCtrl_3_c_P1"             "IsoCtrl_4_a_P1"             "IsoCtrl_4_b_P1"             "IsoCtrl_4_c_P2"            
# [703] "IsoCtrl_5_a_P2"             "IsoCtrl_5_b_P2"             "IsoCtrl_5_c_P2"             "IsoCtrl_6_a_P2"             "IsoCtrl_6_b_P2"             "IsoCtrl_6_c_P2"            
# [709] "IsoCtrl_7_a_P2"             "IsoCtrl_7_b_P2"             "IsoCtrl_7_c_P2"             "PCRwater_A_P2"              "PCRwater_B_P2"              "PCRwater_C_P2"             
# [715] "PCRwater_D_P2"              "PCRwater_a_P1"              "PCRwater_b_P1"              "PCRwater_c_P1"              "PCRwater_d_P1"              "ZymoMock_A_P2"             
# [721] "ZymoMock_B_P2"              "ZymoMock_C_P2"              "ZymoMock_D_P2"              "ZymoMock_a_P1"              "ZymoMock_b_P1"              "ZymoMock_c_P1"             
# [727] "ZymoMock_d_P1" 

# Lots of samples here so lets break down the sample codes.  P1 or P2 has to do with which MiSeq plate.  The number ofter "Adar" is the patient ID.  
# Sample types: IsoCtrl - Isolation Control; PCRwater - PCR grade water; ZymoMock - Mock Community; BLANK - nothing; CBAL - combinded rt and lft BAL; LBAL - lft BAL; 
# NC = Non-Sterile Saline Control; Nasal - Nasal Swab; Oral - Oral wash; PreWash - Scope Pre-Wash; RBAL - rt BAL; Feces - Stool sample; Filter_1 - room air filter;
# Filter_2 - control filter; EBC - Exhaled Breath Condensate; EBC_Blank - EBC Control.  The "a", "b", "c" designations are repeats of the same sample.  If the "a" is "A" 

# it is to distinguish the samples on different plates.

# Now to set up thresholds and remove OTUs below 0.1% of the population
data.tmp<-decostand(data.trim,"total")*100
neg1<-which(data.tmp<0.1)
data.trim[neg1]<-0
pos<-which(colSums(data.trim)>0)
data.good<-data.trim[,pos]
dim(data.good)
# [1]  727 2478  # That's a lot more than usual, but of course there are fecal samples.

# Now the taxonomy info.

tax1 <- read.table("/Users/dunard/Dropbox/Adar_ananlysis/Adar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy", sep = "\t", row.names = 1, header = T)
tax2 <- read.table("/Users/dunard/Dropbox/Adar_ananlysis/Adar.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy", sep = ";", skip = 1)
data.taxonomy <- cbind(tax1[, 1], sub("\\(.*.\\)", "", tax2[, 2]), sub("\\(.*.\\)", "", tax2[, 3]), sub("\\(.*.\\)", "", tax2[, 4]), sub("\\(.*.\\)", "", tax2[, 5]), sub("\\(.*.\\)", "", tax2[, 6]))
rownames(data.taxonomy) <- rownames(tax1)
colnames(data.taxonomy) <- c("Size", "Phylum", "Class", "Order", "Family", "Genus")
data.good.taxonomy <- data.taxonomy[colnames(data.good),]
data.good.taxonomy <- as.data.frame(data.good.taxonomy)


# Let's start with procedural controls:  Mock, Water and Isolation Control.  I guess we can add in blank, but I'm ready for those to look crazy.

mock<-grep("Zymo",rownames(data.good))
iso<-grep("IsoCtrl",rownames(data.good))
water<-grep("PCRwater",rownames(data.good))

data.good.norm<-decostand(data.good,"total")*100 # % Abundance table

mock.order<-names(sort(colMeans(data.good.norm[mock,]),decreasing=T))

#trying to make stacked barplots with ggplot2
library(reshape2)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(ComplexHeatmap)

## This section details how to do stacked barplots with ggplot2.  After graphing things out I concluded that I'd rather a grid of bar graphs.  So it's all commented out.

# mock
# mock.gg<-melt(data.good.norm[mock,mock.order[1:15]]) #reformats the data (see below)
# names(mock.gg)<-c("Name","OTU","Percentage")
# mock.gg[1:10,]
#             Name     OTU Percentage
# 1  ZymoMock_A_P2 Otu0013   17.93411
# 2  ZymoMock_B_P2 Otu0013   16.82218
# 3  ZymoMock_C_P2 Otu0013   19.13741
# 4  ZymoMock_D_P2 Otu0013   16.32108
# 5  ZymoMock_a_P1 Otu0013   22.37149
# 6  ZymoMock_b_P1 Otu0013   20.70973
# 7  ZymoMock_c_P1 Otu0013   18.51426
# 8  ZymoMock_d_P1 Otu0013   20.36144
# 9  ZymoMock_A_P2 Otu0017   18.33127
# 10 ZymoMock_B_P2 Otu0017   20.04041
# 
# The next step is to set the labels so they display nice in the legend
# mock.gg$OTU <- factor(mock.gg$OTU, levels = levels(mock.gg$OTU), labels = paste(data.good.taxonomy[mock.order[1:15],6]," (",levels(mock.gg$OTU),")",sep=""))
# 
# 
# plot_m1<-ggplot(data=mock.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Mock Community Samples", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1))	
# plot_m1
# 
# Isolation Controls
# iso.order<-names(sort(colMeans(data.good.norm[iso,]),decreasing=T))
# temp<-rbind(data.good.norm[iso,iso.order[1:50]],apply(data.good.norm[iso,iso.order[1:50]],2,mean),apply(data.good.norm[iso,iso.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso.gg<-melt(temp) #reformats the data (see below)
# names(iso.gg)<-c("Name","OTU","Percentage")
# iso.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso.gg$OTU <- factor(iso.gg$OTU, levels = levels(iso.gg$OTU), labels = paste(data.good.taxonomy[iso.order[1:50],6]," (",levels(iso.gg$OTU),")",sep=""))
# 
# 
# plot_iso<-ggplot(data=iso.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Controls Samples", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso	
# 
# plot_iso+ scale_fill_brewer(palette = "Set1") will change the colors to brewer palettes but these are limited to the number of colors in the palette.
# the iso controls plots are showing some sample clustering so now going to plot those along with the mean and geo.mean
# 
# isoCtrl1
# iso1.order<-names(sort(colMeans(data.good.norm[iso[1:3],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[1:3],iso1.order[1:50]],apply(data.good.norm[iso[1:3],iso1.order[1:50]],2,mean),apply(data.good.norm[iso[1:3],iso1.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso1.gg<-melt(temp) #reformats the data (see below)
# names(iso1.gg)<-c("Name","OTU","Percentage")
# iso1.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso1.gg$OTU <- factor(iso1.gg$OTU, levels = levels(iso1.gg$OTU), labels = paste(data.good.taxonomy[iso1.order[1:50],6]," (",levels(iso1.gg$OTU),")",sep=""))
# 
# 
# plot_iso1<-ggplot(data=iso1.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso1.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	facet_grid(Name~.) # creates a grid of bargraphs
# 	labs(title="Cross-Comparison of Isolation Control-1 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso1	
# 
# isoCtrl2
# iso2.order<-names(sort(colMeans(data.good.norm[iso[4:6],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[4:6],iso2.order[1:50]],apply(data.good.norm[iso[4:6],iso2.order[1:50]],2,mean),apply(data.good.norm[iso[4:6],iso2.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso2.gg<-melt(temp) #reformats the data (see below)
# names(iso2.gg)<-c("Name","OTU","Percentage")
# iso2.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso2.gg$OTU <- factor(iso2.gg$OTU, levels = levels(iso2.gg$OTU), labels = paste(data.good.taxonomy[iso2.order[1:50],6]," (",levels(iso2.gg$OTU),")",sep=""))
# 
# 
# plot_iso2<-ggplot(data=iso2.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-2 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso2	
# 
# isoCtrl3
# iso3.order<-names(sort(colMeans(data.good.norm[iso[7:9],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[7:9],iso3.order[1:50]],apply(data.good.norm[iso[7:9],iso3.order[1:50]],2,mean),apply(data.good.norm[iso[7:9],iso3.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso3.gg<-melt(temp) #reformats the data (see below)
# names(iso3.gg)<-c("Name","OTU","Percentage")
# iso3.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso3.gg$OTU <- factor(iso3.gg$OTU, levels = levels(iso3.gg$OTU), labels = paste(data.good.taxonomy[iso3.order[1:50],6]," (",levels(iso3.gg$OTU),")",sep=""))
# 
# 
# plot_iso3<-ggplot(data=iso3.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-3 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso3	
# 
# isoCtrl4
# iso4.order<-names(sort(colMeans(data.good.norm[iso[10:12],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[10:12],iso4.order[1:50]],apply(data.good.norm[iso[10:12],iso4.order[1:50]],2,mean),apply(data.good.norm[iso[10:12],iso4.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso4.gg<-melt(temp) #reformats the data (see below)
# names(iso4.gg)<-c("Name","OTU","Percentage")
# iso4.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso4.gg$OTU <- factor(iso4.gg$OTU, levels = levels(iso4.gg$OTU), labels = paste(data.good.taxonomy[iso4.order[1:50],6]," (",levels(iso4.gg$OTU),")",sep=""))
# 
# 
# plot_iso4<-ggplot(data=iso4.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-4 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso4	
# 
# isoCtrl5
# iso5.order<-names(sort(colMeans(data.good.norm[iso[13:15],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[13:15],iso5.order[1:50]],apply(data.good.norm[iso[13:15],iso5.order[1:50]],2,mean),apply(data.good.norm[iso[13:15],iso5.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso5.gg<-melt(temp) #reformats the data (see below)
# names(iso5.gg)<-c("Name","OTU","Percentage")
# iso5.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso5.gg$OTU <- factor(iso5.gg$OTU, levels = levels(iso5.gg$OTU), labels = paste(data.good.taxonomy[iso5.order[1:50],6]," (",levels(iso5.gg$OTU),")",sep=""))
# 
# 
# plot_iso5<-ggplot(data=iso5.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-5 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso5	
# 
# isoCtrl6
# iso6.order<-names(sort(colMeans(data.good.norm[iso[16:18],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[16:18],iso6.order[1:50]],apply(data.good.norm[iso[16:18],iso6.order[1:50]],2,mean),apply(data.good.norm[iso[16:18],iso6.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso6.gg<-melt(temp) #reformats the data (see below)
# names(iso6.gg)<-c("Name","OTU","Percentage")
# iso6.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso6.gg$OTU <- factor(iso6.gg$OTU, levels = levels(iso6.gg$OTU), labels = paste(data.good.taxonomy[iso6.order[1:50],6]," (",levels(iso6.gg$OTU),")",sep=""))
# 
# 
# plot_iso6<-ggplot(data=iso6.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-6 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso6	
# 
# isoCtrl7
# iso7.order<-names(sort(colMeans(data.good.norm[iso[19:21],]),decreasing=T))
# temp<-rbind(data.good.norm[iso[19:21],iso7.order[1:50]],apply(data.good.norm[iso[19:21],iso7.order[1:50]],2,mean),apply(data.good.norm[iso[19:21],iso7.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# iso7.gg<-melt(temp) #reformats the data (see below)
# names(iso7.gg)<-c("Name","OTU","Percentage")
# iso7.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# iso7.gg$OTU <- factor(iso7.gg$OTU, levels = levels(iso7.gg$OTU), labels = paste(data.good.taxonomy[iso7.order[1:50],6]," (",levels(iso7.gg$OTU),")",sep=""))
# 
# 
# plot_iso7<-ggplot(data=iso7.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of Isolation Control-7 Replicates", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_iso7	
# 
# Water samples
# water.order<-names(sort(colMeans(data.good.norm[water,]),decreasing=T))
# temp<-rbind(data.good.norm[water,water.order[1:50]],apply(data.good.norm[water,water.order[1:50]],2,mean),apply(data.good.norm[water,water.order[1:50]],2,geo.mean))
# rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
# temp[which(temp==1)]<-0
# water.gg<-melt(temp) #reformats the data (see below)
# names(water.gg)<-c("Name","OTU","Percentage")
# water.gg[1:10,]
# 
# The next step is to set the labels so they display nice in the legend
# water.gg$OTU <- factor(water.gg$OTU, levels = levels(water.gg$OTU), labels = paste(data.good.taxonomy[water.order[1:50],6]," (",levels(water.gg$OTU),")",sep=""))
# 
# 
# plot_water<-ggplot(data=water.gg,aes(y=Percentage,x=Name,fill=OTU))+
# 	theme_bw()+ #This theme gets rid of the grey background
# 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
# 	labs(title="Cross-Comparison of No Template Controls Samples", x="Sample",y="% Relative Abundance")+
# 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
# 	theme(axis.text.x=element_text(angle=90,hjust=1)) #angles axis text
# plot_water	

# Diagnostic grid of replicates Mean and Geo.Mean.

# Mock Community
mock.order<-names(sort(colMeans(data.good.norm[mock,]),decreasing=T))
temp<-rbind(data.good.norm[mock,mock.order[1:20]],apply(data.good.norm[mock,mock.order[1:20]],2,mean),apply(data.good.norm[mock,mock.order[1:20]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
mock.gg<-melt(temp) #reformats the data (see below)
names(mock.gg)<-c("Sample","OTU","Percentage")
mock.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
mock.gg$OTU <- factor(mock.gg$OTU, levels = levels(mock.gg$OTU), labels = paste(data.good.taxonomy[mock.order[1:20],6]," (",levels(mock.gg$OTU),")",sep=""))


plot.mock<-ggplot(data=mock.gg,aes(y=Percentage,x=OTU,fill=Sample))+ #barplot  #for stacked barplot ggplot(data=mock.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Sample~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Mock Controls Replicates", x="OTUs",y="% Relative Abundance")
plot.mock

#Mock barplots with errorbars
temp.gg<-melt(data.good.norm[mock,mock.order[1:20]])
names(temp.gg)<-c("Sample","OTU","Percentage")
# The next step converts this:
temp.gg[1:10,]
#           Sample     OTU Percentage
# 1  ZymoMock_A_P2 Otu0013   17.93411
# 2  ZymoMock_B_P2 Otu0013   16.82218
# 3  ZymoMock_C_P2 Otu0013   19.13741
# 4  ZymoMock_D_P2 Otu0013   16.32108
# 5  ZymoMock_a_P1 Otu0013   22.37149
# 6  ZymoMock_b_P1 Otu0013   20.70973
# 7  ZymoMock_c_P1 Otu0013   18.51426
# 8  ZymoMock_d_P1 Otu0013   20.36144
# 9  ZymoMock_A_P2 Otu0017   18.33127
# 10 ZymoMock_B_P2 Otu0017   20.04041

temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[mock.order[1:20],6]," (",levels(temp.gg$OTU),")",sep=""))
#... to this:
temp.summary<-summarySE(temp.gg,measurevar="Percentage",groupvars="OTU") # Modified from An R Cookbook
temp.summary
#                                          OTU N Mean.Percentage   geo.mean         sd    geo.sd         se
# 1                   Staphylococcus (Otu0013) 8     19.02146217 18.9254966 2.05114738  1.428729 0.72519011
# 2  Enterobacteriaceae_unclassified (Otu0017) 8     15.43971775 15.2690039 2.53070016  1.686548 0.89473762
# 3                         Bacillus (Otu0040) 8     15.37148804 15.0914823 3.29504585  1.945966 1.16497463
# 4                    Lactobacillus (Otu0056) 8     13.09143931 12.4813618 3.82832243  3.210222 1.35351637
# 5                         Listeria (Otu0057) 8     12.25144315 12.1645599 1.46684123  1.548931 0.51860669
# 6  Enterobacteriaceae_unclassified (Otu0041) 8      9.83748154  9.8128771 0.73268808  1.288362 0.25904435
# 7                     Enterococcus (Otu0042) 8      7.58422025  7.5667352 0.55648900  1.271424 0.19674857
# 8                      Pseudomonas (Otu0067) 8      7.18667866  7.1546612 0.70343015  1.406227 0.24870011
# 9             Bacilli_unclassified (Otu1078) 8      0.11700005  0.4657030 0.14079138 12.710336 0.04977727
# 10         Bacillales_unclassified (Otu1241) 8      0.04087582  0.6247478 0.08207193  5.952982 0.02901681
# 11            Bacilli_unclassified (Otu1296) 8      0.01584571  0.7724586 0.04481843  2.257113 0.01584571
# 12         Bacillales_unclassified (Otu1597) 8      0.01518348  0.7683475 0.04294537  2.257113 0.01518348
# 13            Bacilli_unclassified (Otu1802) 8      0.01414795  0.7615930 0.04001646  2.257113 0.01414795
# 14            Bacilli_unclassified (Otu1797) 8      0.01301612  0.7536964 0.03681514  2.257113 0.01301612
# 15                      Prevotella (Otu0001) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
# 16         Bacillales_unclassified (Otu0002) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
# 17                     Veillonella (Otu0003) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
# 18                      Prevotella (Otu0004) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
# 19    Veillonellaceae_unclassified (Otu0005) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
# 20    Pasteurellaceae_unclassified (Otu0006) 8      0.00000000  0.0000000 0.00000000  1.000000 0.00000000
#            ci
# 1  1.71480213
# 2  2.11571828
# 3  2.75472727
# 4  3.20055764
# 5  1.22630996
# 6  0.61254256
# 7  0.46523645
# 8  0.58808232
# 9  0.11770454
# 10 0.06861385
# 11 0.03746915
# 12 0.03590323
# 13 0.03345459
# 14 0.03077823
# 15 0.00000000
# 16 0.00000000
# 17 0.00000000
# 18 0.00000000
# 19 0.00000000
# 20 0.00000000


#Arithmetic Mean
b<-ggplot(temp.summary,aes(x=OTU,y=Mean.Percentage))+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se),width=0.2)+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Mock Controls Replicates", x="OTUs",y="% Relative Abundance")
b

#Geometric Mean
bb<-ggplot(temp.summary,aes(x=OTU,y=geo.mean))+
	geom_bar(stat="identity")+
	#geom_errorbar(aes(ymin=geo.mean/geo.sd,ymax=geo.mean*geo.sd),width=0.2)+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Geometric Mean of Mock Controls Replicates", x="OTUs",y="% Relative Abundance")
bb

# Isolation Controls All
iso.order<-names(sort(colMeans(data.good.norm[iso,]),decreasing=T))
isoAll.gg<-melt(data.good.norm[iso,iso.order[1:50]])
names(isoAll.gg)<-c("Sample","OTU","Percentage")
# The next step converts this:
isoAll.gg[1:10,]
isoAll.gg$OTU <- factor(isoAll.gg$OTU, levels = levels(isoAll.gg$OTU), labels = paste(data.good.taxonomy[iso.order[1:20],6]," (",levels(isoAll.gg$OTU),")",sep=""))
isoAll.summary<-summarySE(isoAll.gg,measurevar="Percentage",groupvars="OTU") # Modified from An R Cookbook

#Arithmetic Mean
a<-ggplot(isoAll.summary,aes(x=OTU,y=Mean.Percentage))+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se),width=0.2)+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Arithmetic Mean of All Isolation Controls", x="OTUs",y="% Relative Abundance")
a

#Geometric Mean
aa<-ggplot(isoAll.summary,aes(x=OTU,y=geo.mean))+
	geom_bar(stat="identity")+
	#geom_errorbar(aes(ymin=geo.mean/geo.sd,ymax=geo.mean*geo.sd),width=0.2)+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Geometric Mean of All Isolation Controls", x="OTUs",y="% Relative Abundance")
aa

# isoCtrl1
iso1.order<-names(sort(colMeans(data.good.norm[iso[1:3],]),decreasing=T))
temp<-rbind(data.good.norm[iso[1:3],iso1.order[1:50]],apply(data.good.norm[iso[1:3],iso1.order[1:50]],2,mean),apply(data.good.norm[iso[1:3],iso1.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso1.gg<-melt(temp) #reformats the data (see below)
names(iso1.gg)<-c("Name","OTU","Percentage")

#The next step is to set the labels so they display nice in the legend
iso1.gg$OTU <- factor(iso1.gg$OTU, levels = levels(iso1.gg$OTU), labels = paste(data.good.taxonomy[iso1.order[1:50],6]," (",levels(iso1.gg$OTU),")",sep=""))


plot.iso1<-ggplot(data=iso1.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso1.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-1 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso1

# isoCtrl2
iso2.order<-names(sort(colMeans(data.good.norm[iso[4:6],]),decreasing=T))
temp<-rbind(data.good.norm[iso[4:6],iso2.order[1:50]],apply(data.good.norm[iso[4:6],iso2.order[1:50]],2,mean),apply(data.good.norm[iso[4:6],iso2.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso2.gg<-melt(temp) #reformats the data (see below)
names(iso2.gg)<-c("Name","OTU","Percentage")
iso2.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso2.gg$OTU <- factor(iso2.gg$OTU, levels = levels(iso2.gg$OTU), labels = paste(data.good.taxonomy[iso2.order[1:50],6]," (",levels(iso2.gg$OTU),")",sep=""))


plot.iso2<-ggplot(data=iso2.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso2.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-2 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso2

# isoCtrl3
iso3.order<-names(sort(colMeans(data.good.norm[iso[7:9],]),decreasing=T))
temp<-rbind(data.good.norm[iso[7:9],iso3.order[1:50]],apply(data.good.norm[iso[7:9],iso3.order[1:50]],2,mean),apply(data.good.norm[iso[7:9],iso3.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso3.gg<-melt(temp) #reformats the data (see below)
names(iso3.gg)<-c("Name","OTU","Percentage")
iso3.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso3.gg$OTU <- factor(iso3.gg$OTU, levels = levels(iso3.gg$OTU), labels = paste(data.good.taxonomy[iso3.order[1:50],6]," (",levels(iso3.gg$OTU),")",sep=""))


plot.iso3<-ggplot(data=iso3.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso3.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-3 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso3

# isoCtrl4
iso4.order<-names(sort(colMeans(data.good.norm[iso[10:12],]),decreasing=T))
temp<-rbind(data.good.norm[iso[10:12],iso4.order[1:50]],apply(data.good.norm[iso[10:12],iso4.order[1:50]],2,mean),apply(data.good.norm[iso[10:12],iso4.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso4.gg<-melt(temp) #reformats the data (see below)
names(iso4.gg)<-c("Name","OTU","Percentage")
iso4.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso4.gg$OTU <- factor(iso4.gg$OTU, levels = levels(iso4.gg$OTU), labels = paste(data.good.taxonomy[iso4.order[1:50],6]," (",levels(iso4.gg$OTU),")",sep=""))


plot.iso4<-ggplot(data=iso4.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso4.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-4 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso4

# isoCtrl5
iso5.order<-names(sort(colMeans(data.good.norm[iso[13:15],]),decreasing=T))
temp<-rbind(data.good.norm[iso[13:15],iso5.order[1:50]],apply(data.good.norm[iso[13:15],iso5.order[1:50]],2,mean),apply(data.good.norm[iso[13:15],iso5.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso5.gg<-melt(temp) #reformats the data (see below)
names(iso5.gg)<-c("Name","OTU","Percentage")
iso5.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso5.gg$OTU <- factor(iso5.gg$OTU, levels = levels(iso5.gg$OTU), labels = paste(data.good.taxonomy[iso5.order[1:50],6]," (",levels(iso5.gg$OTU),")",sep=""))


plot.iso5<-ggplot(data=iso5.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso5.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-5 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso5

# isoCtrl6
iso6.order<-names(sort(colMeans(data.good.norm[iso[16:18],]),decreasing=T))
temp<-rbind(data.good.norm[iso[16:18],iso6.order[1:50]],apply(data.good.norm[iso[16:18],iso6.order[1:50]],2,mean),apply(data.good.norm[iso[16:18],iso6.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso6.gg<-melt(temp) #reformats the data (see below)
names(iso6.gg)<-c("Name","OTU","Percentage")
iso6.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso6.gg$OTU <- factor(iso6.gg$OTU, levels = levels(iso6.gg$OTU), labels = paste(data.good.taxonomy[iso6.order[1:50],6]," (",levels(iso6.gg$OTU),")",sep=""))


plot.iso6<-ggplot(data=iso6.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso6.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-6 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso6

# isoCtrl7
iso7.order<-names(sort(colMeans(data.good.norm[iso[19:21],]),decreasing=T))
temp<-rbind(data.good.norm[iso[19:21],iso7.order[1:50]],apply(data.good.norm[iso[19:21],iso7.order[1:50]],2,mean),apply(data.good.norm[iso[19:21],iso7.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
iso7.gg<-melt(temp) #reformats the data (see below)
names(iso7.gg)<-c("Name","OTU","Percentage")
iso7.gg[1:10,]

#The next step is to set the labels so they display nice in the legend
iso7.gg$OTU <- factor(iso7.gg$OTU, levels = levels(iso7.gg$OTU), labels = paste(data.good.taxonomy[iso7.order[1:50],6]," (",levels(iso7.gg$OTU),")",sep=""))


plot.iso7<-ggplot(data=iso7.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=iso7.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Isolation Control-7 Replicates", x="OTUs",y="% Relative Abundance")
plot.iso7

# Water
water.order<-names(sort(colMeans(data.good.norm[water,]),decreasing=T))
temp<-rbind(data.good.norm[water,water.order[1:50]],apply(data.good.norm[water,water.order[1:50]],2,mean),apply(data.good.norm[water,water.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
water.gg<-melt(temp) #reformats the data (see below)
names(water.gg)<-c("Name","OTU","Percentage")
water.gg[1:10,]
water.gg$OTU <- factor(water.gg$OTU, levels = levels(water.gg$OTU), labels = paste(data.good.taxonomy[water.order[1:50],6]," (",levels(water.gg$OTU),")",sep=""))


plot.water<-ggplot(data=water.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=water.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of No Template Control Replicates", x="OTUs",y="% Relative Abundance")
plot.water

#Other controls
ae<-grep("AE_",rownames(data.good))
blank<-grep("BLANK|EMPTY",rownames(data.good))

#AE Buffer
ae.order<-names(sort(colMeans(data.good.norm[ae,]),decreasing=T))
temp<-rbind(data.good.norm[ae,ae.order[1:50]],apply(data.good.norm[ae,ae.order[1:50]],2,mean),apply(data.good.norm[ae,ae.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
ae.gg<-melt(temp) #reformats the data (see below)
names(ae.gg)<-c("Name","OTU","Percentage")
ae.gg[1:10,]
ae.gg$OTU <- factor(ae.gg$OTU, levels = levels(ae.gg$OTU), labels = paste(data.good.taxonomy[ae.order[1:50],6]," (",levels(ae.gg$OTU),")",sep=""))


plot.ae<-ggplot(data=ae.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=ae.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of AE Buffer Replicates", x="OTUs",y="% Relative Abundance")
plot.ae

#Blank Wells
blank.order<-names(sort(colMeans(data.good.norm[blank,]),decreasing=T))
temp<-rbind(data.good.norm[blank,blank.order[1:50]],apply(data.good.norm[blank,blank.order[1:50]],2,mean),apply(data.good.norm[blank,blank.order[1:50]],2,geo.mean))
rownames(temp)[c(nrow(temp)-1,nrow(temp))]<-c("Arith.Mean","Geo.Mean")
temp[which(temp==1)]<-0
blank.gg<-melt(temp) #reformats the data (see below)
names(blank.gg)<-c("Name","OTU","Percentage")
blank.gg[1:10,]
blank.gg$OTU <- factor(blank.gg$OTU, levels = levels(blank.gg$OTU), labels = paste(data.good.taxonomy[blank.order[1:50],6]," (",levels(blank.gg$OTU),")",sep=""))


plot.blank<-ggplot(data=blank.gg,aes(y=Percentage,x=OTU,fill=Name))+ #barplot  #for stacked barplot ggplot(data=blank.gg,aes(y=Percentage,x=Name,fill=OTU))+
 	theme_bw()+ #This theme gets rid of the grey background
 	geom_bar(stat="identity", position=position_stack(reverse=T))+ #ggplot2 v 2.2.0 introduced a new position aesthetic which default stacks least to most; this fixes.
 	facet_grid(Name~.)+ # creates a grid of bargraphs
 	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Cross-Comparison of Blank Wells", x="OTUs",y="% Relative Abundance")
plot.blank

### Based on the control results OTU0004, OTU0020, and OTU0008 may be worth eliminating as actual contamination.
### That said they are overall weak so I plan on doing one analysis with them included and repeating it with them excluded.

### Now to think about the sample types.  We have multiple sample types from multiple patients with replicates.  Grouping variables need
### to include Sample, Subject, and Replicate.

### For Sample we have: CBAL, EBC, LBAL, Saline (NC), Nasal, Oral, RBAL, Prewash, Filter-1 (sample),Filter-2 (control), EBC.blank, Nasal.blank, Feces and Feces.blank.
### For Subjects we have: 3057,3116,3172,3183,3194,3219,3264,3301,3356,3378,3389,3482,3493,3518,3552,3585,3600,3611,3644,3769,3840,3851
### for Replicates we have: a,b,and c.

cbal<-grep("CBAL",rownames(data.good))
ebc.ctl<-grep("EBC_Blank|EBC_Control",rownames(data.good))
ebc<-setdiff(grep("EBC",rownames(data.good)),ebc.ctl)
lbal<-grep("LBAL",rownames(data.good))
rbal<-grep("RBAL",rownames(data.good))
saline<-grep("_NC_",rownames(data.good))
nasal.blank<-grep("Nasal_Blank",rownames(data.good))
nasal<-setdiff(grep("Nasal",rownames(data.good)),nasal.blank)
oral<-grep("Oral",rownames(data.good))
prewash<-grep("Pre[Ww]ash",rownames(data.good))
filter.air<-grep("Filter_1_",rownames(data.good))
filter.ctl<-grep("Filter_2_",rownames(data.good))
feces.blank<-grep("Feces_Blank",rownames(data.good))
feces<-setdiff(grep("Feces",rownames(data.good)),feces.blank)

s.3057<-grep("3057",rownames(data.good))
s.3116<-grep("3116",rownames(data.good))
s.3172<-grep("3172",rownames(data.good))
s.3183<-grep("3183",rownames(data.good))
s.3194<-grep("3194",rownames(data.good))
s.3219<-grep("3219",rownames(data.good))
s.3264<-grep("3264",rownames(data.good))
s.3301<-grep("3301",rownames(data.good))
s.3356<-grep("3356",rownames(data.good))
s.3378<-grep("3378",rownames(data.good))
s.3389<-grep("3389",rownames(data.good))
s.3482<-grep("3482",rownames(data.good))
s.3493<-grep("3493",rownames(data.good))
s.3518<-grep("3518",rownames(data.good))
s.3552<-grep("3552",rownames(data.good))
s.3585<-grep("3585",rownames(data.good))
s.3600<-grep("3600",rownames(data.good))
s.3611<-grep("3611",rownames(data.good))
s.3644<-grep("3644",rownames(data.good))
s.3769<-grep("3769",rownames(data.good))
s.3840<-grep("3840",rownames(data.good))
s.3851<-grep("3851",rownames(data.good))



subject<-character()
subject[grep("3057",rownames(data.good))]<-3057
subject[grep("3116",rownames(data.good))]<-3116
subject[grep("3172",rownames(data.good))]<-3172
subject[grep("3183",rownames(data.good))]<-3183
subject[grep("3194",rownames(data.good))]<-3194
subject[grep("3219",rownames(data.good))]<-3219
subject[grep("3264",rownames(data.good))]<-3264
subject[grep("3301",rownames(data.good))]<-3301
subject[grep("3356",rownames(data.good))]<-3356
subject[grep("3378",rownames(data.good))]<-3378
subject[grep("3389",rownames(data.good))]<-3389
subject[grep("3482",rownames(data.good))]<-3482
subject[grep("3493",rownames(data.good))]<-3493
subject[grep("3518",rownames(data.good))]<-3518
subject[grep("3552",rownames(data.good))]<-3552
subject[grep("3585",rownames(data.good))]<-3585
subject[grep("3600",rownames(data.good))]<-3600
subject[grep("3611",rownames(data.good))]<-3611
subject[grep("3644",rownames(data.good))]<-3644
subject[grep("3769",rownames(data.good))]<-3769
subject[grep("3840",rownames(data.good))]<-3840
subject[grep("3851",rownames(data.good))]<-3851
subject[mock]<-"Mock"
subject[c(ae,iso,blank,water)]<-"Control"
subject<-as.factor(subject)

sample<-character()
sample[cbal]<-"CBAL"
sample[lbal]<-"LBAL"
sample[rbal]<-"RBAL"
sample[cbal]<-"CBAL"
sample[ebc.ctl]<-"EBC.Blank"
sample[ebc]<-"EBC"
sample[saline]<-"Saline"
sample[nasal.blank]<-"Nasal.Blank"
sample[nasal]<-"Nasal.Swap"
sample[oral]<-"Oral"
sample[prewash]<-"Pre.Wash"
sample[filter.air]<-"Air.Filter"
sample[filter.ctl]<-"Control.Filter"
sample[feces.blank]<-"Feces.Control"
sample[feces]<-"Feces"
sample[mock]<-"Mock"
sample[iso]<-"Isolation.Controls"
sample[water]<-"No.Template"
sample[ae]<-"AE.Buffer"
sample[blank]<-"Empty.Wells"
sample<-as.factor(sample)

### Now lets see how chaotic this looks in it's entirity
data.good.hel<-decostand(data.good,"hellinger")
data.good.pca<-rda(data.good.hel)
head(summary(data.good.pca))
#0.1943 0.10114
plot(data.good.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (19.4% Explained)",ylab="PC2 (10.1% Explained)",main="Basic PCA of Samples")
points(data.good.pca,pch=19,cex=0.25,display="sites",col=rainbow(20)[as.numeric(sample)])
ordispider(data.good.pca,sample,lty=3,label=T,cex=.5,show.group=sample,col=rainbow(20)[as.numeric(sample)])

###Not surprisingly, it's chaos.  Lets do a bunch of patient sample replicates barplots

#Lets loop it through to create all the graphs we need (Ordered by their group)
i<-1
n<-1
sample.vec<-c("oral","cbal","lbal","rbal","ebc","ebc.ctl","saline","nasal.blank","nasal","prewash","filter.air","filter.ctl","feces","feces.blank")
subject.vec<-c(3057,3116,3172,3183,3194,3219,3264,3301,3356,3378,3389,3482,3493,3518,3552,3585,3600,3611,3644,3769,3840,3851)

for(i in 1:length(sample.vec)){
	temp.sample<-sample.vec[i]
	
	for(n in 1:length(subject.vec)){
		temp.order<-names(sort(colMeans(data.good.norm[get(temp.sample),]),decreasing=T))
		temp.subject<-paste("s.",subject.vec[n],sep="")
		if(length(intersect(get(temp.sample),get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any ",sample.vec[i]," samples",sep="")) 
			next()
				} else if(length(intersect(get(temp.sample),get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 ",sample.vec[i]," sample",sep="")) 
					temp2.order<-names(sort(data.good.norm[intersect(get(temp.sample),get(temp.subject)),],decreasing=T))
					temp.gg<-melt(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp2.order[1:100]])
					temp.gg<-cbind(rownames(temp.gg),temp.gg)
					names(temp.gg)<-c("OTU","Percentage")
					temp.gg$OTU <- factor(temp.gg$OTU, levels = temp.gg$OTU, labels = paste(data.good.taxonomy[temp2.order[1:100],6]," (",temp.gg$OTU,")",sep=""))

					temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU))+ 
 						theme_bw()+ #This theme gets rid of the grey background
 						geom_bar(stat="identity",position=position_stack(reverse=F))+ 
 						scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 						theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
						labs(title=paste(temp.sample," replicates from subject ",temp.subject,sep=""), x="OTUs",y="% Relative Abundance")
					temp.plot
 					ggsave(paste(temp.sample," sample from subject ",subject.vec[n]," barplots.pdf",sep=""),device="pdf",path="/Users/dunard/Dropbox/Adar_ananlysis/Figures/")
					next()
					} else 			
					temp.df<-rbind(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],apply(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],2,mean),apply(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],2,geo.mean))
					rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
					temp.df[which(temp==1)]<-0
					temp.gg<-melt(temp.df) #reformats the data (see below)
					names(temp.gg)<-c("Sample","OTU","Percentage")
					temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[temp.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))


					temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 						theme_bw()+ #This theme gets rid of the grey background
 						geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 						facet_grid(Sample~.)+ # creates a grid of bargraphs
 						scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 						theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
						labs(title=paste(temp.sample," replicates from subject ",temp.subject,sep=""), x="OTUs",y="% Relative Abundance")
					temp.plot
 					ggsave(paste(temp.sample," sample from subject ",subject.vec[n]," barplots.pdf",sep=""),device="pdf",path="/Users/dunard/Dropbox/Adar_ananlysis/Figures/")
			}
		  
		}
	
## Now lets try for Rank ordered
i<-1
n<-1
sample.vec<-c("oral","cbal","lbal","rbal","ebc","ebc.ctl","saline","nasal.blank","nasal","prewash","filter.air","filter.ctl","feces","feces.blank")
subject.vec<-c(3057,3116,3172,3183,3194,3219,3264,3301,3356,3378,3389,3482,3493,3518,3552,3585,3600,3611,3644,3769,3840,3851)

for(i in 1:length(sample.vec)){
	temp.sample<-sample.vec[i]
	
	for(n in 1:length(subject.vec)){
		
		temp.subject<-paste("s.",subject.vec[n],sep="")
		if(length(intersect(get(temp.sample),get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any ",sample.vec[i]," samples",sep="")) 
			next()
				} else if(length(intersect(get(temp.sample),get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 ",sample.vec[i]," sample",sep="")) 
					temp2.order<-names(sort(data.good.norm[intersect(get(temp.sample),get(temp.subject)),],decreasing=T))
					temp.gg<-melt(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp2.order[1:100]])
					temp.gg<-cbind(rownames(temp.gg),temp.gg)
					names(temp.gg)<-c("OTU","Percentage")
					temp.gg$OTU <- factor(temp.gg$OTU, levels = temp.gg$OTU, labels = paste(data.good.taxonomy[temp2.order[1:100],6]," (",temp.gg$OTU,")",sep=""))

					temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU))+ 
 						theme_bw()+ #This theme gets rid of the grey background
 						geom_bar(stat="identity",position=position_stack(reverse=F))+ 
 						scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 						theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
						labs(title=paste("rank abundance ",temp.sample," replicates from subject ",temp.subject,sep=""), x="OTUs",y="% Relative Abundance")
					temp.plot
 					ggsave(paste("Rank abundance ",temp.sample," sample from subject ",subject.vec[n]," barplots.pdf",sep=""),device="pdf",path="~/Dropbox/Adar_ananlysis/Individual Plots with Background/")
					next()
					} else 			
					temp.order<-names(sort(colMeans(data.good.norm[intersect(get(temp.sample),get(temp.subject)),]),decreasing=T))
					temp.df<-rbind(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],apply(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],2,mean),apply(data.good.norm[intersect(get(temp.sample),get(temp.subject)),temp.order[1:100]],2,geo.mean))
					rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
					temp.df[which(temp==1)]<-0
					temp.gg<-melt(temp.df) #reformats the data (see below)
					names(temp.gg)<-c("Sample","OTU","Percentage")
					temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[temp.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))


					temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 						theme_bw()+ #This theme gets rid of the grey background
 						geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 						facet_grid(Sample~.)+ # creates a grid of bargraphs
 						scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 						theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
						labs(title=paste("rank abundance ",temp.sample," replicates from subject ",temp.subject,sep=""), x="OTUs",y="% Relative Abundance")
					temp.plot
 					ggsave(paste("Rank abundance ",temp.sample," sample from subject ",subject.vec[n]," barplots.pdf",sep=""),device="pdf",path="~/Dropbox/Adar_ananlysis/Individual Plots with Background/")
			}
		  
		}
	
# Creating a distance matrix so that when I have the ddPCR data I can look at how BC dist changes depending on [DNA]

data.good.bc<-vegdist(data.good,"bray")
data.good.bc.mat<-as.matrix(data.good.bc)
#data.good.bc.mat[intersect(oral,s.3057),intersect(oral,s.3057)] will get the group comparisons
# I'll come back to this when I have a working numeric keypad.

### Again and again I saw the OTUs that showed up in the controls.  I'm going to create a data set without those.
data.good2<-data.good[-c(mock,water,blank,ae,iso),-c(4,8,20)]
filter.good2<-data.good[c(filter.air,filter.ctl),-c(4,8,20)]

filters<-character()
filters[1:52]<-"Air Filter"
filters[53:102]<-"Filter Control"
filters<-as.factor(filters)
filter.good2.hel<-decostand(filter.good2,"hellinger")
filter.good2.pca<-rda(filter.good2.hel)
head(summary(filter.good2.pca))
#0.1840 0.09167

plot(filter.good2.pca,type="n",font=2,font.lab=2,xlab="PC1 (18.4% Explained)",ylab="PC2 (9.17% Explained)",main = "PCA of Air Filter Samples and Filter Controls")
points(filter.good2.pca,pch=19,cex=0.7,col=as.numeric(filters))
ordispider(filter.good2.pca,filters,show.groups="Air Filter",col="black")
ordispider(filter.good2.pca,filters,show.groups="Filter Control",col="red")
legend("topright",levels(filters),pch=15,col=c("black","red"))
legend("bottomright",c("Adonis p-value = 0.048"))

filter.subjects<-subject[c(filter.air,filter.ctl)]

plot(filter.good2.pca,type="n",font=2,font.lab=2,xlab="PC1 (18.4% Explained)",ylab="PC2 (9.17% Explained)",main = "PCA of Filters With Replicates")
points(filter.good2.pca,pch=19,cex=0.7,col=as.numeric(filters))
ordispider(filter.good2.pca,filter.subjects,label=T)

#Now removing filters from data.good2
data.good2<-data.good[-c(mock,water,blank,ae,iso,filter.air,filter.ctl),-c(4,8,20)]
subject2<-subject[-c(mock,water,blank,ae,iso,filter.air,filter.ctl)]
subject2<-as.factor(as.character(subject2))
sample2<-sample[-c(mock,water,blank,ae,iso,filter.air,filter.ctl)]
sample2<-as.factor(as.character(sample2))
data.good2.norm<-decostand(data.good2,"total")*100

data.good2.hel<-decostand(data.good2,"hellinger")
data.good2.pca<-rda(data.good2.hel)
head(summary(data.good2.pca))
#0.1866 0.03286

sample.cols<-brewer.pal(11,"Spectral")
sample.cols<-c("black",sample.cols)

plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type")
points(data.good2.pca,pch=19,cex=0.7,col=sample.cols[as.numeric(sample2)])
ordispider(data.good2.pca,sample2,show.groups=sample2,col=sample.cols,label=T)
legend("topright",levels(sample2),pch=15,col=sample.cols)

#Now to show off groups
levels(sample2)
#  [1] "CBAL"          "EBC"           "EBC.Blank"     "Feces"         "Feces.Control" "LBAL"         
#  [7] "Nasal.Blank"   "Nasal.Swap"    "Oral"          "Pre.Wash"      "RBAL"          "Saline"  

#CBAL
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:CBAL")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"CBAL",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"CBAL",col="red")
ordispider(data.good2.pca,sample2,show.groups="CBAL",col="red",label=T)
legend("topright",c("CBAL","Others"),pch=15,col=c("red","black"))

#EBC
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:EBC")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"EBC",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"EBC",col="red")
ordispider(data.good2.pca,sample2,show.groups="EBC",col="red",label=T)
legend("topright",c("EBC","Others"),pch=15,col=c("red","black"))

#EBC.Blank
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:EBC.Blank")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"EBC.Blank",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"EBC.Blank",col="red")
ordispider(data.good2.pca,sample2,show.groups="EBC.Blank",col="red",label=T)
legend("topright",c("EBC.Blank","Others"),pch=15,col=c("red","black"))

#Feces
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Feces")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Feces",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Feces",col="red")
ordispider(data.good2.pca,sample2,show.groups="Feces",col="red",label=T)
legend("topright",c("Feces","Others"),pch=15,col=c("red","black"))

#Feces.Control
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Feces.Control")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Feces.Control",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Feces.Control",col="red")
ordispider(data.good2.pca,sample2,show.groups="Feces.Control",col="red",label=T)
legend("topright",c("Feces.Control","Others"),pch=15,col=c("red","black"))

#LBAL
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:LBAL")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"LBAL",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"LBAL",col="red")
ordispider(data.good2.pca,sample2,show.groups="LBAL",col="red",label=T)
legend("topright",c("LBAL","Others"),pch=15,col=c("red","black"))

#Nasal.Blank
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Nasal.Blank")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Nasal.Blank",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Nasal.Blank",col="red")
ordispider(data.good2.pca,sample2,show.groups="Nasal.Blank",col="red",label=T)
legend("topright",c("Nasal.Blank","Others"),pch=15,col=c("red","black"))

#Nasal.Swap
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Nasal.Swap")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Nasal.Swap",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Nasal.Swap",col="red")
ordispider(data.good2.pca,sample2,show.groups="Nasal.Swap",col="red",label=T)
legend("topright",c("Nasal.Swap","Others"),pch=15,col=c("red","black"))

#Oral
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Oral")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Oral",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Oral",col="red")
ordispider(data.good2.pca,sample2,show.groups="Oral",col="red",label=T)
legend("topright",c("Oral","Others"),pch=15,col=c("red","black"))

#Pre.Wash
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Pre.Wash")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Pre.Wash",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Pre.Wash",col="red")
ordispider(data.good2.pca,sample2,show.groups="Pre.Wash",col="red",label=T)
legend("topright",c("Pre.Wash","Others"),pch=15,col=c("red","black"))

#RBAL
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:RBAL")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"RBAL",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"RBAL",col="red")
ordispider(data.good2.pca,sample2,show.groups="RBAL",col="red",label=T)
legend("topright",c("RBAL","Others"),pch=15,col=c("red","black"))

#Saline
plot(data.good2.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Saline")
points(data.good2.pca,pch=19,cex=0.7,select=!sample2%in%"Saline",col="black")
points(data.good2.pca,pch=19,cex=0.7,select=sample2%in%"Saline",col="red")
ordispider(data.good2.pca,sample2,show.groups="Saline",col="red",label=T)
legend("topright",c("Saline","Others"),pch=15,col=c("red","black"))

## Recreating rank ordered plots after removing background
cbal2<-grep("CBAL",rownames(data.good2))
ebc.ctl2<-grep("EBC_Blank|EBC_Control",rownames(data.good2))
ebc2<-setdiff(grep("EBC",rownames(data.good2)),ebc.ctl2)
lbal2<-grep("LBAL",rownames(data.good2))
rbal2<-grep("RBAL",rownames(data.good2))
saline2<-grep("_NC_",rownames(data.good2))
nasal.blank2<-grep("Nasal_Blank",rownames(data.good2))
nasal2<-setdiff(grep("Nasal",rownames(data.good2)),nasal.blank2)
oral2<-grep("Oral",rownames(data.good2))
prewash2<-grep("Pre[Ww]ash",rownames(data.good2))
feces.blank2<-grep("Feces_Blank",rownames(data.good2))
feces2<-setdiff(grep("Feces",rownames(data.good2)),feces.blank2)

s.3057.2<-grep("3057",rownames(data.good2))
s.3116.2<-grep("3116",rownames(data.good2))
s.3172.2<-grep("3172",rownames(data.good2))
s.3183.2<-grep("3183",rownames(data.good2))
s.3194.2<-grep("3194",rownames(data.good2))
s.3219.2<-grep("3219",rownames(data.good2))
s.3264.2<-grep("3264",rownames(data.good2))
s.3301.2<-grep("3301",rownames(data.good2))
s.3356.2<-grep("3356",rownames(data.good2))
s.3378.2<-grep("3378",rownames(data.good2))
s.3389.2<-grep("3389",rownames(data.good2))
s.3482.2<-grep("3482",rownames(data.good2))
s.3493.2<-grep("3493",rownames(data.good2))
s.3518.2<-grep("3518",rownames(data.good2))
s.3552.2<-grep("3552",rownames(data.good2))
s.3585.2<-grep("3585",rownames(data.good2))
s.3600.2<-grep("3600",rownames(data.good2))
s.3611.2<-grep("3611",rownames(data.good2))
s.3644.2<-grep("3644",rownames(data.good2))
s.3769.2<-grep("3769",rownames(data.good2))
s.3840.2<-grep("3840",rownames(data.good2))
s.3851.2<-grep("3851",rownames(data.good2))

i<-1
n<-1
sample.vec<-c("oral2","cbal2","lbal2","rbal2","ebc2","ebc.ctl2","saline2","nasal.blank2","nasal2","prewash2","feces2","feces.blank2")
subject.vec<-c(3057,3116,3172,3183,3194,3219,3264,3301,3356,3378,3389,3482,3493,3518,3552,3585,3600,3611,3644,3769,3840,3851)


# Keyboard is here so lets tackle the distance matrix
data.good2.bc<-vegdist(data.good2,"bray")
data.good2.bc.mat<-as.matrix(data.good2.bc)	
data.good2.bc.mat[intersect(oral2,s.3057.2),intersect(oral2,s.3057.2)]

##Oral
data.good2.bc.mat[intersect(oral2,s.3057.2),intersect(oral2,s.3057.2)]
#                     Adar_3057_Oral_a_P1 Adar_3057_Oral_b_P1 Adar_3057_Oral_c_P1
# Adar_3057_Oral_a_P1          0.00000000          0.05550054          0.11467789
# Adar_3057_Oral_b_P1          0.05550054          0.00000000          0.06613876
# Adar_3057_Oral_c_P1          0.11467789          0.06613876          0.00000000
data.good2.bc.mat[intersect(oral2,s.3116.2),intersect(oral2,s.3116.2)]
#                     Adar_3116_Oral_a_P1 Adar_3116_Oral_b_P1 Adar_3116_Oral_c_P1
# Adar_3116_Oral_a_P1           0.0000000           0.1166815           0.2303657
# Adar_3116_Oral_b_P1           0.1166815           0.0000000           0.1252564
# Adar_3116_Oral_c_P1           0.2303657           0.1252564           0.0000000
data.good2.bc.mat[intersect(oral2,s.3172.2),intersect(oral2,s.3172.2)]
#                     Adar_3172_Oral_a_P1 Adar_3172_Oral_b_P1 Adar_3172_Oral_c_P1
# Adar_3172_Oral_a_P1           0.0000000           0.1087312           0.1200913
# Adar_3172_Oral_b_P1           0.1087312           0.0000000           0.0402582
# Adar_3172_Oral_c_P1           0.1200913           0.0402582           0.0000000
data.good2.bc.mat[intersect(oral2,s.3183.2),intersect(oral2,s.3183.2)]
#                     Adar_3183_Oral_a_P1 Adar_3183_Oral_b_P1 Adar_3183_Oral_c_P1
# Adar_3183_Oral_a_P1          0.00000000          0.09109200          0.05853573
# Adar_3183_Oral_b_P1          0.09109200          0.00000000          0.04715937
# Adar_3183_Oral_c_P1          0.05853573          0.04715937          0.00000000
data.good2.bc.mat[intersect(oral2,s.3194.2),intersect(oral2,s.3194.2)]
#                     Adar_3194_Oral_a_P1 Adar_3194_Oral_b_P1 Adar_3194_Oral_c_P1
# Adar_3194_Oral_a_P1           0.0000000          0.12131591          0.11198881
# Adar_3194_Oral_b_P1           0.1213159          0.00000000          0.02316507
# Adar_3194_Oral_c_P1           0.1119888          0.02316507          0.00000000
data.good2.bc.mat[intersect(oral2,s.3219.2),intersect(oral2,s.3219.2)]
#                     Adar_3219_Oral_a_P1 Adar_3219_Oral_b_P1 Adar_3219_Oral_c_P1
# Adar_3219_Oral_a_P1          0.00000000          0.04870848          0.04779412
# Adar_3219_Oral_b_P1          0.04870848          0.00000000          0.04403023
# Adar_3219_Oral_c_P1          0.04779412          0.04403023          0.00000000
data.good2.bc.mat[intersect(oral2,s.3264.2),intersect(oral2,s.3264.2)]
#                     Adar_3264_Oral_a_P1 Adar_3264_Oral_b_P1 Adar_3264_Oral_c_P1
# Adar_3264_Oral_a_P1          0.00000000          0.04142437          0.07058987
# Adar_3264_Oral_b_P1          0.04142437          0.00000000          0.04127131
# Adar_3264_Oral_c_P1          0.07058987          0.04127131          0.00000000
data.good2.bc.mat[intersect(oral2,s.3301.2),intersect(oral2,s.3301.2)]
#                     Adar_3301_Oral_a_P1 Adar_3301_Oral_b_P1 Adar_3301_Oral_c_P1
# Adar_3301_Oral_a_P1                   0                   1                   1
# Adar_3301_Oral_b_P1                   1                   0                   1
# Adar_3301_Oral_c_P1                   1                   1                   0
data.good2.bc.mat[intersect(oral2,s.3356.2),intersect(oral2,s.3356.2)]
#                     Adar_3356_Oral_a_P1 Adar_3356_Oral_b_P1 Adar_3356_Oral_c_P1
# Adar_3356_Oral_a_P1          0.00000000          0.03729151          0.05378115
# Adar_3356_Oral_b_P1          0.03729151          0.00000000          0.04209037
# Adar_3356_Oral_c_P1          0.05378115          0.04209037          0.00000000
data.good2.bc.mat[intersect(oral2,s.3378.2),intersect(oral2,s.3378.2)]
#                     Adar_3378_Oral_a_P1 Adar_3378_Oral_b_P1 Adar_3378_Oral_c_P2
# Adar_3378_Oral_a_P1          0.00000000          0.04707016           0.1338270
# Adar_3378_Oral_b_P1          0.04707016          0.00000000           0.1003222
# Adar_3378_Oral_c_P2          0.13382701          0.10032223           0.0000000
data.good2.bc.mat[intersect(oral2,s.3389.2),intersect(oral2,s.3389.2)]
#                     Adar_3389_Oral_a_P2 Adar_3389_Oral_b_P2 Adar_3389_Oral_c_P2
# Adar_3389_Oral_a_P2          0.00000000          0.05450060          0.08267663
# Adar_3389_Oral_b_P2          0.05450060          0.00000000          0.08183725
# Adar_3389_Oral_c_P2          0.08267663          0.08183725          0.00000000
data.good2.bc.mat[intersect(oral2,s.3482.2),intersect(oral2,s.3482.2)]
#                     Adar_3482_Oral_a_P2 Adar_3482_Oral_b_P2 Adar_3482_Oral_c_P2
# Adar_3482_Oral_a_P2           0.0000000                   1           0.8787627
# Adar_3482_Oral_b_P2           1.0000000                   0           1.0000000
# Adar_3482_Oral_c_P2           0.8787627                   1           0.0000000
data.good2.bc.mat[intersect(oral2,s.3493.2),intersect(oral2,s.3493.2)]
#                     Adar_3493_Oral_a_P2 Adar_3493_Oral_b_P2 Adar_3493_Oral_c_P2
# Adar_3493_Oral_a_P2           0.0000000           0.1599834           0.1014600
# Adar_3493_Oral_b_P2           0.1599834           0.0000000           0.1182792
# Adar_3493_Oral_c_P2           0.1014600           0.1182792           0.0000000
data.good2.bc.mat[intersect(oral2,s.3518.2),intersect(oral2,s.3518.2)]
#                     Adar_3518_Oral_a_P2 Adar_3518_Oral_b_P2 Adar_3518_Oral_c_P2
# Adar_3518_Oral_a_P2           0.0000000           0.1178867           0.2387174
# Adar_3518_Oral_b_P2           0.1178867           0.0000000           0.1390414
# Adar_3518_Oral_c_P2           0.2387174           0.1390414           0.0000000
data.good2.bc.mat[intersect(oral2,s.3552.2),intersect(oral2,s.3552.2)]
#                     Adar_3552_Oral_a_P2 Adar_3552_Oral_b_P2 Adar_3552_Oral_c_P2
# Adar_3552_Oral_a_P2           0.0000000           0.2010405           0.3320246
# Adar_3552_Oral_b_P2           0.2010405           0.0000000           0.1467522
# Adar_3552_Oral_c_P2           0.3320246           0.1467522           0.0000000
data.good2.bc.mat[intersect(oral2,s.3585.2),intersect(oral2,s.3585.2)]
#                     Adar_3585_Oral_a_P2 Adar_3585_Oral_b_P2 Adar_3585_Oral_c_P2
# Adar_3585_Oral_a_P2           0.0000000          0.15980338          0.10349294
# Adar_3585_Oral_b_P2           0.1598034          0.00000000          0.07822056
# Adar_3585_Oral_c_P2           0.1034929          0.07822056          0.00000000
data.good2.bc.mat[intersect(oral2,s.3600.2),intersect(oral2,s.3600.2)]
#                     Adar_3600_Oral_a_P2 Adar_3600_Oral_b_P2 Adar_3600_Oral_c_P2
# Adar_3600_Oral_a_P2           0.0000000           0.2469319           0.1207960
# Adar_3600_Oral_b_P2           0.2469319           0.0000000           0.2800896
# Adar_3600_Oral_c_P2           0.1207960           0.2800896           0.0000000
data.good2.bc.mat[intersect(oral2,s.3611.2),intersect(oral2,s.3611.2)]
#                     Adar_3611_Oral_a_P2 Adar_3611_Oral_b_P2 Adar_3611_Oral_c_P2
# Adar_3611_Oral_a_P2           0.0000000           0.1150711           0.1141943
# Adar_3611_Oral_b_P2           0.1150711           0.0000000           0.1004470
# Adar_3611_Oral_c_P2           0.1141943           0.1004470           0.0000000
data.good2.bc.mat[intersect(oral2,s.3644.2),intersect(oral2,s.3644.2)]
#                     Adar_3644_Oral_a_P2 Adar_3644_Oral_b_P2 Adar_3644_Oral_c_P2
# Adar_3644_Oral_a_P2          0.00000000          0.09375936           0.2610470
# Adar_3644_Oral_b_P2          0.09375936          0.00000000           0.2653187
# Adar_3644_Oral_c_P2          0.26104701          0.26531866           0.0000000
data.good2.bc.mat[intersect(oral2,s.3769.2),intersect(oral2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(oral2,s.3840.2),intersect(oral2,s.3840.2)]
#                     Adar_3840_Oral_a_P2 Adar_3840_Oral_b_P2 Adar_3840_Oral_c_P2
# Adar_3840_Oral_a_P2          0.00000000           0.1142393          0.06541703
# Adar_3840_Oral_b_P2          0.11423935           0.0000000          0.13518882
# Adar_3840_Oral_c_P2          0.06541703           0.1351888          0.00000000
data.good2.bc.mat[intersect(oral2,s.3851.2),intersect(oral2,s.3851.2)]
#                     Adar_3851_Oral_a_P2 Adar_3851_Oral_b_P2 Adar_3851_Oral_c_P2
# Adar_3851_Oral_a_P2           0.0000000           0.0760316           0.1626726
# Adar_3851_Oral_b_P2           0.0760316           0.0000000           0.2116377
# Adar_3851_Oral_c_P2           0.1626726           0.2116377           0.0000000

##CBAL
data.good2.bc.mat[intersect(cbal2,s.3057.2),intersect(cbal2,s.3057.2)]
#                       Adar_3057_CBAL_1_a_P1 Adar_3057_CBAL_1_b_P1 Adar_3057_CBAL_1_c_P1
# Adar_3057_CBAL_1_a_P1             0.0000000             0.8064621             0.6407422
# Adar_3057_CBAL_1_b_P1             0.8064621             0.0000000             0.7913482
# Adar_3057_CBAL_1_c_P1             0.6407422             0.7913482             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3116.2),intersect(cbal2,s.3116.2)]
#                       Adar_3116_CBAL_1_a_P1 Adar_3116_CBAL_1_b_P1 Adar_3116_CBAL_1_c_P1 Adar_3116_CBAL_2_a_P1
# Adar_3116_CBAL_1_a_P1             0.0000000              0.994618             1.0000000             1.0000000
# Adar_3116_CBAL_1_b_P1             0.9946180              0.000000             1.0000000             1.0000000
# Adar_3116_CBAL_1_c_P1             1.0000000              1.000000             0.0000000             0.9980961
# Adar_3116_CBAL_2_a_P1             1.0000000              1.000000             0.9980961             0.0000000
# Adar_3116_CBAL_2_b_P1             1.0000000              1.000000             1.0000000             0.9949401
# Adar_3116_CBAL_2_c_P1             0.9875449              1.000000             0.9760565             0.9458492
#                       Adar_3116_CBAL_2_b_P1 Adar_3116_CBAL_2_c_P1
# Adar_3116_CBAL_1_a_P1             1.0000000             0.9875449
# Adar_3116_CBAL_1_b_P1             1.0000000             1.0000000
# Adar_3116_CBAL_1_c_P1             1.0000000             0.9760565
# Adar_3116_CBAL_2_a_P1             0.9949401             0.9458492
# Adar_3116_CBAL_2_b_P1             0.0000000             1.0000000
# Adar_3116_CBAL_2_c_P1             1.0000000             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3172.2),intersect(cbal2,s.3172.2)]
#                       Adar_3172_CBAL_1_a_P1 Adar_3172_CBAL_1_b_P1 Adar_3172_CBAL_1_c_P1 Adar_3172_CBAL_2_a_P1
# Adar_3172_CBAL_1_a_P1             0.0000000             0.3790253             0.2781217             0.2591047
# Adar_3172_CBAL_1_b_P1             0.3790253             0.0000000             0.2850024             0.4826603
# Adar_3172_CBAL_1_c_P1             0.2781217             0.2850024             0.0000000             0.3613434
# Adar_3172_CBAL_2_a_P1             0.2591047             0.4826603             0.3613434             0.0000000
# Adar_3172_CBAL_2_b_P1             0.3996787             0.4244927             0.3913022             0.3347381
# Adar_3172_CBAL_2_c_P1             0.2131140             0.3939941             0.3048657             0.2616203
#                       Adar_3172_CBAL_2_b_P1 Adar_3172_CBAL_2_c_P1
# Adar_3172_CBAL_1_a_P1             0.3996787             0.2131140
# Adar_3172_CBAL_1_b_P1             0.4244927             0.3939941
# Adar_3172_CBAL_1_c_P1             0.3913022             0.3048657
# Adar_3172_CBAL_2_a_P1             0.3347381             0.2616203
# Adar_3172_CBAL_2_b_P1             0.0000000             0.3156301
# Adar_3172_CBAL_2_c_P1             0.3156301             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3183.2),intersect(cbal2,s.3183.2)]
#                       Adar_3183_CBAL_1_a_P1 Adar_3183_CBAL_1_b_P1 Adar_3183_CBAL_1_c_P1 Adar_3183_CBAL_2_a_P1
# Adar_3183_CBAL_1_a_P1             0.0000000             0.1611712             0.1744301             0.2518341
# Adar_3183_CBAL_1_b_P1             0.1611712             0.0000000             0.1263878             0.1604086
# Adar_3183_CBAL_1_c_P1             0.1744301             0.1263878             0.0000000             0.1681120
# Adar_3183_CBAL_2_a_P1             0.2518341             0.1604086             0.1681120             0.0000000
# Adar_3183_CBAL_2_b_P1             0.2035496             0.1797169             0.1753541             0.1775426
# Adar_3183_CBAL_2_c_P1             0.2265098             0.2546781             0.2593109             0.2582411
#                       Adar_3183_CBAL_2_b_P1 Adar_3183_CBAL_2_c_P1
# Adar_3183_CBAL_1_a_P1             0.2035496             0.2265098
# Adar_3183_CBAL_1_b_P1             0.1797169             0.2546781
# Adar_3183_CBAL_1_c_P1             0.1753541             0.2593109
# Adar_3183_CBAL_2_a_P1             0.1775426             0.2582411
# Adar_3183_CBAL_2_b_P1             0.0000000             0.1590394
# Adar_3183_CBAL_2_c_P1             0.1590394             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3194.2),intersect(cbal2,s.3194.2)]
#                       Adar_3194_CBAL_1_a_P1 Adar_3194_CBAL_1_b_P1 Adar_3194_CBAL_1_c_P1 Adar_3194_CBAL_2_a_P1
# Adar_3194_CBAL_1_a_P1             0.0000000             0.1478271             0.1152035             0.2030179
# Adar_3194_CBAL_1_b_P1             0.1478271             0.0000000             0.1369994             0.2454901
# Adar_3194_CBAL_1_c_P1             0.1152035             0.1369994             0.0000000             0.2063373
# Adar_3194_CBAL_2_a_P1             0.2030179             0.2454901             0.2063373             0.0000000
# Adar_3194_CBAL_2_b_P1             0.1447417             0.2068187             0.1760953             0.1738322
# Adar_3194_CBAL_2_c_P1             0.1538473             0.1842896             0.1744247             0.1962483
#                       Adar_3194_CBAL_2_b_P1 Adar_3194_CBAL_2_c_P1
# Adar_3194_CBAL_1_a_P1             0.1447417             0.1538473
# Adar_3194_CBAL_1_b_P1             0.2068187             0.1842896
# Adar_3194_CBAL_1_c_P1             0.1760953             0.1744247
# Adar_3194_CBAL_2_a_P1             0.1738322             0.1962483
# Adar_3194_CBAL_2_b_P1             0.0000000             0.1523714
# Adar_3194_CBAL_2_c_P1             0.1523714             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3219.2),intersect(cbal2,s.3219.2)]
#                       Adar_3219_CBAL_1_a_P1 Adar_3219_CBAL_1_b_P1 Adar_3219_CBAL_1_c_P1 Adar_3219_CBAL_2_a_P1
# Adar_3219_CBAL_1_a_P1             0.0000000             0.9963805             1.0000000             0.9811923
# Adar_3219_CBAL_1_b_P1             0.9963805             0.0000000             0.9555110             1.0000000
# Adar_3219_CBAL_1_c_P1             1.0000000             0.9555110             0.0000000             1.0000000
# Adar_3219_CBAL_2_a_P1             0.9811923             1.0000000             1.0000000             0.0000000
# Adar_3219_CBAL_2_b_P1             0.9949996             0.9701429             0.9260808             0.9077918
# Adar_3219_CBAL_2_c_P1             0.9923644             0.9973859             1.0000000             0.9838024
#                       Adar_3219_CBAL_2_b_P1 Adar_3219_CBAL_2_c_P1
# Adar_3219_CBAL_1_a_P1             0.9949996             0.9923644
# Adar_3219_CBAL_1_b_P1             0.9701429             0.9973859
# Adar_3219_CBAL_1_c_P1             0.9260808             1.0000000
# Adar_3219_CBAL_2_a_P1             0.9077918             0.9838024
# Adar_3219_CBAL_2_b_P1             0.0000000             0.9910798
# Adar_3219_CBAL_2_c_P1             0.9910798             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3264.2),intersect(cbal2,s.3264.2)]
#                       Adar_3264_CBAL_1_a_P1 Adar_3264_CBAL_1_b_P1 Adar_3264_CBAL_1_c_P1 Adar_3264_CBAL_2_a_P1
# Adar_3264_CBAL_1_a_P1             0.0000000             0.4036091             0.5989357             0.9278277
# Adar_3264_CBAL_1_b_P1             0.4036091             0.0000000             0.6282742             0.9480495
# Adar_3264_CBAL_1_c_P1             0.5989357             0.6282742             0.0000000             0.9202648
# Adar_3264_CBAL_2_a_P1             0.9278277             0.9480495             0.9202648             0.0000000
# Adar_3264_CBAL_2_b_P1             0.8136610             0.8148945             0.7509174             0.8368637
# Adar_3264_CBAL_2_c_P1             0.9573085             0.9861788             0.8630421             0.9124011
#                       Adar_3264_CBAL_2_b_P1 Adar_3264_CBAL_2_c_P1
# Adar_3264_CBAL_1_a_P1             0.8136610             0.9573085
# Adar_3264_CBAL_1_b_P1             0.8148945             0.9861788
# Adar_3264_CBAL_1_c_P1             0.7509174             0.8630421
# Adar_3264_CBAL_2_a_P1             0.8368637             0.9124011
# Adar_3264_CBAL_2_b_P1             0.0000000             0.8197970
# Adar_3264_CBAL_2_c_P1             0.8197970             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3301.2),intersect(cbal2,s.3301.2)]
#                       Adar_3301_CBAL_1_a_P1 Adar_3301_CBAL_1_b_P1 Adar_3301_CBAL_1_c_P1 Adar_3301_CBAL_2_a_P1
# Adar_3301_CBAL_1_a_P1             0.0000000             0.3810463             0.3341792             0.4224708
# Adar_3301_CBAL_1_b_P1             0.3810463             0.0000000             0.3176016             0.3309989
# Adar_3301_CBAL_1_c_P1             0.3341792             0.3176016             0.0000000             0.3902727
# Adar_3301_CBAL_2_a_P1             0.4224708             0.3309989             0.3902727             0.0000000
# Adar_3301_CBAL_2_b_P1             0.4788310             0.3946747             0.4149891             0.1859369
# Adar_3301_CBAL_2_c_P1             0.3914667             0.2981180             0.3753274             0.1868852
#                       Adar_3301_CBAL_2_b_P1 Adar_3301_CBAL_2_c_P1
# Adar_3301_CBAL_1_a_P1             0.4788310             0.3914667
# Adar_3301_CBAL_1_b_P1             0.3946747             0.2981180
# Adar_3301_CBAL_1_c_P1             0.4149891             0.3753274
# Adar_3301_CBAL_2_a_P1             0.1859369             0.1868852
# Adar_3301_CBAL_2_b_P1             0.0000000             0.2106058
# Adar_3301_CBAL_2_c_P1             0.2106058             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3356.2),intersect(cbal2,s.3356.2)]
#                       Adar_3356_CBAL_1_a_P1 Adar_3356_CBAL_1_b_P1 Adar_3356_CBAL_1_c_P1 Adar_3356_CBAL_2_a_P1
# Adar_3356_CBAL_1_a_P1             0.0000000             0.3237813             0.2517844             0.2615326
# Adar_3356_CBAL_1_b_P1             0.3237813             0.0000000             0.2130616             0.3256956
# Adar_3356_CBAL_1_c_P1             0.2517844             0.2130616             0.0000000             0.3212422
# Adar_3356_CBAL_2_a_P1             0.2615326             0.3256956             0.3212422             0.0000000
# Adar_3356_CBAL_2_b_P1             0.3160269             0.3174930             0.3370067             0.1695520
# Adar_3356_CBAL_2_c_P1             0.3150886             0.3231923             0.3003140             0.2827488
#                       Adar_3356_CBAL_2_b_P1 Adar_3356_CBAL_2_c_P1
# Adar_3356_CBAL_1_a_P1             0.3160269             0.3150886
# Adar_3356_CBAL_1_b_P1             0.3174930             0.3231923
# Adar_3356_CBAL_1_c_P1             0.3370067             0.3003140
# Adar_3356_CBAL_2_a_P1             0.1695520             0.2827488
# Adar_3356_CBAL_2_b_P1             0.0000000             0.2895290
# Adar_3356_CBAL_2_c_P1             0.2895290             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3378.2),intersect(cbal2,s.3378.2)]
#                       Adar_3378_CBAL_1_a_P1 Adar_3378_CBAL_1_b_P1 Adar_3378_CBAL_1_c_P1 Adar_3378_CBAL_2_a_P1
# Adar_3378_CBAL_1_a_P1             0.0000000             0.3419265             0.5929328             0.6131082
# Adar_3378_CBAL_1_b_P1             0.3419265             0.0000000             0.5855292             0.6305131
# Adar_3378_CBAL_1_c_P1             0.5929328             0.5855292             0.0000000             0.8304928
# Adar_3378_CBAL_2_a_P1             0.6131082             0.6305131             0.8304928             0.0000000
# Adar_3378_CBAL_2_b_P1             0.4681813             0.4847438             0.6910051             0.4660697
# Adar_3378_CBAL_2_c_P1             0.5253748             0.5598291             0.7610995             0.4964489
#                       Adar_3378_CBAL_2_b_P1 Adar_3378_CBAL_2_c_P1
# Adar_3378_CBAL_1_a_P1             0.4681813             0.5253748
# Adar_3378_CBAL_1_b_P1             0.4847438             0.5598291
# Adar_3378_CBAL_1_c_P1             0.6910051             0.7610995
# Adar_3378_CBAL_2_a_P1             0.4660697             0.4964489
# Adar_3378_CBAL_2_b_P1             0.0000000             0.3210180
# Adar_3378_CBAL_2_c_P1             0.3210180             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3389.2),intersect(cbal2,s.3389.2)]
#                       Adar_3389_CBAL_1_a_P2 Adar_3389_CBAL_1_b_P2 Adar_3389_CBAL_1_c_P2
# Adar_3389_CBAL_1_a_P2             0.0000000             0.9317523             0.9261380
# Adar_3389_CBAL_1_b_P2             0.9317523             0.0000000             0.9719414
# Adar_3389_CBAL_1_c_P2             0.9261380             0.9719414             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3482.2),intersect(cbal2,s.3482.2)]
#                       Adar_3482_CBAL_1_a_P2 Adar_3482_CBAL_1_b_P2 Adar_3482_CBAL_1_c_P2
# Adar_3482_CBAL_1_a_P2             0.0000000             0.8267712             0.7246235
# Adar_3482_CBAL_1_b_P2             0.8267712             0.0000000             0.7783676
# Adar_3482_CBAL_1_c_P2             0.7246235             0.7783676             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3493.2),intersect(cbal2,s.3493.2)]
#                       Adar_3493_CBAL_1_a_P2 Adar_3493_CBAL_1_b_P2 Adar_3493_CBAL_1_c_P2
# Adar_3493_CBAL_1_a_P2             0.0000000             0.3013424             0.1960322
# Adar_3493_CBAL_1_b_P2             0.3013424             0.0000000             0.3235684
# Adar_3493_CBAL_1_c_P2             0.1960322             0.3235684             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3518.2),intersect(cbal2,s.3518.2)]
#                       Adar_3518_CBAL_1_a_P2 Adar_3518_CBAL_1_b_P2 Adar_3518_CBAL_1_c_P2
# Adar_3518_CBAL_1_a_P2             0.0000000             0.1390393             0.1051614
# Adar_3518_CBAL_1_b_P2             0.1390393             0.0000000             0.1049056
# Adar_3518_CBAL_1_c_P2             0.1051614             0.1049056             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3552.2),intersect(cbal2,s.3552.2)]
#                       Adar_3552_CBAL_1_a_P2 Adar_3552_CBAL_1_b_P2 Adar_3552_CBAL_1_c_P2 Adar_3552_CBAL_2_a_P2
# Adar_3552_CBAL_1_a_P2             0.0000000             0.1672821             0.1956788             0.2275939
# Adar_3552_CBAL_1_b_P2             0.1672821             0.0000000             0.2182796             0.2315574
# Adar_3552_CBAL_1_c_P2             0.1956788             0.2182796             0.0000000             0.2332461
# Adar_3552_CBAL_2_a_P2             0.2275939             0.2315574             0.2332461             0.0000000
# Adar_3552_CBAL_2_b_P2             0.2007065             0.2221753             0.2240620             0.1278313
# Adar_3552_CBAL_2_c_P2             0.1951886             0.2121120             0.2007816             0.1647920
#                       Adar_3552_CBAL_2_b_P2 Adar_3552_CBAL_2_c_P2
# Adar_3552_CBAL_1_a_P2             0.2007065             0.1951886
# Adar_3552_CBAL_1_b_P2             0.2221753             0.2121120
# Adar_3552_CBAL_1_c_P2             0.2240620             0.2007816
# Adar_3552_CBAL_2_a_P2             0.1278313             0.1647920
# Adar_3552_CBAL_2_b_P2             0.0000000             0.1534376
# Adar_3552_CBAL_2_c_P2             0.1534376             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3585.2),intersect(cbal2,s.3585.2)]
#                       Adar_3585_CBAL_1_a_P2 Adar_3585_CBAL_1_b_P2 Adar_3585_CBAL_1_c_P2 Adar_3585_CBAL_2_a_P2
# Adar_3585_CBAL_1_a_P2             0.0000000             0.4269771             0.4613610             0.4974292
# Adar_3585_CBAL_1_b_P2             0.4269771             0.0000000             0.4209146             0.6007391
# Adar_3585_CBAL_1_c_P2             0.4613610             0.4209146             0.0000000             0.6110175
# Adar_3585_CBAL_2_a_P2             0.4974292             0.6007391             0.6110175             0.0000000
# Adar_3585_CBAL_2_b_P2             0.5385666             0.6574864             0.6867105             0.3503161
# Adar_3585_CBAL_2_c_P2             0.6074061             0.6418348             0.6789347             0.4602150
#                       Adar_3585_CBAL_2_b_P2 Adar_3585_CBAL_2_c_P2
# Adar_3585_CBAL_1_a_P2             0.5385666             0.6074061
# Adar_3585_CBAL_1_b_P2             0.6574864             0.6418348
# Adar_3585_CBAL_1_c_P2             0.6867105             0.6789347
# Adar_3585_CBAL_2_a_P2             0.3503161             0.4602150
# Adar_3585_CBAL_2_b_P2             0.0000000             0.4621429
# Adar_3585_CBAL_2_c_P2             0.4621429             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3600.2),intersect(cbal2,s.3600.2)]
#                       Adar_3600_CBAL_1_a_P2 Adar_3600_CBAL_1_b_P2 Adar_3600_CBAL_1_c_P2
# Adar_3600_CBAL_1_a_P2             0.0000000             0.9990254             0.9992685
# Adar_3600_CBAL_1_b_P2             0.9990254             0.0000000             0.1658428
# Adar_3600_CBAL_1_c_P2             0.9992685             0.1658428             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3611.2),intersect(cbal2,s.3611.2)]
#                       Adar_3611_CBAL_1_a_P2 Adar_3611_CBAL_1_b_P2 Adar_3611_CBAL_1_c_P2 Adar_3611_CBAL_2_a_P2
# Adar_3611_CBAL_1_a_P2             0.0000000             0.3597794            0.11732789            0.13317955
# Adar_3611_CBAL_1_b_P2             0.3597794             0.0000000            0.43912449            0.44477988
# Adar_3611_CBAL_1_c_P2             0.1173279             0.4391245            0.00000000            0.08229879
# Adar_3611_CBAL_2_a_P2             0.1331795             0.4447799            0.08229879            0.00000000
# Adar_3611_CBAL_2_b_P2             0.1046182             0.3484837            0.18174125            0.15807293
# Adar_3611_CBAL_2_c_P2             0.3145462             0.5817880            0.23009143            0.24133783
#                       Adar_3611_CBAL_2_b_P2 Adar_3611_CBAL_2_c_P2
# Adar_3611_CBAL_1_a_P2             0.1046182             0.3145462
# Adar_3611_CBAL_1_b_P2             0.3484837             0.5817880
# Adar_3611_CBAL_1_c_P2             0.1817412             0.2300914
# Adar_3611_CBAL_2_a_P2             0.1580729             0.2413378
# Adar_3611_CBAL_2_b_P2             0.0000000             0.3207887
# Adar_3611_CBAL_2_c_P2             0.3207887             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3644.2),intersect(cbal2,s.3644.2)]
#                       Adar_3644_CBAL_1_a_P2 Adar_3644_CBAL_1_b_P2 Adar_3644_CBAL_1_c_P2 Adar_3644_CBAL_2_a_P2
# Adar_3644_CBAL_1_a_P2             0.0000000             0.1032269             0.1117287             0.2408272
# Adar_3644_CBAL_1_b_P2             0.1032269             0.0000000             0.1142744             0.2857423
# Adar_3644_CBAL_1_c_P2             0.1117287             0.1142744             0.0000000             0.2418487
# Adar_3644_CBAL_2_a_P2             0.2408272             0.2857423             0.2418487             0.0000000
# Adar_3644_CBAL_2_b_P2             0.3745650             0.4142787             0.4101090             0.2195013
# Adar_3644_CBAL_2_c_P2             0.2944322             0.3396929             0.3124991             0.1199112
#                       Adar_3644_CBAL_2_b_P2 Adar_3644_CBAL_2_c_P2
# Adar_3644_CBAL_1_a_P2             0.3745650             0.2944322
# Adar_3644_CBAL_1_b_P2             0.4142787             0.3396929
# Adar_3644_CBAL_1_c_P2             0.4101090             0.3124991
# Adar_3644_CBAL_2_a_P2             0.2195013             0.1199112
# Adar_3644_CBAL_2_b_P2             0.0000000             0.1468187
# Adar_3644_CBAL_2_c_P2             0.1468187             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3769.2),intersect(cbal2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(cbal2,s.3840.2),intersect(cbal2,s.3840.2)]
#                       Adar_3840_CBAL_1_a_P2 Adar_3840_CBAL_1_b_P2 Adar_3840_CBAL_1_c_P2
# Adar_3840_CBAL_1_a_P2             0.0000000             0.4651044             0.3791271
# Adar_3840_CBAL_1_b_P2             0.4651044             0.0000000             0.2690041
# Adar_3840_CBAL_1_c_P2             0.3791271             0.2690041             0.0000000

data.good2.bc.mat[intersect(cbal2,s.3851.2),intersect(cbal2,s.3851.2)]
#                       Adar_3851_CBAL_1_a_P2 Adar_3851_CBAL_1_b_P2 Adar_3851_CBAL_1_c_P2 Adar_3851_CBAL_2_a_P2
# Adar_3851_CBAL_1_a_P2             0.0000000             0.6746955             0.7068358             0.9672109
# Adar_3851_CBAL_1_b_P2             0.6746955             0.0000000             0.6364131             0.9574018
# Adar_3851_CBAL_1_c_P2             0.7068358             0.6364131             0.0000000             0.9658162
# Adar_3851_CBAL_2_a_P2             0.9672109             0.9574018             0.9658162             0.0000000
# Adar_3851_CBAL_2_b_P2             0.9973189             0.9956844             1.0000000             0.2484471
# Adar_3851_CBAL_2_c_P2             0.9884591             0.9836571             0.9939702             0.2126547
#                       Adar_3851_CBAL_2_b_P2 Adar_3851_CBAL_2_c_P2
# Adar_3851_CBAL_1_a_P2             0.9973189             0.9884591
# Adar_3851_CBAL_1_b_P2             0.9956844             0.9836571
# Adar_3851_CBAL_1_c_P2             1.0000000             0.9939702
# Adar_3851_CBAL_2_a_P2             0.2484471             0.2126547
# Adar_3851_CBAL_2_b_P2             0.0000000             0.1033426
# Adar_3851_CBAL_2_c_P2             0.1033426             0.0000000

##LBAL
data.good2.bc.mat[intersect(lbal2,s.3057.2),intersect(lbal2,s.3057.2)]
#                     Adar_3057_LBAL_a_P1 Adar_3057_LBAL_b_P1 Adar_3057_LBAL_c_P1
# Adar_3057_LBAL_a_P1           0.0000000           0.8691311           0.7546030
# Adar_3057_LBAL_b_P1           0.8691311           0.0000000           0.8539877
# Adar_3057_LBAL_c_P1           0.7546030           0.8539877           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3116.2),intersect(lbal2,s.3116.2)]
#                     Adar_3116_LBAL_a_P1 Adar_3116_LBAL_b_P1 Adar_3116_LBAL_c_P1
# Adar_3116_LBAL_a_P1           0.0000000           0.9854922           0.9784786
# Adar_3116_LBAL_b_P1           0.9854922           0.0000000           0.9798791
# Adar_3116_LBAL_c_P1           0.9784786           0.9798791           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3172.2),intersect(lbal2,s.3172.2)]
#                     Adar_3172_LBAL_a_P1 Adar_3172_LBAL_b_P1 Adar_3172_LBAL_c_P1
# Adar_3172_LBAL_a_P1           0.0000000           0.4068717           0.3867899
# Adar_3172_LBAL_b_P1           0.4068717           0.0000000           0.3247831
# Adar_3172_LBAL_c_P1           0.3867899           0.3247831           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3183.2),intersect(lbal2,s.3183.2)]
#                     Adar_3183_LBAL_a_P1 Adar_3183_LBAL_b_P1 Adar_3183_LBAL_c_P1
# Adar_3183_LBAL_a_P1           0.0000000           0.1547821           0.1757425
# Adar_3183_LBAL_b_P1           0.1547821           0.0000000           0.1417626
# Adar_3183_LBAL_c_P1           0.1757425           0.1417626           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3194.2),intersect(lbal2,s.3194.2)]
#                     Adar_3194_LBAL_a_P1 Adar_3194_LBAL_b_P1 Adar_3194_LBAL_c_P1
# Adar_3194_LBAL_a_P1           0.0000000           0.2496043           0.3428519
# Adar_3194_LBAL_b_P1           0.2496043           0.0000000           0.2064476
# Adar_3194_LBAL_c_P1           0.3428519           0.2064476           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3219.2),intersect(lbal2,s.3219.2)]
#                     Adar_3219_LBAL_a_P1 Adar_3219_LBAL_b_P1 Adar_3219_LBAL_c_P1
# Adar_3219_LBAL_a_P1           0.0000000                   1           0.8659698
# Adar_3219_LBAL_b_P1           1.0000000                   0           1.0000000
# Adar_3219_LBAL_c_P1           0.8659698                   1           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3264.2),intersect(lbal2,s.3264.2)]
#                     Adar_3264_LBAL_a_P1 Adar_3264_LBAL_b_P1 Adar_3264_LBAL_c_P1
# Adar_3264_LBAL_a_P1           0.0000000           0.4026169           0.4055817
# Adar_3264_LBAL_b_P1           0.4026169           0.0000000           0.3891938
# Adar_3264_LBAL_c_P1           0.4055817           0.3891938           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3301.2),intersect(lbal2,s.3301.2)]
#                     Adar_3301_LBAL_a_P1 Adar_3301_LBAL_b_P1 Adar_3301_LBAL_c_P1
# Adar_3301_LBAL_a_P1           0.0000000           0.3057190           0.3638998
# Adar_3301_LBAL_b_P1           0.3057190           0.0000000           0.3235398
# Adar_3301_LBAL_c_P1           0.3638998           0.3235398           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3356.2),intersect(lbal2,s.3356.2)]
#                     Adar_3356_LBAL_a_P1 Adar_3356_LBAL_b_P1 Adar_3356_LBAL_c_P1
# Adar_3356_LBAL_a_P1           0.0000000           0.3447412           0.3412773
# Adar_3356_LBAL_b_P1           0.3447412           0.0000000           0.3364013
# Adar_3356_LBAL_c_P1           0.3412773           0.3364013           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3378.2),intersect(lbal2,s.3378.2)]
#                     Adar_3378_LBAL_a_P1 Adar_3378_LBAL_b_P1 Adar_3378_LBAL_c_P1
# Adar_3378_LBAL_a_P1           0.0000000              1.0000           0.8783723
# Adar_3378_LBAL_b_P1           1.0000000              0.0000           0.9858000
# Adar_3378_LBAL_c_P1           0.8783723              0.9858           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3389.2),intersect(lbal2,s.3389.2)]
#                     Adar_3389_LBAL_a_P2 Adar_3389_LBAL_b_P2 Adar_3389_LBAL_c_P2
# Adar_3389_LBAL_a_P2           0.0000000           0.9233794           0.8354468
# Adar_3389_LBAL_b_P2           0.9233794           0.0000000           0.9036523
# Adar_3389_LBAL_c_P2           0.8354468           0.9036523           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3482.2),intersect(lbal2,s.3482.2)]
#                     Adar_3482_LBAL_a_P2 Adar_3482_LBAL_b_P2 Adar_3482_LBAL_c_P2
# Adar_3482_LBAL_a_P2           0.0000000           0.9071497           0.4329755
# Adar_3482_LBAL_b_P2           0.9071497           0.0000000           0.9143814
# Adar_3482_LBAL_c_P2           0.4329755           0.9143814           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3493.2),intersect(lbal2,s.3493.2)]
#                     Adar_3493_LBAL_a_P2 Adar_3493_LBAL_b_P2 Adar_3493_LBAL_c_P2
# Adar_3493_LBAL_a_P2           0.0000000           0.6634216           0.6282972
# Adar_3493_LBAL_b_P2           0.6634216           0.0000000           0.6101092
# Adar_3493_LBAL_c_P2           0.6282972           0.6101092           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3518.2),intersect(lbal2,s.3518.2)]
#                     Adar_3518_LBAL_a_P2 Adar_3518_LBAL_b_P2 Adar_3518_LBAL_c_P2
# Adar_3518_LBAL_a_P2           0.0000000           0.1533809           0.1193011
# Adar_3518_LBAL_b_P2           0.1533809           0.0000000           0.2046490
# Adar_3518_LBAL_c_P2           0.1193011           0.2046490           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3552.2),intersect(lbal2,s.3552.2)]
#                     Adar_3552_LBAL_a_P2 Adar_3552_LBAL_b_P2 Adar_3552_LBAL_c_P2
# Adar_3552_LBAL_a_P2           0.0000000           0.4071888           0.2896017
# Adar_3552_LBAL_b_P2           0.4071888           0.0000000           0.2852884
# Adar_3552_LBAL_c_P2           0.2896017           0.2852884           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3585.2),intersect(lbal2,s.3585.2)]
#                     Adar_3585_LBAL_a_P2 Adar_3585_LBAL_b_P2 Adar_3585_LBAL_c_P2
# Adar_3585_LBAL_a_P2           0.0000000           0.6866952           0.5745741
# Adar_3585_LBAL_b_P2           0.6866952           0.0000000           0.4704237
# Adar_3585_LBAL_c_P2           0.5745741           0.4704237           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3600.2),intersect(lbal2,s.3600.2)]
#                     Adar_3600_LBAL_a_P2 Adar_3600_LBAL_b_P2 Adar_3600_LBAL_c_P2
# Adar_3600_LBAL_a_P2           0.0000000           0.3359147           0.2405050
# Adar_3600_LBAL_b_P2           0.3359147           0.0000000           0.1982567
# Adar_3600_LBAL_c_P2           0.2405050           0.1982567           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3611.2),intersect(lbal2,s.3611.2)]
#                     Adar_3611_LBAL_a_P2 Adar_3611_LBAL_b_P2 Adar_3611_LBAL_c_P2
# Adar_3611_LBAL_a_P2          0.00000000          0.11535590          0.08044559
# Adar_3611_LBAL_b_P2          0.11535590          0.00000000          0.09763604
# Adar_3611_LBAL_c_P2          0.08044559          0.09763604          0.00000000

data.good2.bc.mat[intersect(lbal2,s.3644.2),intersect(lbal2,s.3644.2)]
#                     Adar_3644_LBAL_a_P2 Adar_3644_LBAL_b_P2 Adar_3644_LBAL_c_P2
# Adar_3644_LBAL_a_P2           0.0000000           0.5342316           0.4653982
# Adar_3644_LBAL_b_P2           0.5342316           0.0000000           0.5443947
# Adar_3644_LBAL_c_P2           0.4653982           0.5443947           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3769.2),intersect(lbal2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(lbal2,s.3840.2),intersect(lbal2,s.3840.2)]
#                     Adar_3840_LBAL_a_P2 Adar_3840_LBAL_b_P2 Adar_3840_LBAL_c_P2
# Adar_3840_LBAL_a_P2           0.0000000           0.3001038           0.3255986
# Adar_3840_LBAL_b_P2           0.3001038           0.0000000           0.3372318
# Adar_3840_LBAL_c_P2           0.3255986           0.3372318           0.0000000

data.good2.bc.mat[intersect(lbal2,s.3851.2),intersect(lbal2,s.3851.2)]
#                     Adar_3851_LBAL_a_P2 Adar_3851_LBAL_b_P2 Adar_3851_LBAL_c_P2
# Adar_3851_LBAL_a_P2           0.0000000           0.4719286           0.4219527
# Adar_3851_LBAL_b_P2           0.4719286           0.0000000           0.5938986
# Adar_3851_LBAL_c_P2           0.4219527           0.5938986           0.0000000

##RBAL
data.good2.bc.mat[intersect(rbal2,s.3057.2),intersect(rbal2,s.3057.2)]
#                     Adar_3057_RBAL_a_P1 Adar_3057_RBAL_b_P1 Adar_3057_RBAL_c_P1
# Adar_3057_RBAL_a_P1           0.0000000           0.4573208           0.4040763
# Adar_3057_RBAL_b_P1           0.4573208           0.0000000           0.4020968
# Adar_3057_RBAL_c_P1           0.4040763           0.4020968           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3116.2),intersect(rbal2,s.3116.2)]
#                     Adar_3116_RBAL_a_P1 Adar_3116_RBAL_b_P1 Adar_3116_RBAL_c_P1
# Adar_3116_RBAL_a_P1           0.0000000           0.9647271           0.9796171
# Adar_3116_RBAL_b_P1           0.9647271           0.0000000           0.9683415
# Adar_3116_RBAL_c_P1           0.9796171           0.9683415           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3172.2),intersect(rbal2,s.3172.2)]
#                     Adar_3172_RBAL_a_P1 Adar_3172_RBAL_b_P1 Adar_3172_RBAL_c_P1
# Adar_3172_RBAL_a_P1           0.0000000           0.3264653           0.2921797
# Adar_3172_RBAL_b_P1           0.3264653           0.0000000           0.1883176
# Adar_3172_RBAL_c_P1           0.2921797           0.1883176           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3183.2),intersect(rbal2,s.3183.2)]
#                     Adar_3183_RBAL_a_P1 Adar_3183_RBAL_b_P1 Adar_3183_RBAL_c_P1
# Adar_3183_RBAL_a_P1           0.0000000           0.1932947           0.2905169
# Adar_3183_RBAL_b_P1           0.1932947           0.0000000           0.2305535
# Adar_3183_RBAL_c_P1           0.2905169           0.2305535           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3194.2),intersect(rbal2,s.3194.2)]
#                     Adar_3194_RBAL_a_P1 Adar_3194_RBAL_b_P1 Adar_3194_RBAL_c_P1
# Adar_3194_RBAL_a_P1          0.00000000           0.1137645          0.09300278
# Adar_3194_RBAL_b_P1          0.11376447           0.0000000          0.10993881
# Adar_3194_RBAL_c_P1          0.09300278           0.1099388          0.00000000

data.good2.bc.mat[intersect(rbal2,s.3219.2),intersect(rbal2,s.3219.2)]
#                     Adar_3219_RBAL_a_P1 Adar_3219_RBAL_b_P1 Adar_3219_RBAL_c_P1
# Adar_3219_RBAL_a_P1           0.0000000           0.5444066           0.6860748
# Adar_3219_RBAL_b_P1           0.5444066           0.0000000           0.7036382
# Adar_3219_RBAL_c_P1           0.6860748           0.7036382           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3264.2),intersect(rbal2,s.3264.2)]
#                     Adar_3264_RBAL_a_P1 Adar_3264_RBAL_b_P1 Adar_3264_RBAL_c_P1
# Adar_3264_RBAL_a_P1           0.0000000           0.9858156                   1
# Adar_3264_RBAL_b_P1           0.9858156           0.0000000                   1
# Adar_3264_RBAL_c_P1           1.0000000           1.0000000                   0

data.good2.bc.mat[intersect(rbal2,s.3301.2),intersect(rbal2,s.3301.2)]
#                     Adar_3301_RBAL_a_P1 Adar_3301_RBAL_b_P1 Adar_3301_RBAL_c_P1
# Adar_3301_RBAL_a_P1           0.0000000           0.8272582           0.4839644
# Adar_3301_RBAL_b_P1           0.8272582           0.0000000           0.7186628
# Adar_3301_RBAL_c_P1           0.4839644           0.7186628           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3356.2),intersect(rbal2,s.3356.2)]
#                     Adar_3356_RBAL_a_P1 Adar_3356_RBAL_b_P1 Adar_3356_RBAL_c_P1
# Adar_3356_RBAL_a_P1           0.0000000           0.2541297           0.2602833
# Adar_3356_RBAL_b_P1           0.2541297           0.0000000           0.2656791
# Adar_3356_RBAL_c_P1           0.2602833           0.2656791           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3378.2),intersect(rbal2,s.3378.2)]
#                     Adar_3378_RBAL_a_P1 Adar_3378_RBAL_b_P1 Adar_3378_RBAL_c_P1
# Adar_3378_RBAL_a_P1           0.0000000           0.3557384           0.3041321
# Adar_3378_RBAL_b_P1           0.3557384           0.0000000           0.2760292
# Adar_3378_RBAL_c_P1           0.3041321           0.2760292           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3389.2),intersect(rbal2,s.3389.2)]
#                     Adar_3389_RBAL_a_P2 Adar_3389_RBAL_b_P2 Adar_3389_RBAL_c_P2
# Adar_3389_RBAL_a_P2           0.0000000           0.9983065                   1
# Adar_3389_RBAL_b_P2           0.9983065           0.0000000                   1
# Adar_3389_RBAL_c_P2           1.0000000           1.0000000                   0

data.good2.bc.mat[intersect(rbal2,s.3482.2),intersect(rbal2,s.3482.2)]
#                     Adar_3482_RBAL_a_P2 Adar_3482_RBAL_b_P2 Adar_3482_RBAL_c_P2
# Adar_3482_RBAL_a_P2           0.0000000           0.7843000           0.8588215
# Adar_3482_RBAL_b_P2           0.7843000           0.0000000           0.7014932
# Adar_3482_RBAL_c_P2           0.8588215           0.7014932           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3493.2),intersect(rbal2,s.3493.2)]
#                     Adar_3493_RBAL_a_P2 Adar_3493_RBAL_b_P2 Adar_3493_RBAL_c_P2
# Adar_3493_RBAL_a_P2           0.0000000           0.2402324           0.2330744
# Adar_3493_RBAL_b_P2           0.2402324           0.0000000           0.2477377
# Adar_3493_RBAL_c_P2           0.2330744           0.2477377           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3518.2),intersect(rbal2,s.3518.2)]
#                     Adar_3518_RBAL_a_P2 Adar_3518_RBAL_b_P2 Adar_3518_RBAL_c_P2
# Adar_3518_RBAL_a_P2           0.0000000           0.3103704           0.3435374
# Adar_3518_RBAL_b_P2           0.3103704           0.0000000           0.3048765
# Adar_3518_RBAL_c_P2           0.3435374           0.3048765           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3552.2),intersect(rbal2,s.3552.2)]
#                     Adar_3552_RBAL_a_P2 Adar_3552_RBAL_b_P2 Adar_3552_RBAL_c_P2
# Adar_3552_RBAL_a_P2           0.0000000           0.1724486           0.1552146
# Adar_3552_RBAL_b_P2           0.1724486           0.0000000           0.1969753
# Adar_3552_RBAL_c_P2           0.1552146           0.1969753           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3585.2),intersect(rbal2,s.3585.2)]
#                     Adar_3585_RBAL_a_P2 Adar_3585_RBAL_b_P2 Adar_3585_RBAL_c_P2
# Adar_3585_RBAL_a_P2           0.0000000           0.5350912           0.4884355
# Adar_3585_RBAL_b_P2           0.5350912           0.0000000           0.5130556
# Adar_3585_RBAL_c_P2           0.4884355           0.5130556           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3600.2),intersect(rbal2,s.3600.2)]
#                     Adar_3600_RBAL_a_P2 Adar_3600_RBAL_b_P2 Adar_3600_RBAL_c_P2
# Adar_3600_RBAL_a_P2           0.0000000           0.9833900           0.9853464
# Adar_3600_RBAL_b_P2           0.9833900           0.0000000           0.2701999
# Adar_3600_RBAL_c_P2           0.9853464           0.2701999           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3611.2),intersect(rbal2,s.3611.2)]
#                     Adar_3611_RBAL_a_P2 Adar_3611_RBAL_b_P2 Adar_3611_RBAL_c_P2
# Adar_3611_RBAL_a_P2          0.00000000          0.06603241           0.1587704
# Adar_3611_RBAL_b_P2          0.06603241          0.00000000           0.1537180
# Adar_3611_RBAL_c_P2          0.15877044          0.15371797           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3644.2),intersect(rbal2,s.3644.2)]
#                     Adar_3644_RBAL_a_P2 Adar_3644_RBAL_b_P2 Adar_3644_RBAL_c_P2
# Adar_3644_RBAL_a_P2          0.00000000          0.05147993          0.12026283
# Adar_3644_RBAL_b_P2          0.05147993          0.00000000          0.09995583
# Adar_3644_RBAL_c_P2          0.12026283          0.09995583          0.00000000

data.good2.bc.mat[intersect(rbal2,s.3769.2),intersect(rbal2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(rbal2,s.3840.2),intersect(rbal2,s.3840.2)]
#                     Adar_3840_RBAL_a_P2 Adar_3840_RBAL_b_P2 Adar_3840_RBAL_c_P2
# Adar_3840_RBAL_a_P2           0.0000000           0.2299599           0.1718772
# Adar_3840_RBAL_b_P2           0.2299599           0.0000000           0.1746311
# Adar_3840_RBAL_c_P2           0.1718772           0.1746311           0.0000000

data.good2.bc.mat[intersect(rbal2,s.3851.2),intersect(rbal2,s.3851.2)]
#                     Adar_3851_RBAL_a_P2 Adar_3851_RBAL_b_P2 Adar_3851_RBAL_c_P2
# Adar_3851_RBAL_a_P2           0.0000000           0.8238861           0.7732364
# Adar_3851_RBAL_b_P2           0.8238861           0.0000000           0.7264965
# Adar_3851_RBAL_c_P2           0.7732364           0.7264965           0.0000000

##EBC
data.good2.bc.mat[intersect(ebc2,s.3057.2),intersect(ebc2,s.3057.2)]
#                    Adar_3057_EBC_a_P1 Adar_3057_EBC_b_P1 Adar_3057_EBC_c_P1
# Adar_3057_EBC_a_P1          0.0000000          0.8765432               0.95
# Adar_3057_EBC_b_P1          0.8765432          0.0000000               0.68
# Adar_3057_EBC_c_P1          0.9500000          0.6800000               0.00

data.good2.bc.mat[intersect(ebc2,s.3116.2),intersect(ebc2,s.3116.2)]
#                    Adar_3116_EBC_a_P1 Adar_3116_EBC_b_P1 Adar_3116_EBC_c_P1
# Adar_3116_EBC_a_P1          0.0000000                  1          0.9991857
# Adar_3116_EBC_b_P1          1.0000000                  0          1.0000000
# Adar_3116_EBC_c_P1          0.9991857                  1          0.0000000

data.good2.bc.mat[intersect(ebc2,s.3172.2),intersect(ebc2,s.3172.2)]
 #                   Adar_3172_EBC_a_P1 Adar_3172_EBC_b_P1 Adar_3172_EBC_c_P1
# Adar_3172_EBC_a_P1          0.0000000          0.9976098          1.0000000
# Adar_3172_EBC_b_P1          0.9976098          0.0000000          0.9998828
# Adar_3172_EBC_c_P1          1.0000000          0.9998828          0.0000000

data.good2.bc.mat[intersect(ebc2,s.3183.2),intersect(ebc2,s.3183.2)]
#                    Adar_3183_EBC_a_P1 Adar_3183_EBC_b_P1 Adar_3183_EBC_c_P1
# Adar_3183_EBC_a_P1          0.0000000          0.9526033                  1
# Adar_3183_EBC_b_P1          0.9526033          0.0000000                  1
# Adar_3183_EBC_c_P1          1.0000000          1.0000000                  0

data.good2.bc.mat[intersect(ebc2,s.3194.2),intersect(ebc2,s.3194.2)]
#                    Adar_3194_EBC_a_P1 Adar_3194_EBC_b_P1 Adar_3194_EBC_c_P1
# Adar_3194_EBC_a_P1                  0                  1                  1
# Adar_3194_EBC_b_P1                  1                  0                  1
# Adar_3194_EBC_c_P1                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3219.2),intersect(ebc2,s.3219.2)]  #not sure why this one caught the control as well
#                    Adar_3219_EBC_a_P1 Adar_3219_EBC_b_P1 Adar_3219_EBC_c_P1
# Adar_3219_EBC_a_P1          0.0000000          0.9836334                  1
# Adar_3219_EBC_b_P1          0.9836334          0.0000000                  1
# Adar_3219_EBC_c_P1          1.0000000          1.0000000                  0

data.good2.bc.mat[intersect(ebc2,s.3264.2),intersect(ebc2,s.3264.2)]
#                    Adar_3264_EBC_a_P1 Adar_3264_EBC_b_P1 Adar_3264_EBC_c_P1
# Adar_3264_EBC_a_P1                  0                  1                  1
# Adar_3264_EBC_b_P1                  1                  0                  1
# Adar_3264_EBC_c_P1                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3301.2),intersect(ebc2,s.3301.2)]
#                    Adar_3301_EBC_a_P1 Adar_3301_EBC_b_P1 Adar_3301_EBC_c_P1
# Adar_3301_EBC_a_P1                  0                  1                  1
# Adar_3301_EBC_b_P1                  1                  0                  1
# Adar_3301_EBC_c_P1                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3356.2),intersect(ebc2,s.3356.2)]
#                    Adar_3356_EBC_a_P1 Adar_3356_EBC_b_P1 Adar_3356_EBC_c_P1
# Adar_3356_EBC_a_P1                  0                  1                  1
# Adar_3356_EBC_b_P1                  1                  0                  1
# Adar_3356_EBC_c_P1                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3378.2),intersect(ebc2,s.3378.2)]
#                    Adar_3378_EBC_a_P1 Adar_3378_EBC_b_P1 Adar_3378_EBC_c_P2
# Adar_3378_EBC_a_P1          0.0000000                  1          0.9808043
# Adar_3378_EBC_b_P1          1.0000000                  0          1.0000000
# Adar_3378_EBC_c_P2          0.9808043                  1          0.0000000

data.good2.bc.mat[intersect(ebc2,s.3389.2),intersect(ebc2,s.3389.2)]
#                    Adar_3389_EBC_a_P2 Adar_3389_EBC_b_P2 Adar_3389_EBC_c_P2
# Adar_3389_EBC_a_P2          0.0000000                  1          0.6580978
# Adar_3389_EBC_b_P2          1.0000000                  0          1.0000000
# Adar_3389_EBC_c_P2          0.6580978                  1          0.0000000

data.good2.bc.mat[intersect(ebc2,s.3482.2),intersect(ebc2,s.3482.2)]
#                    Adar_3482_EBC_a_P2 Adar_3482_EBC_b_P2 Adar_3482_EBC_c_P2
# Adar_3482_EBC_a_P2                  0                  1                  1
# Adar_3482_EBC_b_P2                  1                  0                  1
# Adar_3482_EBC_c_P2                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3493.2),intersect(ebc2,s.3493.2)]
#                    Adar_3493_EBC_a_P2 Adar_3493_EBC_b_P2 Adar_3493_EBC_c_P2
# Adar_3493_EBC_a_P2                  0                  1                  1
# Adar_3493_EBC_b_P2                  1                  0                  1
# Adar_3493_EBC_c_P2                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3518.2),intersect(ebc2,s.3518.2)]
#                    Adar_3518_EBC_a_P2 Adar_3518_EBC_b_P2 Adar_3518_EBC_c_P2
# Adar_3518_EBC_a_P2                  0                  1                  1
# Adar_3518_EBC_b_P2                  1                  0                  1
# Adar_3518_EBC_c_P2                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3552.2),intersect(ebc2,s.3552.2)]
#                      Adar_3552_EBC_1_a_P2 Adar_3552_EBC_1_b_P2 Adar_3552_EBC_1_c_P2 Adar_3552_EBC_2_a_P2
# Adar_3552_EBC_1_a_P2            0.0000000                    1            0.9831503             1.000000
# Adar_3552_EBC_1_b_P2            1.0000000                    0            1.0000000             1.000000
# Adar_3552_EBC_1_c_P2            0.9831503                    1            0.0000000             1.000000
# Adar_3552_EBC_2_a_P2            1.0000000                    1            1.0000000             0.000000
# Adar_3552_EBC_2_b_P2            1.0000000                    1            1.0000000             0.946355
# Adar_3552_EBC_2_c_P2            1.0000000                    1            1.0000000             1.000000
#                      Adar_3552_EBC_2_b_P2 Adar_3552_EBC_2_c_P2
# Adar_3552_EBC_1_a_P2             1.000000                    1
# Adar_3552_EBC_1_b_P2             1.000000                    1
# Adar_3552_EBC_1_c_P2             1.000000                    1
# Adar_3552_EBC_2_a_P2             0.946355                    1
# Adar_3552_EBC_2_b_P2             0.000000                    1
# Adar_3552_EBC_2_c_P2             1.000000                    0

data.good2.bc.mat[intersect(ebc2,s.3585.2),intersect(ebc2,s.3585.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(ebc2,s.3600.2),intersect(ebc2,s.3600.2)]
#                    Adar_3600_EBC_a_P2 Adar_3600_EBC_b_P2 Adar_3600_EBC_c_P2
# Adar_3600_EBC_a_P2          0.0000000          0.9273973                  1
# Adar_3600_EBC_b_P2          0.9273973          0.0000000                  1
# Adar_3600_EBC_c_P2          1.0000000          1.0000000                  0

data.good2.bc.mat[intersect(ebc2,s.3611.2),intersect(ebc2,s.3611.2)]
#                    Adar_3611_EBC_a_P2 Adar_3611_EBC_b_P2 Adar_3611_EBC_c_P2
# Adar_3611_EBC_a_P2          0.0000000          0.7901711                  1
# Adar_3611_EBC_b_P2          0.7901711          0.0000000                  1
# Adar_3611_EBC_c_P2          1.0000000          1.0000000                  0

data.good2.bc.mat[intersect(ebc2,s.3644.2),intersect(ebc2,s.3644.2)]
#                    Adar_3644_EBC_a_P2 Adar_3644_EBC_b_P2 Adar_3644_EBC_c_P2
# Adar_3644_EBC_a_P2                  0                  1                  1
# Adar_3644_EBC_b_P2                  1                  0                  1
# Adar_3644_EBC_c_P2                  1                  1                  0

data.good2.bc.mat[intersect(ebc2,s.3769.2),intersect(ebc2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(ebc2,s.3840.2),intersect(ebc2,s.3840.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(ebc2,s.3851.2),intersect(ebc2,s.3851.2)]
# <0 x 0 matrix>


##EBC Controls
data.good2.bc.mat[intersect(ebc.ctl2,s.3057.2),intersect(ebc.ctl2,s.3057.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3116.2),intersect(ebc.ctl2,s.3116.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3172.2),intersect(ebc.ctl2,s.3172.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3183.2),intersect(ebc.ctl2,s.3183.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3194.2),intersect(ebc.ctl2,s.3194.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3219.2),intersect(ebc.ctl2,s.3219.2)]
#                          Adar_3219_EBC_Blank_a_P1 Adar_3219_EBC_Blank_b_P1 Adar_3219_EBC_Blank_c_P1
# Adar_3219_EBC_Blank_a_P1                        0                        1                        1
# Adar_3219_EBC_Blank_b_P1                        1                        0                        1
# Adar_3219_EBC_Blank_c_P1                        1                        1                        0

data.good2.bc.mat[intersect(ebc.ctl2,s.3264.2),intersect(ebc.ctl2,s.3264.2)]
#                          Adar_3264_EBC_Blank_a_P1 Adar_3264_EBC_Blank_b_P1 Adar_3264_EBC_Blank_c_P1
# Adar_3264_EBC_Blank_a_P1                0.0000000                0.9889291                0.8782112
# Adar_3264_EBC_Blank_b_P1                0.9889291                0.0000000                0.9902931
# Adar_3264_EBC_Blank_c_P1                0.8782112                0.9902931                0.0000000

data.good2.bc.mat[intersect(ebc.ctl2,s.3301.2),intersect(ebc.ctl2,s.3301.2)]
#                          Adar_3301_EBC_Blank_a_P1 Adar_3301_EBC_Blank_b_P1 Adar_3301_EBC_Blank_c_P1
# Adar_3301_EBC_Blank_a_P1                        0                        1                        1
# Adar_3301_EBC_Blank_b_P1                        1                        0                        1
# Adar_3301_EBC_Blank_c_P1                        1                        1                        0

data.good2.bc.mat[intersect(ebc.ctl2,s.3356.2),intersect(ebc.ctl2,s.3356.2)]
#                          Adar_3356_EBC_Blank_a_P1 Adar_3356_EBC_Blank_b_P1 Adar_3356_EBC_Blank_c_P1
# Adar_3356_EBC_Blank_a_P1                0.0000000                0.6468393                0.7395849
# Adar_3356_EBC_Blank_b_P1                0.6468393                0.0000000                0.1408568
# Adar_3356_EBC_Blank_c_P1                0.7395849                0.1408568                0.0000000

data.good2.bc.mat[intersect(ebc.ctl2,s.3378.2),intersect(ebc.ctl2,s.3378.2)]
#                          Adar_3378_EBC_Blank_a_P1 Adar_3378_EBC_Blank_b_P1 Adar_3378_EBC_Blank_c_P1
# Adar_3378_EBC_Blank_a_P1                0.0000000                        1                0.9820139
# Adar_3378_EBC_Blank_b_P1                1.0000000                        0                1.0000000
# Adar_3378_EBC_Blank_c_P1                0.9820139                        1                0.0000000

data.good2.bc.mat[intersect(ebc.ctl2,s.3389.2),intersect(ebc.ctl2,s.3389.2)]
#                            Adar_3389_EBC_Blank_a_P2 Adar_3389_EBC_Blank_b_P2 Adar_3389_EBC_Blank_c_P2
# Adar_3389_EBC_Blank_a_P2                  0.0000000                0.9570451                1.0000000
# Adar_3389_EBC_Blank_b_P2                  0.9570451                0.0000000                1.0000000
# Adar_3389_EBC_Blank_c_P2                  1.0000000                1.0000000                0.0000000
# Adar_3389_EBC_Control_a_P2                1.0000000                0.9988532                0.9937811
# Adar_3389_EBC_Control_b_P2                0.9856238                0.9987915                0.9980551
# Adar_3389_EBC_Control_c_P2                1.0000000                0.9990802                0.9967598
#                            Adar_3389_EBC_Control_a_P2 Adar_3389_EBC_Control_b_P2 Adar_3389_EBC_Control_c_P2
# Adar_3389_EBC_Blank_a_P2                    1.0000000                  0.9856238                  1.0000000
# Adar_3389_EBC_Blank_b_P2                    0.9988532                  0.9987915                  0.9990802
# Adar_3389_EBC_Blank_c_P2                    0.9937811                  0.9980551                  0.9967598
# Adar_3389_EBC_Control_a_P2                  0.0000000                  0.9901039                  0.9829181
# Adar_3389_EBC_Control_b_P2                  0.9901039                  0.0000000                  0.9916724
# Adar_3389_EBC_Control_c_P2                  0.9829181                  0.9916724                  0.0000000

data.good2.bc.mat[intersect(ebc.ctl2,s.3482.2),intersect(ebc.ctl2,s.3482.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3493.2),intersect(ebc.ctl2,s.3493.2)]
#                          Adar_3493_EBC_Blank_a_P2 Adar_3493_EBC_Blank_b_P2 Adar_3493_EBC_Blank_c_P2
# Adar_3493_EBC_Blank_a_P2                0.0000000                0.9268329                        1
# Adar_3493_EBC_Blank_b_P2                0.9268329                0.0000000                        1
# Adar_3493_EBC_Blank_c_P2                1.0000000                1.0000000                        0

data.good2.bc.mat[intersect(ebc.ctl2,s.3518.2),intersect(ebc.ctl2,s.3518.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3552.2),intersect(ebc.ctl2,s.3552.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3585.2),intersect(ebc.ctl2,s.3585.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3600.2),intersect(ebc.ctl2,s.3600.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3611.2),intersect(ebc.ctl2,s.3611.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3644.2),intersect(ebc.ctl2,s.3644.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3769.2),intersect(ebc.ctl2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(ebc.ctl2,s.3840.2),intersect(ebc.ctl2,s.3840.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(ebc.ctl2,s.3851.2),intersect(ebc.ctl2,s.3851.2)]
# <0 x 0 matrix>


##Nasal Blank 
data.good2.bc.mat[intersect(nasal.blank2,s.3057.2),intersect(nasal.blank2,s.3057.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3116.2),intersect(nasal.blank2,s.3116.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3172.2),intersect(nasal.blank2,s.3172.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3183.2),intersect(nasal.blank2,s.3183.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3194.2),intersect(nasal.blank2,s.3194.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3219.2),intersect(nasal.blank2,s.3219.2)]
#                            Adar_3219_Nasal_Blank_a_P1 Adar_3219_Nasal_Blank_b_P1 Adar_3219_Nasal_Blank_c_P1
# Adar_3219_Nasal_Blank_a_P1                  0.0000000                  0.9978151                  0.9966865
# Adar_3219_Nasal_Blank_b_P1                  0.9978151                  0.0000000                  0.9771791
# Adar_3219_Nasal_Blank_c_P1                  0.9966865                  0.9771791                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3264.2),intersect(nasal.blank2,s.3264.2)]
#                            Adar_3264_Nasal_Blank_a_P1 Adar_3264_Nasal_Blank_b_P1 Adar_3264_Nasal_Blank_c_P1
# Adar_3264_Nasal_Blank_a_P1                  0.0000000                  0.9848255                  0.8201973
# Adar_3264_Nasal_Blank_b_P1                  0.9848255                  0.0000000                  1.0000000
# Adar_3264_Nasal_Blank_c_P1                  0.8201973                  1.0000000                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3301.2),intersect(nasal.blank2,s.3301.2)]
#                            Adar_3301_Nasal_Blank_a_P1 Adar_3301_Nasal_Blank_b_P1 Adar_3301_Nasal_Blank_c_P1
# Adar_3301_Nasal_Blank_a_P1                          0                          1                          1
# Adar_3301_Nasal_Blank_b_P1                          1                          0                          1
# Adar_3301_Nasal_Blank_c_P1                          1                          1                          0

data.good2.bc.mat[intersect(nasal.blank2,s.3356.2),intersect(nasal.blank2,s.3356.2)]
#                            Adar_3356_Nasal_Blank_a_P1 Adar_3356_Nasal_Blank_b_P1 Adar_3356_Nasal_Blank_c_P1
# Adar_3356_Nasal_Blank_a_P1                  0.0000000                  0.9957189                  0.9945392
# Adar_3356_Nasal_Blank_b_P1                  0.9957189                  0.0000000                  0.9984140
# Adar_3356_Nasal_Blank_c_P1                  0.9945392                  0.9984140                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3378.2),intersect(nasal.blank2,s.3378.2)]
#                            Adar_3378_Nasal_Blank_a_P1 Adar_3378_Nasal_Blank_b_P1 Adar_3378_Nasal_Blank_c_P2
# Adar_3378_Nasal_Blank_a_P1                  0.0000000                  0.5732790                  0.7296556
# Adar_3378_Nasal_Blank_b_P1                  0.5732790                  0.0000000                  0.5564824
# Adar_3378_Nasal_Blank_c_P2                  0.7296556                  0.5564824                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3389.2),intersect(nasal.blank2,s.3389.2)]
#                            Adar_3389_Nasal_Blank_a_P2 Adar_3389_Nasal_Blank_b_P2 Adar_3389_Nasal_Blank_c_P2
# Adar_3389_Nasal_Blank_a_P2                  0.0000000                  1.0000000                  0.9977442
# Adar_3389_Nasal_Blank_b_P2                  1.0000000                  0.0000000                  0.7786589
# Adar_3389_Nasal_Blank_c_P2                  0.9977442                  0.7786589                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3482.2),intersect(nasal.blank2,s.3482.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3493.2),intersect(nasal.blank2,s.3493.2)]
#                            Adar_3493_Nasal_Blank_a_P2 Adar_3493_Nasal_Blank_b_P2
# Adar_3493_Nasal_Blank_a_P2                          0                  1.0000000
# Adar_3493_Nasal_Blank_b_P2                          1                  0.0000000
# Adar_3493_Nasal_Blank_c_P2                          1                  0.9647811
#                            Adar_3493_Nasal_Blank_c_P2
# Adar_3493_Nasal_Blank_a_P2                  1.0000000
# Adar_3493_Nasal_Blank_b_P2                  0.9647811
# Adar_3493_Nasal_Blank_c_P2                  0.0000000

data.good2.bc.mat[intersect(nasal.blank2,s.3518.2),intersect(nasal.blank2,s.3518.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3552.2),intersect(nasal.blank2,s.3552.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3585.2),intersect(nasal.blank2,s.3585.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3600.2),intersect(nasal.blank2,s.3600.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3611.2),intersect(nasal.blank2,s.3611.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3644.2),intersect(nasal.blank2,s.3644.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3769.2),intersect(nasal.blank2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(nasal.blank2,s.3840.2),intersect(nasal.blank2,s.3840.2)]
# <0 x 0 matrix>

data.good2.bc.mat[intersect(nasal.blank2,s.3851.2),intersect(nasal.blank2,s.3851.2)]
# <0 x 0 matrix>


##Nasal Swab
data.good2.bc.mat[intersect(nasal2,s.3057.2),intersect(nasal2,s.3057.2)]
#                      Adar_3057_Nasal_a_P1 Adar_3057_Nasal_b_P1 Adar_3057_Nasal_c_P1
# Adar_3057_Nasal_a_P1            0.0000000            0.5134686            0.5169633
# Adar_3057_Nasal_b_P1            0.5134686            0.0000000            0.3606686
# Adar_3057_Nasal_c_P1            0.5169633            0.3606686            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3116.2),intersect(nasal2,s.3116.2)]
#                      Adar_3116_Nasal_a_P1 Adar_3116_Nasal_b_P1 Adar_3116_Nasal_c_P1
# Adar_3116_Nasal_a_P1           0.00000000            0.1171518           0.05262423
# Adar_3116_Nasal_b_P1           0.11715181            0.0000000           0.12273734
# Adar_3116_Nasal_c_P1           0.05262423            0.1227373           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3172.2),intersect(nasal2,s.3172.2)]
#                      Adar_3172_Nasal_a_P1 Adar_3172_Nasal_b_P1 Adar_3172_Nasal_c_P1
# Adar_3172_Nasal_a_P1            0.0000000            0.2662328            0.2876449
# Adar_3172_Nasal_b_P1            0.2662328            0.0000000            0.2061846
# Adar_3172_Nasal_c_P1            0.2876449            0.2061846            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3183.2),intersect(nasal2,s.3183.2)]
#                      Adar_3183_Nasal_a_P1 Adar_3183_Nasal_b_P1 Adar_3183_Nasal_c_P1
# Adar_3183_Nasal_a_P1            0.0000000           0.14009053           0.14304719
# Adar_3183_Nasal_b_P1            0.1400905           0.00000000           0.03777693
# Adar_3183_Nasal_c_P1            0.1430472           0.03777693           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3194.2),intersect(nasal2,s.3194.2)]
#                      Adar_3194_Nasal_a_P1 Adar_3194_Nasal_b_P1 Adar_3194_Nasal_c_P1
# Adar_3194_Nasal_a_P1            0.0000000            0.5184225            0.4804189
# Adar_3194_Nasal_b_P1            0.5184225            0.0000000            0.5445393
# Adar_3194_Nasal_c_P1            0.4804189            0.5445393            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3219.2),intersect(nasal2,s.3219.2)]
#                      Adar_3219_Nasal_a_P1 Adar_3219_Nasal_b_P1 Adar_3219_Nasal_c_P1
# Adar_3219_Nasal_a_P1            0.0000000            0.7079300            0.6914929
# Adar_3219_Nasal_b_P1            0.7079300            0.0000000            0.7847128
# Adar_3219_Nasal_c_P1            0.6914929            0.7847128            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3264.2),intersect(nasal2,s.3264.2)]
#                      Adar_3264_Nasal_a_P1 Adar_3264_Nasal_b_P1 Adar_3264_Nasal_c_P1
# Adar_3264_Nasal_a_P1            0.0000000            0.3596267            0.1525113
# Adar_3264_Nasal_b_P1            0.3596267            0.0000000            0.3073843
# Adar_3264_Nasal_c_P1            0.1525113            0.3073843            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3301.2),intersect(nasal2,s.3301.2)]
#                      Adar_3301_Nasal_a_P1 Adar_3301_Nasal_b_P1 Adar_3301_Nasal_c_P1
# Adar_3301_Nasal_a_P1           0.00000000           0.07899643           0.08024264
# Adar_3301_Nasal_b_P1           0.07899643           0.00000000           0.05179987
# Adar_3301_Nasal_c_P1           0.08024264           0.05179987           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3356.2),intersect(nasal2,s.3356.2)]
#                      Adar_3356_Nasal_a_P1 Adar_3356_Nasal_b_P1 Adar_3356_Nasal_c_P1
# Adar_3356_Nasal_a_P1           0.00000000           0.09763466           0.04766589
# Adar_3356_Nasal_b_P1           0.09763466           0.00000000           0.07280102
# Adar_3356_Nasal_c_P1           0.04766589           0.07280102           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3378.2),intersect(nasal2,s.3378.2)]
#                      Adar_3378_Nasal_a_P1 Adar_3378_Nasal_b_P1 Adar_3378_Nasal_c_P2
# Adar_3378_Nasal_a_P1            0.0000000            0.4897200            0.4864518
# Adar_3378_Nasal_b_P1            0.4897200            0.0000000            0.5606495
# Adar_3378_Nasal_c_P2            0.4864518            0.5606495            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3389.2),intersect(nasal2,s.3389.2)]
#                      Adar_3389_Nasal_a_P2 Adar_3389_Nasal_b_P2 Adar_3389_Nasal_c_P2
# Adar_3389_Nasal_a_P2            0.0000000           0.11035452           0.10197664
# Adar_3389_Nasal_b_P2            0.1103545           0.00000000           0.07853403
# Adar_3389_Nasal_c_P2            0.1019766           0.07853403           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3482.2),intersect(nasal2,s.3482.2)]
#                      Adar_3482_Nasal_a_P2 Adar_3482_Nasal_b_P2 Adar_3482_Nasal_c_P2
# Adar_3482_Nasal_a_P2            0.0000000            0.4030123             0.525641
# Adar_3482_Nasal_b_P2            0.4030123            0.0000000             0.373814
# Adar_3482_Nasal_c_P2            0.5256410            0.3738140             0.000000

data.good2.bc.mat[intersect(nasal2,s.3493.2),intersect(nasal2,s.3493.2)]
#                      Adar_3493_Nasal_a_P2 Adar_3493_Nasal_b_P2 Adar_3493_Nasal_c_P2
# Adar_3493_Nasal_a_P2            0.0000000            0.2508075            0.3418505
# Adar_3493_Nasal_b_P2            0.2508075            0.0000000            0.2748743
# Adar_3493_Nasal_c_P2            0.3418505            0.2748743            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3518.2),intersect(nasal2,s.3518.2)]
#                      Adar_3518_Nasal_a_P2 Adar_3518_Nasal_b_P2 Adar_3518_Nasal_c_P2
# Adar_3518_Nasal_a_P2           0.00000000           0.04812407           0.09024684
# Adar_3518_Nasal_b_P2           0.04812407           0.00000000           0.08531247
# Adar_3518_Nasal_c_P2           0.09024684           0.08531247           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3552.2),intersect(nasal2,s.3552.2)]
#                      Adar_3552_Nasal_a_P2 Adar_3552_Nasal_b_P2 Adar_3552_Nasal_c_P2
# Adar_3552_Nasal_a_P2           0.00000000           0.06497978            0.2305534
# Adar_3552_Nasal_b_P2           0.06497978           0.00000000            0.2663283
# Adar_3552_Nasal_c_P2           0.23055340           0.26632826            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3585.2),intersect(nasal2,s.3585.2)]
#                      Adar_3585_Nasal_a_P2 Adar_3585_Nasal_b_P2 Adar_3585_Nasal_c_P2
# Adar_3585_Nasal_a_P2            0.0000000           0.11719394           0.13991828
# Adar_3585_Nasal_b_P2            0.1171939           0.00000000           0.06424706
# Adar_3585_Nasal_c_P2            0.1399183           0.06424706           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3600.2),intersect(nasal2,s.3600.2)]
#                      Adar_3600_Nasal_a_P2 Adar_3600_Nasal_b_P2 Adar_3600_Nasal_c_P2
# Adar_3600_Nasal_a_P2            0.0000000           0.14408563           0.14323655
# Adar_3600_Nasal_b_P2            0.1440856           0.00000000           0.06048082
# Adar_3600_Nasal_c_P2            0.1432366           0.06048082           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3611.2),intersect(nasal2,s.3611.2)]
#                      Adar_3611_Nasal_a_P2 Adar_3611_Nasal_b_P2 Adar_3611_Nasal_c_P2
# Adar_3611_Nasal_a_P2            0.0000000            0.1618252            0.1873131
# Adar_3611_Nasal_b_P2            0.1618252            0.0000000            0.1767191
# Adar_3611_Nasal_c_P2            0.1873131            0.1767191            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3644.2),intersect(nasal2,s.3644.2)]
#                      Adar_3644_Nasal_a_P2 Adar_3644_Nasal_b_P2 Adar_3644_Nasal_c_P2
# Adar_3644_Nasal_a_P2           0.00000000           0.08462419           0.05066879
# Adar_3644_Nasal_b_P2           0.08462419           0.00000000           0.06499814
# Adar_3644_Nasal_c_P2           0.05066879           0.06499814           0.00000000

data.good2.bc.mat[intersect(nasal2,s.3769.2),intersect(nasal2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(nasal2,s.3840.2),intersect(nasal2,s.3840.2)]
#                      Adar_3840_Nasal_a_P2 Adar_3840_Nasal_b_P2 Adar_3840_Nasal_c_P2
# Adar_3840_Nasal_a_P2            0.0000000            0.4160265            0.1089219
# Adar_3840_Nasal_b_P2            0.4160265            0.0000000            0.3335227
# Adar_3840_Nasal_c_P2            0.1089219            0.3335227            0.0000000

data.good2.bc.mat[intersect(nasal2,s.3851.2),intersect(nasal2,s.3851.2)]
#                      Adar_3851_Nasal_a_P2 Adar_3851_Nasal_b_P2 Adar_3851_Nasal_c_P2
# Adar_3851_Nasal_a_P2            0.0000000            0.5052213            0.4188708
# Adar_3851_Nasal_b_P2            0.5052213            0.0000000            0.5259852
# Adar_3851_Nasal_c_P2            0.4188708            0.5259852            0.0000000


##Prewash
data.good2.bc.mat[intersect(prewash2,s.3057.2),intersect(prewash2,s.3057.2)]
#                        Adar_3057_PreWash_a_P1 Adar_3057_PreWash_b_P1 Adar_3057_PreWash_c_P1
# Adar_3057_PreWash_a_P1                      0                      1                      1
# Adar_3057_PreWash_b_P1                      1                      0                      1
# Adar_3057_PreWash_c_P1                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3116.2),intersect(prewash2,s.3116.2)]
#                        Adar_3116_PreWash_a_P1 Adar_3116_PreWash_b_P1 Adar_3116_PreWash_c_P1
# Adar_3116_PreWash_a_P1                      0                      1                      1
# Adar_3116_PreWash_b_P1                      1                      0                      1
# Adar_3116_PreWash_c_P1                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3172.2),intersect(prewash2,s.3172.2)]
#                        Adar_3172_PreWash_a_P1 Adar_3172_PreWash_b_P1 Adar_3172_PreWash_c_P1
# Adar_3172_PreWash_a_P1              0.0000000                      1              0.9861859
# Adar_3172_PreWash_b_P1              1.0000000                      0              1.0000000
# Adar_3172_PreWash_c_P1              0.9861859                      1              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3183.2),intersect(prewash2,s.3183.2)]
#                        Adar_3183_PreWash_a_P1 Adar_3183_PreWash_b_P1 Adar_3183_PreWash_c_P1
# Adar_3183_PreWash_a_P1              0.0000000                      1              0.9898136
# Adar_3183_PreWash_b_P1              1.0000000                      0              1.0000000
# Adar_3183_PreWash_c_P1              0.9898136                      1              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3194.2),intersect(prewash2,s.3194.2)]
#                        Adar_3194_PreWash_a_P1 Adar_3194_PreWash_b_P1 Adar_3194_PreWash_c_P1
# Adar_3194_PreWash_a_P1                      0                      1                      1
# Adar_3194_PreWash_b_P1                      1                      0                      1
# Adar_3194_PreWash_c_P1                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3219.2),intersect(prewash2,s.3219.2)]
#                        Adar_3219_PreWash_a_P1 Adar_3219_PreWash_b_P1 Adar_3219_PreWash_c_P1
# Adar_3219_PreWash_a_P1                      0                      1                      1
# Adar_3219_PreWash_b_P1                      1                      0                      1
# Adar_3219_PreWash_c_P1                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3264.2),intersect(prewash2,s.3264.2)]
#                        Adar_3264_PreWash_a_P1 Adar_3264_PreWash_b_P1 Adar_3264_PreWash_c_P1
# Adar_3264_PreWash_a_P1                      0              1.0000000              1.0000000
# Adar_3264_PreWash_b_P1                      1              0.0000000              0.9961718
# Adar_3264_PreWash_c_P1                      1              0.9961718              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3301.2),intersect(prewash2,s.3301.2)]
#                        Adar_3301_PreWash_a_P1 Adar_3301_PreWash_b_P1 Adar_3301_PreWash_c_P1
# Adar_3301_PreWash_a_P1              0.0000000              0.9652388              1.0000000
# Adar_3301_PreWash_b_P1              0.9652388              0.0000000              0.9978448
# Adar_3301_PreWash_c_P1              1.0000000              0.9978448              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3356.2),intersect(prewash2,s.3356.2)]
#                        Adar_3356_PreWash_a_P1 Adar_3356_PreWash_b_P1 Adar_3356_PreWash_c_P1
# Adar_3356_PreWash_a_P1              0.0000000              0.9972556              0.9995935
# Adar_3356_PreWash_b_P1              0.9972556              0.0000000              0.9900000
# Adar_3356_PreWash_c_P1              0.9995935              0.9900000              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3378.2),intersect(prewash2,s.3378.2)]
#                        Adar_3378_PreWash_a_P1 Adar_3378_PreWash_b_P1 Adar_3378_PreWash_c_P2
# Adar_3378_PreWash_a_P1              0.0000000              0.9395833              0.6988173
# Adar_3378_PreWash_b_P1              0.9395833              0.0000000              0.8860697
# Adar_3378_PreWash_c_P2              0.6988173              0.8860697              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3389.2),intersect(prewash2,s.3389.2)]
#                        Adar_3389_PreWash_a_P2 Adar_3389_PreWash_b_P2 Adar_3389_PreWash_c_P2
# Adar_3389_PreWash_a_P2              0.0000000              1.0000000              0.9873276
# Adar_3389_PreWash_b_P2              1.0000000              0.0000000              0.9976082
# Adar_3389_PreWash_c_P2              0.9873276              0.9976082              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3482.2),intersect(prewash2,s.3482.2)]
#                        Adar_3482_PreWash_a_P2 Adar_3482_PreWash_b_P2 Adar_3482_PreWash_c_P2
# Adar_3482_PreWash_a_P2              0.0000000              0.7992771              0.9233183
# Adar_3482_PreWash_b_P2              0.7992771              0.0000000              0.7739269
# Adar_3482_PreWash_c_P2              0.9233183              0.7739269              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3493.2),intersect(prewash2,s.3493.2)]
#                        Adar_3493_PreWash_a_P2 Adar_3493_PreWash_b_P2 Adar_3493_PreWash_c_P2
# Adar_3493_PreWash_a_P2                      0                      1                      1
# Adar_3493_PreWash_b_P2                      1                      0                      1
# Adar_3493_PreWash_c_P2                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3518.2),intersect(prewash2,s.3518.2)]
#                        Adar_3518_PreWash_a_P2 Adar_3518_PreWash_b_P2 Adar_3518_PreWash_c_P2
# Adar_3518_PreWash_a_P2              0.0000000              0.9944873              1.0000000
# Adar_3518_PreWash_b_P2              0.9944873              0.0000000              0.9978142
# Adar_3518_PreWash_c_P2              1.0000000              0.9978142              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3552.2),intersect(prewash2,s.3552.2)]
#                        Adar_3552_PreWash_a_P2 Adar_3552_PreWash_b_P2 Adar_3552_PreWash_c_P2
# Adar_3552_PreWash_a_P2               0.000000               1.000000               0.991338
# Adar_3552_PreWash_b_P2               1.000000               0.000000               0.973801
# Adar_3552_PreWash_c_P2               0.991338               0.973801               0.000000

data.good2.bc.mat[intersect(prewash2,s.3585.2),intersect(prewash2,s.3585.2)]
#                        Adar_3585_PreWash_a_P2 Adar_3585_PreWash_b_P2 Adar_3585_PreWash_c_P2
# Adar_3585_PreWash_a_P2                      0              1.0000000              1.0000000
# Adar_3585_PreWash_b_P2                      1              0.0000000              0.9742608
# Adar_3585_PreWash_c_P2                      1              0.9742608              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3600.2),intersect(prewash2,s.3600.2)]
#                        Adar_3600_PreWash_a_P2 Adar_3600_PreWash_b_P2 Adar_3600_PreWash_c_P2
# Adar_3600_PreWash_a_P2              0.0000000              1.0000000              0.8843844
# Adar_3600_PreWash_b_P2              1.0000000              0.0000000              0.9239748
# Adar_3600_PreWash_c_P2              0.8843844              0.9239748              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3611.2),intersect(prewash2,s.3611.2)]
#                        Adar_3611_PreWash_a_P2 Adar_3611_PreWash_b_P2 Adar_3611_PreWash_c_P2
# Adar_3611_PreWash_a_P2              0.0000000              0.9165886              0.9078810
# Adar_3611_PreWash_b_P2              0.9165886              0.0000000              0.8618297
# Adar_3611_PreWash_c_P2              0.9078810              0.8618297              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3644.2),intersect(prewash2,s.3644.2)]
#                        Adar_3644_PreWash_a_P2 Adar_3644_PreWash_b_P2 Adar_3644_PreWash_c_P2
# Adar_3644_PreWash_a_P2                      0                      1                      1
# Adar_3644_PreWash_b_P2                      1                      0                      1
# Adar_3644_PreWash_c_P2                      1                      1                      0

data.good2.bc.mat[intersect(prewash2,s.3769.2),intersect(prewash2,s.3769.2)]
# <0 x 0 matrix>
data.good2.bc.mat[intersect(prewash2,s.3840.2),intersect(prewash2,s.3840.2)]
#                        Adar_3840_PreWash_a_P2 Adar_3840_PreWash_b_P2 Adar_3840_PreWash_c_P2
# Adar_3840_PreWash_a_P2                      0              1.0000000              1.0000000
# Adar_3840_PreWash_b_P2                      1              0.0000000              0.9949734
# Adar_3840_PreWash_c_P2                      1              0.9949734              0.0000000

data.good2.bc.mat[intersect(prewash2,s.3851.2),intersect(prewash2,s.3851.2)]
#                        Adar_3851_PreWash_a_P2 Adar_3851_PreWash_b_P2 Adar_3851_PreWash_c_P2
# Adar_3851_PreWash_a_P2                      0                      1                      1
# Adar_3851_PreWash_b_P2                      1                      0                      1
# Adar_3851_PreWash_c_P2                      1                      1                      0

## There has got to be an easier way to do this...
subject2.vec<-c("s.3057.2","s.3116.2","s.3172.2","s.3183.2","s.3194.2","s.3219.2","s.3264.2","s.3301.2","s.3356.2","s.3378.2","s.3389.2","s.3482.2","s.3493.2","s.3518.2","s.3552.2","s.3585.2","s.3600.2","s.3611.2","s.3644.2","s.3769.2","s.3840.2","s.3851.2")


dist.table<-list("list",length(subject.vec))
i<-1
for(i in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[i],".2",sep="")
	dist.table[[i]]<-unique(melt(data.good2.bc.mat[intersect(oral2,get(temp.subject)),intersect(oral2,get(temp.subject))])[,3])
}

## Need to recreate the PCAs after removing background and Feces because I don't think that I'm going to use those samples.
#Now removing filters and feces from data.good
data.good3<-data.good[-c(feces,feces.blank,mock,water,blank,ae,iso,filter.air,filter.ctl),-c(4,8,20)]
subject3<-subject[-c(feces,feces.blank,mock,water,blank,ae,iso,filter.air,filter.ctl)]
subject3<-as.factor(as.character(subject3))
sample3<-sample[-c(feces,feces.blank,mock,water,blank,ae,iso,filter.air,filter.ctl)]
sample3<-as.factor(as.character(sample3))
data.good3.norm<-decostand(data.good3,"total")*100

data.good3.hel<-decostand(data.good3,"hellinger")
data.good3.pca<-rda(data.good3.hel)
head(summary(data.good3.pca))
#0.1905 0.03453

sample.cols<-brewer.pal(11,"Spectral")
sample.cols<-c("black",sample.cols)

plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (19.05% Explained)",ylab="PC2 (3.45% Explained)",main = "PCA of Samples by Sample Type")
points(data.good3.pca,pch=19,cex=0.7,col=sample.cols[as.numeric(sample3)])
ordispider(data.good3.pca,sample3,show.groups=sample3,col=sample.cols,label=T)
legend("topright",levels(sample3),pch=15,col=sample.cols)

#Now to show off groups
levels(sample3)
#  [1] "CBAL"        "EBC"         "EBC.Blank"   "LBAL"        "Nasal.Blank" "Nasal.Swap"  "Oral"       
#  [8] "Pre.Wash"    "RBAL"        "Saline"     

#CBAL
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:CBAL")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"CBAL",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"CBAL",col="red")
ordispider(data.good3.pca,sample3,show.groups="CBAL",col="red",label=T)
legend("topright",c("CBAL","Others"),pch=15,col=c("red","black"))

#EBC
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:EBC")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"EBC",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"EBC",col="red")
ordispider(data.good3.pca,sample3,show.groups="EBC",col="red",label=T)
legend("topright",c("EBC","Others"),pch=15,col=c("red","black"))

#EBC.Blank
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:EBC.Blank")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"EBC.Blank",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"EBC.Blank",col="red")
ordispider(data.good3.pca,sample3,show.groups="EBC.Blank",col="red",label=T)
legend("topright",c("EBC.Blank","Others"),pch=15,col=c("red","black"))

#LBAL
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:LBAL")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"LBAL",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"LBAL",col="red")
ordispider(data.good3.pca,sample3,show.groups="LBAL",col="red",label=T)
legend("topright",c("LBAL","Others"),pch=15,col=c("red","black"))

#Nasal.Blank
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Nasal.Blank")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"Nasal.Blank",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"Nasal.Blank",col="red")
ordispider(data.good3.pca,sample3,show.groups="Nasal.Blank",col="red",label=T)
legend("topright",c("Nasal.Blank","Others"),pch=15,col=c("red","black"))

#Nasal.Swap
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Nasal.Swap")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"Nasal.Swap",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"Nasal.Swap",col="red")
ordispider(data.good3.pca,sample3,show.groups="Nasal.Swap",col="red",label=T)
legend("topright",c("Nasal.Swap","Others"),pch=15,col=c("red","black"))

#Oral
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Oral")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"Oral",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"Oral",col="red")
ordispider(data.good3.pca,sample3,show.groups="Oral",col="red",label=T)
legend("topright",c("Oral","Others"),pch=15,col=c("red","black"))

#Pre.Wash
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Pre.Wash")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"Pre.Wash",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"Pre.Wash",col="red")
ordispider(data.good3.pca,sample3,show.groups="Pre.Wash",col="red",label=T)
legend("topright",c("Pre.Wash","Others"),pch=15,col=c("red","black"))

#RBAL
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:RBAL")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"RBAL",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"RBAL",col="red")
ordispider(data.good3.pca,sample3,show.groups="RBAL",col="red",label=T)
legend("topright",c("RBAL","Others"),pch=15,col=c("red","black"))

#Saline
plot(data.good3.pca,display="sites",type="n",font=2,font.lab=2,xlab="PC1 (18.7% Explained)",ylab="PC2 (3.3% Explained)",main = "PCA of Samples by Sample Type:Saline")
points(data.good3.pca,pch=19,cex=0.7,select=!sample3%in%"Saline",col="black")
points(data.good3.pca,pch=19,cex=0.7,select=sample3%in%"Saline",col="red")
ordispider(data.good3.pca,sample3,show.groups="Saline",col="red",label=T)
legend("topright",c("Saline","Others"),pch=15,col=c("red","black"))


## Some adonis to demonstrate samples different from controls
cbal3<-grep("CBAL",rownames(data.good3))
ebc.ctl3<-grep("EBC_Blank|EBC_Control",rownames(data.good3))
ebc3<-setdiff(grep("EBC",rownames(data.good3)),ebc.ctl3)
lbal3<-grep("LBAL",rownames(data.good3))
rbal3<-grep("RBAL",rownames(data.good3))
saline3<-grep("_NC_",rownames(data.good3))
nasal.blank3<-grep("Nasal_Blank",rownames(data.good3))
nasal3<-setdiff(grep("Nasal",rownames(data.good3)),nasal.blank3)
oral3<-grep("Oral",rownames(data.good3))
prewash3<-grep("Pre[Ww]ash",rownames(data.good3))

#EBC samples
adonis(data.good3.hel[c(ebc3,ebc.ctl3),]~sample3[c(ebc3,ebc.ctl3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 							   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# 	sample3[c(ebc3, ebc.ctl3)]  1     1.281 1.28072  1.3288 0.01654  0.003 **
# 	Residuals                  79    76.144 0.96385         0.98346          
# 	Total                      80    77.425                 1.00000          
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#BALs vs. Pre-Wash
adonis(data.good3.hel[c(rbal3,prewash3),]~sample3[c(rbal3,prewash3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 								 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# 	sample3[c(rbal3, prewash3)]   1    10.243 10.2428  12.621 0.09238  0.001 ***
# 	Residuals                   124   100.637  0.8116         0.90762           
# 	Total                       125   110.880                 1.00000           
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(data.good3.hel[c(rbal3,prewash3),]~sample3[c(rbal3,prewash3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 								 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# 	sample3[c(rbal3, prewash3)]   1    10.243 10.2428  12.621 0.09238  0.001 ***
# 	Residuals                   124   100.637  0.8116         0.90762           
# 	Total                       125   110.880                 1.00000           
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(data.good3.hel[c(lbal3,prewash3),]~sample3[c(lbal3,prewash3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 								 Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# 	sample3[c(lbal3, prewash3)]   1      9.18  9.1801  11.091 0.0821  0.001 ***
# 	Residuals                   124    102.64  0.8277         0.9179           
# 	Total                       125    111.81                 1.0000           
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Nasal vs Control
adonis(data.good3.hel[c(nasal3,nasal.blank3),]~sample3[c(nasal3,nasal.blank3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 									 Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# 	sample3[c(nasal3, nasal.blank3)]  1     4.785  4.7848  6.2377 0.07069  0.001 ***
# 	Residuals                        82    62.900  0.7671         0.92931           
# 	Total                            83    67.684                 1.00000           
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Oral vs. Saline
adonis(data.good3.hel[c(oral3,saline3),]~sample3[c(oral3,saline3)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 							   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# 	sample3[c(oral3, saline3)]  1     8.548  8.5480  16.872 0.16894  0.001 ***
# 	Residuals                  83    42.051  0.5066         0.83106           
# 	Total                      84    50.599                 1.00000           
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Filters for S&G
adonis(data.good.hel[c(filter.air,filter.ctl),]~sample[c(filter.air,filter.ctl)],method="euclidean")
# 	Permutation: free
# 	Number of permutations: 999
# 
# 	Terms added sequentially (first to last)
# 
# 									   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# 	sample[c(filter.air, filter.ctl)]   1    0.2766 0.27656  1.5317 0.01509  0.054 .
# 	Residuals                         100   18.0554 0.18055         0.98491         
# 	Total                             101   18.3320                 1.00000         
# 	---
# 	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### I want to create some Grouped Plots with Error faceted to get Control Mean, Control GeoMean, Sample Mean, Sample GeoMean.
data.good3.norm<-decostand(data.good3,"total")*100

##EBC Compare
ebc.compare.order<-names(sort(colMeans(data.good3.norm[ebc3,]), decreasing=T))

#This next block is to get some stats for the control and sample groups.  First Controls.
temp.gg<-melt(data.good3.norm[ebc.ctl3,ebc.compare.order[1:100]])
names(temp.gg)<-c("Sample","OTU","Percentage")
temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[ebc.compare.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))
temp.summary<-summarySE(temp.gg, measurevar="Percentage", groupvars="OTU")
head(temp.summary)
ebc.summary<-temp.summary[,c(1,3,4,6)]
colnames(ebc.summary)[2:4]<-c("Ctl.Mean","Ctl.GeoMean","Ctl.SEM")
#next samples
temp.gg<-melt(data.good3.norm[ebc3,ebc.compare.order[1:100]])
names(temp.gg)<-c("Sample","OTU","Percentage")
temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[ebc.compare.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))
temp.summary<-summarySE(temp.gg, measurevar="Percentage", groupvars="OTU")
ebc.summary<-cbind(ebc.summary,temp.summary[,c(3,4,6)])
colnames(ebc.summary)[5:7]<-c("EBC.Mean","EBC.GeoMean","EBC.SEM")
head(ebc.summary)

ebc.tmp<-rbind(ebc.summary[,2],ebc.summary[,3],ebc.summary[,5],ebc.summary[,6])
rownames(ebc.tmp)<-c("Ctl.Mean","Ctl.GeoMean","EBC.Mean","EBC.GeoMean")
colnames(ebc.tmp)<-ebc.summary[,1]
ebc.gg<-melt(ebc.tmp)
colnames(ebc.gg)<-c("Name","OTU","Percentage")
#the next line interleaves the SEM or 0 in a way that makes sense for ebc.gg
ebc.gg<-cbind(ebc.gg,as.vector(rbind(ebc.summary$Ctl.SEM,rep(0,100),ebc.summary$EBC.SEM,rep(0,100))))
colnames(ebc.gg)[4]<-"SEM"

#This plots Mean with error and geo mean for controls and samples
z<-ggplot(ebc.gg,aes(x=OTU,y=Percentage,fill=rep(data.good.taxonomy[ebc.compare.order[1:100],2],4)))+
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Name~.)+
	geom_errorbar(aes(ymin=Percentage-SEM,ymax=Percentage+SEM,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across EBC Controls and Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z
##

##Nasal Compare
nasal.compare.order<-names(sort(colMeans(data.good3.norm[nasal3,]), decreasing=T))

#This next block is to get some stats for the control and sample groups.  First Controls.
temp.gg<-melt(data.good3.norm[nasal.blank3,nasal.compare.order[1:100]])
names(temp.gg)<-c("Sample","OTU","Percentage")
temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[nasal.compare.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))
temp.summary<-summarySE(temp.gg, measurevar="Percentage", groupvars="OTU")
head(temp.summary)
nasal.summary<-temp.summary[,c(1,3,4,6)]
colnames(nasal.summary)[2:4]<-c("Ctl.Mean","Ctl.GeoMean","Ctl.SEM")
#next samples
temp.gg<-melt(data.good3.norm[nasal3,nasal.compare.order[1:100]])
names(temp.gg)<-c("Sample","OTU","Percentage")
temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU), labels = paste(data.good.taxonomy[nasal.compare.order[1:100],6]," (",levels(temp.gg$OTU),")",sep=""))
temp.summary<-summarySE(temp.gg, measurevar="Percentage", groupvars="OTU")
nasal.summary<-cbind(nasal.summary,temp.summary[,c(3,4,6)])
colnames(nasal.summary)[5:7]<-c("Nasal.Mean","Nasal.GeoMean","Nasal.SEM")
head(nasal.summary)

nasal.tmp<-rbind(nasal.summary[,2],nasal.summary[,3],nasal.summary[,5],nasal.summary[,6])
rownames(nasal.tmp)<-c("Ctl.Mean","Ctl.GeoMean","Nasal.Mean","Nasal.GeoMean")
colnames(nasal.tmp)<-nasal.summary[,1]
nasal.gg<-melt(nasal.tmp)
colnames(nasal.gg)<-c("Name","OTU","Percentage")
#the next line interleaves the SEM or 0 in a way that makes sense for nasal.gg
nasal.gg<-cbind(nasal.gg,as.vector(rbind(nasal.summary$Ctl.SEM,rep(0,100),nasal.summary$Nasal.SEM,rep(0,100))))
colnames(nasal.gg)[4]<-"SEM"

#This plots Mean with error and geo mean for controls and samples
z<-ggplot(nasal.gg,aes(x=OTU,y=Percentage,fill=rep(data.good.taxonomy[nasal.compare.order[1:100],2],4)))+
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Name~.)+
	geom_errorbar(aes(ymin=Percentage-SEM,ymax=Percentage+SEM,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Nasal Controls and Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z
##
#### I think I'm going to have to take the mean of the geo mean values to properly show a positive signal

temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[oral2,]),decreasing=T))

	if(length(intersect(oral2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(oral2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(oral2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(oral2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
oral.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
oral.summary<-cbind(oral.summary,c(rep("Oral Mean",100),rep("Oral Geo.Mean",100)))
colnames(oral.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[saline2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(saline2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(saline2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(saline2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(saline2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
saline.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
saline.summary<-cbind(saline.summary,c(rep(" Saline Mean",100),rep(" Saline Geo.Mean",100)))
colnames(saline.summary)[4]<-"Data"

oral.combined<-rbind(saline.summary,oral.summary)
### Saline was not replicated so it's a bad example to work with
###Lets try Nasal samples and controls
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[nasal2,]),decreasing=T))

	if(length(intersect(nasal2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(nasal2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(nasal2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(nasal2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
nasal.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
nasal.summary<-cbind(nasal.summary,c(rep("Nasal Mean",100),rep("Nasal Geo.Mean",100)))
colnames(nasal.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[nasal.blank2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(nasal.blank2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(nasal.blank2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(nasal.blank2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(nasal.blank2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
nasal.blank.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
nasal.blank.summary<-cbind(nasal.blank.summary,c(rep(" nasal.blank Mean",100),rep(" nasal.blank Geo.Mean",100)))
colnames(nasal.blank.summary)[4]<-"Data"

nasal.combined<-rbind(nasal.blank.summary,nasal.summary)

# Now to plot it
# Nasal and controls Mean of Means
z<-ggplot(nasal.combined,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],4)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Nasal Samples and Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

# Now Nasal alone
z<-ggplot(nasal.summary,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Nasal Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

#Now Controls but we have to re-run the above to get it's own order
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[nasal.blank2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(nasal.blank2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(nasal.blank2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(nasal.blank2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(nasal.blank2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
nasal.blank.summary2<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
nasal.blank.summary2<-cbind(nasal.blank.summary2,c(rep(" nasal.blank Mean",100),rep(" nasal.blank Geo.Mean",100)))
colnames(nasal.blank.summary2)[4]<-"Data"

z<-ggplot(nasal.blank.summary2,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Nasal Control Swabs", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z


## Now lets do EBC and controls
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[ebc2,]),decreasing=T))

	if(length(intersect(ebc2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(ebc2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(ebc2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(ebc2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
ebc.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
ebc.summary<-cbind(ebc.summary,c(rep("EBC Mean",100),rep("EBC Geo.Mean",100)))
colnames(ebc.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[ebc.ctl2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(ebc.ctl2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(ebc.ctl2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(ebc.ctl2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(ebc.ctl2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
ebc.ctl.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
ebc.ctl.summary<-cbind(ebc.ctl.summary,c(rep(" EBC Control Mean",100),rep(" EBC Control Geo.Mean",100)))
colnames(ebc.ctl.summary)[4]<-"Data"

ebc.combined<-rbind(ebc.ctl.summary,ebc.summary)

# Now to plot it
# ebc and controls Mean of Means
z<-ggplot(ebc.combined,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],4)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Exhaled Breath Condensate Samples and Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

# Now ebc alone
z<-ggplot(ebc.summary,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across EBC Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

#Now Controls but we have to re-run the above to get it's own order
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[ebc.ctl2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(ebc.ctl2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(ebc.ctl2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(ebc.ctl2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(ebc.ctl2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
ebc.ctl.summary2<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
ebc.ctl.summary2<-cbind(ebc.ctl.summary2,c(rep(" ebc.ctl Mean",100),rep(" ebc.ctl Geo.Mean",100)))
colnames(ebc.ctl.summary2)[4]<-"Data"

z<-ggplot(ebc.ctl.summary2,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across EBC Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

##LBALs and Controls

temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[lbal2,]),decreasing=T))

	if(length(intersect(lbal2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(lbal2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(lbal2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(lbal2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
lbal.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
lbal.summary<-cbind(lbal.summary,c(rep("LBAL Mean",100),rep("LBAL Geo.Mean",100)))
colnames(lbal.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[prewash2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(prewash2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(prewash2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
prewash.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
prewash.summary<-cbind(prewash.summary,c(rep(" Scope Wash Mean",100),rep(" Scope Wash Geo.Mean",100)))
colnames(prewash.summary)[4]<-"Data"

lbal.combined<-rbind(prewash.summary,lbal.summary)

# Now to plot it
# lbal and controls Mean of Means
z<-ggplot(lbal.combined,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],4)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across LBAL Samples and Scope Wash Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

# Now lbal alone
z<-ggplot(lbal.summary,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across LBAL Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

#Now Controls but we have to re-run the above to get it's own order
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[prewash2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(prewash2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(prewash2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
prewash.summary2<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
prewash.summary2<-cbind(prewash.summary2,c(rep(" Scope Wash Mean",100),rep(" Scope Wash Geo.Mean",100)))
colnames(prewash.summary2)[4]<-"Data"

z<-ggplot(prewash.summary2,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Scope Wash Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z


##RBAL and Controls

temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[rbal2,]),decreasing=T))

	if(length(intersect(rbal2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(rbal2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(rbal2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(rbal2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
rbal.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
rbal.summary<-cbind(rbal.summary,c(rep("RBAL Mean",100),rep("RBAL Geo.Mean",100)))
colnames(rbal.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[prewash2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(prewash2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(prewash2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
prewash.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
prewash.summary<-cbind(prewash.summary,c(rep(" Scope Wash Mean",100),rep(" Scope Wash Geo.Mean",100)))
colnames(prewash.summary)[4]<-"Data"

rbal.combined<-rbind(prewash.summary,rbal.summary)

# Now to plot it
# rbal and controls Mean of Means
z<-ggplot(rbal.combined,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],4)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across RBAL Samples and Scope Wash Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

# Now rbal alone
z<-ggplot(rbal.summary,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across RBAL Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

#We don't have to run Pre wash alone again


##CBALs and Controls

temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
#There is probably an easier way of doing this, but the purpose is to get a dataframe with the means of replicates for each subject
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	temp.order<-names(sort(colMeans(data.good2.norm[cbal2,]),decreasing=T))

	if(length(intersect(cbal2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(cbal2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(cbal2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(cbal2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
cbal.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
cbal.summary<-cbind(cbal.summary,c(rep("Combined BAL Mean",100),rep("Combined BAL Geo.Mean",100)))
colnames(cbal.summary)[4]<-"Data"

##Now lets do the above for the controls and combine the datasets
temp.means<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
temp.geomeans<-matrix(0,nrow=length(subject.vec),ncol=ncol(data.good2))
n<-1
for(n in 1:length(subject.vec)){
	temp.subject<-paste("s.",subject.vec[n],".2",sep="")
	#temp.order<-names(sort(colMeans(data.good2.norm[prewash2,]),decreasing=T)) #Removing in control samples to keep order the same

	if(length(intersect(prewash2,get(temp.subject)))==0){
			print(paste("Subject ",subject.vec[n]," does not have any samples",sep="")) 
			next()
			}else if(length(intersect(prewash2,get(temp.subject)))==1) {
					print(paste("Subject ",subject.vec[n]," has 1 sample",sep="")) 
			next()
			}else
			temp.means[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,mean)
			temp.geomeans[n,]<-apply(data.good2.norm[intersect(prewash2,get(temp.subject)),temp.order],2,geo.mean)
			
			}
colnames(temp.means)<-temp.order
colnames(temp.geomeans)<-temp.order
#first temp.means
temp.mean.gg<-melt(temp.means)
names(temp.mean.gg)<-c("Subject","OTU","Percentage")
temp.mean.gg$OTU <- factor(temp.mean.gg$OTU, levels = levels(temp.mean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.mean.gg$OTU),")",sep=""))
temp.mean.summary<-summarySE(temp.mean.gg, measurevar="Percentage", groupvars="OTU")
#next temp.geomeans
temp.geomean.gg<-melt(temp.geomeans)
names(temp.geomean.gg)<-c("Subject","OTU","Percentage")
temp.geomean.gg$OTU <- factor(temp.geomean.gg$OTU, levels = levels(temp.geomean.gg$OTU), labels = paste(data.good.taxonomy[temp.order,6]," (",levels(temp.geomean.gg$OTU),")",sep=""))
temp.geomean.summary<-summarySE(temp.geomean.gg, measurevar="Percentage", groupvars="OTU")
# Put into single dataframe
prewash.summary<-rbind(temp.mean.summary[1:100,c(1,3,6)],temp.geomean.summary[1:100,c(1,3,6)])
prewash.summary<-cbind(prewash.summary,c(rep(" Scope Wash Mean",100),rep(" Scope Wash Geo.Mean",100)))
colnames(prewash.summary)[4]<-"Data"

cbal.combined<-rbind(prewash.summary,cbal.summary)

# Now to plot it
# cbal and controls Mean of Means
z<-ggplot(cbal.combined,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],4)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Combined BAL Samples and Scope Wash Controls", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

# Now cbal alone
z<-ggplot(cbal.summary,aes(x=OTU,y=Mean.Percentage,fill=rep(data.good.taxonomy[temp.order[1:100],2],2)))+ #remember to change the repeat number to number of facets
	theme_bw()+
	geom_bar(stat="identity")+
	facet_grid(Data~.)+
	geom_errorbar(aes(ymin=Mean.Percentage-se,ymax=Mean.Percentage+se,width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Comparison of Arithmetic and Geometric Means Across Combined BAL Samples", x="OTUs",y="% Relative Abundance",fill="Phylum") #fill is the command for legend title
z

#We don't need to rerun Prewash alone

## Pulling the mock dilutions from a different data set (Replication and Cocentration2.rData)
# Neat Mock
mock.neat.order<-names(sort(colMeans(otu.good.norm[mock.neat,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.neat,mock.neat.order[1:100]],apply(otu.good.norm[mock.neat,mock.neat.order[1:100]],2,mean),apply(otu.good.norm[mock.neat,mock.neat.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community", x="OTUs",y="% Relative Abundance")
		temp.plot

#Mock 1:10
mock.d1.order<-names(sort(colMeans(otu.good.norm[mock.d1,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.d1,mock.d1.order[1:100]],apply(otu.good.norm[mock.d1,mock.d1.order[1:100]],2,mean),apply(otu.good.norm[mock.d1,mock.d1.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community 1:10 Dilution", x="OTUs",y="% Relative Abundance")
		temp.plot

# Mock 1:100
mock.d2.order<-names(sort(colMeans(otu.good.norm[mock.d2,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.d2,mock.d2.order[1:100]],apply(otu.good.norm[mock.d2,mock.d2.order[1:100]],2,mean),apply(otu.good.norm[mock.d2,mock.d2.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community 1:100 Dilution", x="OTUs",y="% Relative Abundance")
		temp.plot

#Mock 1:1000
mock.d3.order<-names(sort(colMeans(otu.good.norm[mock.d3,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.d3,mock.d3.order[1:100]],apply(otu.good.norm[mock.d3,mock.d3.order[1:100]],2,mean),apply(otu.good.norm[mock.d3,mock.d3.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community 1:1000 Dilution", x="OTUs",y="% Relative Abundance")
		temp.plot

#Mock 1:10,000
mock.d4.order<-names(sort(colMeans(otu.good.norm[mock.d4,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.d4,mock.d4.order[1:100]],apply(otu.good.norm[mock.d4,mock.d4.order[1:100]],2,mean),apply(otu.good.norm[mock.d4,mock.d4.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community 1:10000 Dilution", x="OTUs",y="% Relative Abundance")
		temp.plot

#Mock 1:100,000
mock.d5.order<-names(sort(colMeans(otu.good.norm[mock.d5,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[mock.d5,mock.d5.order[1:100]],apply(otu.good.norm[mock.d5,mock.d5.order[1:100]],2,mean),apply(otu.good.norm[mock.d5,mock.d5.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Mock Community 1:100000 Dilution", x="OTUs",y="% Relative Abundance")
		temp.plot

#Water connected with that data set

water.order<-names(sort(colMeans(otu.good.norm[water,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[water,water.order[1:100]],apply(otu.good.norm[water,water.order[1:100]],2,mean),apply(otu.good.norm[water,water.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Non-Template Controls", x="OTUs",y="% Relative Abundance")
		temp.plot

# AE connected with data
ae.order<-names(sort(colMeans(otu.good.norm[ae,]),decreasing=T))
		temp.df<-rbind(otu.good.norm[ae,ae.order[1:100]],apply(otu.good.norm[ae,ae.order[1:100]],2,mean),apply(otu.good.norm[ae,ae.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="AE Buffer Controls", x="OTUs",y="% Relative Abundance")
		temp.plot

#How about water ordered by mock 10^5 dilution
		temp.df<-rbind(otu.good.norm[water,mock.neat.order[1:100]],apply(otu.good.norm[water,mock.neat.order[1:100]],2,mean),apply(otu.good.norm[water,mock.neat.order[1:100]],2,geo.mean))
		rownames(temp.df)[c(nrow(temp.df)-1,nrow(temp.df))]<-c("Arith.Mean","Geo.Mean")
		#temp.df[which(temp==1)]<-0
		temp.gg<-melt(temp.df) #reformats the data (see below)
		names(temp.gg)<-c("Sample","OTU","Percentage")
		temp.gg$OTU <- factor(temp.gg$OTU, levels = levels(temp.gg$OTU))


		temp.plot<-ggplot(data=temp.gg,aes(y=Percentage,x=OTU,fill=Sample))+ 
 			theme_bw()+ #This theme gets rid of the grey background
 			geom_bar(stat="identity", position=position_stack(reverse=T))+ 
 			facet_grid(Sample~.)+ # creates a grid of bargraphs
 			scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 			theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
			labs(title="Non-Template Controls Ordered by 10^5 Dilution of Mock", x="OTUs",y="% Relative Abundance")
		temp.plot

### Lets see about Heatmap display of the bargraph data.  Started fresh R workspace.
# Make ComplexHeatmap for Adar Project
# for P.aeruginosa dilutions
# This is for the data in otu.good.norm
load('~/Dropbox/Replication and Concentration Study/Replication Study Workspace2.rData')
library(reshape2)
library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(vegan)
library(mvabund)
sem<-function(x){sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))}
geo.mean<-function(x, na.rm=TRUE){
   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
 }
library(RColorBrewer)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                       conf.interval=.95, .drop=TRUE) {
     library(plyr)
 
     # New version of length which can handle NA's: if na.rm==T, don't count them
     length2 <- function (x, na.rm=FALSE) {
         if (na.rm) sum(!is.na(x))
         else       length(x)
     }
 
    #Geometric Mean functiongeo.mean<-function(x, na.rm=TRUE){
   geo.mean<-function(x, na.rm=TRUE){
   	exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
 	}
       
     # This does the summary. For each group's data frame, return a vector with
     # N, mean, and sd
     datac <- ddply(data, groupvars, .drop=.drop,
       .fun = function(xx, col) {
         c(N    = length2(xx[[col]], na.rm=na.rm),
           mean = mean   (xx[[col]], na.rm=na.rm),
          geo.mean=geo.mean(xx[[col]],na.rm=na.rm),
           sd   = sd     (xx[[col]], na.rm=na.rm)       
         #geo.sd=10^(sd(decostand(xx[[col]],"log"),na.rm=na.rm))
         )
       },
       measurevar
     )
 
     # Rename the "mean" column    
     datac <- rename(datac, c("mean" = paste("Mean.",measurevar,sep="")))
 
     datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
 	datac$geo.mean[which(datac$geo.mean==1)]<-0 #Geometric Mean will take 0 abundance and transform to exactly 1.  This looks for values that are exactly 1 and sets them back to 0.
 	
     # Confidence interval multiplier for standard error
     # Calculate t-statistic for confidence interval: 
     # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
     ciMult <- qt(conf.interval/2 + .5, datac$N-1)
     datac$ci <- datac$se * ciMult
 	
     return(datac)
 }
library(Rtsne)
library(corrplot)
library(ComplexHeatmap)
library(vegan3d)
library(rgl)
library(ghibli)
library(circlize)

#Identify groups for Complex Heatmap
pa.dilutions<-c(rep("1:10",12),rep("1:100",12),rep("1:1000",12),rep("1:10000",12),rep("1:100000",12),rep(" Neat",12))
pa.dilutions<-as.factor(pa.dilutions)
pa.order<-names(sort(colMeans(otu.good.norm[pa,]),decreasing=T))
# to keep track of kits
pa.kitvec<-c(rep("Kit1",4),rep("Kit2",4),rep("Kit3",4),rep("Kit1",4),rep("Kit2",4),rep("Kit3",4),rep("Kit1",4),rep("Kit2",4),rep("Kit3",4),rep("Kit1",4),rep("Kit2",4),rep("Kit3",4),rep("Kit1",4),rep("Kit2",3),rep("Kit3",4),"Kit2",rep("Kit1",4),rep("Kit2",4),rep("Kit3",4))
pa.kitvec<-as.factor(pa.kitvec)
kit.cols<-brewer.pal(9,"Set1")
kit.cols<-kit.cols[c(3,4,5)]

#This generates a heatmap separated by dilution with a bar on the right identifying what kit they belong to.
Heatmap(otu.good.norm[pa,pa.order[1:75]],cluster_columns=F,cluster_rows=T,split=pa.dilutions,show_row_names=F,col=colorRamp2(c(0, 20,70,100), c("white", "skyblue", "cornflowerblue","blue")))+
Heatmap(pa.kitvec,name="Kit",show_row_names=F,col=kit.cols)


# Now for mock heatmap.  I'll try to move the neat.
mock.dilutions<-c(rep(" Neat",4),rep("1:10",4),rep("1:100",4),rep("1:1000",3),rep("1:10000",4),rep("1:100000",4))
mock.dilutions<-as factor(mock.dilutions)
mock.order<-names(sort(colMeans(otu.good.norm[mock,]),decreasing=T))

Heatmap(otu.good.norm[c(mock.neat,mock.d1,mock.d2,mock.d3,mock.d4,mock.d5),mock.order[1:75]],cluster_columns=F,cluster_rows=F,col=mapcols(100),row_split=mock.dilutions)
# Loop to create a mean of replicate medians.
i<-0
n<-1:4
mock.medians<-matrix(0,nrow=7,ncol=75)
for(i in 1:5){
		mock.medians[i,]<-apply(otu.good.norm[c(mock.neat,mock.d1,mock.d2,mock.d4,mock.d5)[n],mock.order[1:75]],2,median)
		n<-n+4
		}
mock.medians<-rbind(mock.medians[1:3,],apply(otu.good.norm[mock.d3,mock.order[1:75]],2,median),mock.medians[4:5,])		
colnames(mock.medians)<-mock.order[1:75]

#ha1 = HeatmapAnnotation(Avg.Median.Abundance = anno_barplot(apply(mock.medians,2,mean),height=unit(2,"cm")))
Heatmap(otu.good.norm[c(mock.neat,mock.d1,mock.d2,mock.d3,mock.d4,mock.d5),mock.order[1:75]],name="%Abundance",cluster_columns=F,cluster_rows=F,show_row_names=F,row_split=mock.dilutions,column_title="Abundance of OTUs Found in Mock Community Dilutions",col=colorRamp2(c(0, 20,70,100), c("white", "skyblue", "cornflowerblue","blue")))
#Heatmap(mock.dilutions,name="Control Replicate",col=ghibli_palette("PonyoMedium"),show_row_names=T)

# Lets try a 3D NMDS where dilutions are the 3rd axes
pa.nmds<-metaMDS(otu.good[pa,], trymax=1000)
pa.nmds.sc<-scores(pa.nmds,display="sites")
x <- pa.nmds.sc[,1]
y<-pa.nmds.sc[,2]
plot3d(0,0,0,xlim=c(min(x),max(x)),ylim=c(min(y),max(y)),z=c(0,5),type="n",xlab="Axis-1",ylab="Axis-2",zlab="Dilution")
spheres3d(x[grep("^N_",rownames(pa.nmds.sc))],y[grep("^N_",rownames(pa.nmds.sc))],0,type="s",radius=0.1,col=kit.cols[pa.kitvec[61:72]])
spheres3d(x[grep("^1_",rownames(pa.nmds.sc))],y[grep("^1_",rownames(pa.nmds.sc))],1,type="s",radius=0.1,col=kit.cols[pa.kitvec[1:12]])
spheres3d(x[grep("^2_",rownames(pa.nmds.sc))],y[grep("^2_",rownames(pa.nmds.sc))],2,type="s",radius=0.1,col=kit.cols[pa.kitvec[13:24]])
spheres3d(x[grep("^3_",rownames(pa.nmds.sc))],y[grep("^3_",rownames(pa.nmds.sc))],3,type="s",radius=0.1,col=kit.cols[pa.kitvec[25:36]])
spheres3d(x[grep("^4_",rownames(pa.nmds.sc))],y[grep("^4_",rownames(pa.nmds.sc))],4,type="s",radius=0.1,col=kit.cols[pa.kitvec[37:48]])
spheres3d(x[c(grep("^5_",rownames(pa.nmds.sc)),60)],y[c(grep("^5_",rownames(pa.nmds.sc)),60)],5,type="s",radius=0.1,col=kit.cols[pa.kitvec[49:60]])
# Add planes for clarity
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(0,0,0,0),alpha=0.1,col="grey")
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(1,1,1,1),alpha=0.1,col="grey")
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(2,2,2,2),alpha=0.1,col="grey")
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(3,3,3,3),alpha=0.1,col="grey")
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(4,4,4,4),alpha=0.1,col="grey")
rgl.quads(c(min(x),min(x),max(x),max(x)),c(min(y),max(y),max(y),min(y)),c(5,5,5,5),alpha=0.1,col="grey")

# Heatmap for the isolation controls in Adar samples
iso.rep.vec<-rep(1:7,each=3)
iso.rep.vec<-as.factor(iso.rep.vec)

# Loop to create a median of replicate values.
i<-0
n<-1:3
iso.medians<-matrix(0,nrow=7,ncol=75)
for(i in 1:7){
		iso.medians[i,]<-apply(data.good.norm[iso[n],iso.order[1:75]],2,median)
		n<-n+3
		}
colnames(iso.medians)<-iso.order[1:75]

ha1 = HeatmapAnnotation(Avg.Median.Abundance = anno_barplot(apply(iso.medians,2,mean),height=unit(2,"cm")))
Heatmap(data.good.norm[iso,iso.order[1:75]],name="%Abundance",cluster_columns=F,cluster_rows=T,show_row_names=F,top_annotation=ha1,column_title="Abundance of OTUs Found in Isolation Controls",col=colorRamp2(c(0, 20,70,100), c("white", "skyblue", "cornflowerblue","blue")))+
Heatmap(iso.rep.vec,name="Control Replicate",col=ghibli_palette("PonyoMedium"),show_row_names=F)

# Loop to create a mean of replicate geometric means.
# i<-0
# n<-1:3
# iso.geo.means<-matrix(0,nrow=7,ncol=75)
# for(i in 1:7){
# 		iso.geo.means[i,]<-apply(data.good.norm[iso[n],iso.order[1:75]],2,geo.mean)
# 		n<-n+3
# 		}
# 		
# iso.geo.means[which(iso.geo.means==1.0000000)]<-0 #changing all the 1's created by geo.mean back to 0
# colnames(iso.geo.means)<-iso.order[1:75]
# iso.mean.geomean<-apply(iso.geo.means,2,mean)
# iso.sem.geomean<-apply(iso.geo.means,2,sem)


# This is looking into a MeanLogTx
# i<-0
# n<-1:3
# data.good.log<-decostand(data.good,"log",logbase=10)
# iso.log.means<-matrix(0,nrow=7,ncol=75)
# for(i in 1:7){
# 		iso.log.means[i,]<-apply(data.good.log[iso[n],iso.order[1:75]],2,mean)
# 		n<-n+3
# 		}
# 
# iso.log.means[which(iso.log.means==1.0000000)]<-0 #changing all the 1's created by log.mean back to 0
# colnames(iso.log.means)<-iso.order[1:75]
# iso.mean.logmean<-apply(iso.log.means,2,mean)
# iso.sem.logmean<-apply(iso.log.means,2,sem)
# 
# ha1 = HeatmapAnnotation(dist1 = anno_barplot(iso.mean.logmean, bar_width = 1, gp = gpar(col = NA, fill = "#000000"), border = FALSE, axis = TRUE))
# Heatmap(data.good.norm[iso,iso.order[1:75]],name="%Abundance",cluster_col=F,cluster_rows=T,show_row_names=F,top_annotation=ha1,top_annotation_height=unit(2,"cm"),column_title="Abundance of OTUs Found in Isolation Controls")+
# Heatmap(iso.rep.vec,name="Control Replicate",col=ghibli_palette("PonyoMedium"),show_row_names=F)


##Comparing EBC with geo.mean, and meanlogtx
# Heatmap for the ebclation controls in Adar samples
ebc.rep.vec<-rep(1:19,each=3)
ebc.rep.vec<-as.factor(ebc.rep.vec)
ebc2.order<-names(sort(colMeans(data.good2.norm[ebc2,]),decreasing=T))

# Loop to create a mean of replicate geometric means.
i<-0
n<-1:3
ebc.geo.means<-matrix(0,nrow=19,ncol=75)
for(i in 1:19){
		ebc.geo.means[i,]<-apply(data.good2.norm[ebc2[n],ebc2.order[1:75]],2,geo.mean)
		n<-n+3
		}
		
ebc.geo.means[which(ebc.geo.means==1.0000000)]<-0 #changing all the 1's created by geo.mean back to 0
colnames(ebc.geo.means)<-ebc2.order[1:75]
ebc.mean.geomean<-apply(ebc.geo.means,2,mean)
ebc.sem.geomean<-apply(ebc.geo.means,2,sem)

ha1 = HeatmapAnnotation(dist1 = anno_barplot(ebc.mean.geomean, bar_width = 1, gp = gpar(col = NA, fill = "#000000"), border = FALSE, axis = TRUE))
Heatmap(data.good2.norm[ebc2,ebc2.order[1:75]],name="%Abundance",cluster_col=F,cluster_rows=T,show_row_names=F,top_annotation=ha1,top_annotation_height=unit(2,"cm"),column_title="Abundance of OTUs Found in Isolation Controls")+
Heatmap(ebc.rep.vec,name="EBC Replicate",col=group.cols[1:19],show_row_names=F)


# This is looking into a MeanLogTx
i<-0
n<-1:3
data.good2.log<-decostand(data.good2,"log",logbase=10)
ebc.log.means<-matrix(0,nrow=19,ncol=75)
for(i in 1:19){
		ebc.log.means[i,]<-apply(data.good2.log[ebc2[n],ebc2.order[1:75]],2,mean)
		n<-n+3
		}

colnames(ebc.log.means)<-ebc2.order[1:75]
ebc.mean.logmean<-apply(ebc.log.means,2,mean)
ebc.sem.logmean<-apply(ebc.log.means,2,sem)

ha1 = HeatmapAnnotation(dist1 = anno_barplot(ebc.mean.logmean, bar_width = 1, gp = gpar(col = NA, fill = "#000000"), border = FALSE, axis = TRUE))
Heatmap(data.good2.norm[ebc2,ebc2.order[1:75]],name="%Abundance",cluster_col=F,cluster_rows=T,show_row_names=F,top_annotation=ha1,top_annotation_height=unit(2,"cm"),column_title="Abundance of OTUs Found in Isolation Controls")+
Heatmap(ebc.rep.vec,name="EBC Replicate",col=group.cols[1:19],show_row_names=F)

## Barplots with error to compare geo.mean and log.mean
ebc.compare<-cbind(ebc.mean.geomean,ebc.sem.geomean,ebc.mean.logmean,ebc.sem.logmean)
rownames(ebc.compare)<-ebc2.order[1:75]
colnames(ebc.compare)<-c("Geo.mean","Geo.SEM","Log.mean","Log.mean.SEM")
ebc.compare<-as.data.frame(ebc.compare)

p1<-ggplot(ebc.compare,aes(x=rownames(ebc.compare),y=Geo.mean))+
	theme_bw()+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Geo.mean-Geo.SEM,ymax=Geo.mean+Geo.SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Geometric Means Across EBC Samples", x="OTUs",y="% Relative Abundance")
	
p2<-ggplot(ebc.compare,aes(x=rownames(ebc.compare),y=Log.mean))+
	theme_bw()+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Log.mean-Log.mean.SEM,ymax=Log.mean+Log.mean.SEM, width=0.2))+
	scale_y_continuous()+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="LogMeans Across EBC Samples", x="OTUs",y="Mean of Log (Counts)")

grid.arrange(p1,p2,nrow=2)

###= Same as above with oral
oral.rep.vec<-rep(1:21,each=3)
oral.rep.vec<-as.factor(oral.rep.vec)
oral2.order<-names(sort(colMeans(data.good2.norm[oral2,]),decreasing=T))

#Geo.mean
i<-0
n<-1:3
oral.geo.means<-matrix(0,nrow=21,ncol=75)
for(i in 1:21){
		oral.geo.means[i,]<-apply(data.good2.norm[oral2[n],oral2.order[1:75]],2,geo.mean)
		n<-n+3
		}
		
oral.geo.means[which(oral.geo.means==1.0000000)]<-0 #changing all the 1's created by geo.mean back to 0
colnames(oral.geo.means)<-oral2.order[1:75]
oral.mean.geomean<-apply(oral.geo.means,2,mean)
oral.sem.geomean<-apply(oral.geo.means,2,sem)

ha1 = HeatmapAnnotation(dist1 = anno_barplot(oral.mean.geomean, bar_width = 1, gp = gpar(col = NA, fill = "#000000"), border = FALSE, axis = TRUE))
Heatmap(data.good2.norm[oral2,oral2.order[1:75]],name="%Abundance",cluster_col=F,cluster_rows=T,show_row_names=F,top_annotation=ha1,top_annotation_height=unit(2,"cm"),column_title="Abundance of OTUs Found in Isolation Controls")+
Heatmap(oral.rep.vec,name="Oral Replicates",col=group.cols[1:21],show_row_names=F)

# Log.mean
i<-0
n<-1:3
data.good2.log<-decostand(data.good2,"log",logbase=10)
oral.log.means<-matrix(0,nrow=21,ncol=75)
for(i in 1:21){
		oral.log.means[i,]<-apply(data.good2.log[oral2[n],oral2.order[1:75]],2,mean)
		n<-n+3
		}

colnames(oral.log.means)<-oral2.order[1:75]
oral.mean.logmean<-apply(oral.log.means,2,mean)
oral.sem.logmean<-apply(oral.log.means,2,sem)

## Barplots with error to compare geo.mean and log.mean
ebc.compare<-cbind(ebc.mean.geomean,ebc.sem.geomean,ebc.mean.logmean,ebc.sem.logmean)
rownames(ebc.compare)<-ebc2.order[1:75]
colnames(ebc.compare)<-c("Geo.mean","Geo.SEM","Log.mean","Log.mean.SEM")
ebc.compare<-as.data.frame(ebc.compare)

p1<-ggplot(ebc.compare,aes(x=rownames(ebc.compare),y=Geo.mean))+
	theme_bw()+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Geo.mean-Geo.SEM,ymax=Geo.mean+Geo.SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Geometric Means Across EBC Samples", x="OTUs",y="% Relative Abundance")
	
p2<-ggplot(ebc.compare,aes(x=rownames(ebc.compare),y=Log.mean))+
	theme_bw()+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Log.mean-Log.mean.SEM,ymax=Log.mean+Log.mean.SEM, width=0.2))+
	scale_y_continuous()+
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="LogMeans Across EBC Samples", x="OTUs",y="Mean of Log (Counts)")

grid.arrange(p1,p2,nrow=2)

##Density plots.  I pulled the BC distances for within and between kits and put them into vectors *.intraBC or *.interBC

gg<-ggplot()+
	geom_density(aes(x=neat.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=neat.intraBC),fill="red",alpha=0.4,bw=0.05)+
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg

gg<-ggplot()+
	geom_density(aes(x=d1.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=d1.intraBC),fill="red",alpha=0.4,bw=0.05)+
	
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg

gg<-ggplot()+
	geom_density(aes(x=d2.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=d2.intraBC),fill="red",alpha=0.4,bw=0.05)+
	
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg

gg<-ggplot()+
	geom_density(aes(x=d3.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=d3.intraBC),fill="red",alpha=0.4,bw=0.05)+
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg

gg<-ggplot()+
	geom_density(aes(x=d4.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=d4.intraBC),fill="red",alpha=0.4,bw=0.05)+
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg

gg<-ggplot()+
	geom_density(aes(x=d5.interBC),fill="blue",alpha=0.4,bw=0.05)+
	geom_density(aes(x=d5.intraBC),fill="red",alpha=0.4,bw=0.05)+
	theme_bw()+
	scale_x_continuous(limits=c(0,1))
gg


###  Barplots of Means of replicate medians
cbal3.order<-names(sort(colMeans(data.good3.norm[cbal3,]),decreasing=T))[1:100]
data.good3.taxonomy<-data.taxonomy[colnames(data.good3.norm),]
data.good3.taxonomy<-as.data.frame(data.good3.taxonomy)
temp.otu.table<-data.good3.norm[,cbal3.order]
colnames(temp.otu.table)<-paste(data.good3.taxonomy[cbal3.order,6],"_",cbal3.order,sep="")

bal.compare.tbl<-gather(data.frame(Sample_Name=rownames(temp.otu.table),Sample_Type=sample3,Subject=subject3,temp.otu.table),key="OTU",value="Percentage",-c(Sample_Name,Sample_Type,Subject),factor_key=T)

### When I tried to get everything to plot together something was out of phase that doesn't happen when plotted singly.  I'll assemble in photohop.

#CBAL
temp.tbl<-filter(bal.compare.tbl,Sample_Type=="CBAL")%>%
    group_by(Sample_Type,Subject,OTU)%>%
    summarize(Rep.Median=median(Percentage))%>%
    group_by(Sample_Type,OTU)%>%
    summarize(Mean=mean(Rep.Median),SEM=sem(Rep.Median))

ggplot(temp.tbl, aes(x=OTU,y=Mean,fill=data.good3.taxonomy[cbal3.order,2]))+
	theme_bw()+
	#facet_grid(Sample_Type~.)+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean-SEM,ymax=Mean+SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""),limits=c(0,40))+ #requires scales library
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Rank Abundance of CBAL Means of Replicate Medians (Ordered by CBAL)", x="OTUs",y="% Relative Abundance")

#EBC
temp.tbl<-filter(bal.compare.tbl,Sample_Type=="EBC")%>%
    group_by(Sample_Type,Subject,OTU)%>%
    summarize(Rep.Median=median(Percentage))%>%
    group_by(Sample_Type,OTU)%>%
    summarize(Mean=mean(Rep.Median),SEM=sem(Rep.Median))

ggplot(temp.tbl, aes(x=OTU,y=Mean,fill=data.good3.taxonomy[cbal3.order,2]))+
	theme_bw()+
	#facet_grid(Sample_Type~.)+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean-SEM,ymax=Mean+SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""),limits=c(0,40))+ #requires scales library
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Rank Abundance of EBC Means of Replicate Medians (Ordered by CBAL)", x="OTUs",y="% Relative Abundance")

#EBC.Blank
temp.tbl<-filter(bal.compare.tbl,Sample_Type=="EBC.Blank")%>%
    group_by(Sample_Type,Subject,OTU)%>%
    summarize(Rep.Median=median(Percentage))%>%
    group_by(Sample_Type,OTU)%>%
    summarize(Mean=mean(Rep.Median),SEM=sem(Rep.Median))

ggplot(temp.tbl, aes(x=OTU,y=Mean,fill=data.good3.taxonomy[cbal3.order,2]))+
	theme_bw()+
	#facet_grid(Sample_Type~.)+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean-SEM,ymax=Mean+SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""),limits=c(0,40))+ #requires scales library
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Rank Abundance of EBC Control Means of Replicate Medians (Ordered by CBAL)", x="OTUs",y="% Relative Abundance")

#Pre.Wash
temp.tbl<-filter(bal.compare.tbl,Sample_Type=="Pre.Wash")%>%
    group_by(Sample_Type,Subject,OTU)%>%
    summarize(Rep.Median=median(Percentage))%>%
    group_by(Sample_Type,OTU)%>%
    summarize(Mean=mean(Rep.Median),SEM=sem(Rep.Median))

ggplot(temp.tbl, aes(x=OTU,y=Mean,fill=data.good3.taxonomy[cbal3.order,2]))+
	theme_bw()+
	#facet_grid(Sample_Type~.)+
	geom_bar(stat="identity")+
	geom_errorbar(aes(ymin=Mean-SEM,ymax=Mean+SEM, width=0.2))+
	scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""),limits=c(0,40))+ #requires scales library
 	theme(axis.text.x=element_text(angle=90,hjust=1))+ #angles axis text
	labs(title="Rank Abundance of Scope Control Means of Replicate Medians (Ordered by CBAL)", x="OTUs",y="% Relative Abundance")













