######################     GxETest.R program ####################################
##                              Overview                                                      #
# This program reads output created by GxEScan.exe.   
# Input files must be in the format created by the --gxeout option of GxEScan.  
# If GxEScan has been used to analyze multiple files of SNPs 
# (e.g. one file per chromosome), please use GxEMerge.exe to combine the files
# prior to executing GxETest.R               
######################################################################################


##                                        TimeLine                                                      ##
# Sep 5, 2014:  Create the functions 
# Sep 18, 2015: fix the problem by sorting the snps first by raw p-values then by rounded p-values 
#               Add one line to excel describing blanks under Ranks
#               Change the reading routine to read in .snpinfo file (SNPID now is created by gxescan)
#               Change .gxescan to gxeout for gxescan output files                
###########################################################################################################


########################## Required R packages #############################
# (1) data.table
# (2) XLConnect
# (3) rJava
# You can use R function "install.packages()" to install required R packages.
#
#############################################################################

############################################################
#####                  Paramter Input                  #####
rm(list=ls(all=TRUE))
library("data.table")
options(java.parameters = "-Xmx1024m")
suppressMessages(library(XLConnect))
library(rJava)

GWIS<-function(){
#####################################################################
#args <- commandArgs(trailingOnly = TRUE)
args<-c( "--home-dir","C:/Users/lopasy/Documents/UKB/cc",  	   # Path to GxEScan output files 
         "--alpha","0.05",                 # Family-wise error rate
		 "--gxescan","cc10",       # Base filename of GxEScan/GxEMerge output
		 "--marginal-g",              # Output Marginal G scan, Section B.2.1
		 "--case-control",            # Output CC GxE test, Section B.2.2
		 "--case-only",               # Output Case only analysis, Section B.2.3
		 "--control-only",            # Output Cntl-only analysis, Section B.2.3
		 "--2df",                     # Output 2df test, Section B.2.4
		 "--3df",                     # Output 3df test, Section B.2.5
## Each of the following four lines specifies the available 2-step methods (see Section B.2.6 of the documentation).  The first value (e.g. 0.05) is the step-1 significance threshold required for subset testing, and the second value (e.g. 5) is the initial bin size for weighted hypothesis testing (Section B.2.7)
		 "--EDGxE","0.05","5",        # EDGxE analysis, Section B.2.6.3
		 "--DGGxE","0.05","5",        # DG|GxE, Section B.2.6.1
		 "--GEGxE","0.05","5",        # GE|GxE, Section B.2.6.2
		 "--GE2df","0.05","5",        # GE|2df, Section B.2.6.4
## Output options
		 "--include-ranks",           # Output ranks for subset testing
		 "--top","50")                # Output top 50 SNPs from each method}
######################################################################################


n.args<- length(args)

for( i in 1:n.args){

if(grepl("--",args[i])){
CurrentArgName<- args[i]
j<-i+1
k<-j+1
l<-k+1
if(CurrentArgName=="--home-dir"){
homedir<- args[j]
} else if(CurrentArgName=="--alpha"){
alpha<- as.numeric(args[j])
} else if(CurrentArgName=="--gxescan"){
mytest<- args[j]
} else if(CurrentArgName=="--marginal-g"){
RunDG<- 1
} else if(CurrentArgName=="--case-control"){
RunCC<- 1
} else if(CurrentArgName=="--case-only"){
RunCO<- 1
} else if(CurrentArgName=="--control-only"){
RunCNTL<- 1
} else if(CurrentArgName=="--empirical-bayes"){
RunEB<- 1
} else if(CurrentArgName=="--2df"){
Run2df<- 1
} else if(CurrentArgName=="--3df"){
Run3df<- 1
} else if(CurrentArgName=="--EDGxE"){
RunEDGxE<-1
EDGxEscreen<- as.numeric(args[j])
EDGxE.bin<- as.integer(args[k])
} else if(CurrentArgName=="--DGGxE"){
RunDGGxE<-1
DGscreen<- as.numeric(args[j])
DGGxE.bin<- as.integer(args[k])
} else if(CurrentArgName=="--DGEB"){
RunDGEB<-1
DGEBscreen<- as.numeric(args[j])
DGEB.bin<- as.integer(args[k])
} else if(CurrentArgName=="--GEGxE"){
RunEGGxE<-1
EGscreen<- as.numeric(args[j])
EGGxE.bin<- as.integer(args[k])
} else if(CurrentArgName=="--GE2df"){
RunEG2df<-1
EG2dfscreen<- as.numeric(args[j])
EG2df.bin<- as.integer(args[k])
} else if(CurrentArgName=="--cocktail"){
RunCocktail<-1
Cocktailscreen<- as.numeric(args[j])
Cocktail.bin<- as.integer(args[k])
cocktail.c<- as.numeric(args[l])
} else if(CurrentArgName=="--include-ranks"){
IncludeRanks<- 1
} else if(CurrentArgName=="--top"){
topnum<- as.integer(args[j])
} 

}
}

if(!exists("IncludeRanks")){IncludeRanks=0}
if(!exists("RunDG")){RunDG=0}
if(!exists("RunCC")){RunCC=0}
if(!exists("RunCO")){RunCO=0}
if(!exists("RunCNTL")){RunCNTL=0}
if(!exists("RunEB")){RunEB=0}
if(!exists("Run2df")){Run2df=0}
if(!exists("Run3df")){Run3df=0}
if(!exists("RunEDGxE")){RunEDGxE=0}
if(!exists("RunDGGxE")){RunDGGxE=0}
if(!exists("RunDGEB")){RunDGEB=0}
if(!exists("RunEGGxE")){RunEGGxE=0}
if(!exists("RunEG2df")){RunEG2df=0}
if(!exists("RunCocktail")){RunCocktail=0}
if(!exists("topnum")){topnum=50}




##########################################################


### Set Working Directory to Location of gxescan output Files
setwd(homedir);

########################### Functions ###############################

### Function 1: Read in gxescan output files and remove NAs ##

ReadMplink<- function(filename){

if(!exists("mytest")){
res <-fread(paste("gxescan_",filename,sep=""), header=T, sep='\t',stringsAsFactors=F)
} else {
res <-fread(paste(mytest,"_",filename,sep=""), header=T, sep='\t',stringsAsFactors=F)
}


# Remove missing values
res.sub <- res[!is.na(P)]

return(res.sub)

}

### Function 2: Save image to excel sheets ##

saveimage <- function(wb,fname, sheetname,picname,row,col) {
sheet=sheetname
createName(wb, name = fname,formula = paste(sheet, idx2cref(c(row, col)), sep = "!"),overwrite = TRUE)
# Note: idx2cref converts indices (row, col) to Excel cell references
# Put the image created above at the corresponding location
addImage(wb, filename = picname, name = fname,originalSize = TRUE)

}

### Function 3: QQ Plot ## 

QQ 	<- function(p, t,max,PlotAll){

## Sort by P-value
setnames(p,"p")
setkey(p,p)
num.snp<-nrow(p)

if(PlotAll==1){
qqplot(-log10(seq(1,num.snp)/(num.snp+1)), -log10(p[,p]), xlim=c(0, max), ylim=c(0, max),main=t, ylab='Observed', xlab='Expected')
abline(0,1)
} else {

exp<- -log10(seq(1,num.snp)/(num.snp+1))
grpNum<-floor(log(num.snp))
pSub<-list()
expSub<-list()
pSub[[1]]<-p[1:5000]
expSub[[1]]<-exp[1:5000]
x<-5000
y<-15000
i=2
while(i<=grpNum & y<=num.snp){
pSub[[i]]<-p[seq((x+1),y,2^(i-1))]
expSub[[i]]<-exp[seq((x+1),y,2^(i-1))]
x<-y
y<-y+(1e4)*2^(i-1)
i=i+1
}
pSub[[i]]<-p[seq((x+1),num.snp,2^(i-1))]
expSub[[i]]<-exp[seq((x+1),num.snp,2^(i-1))]
p<-c(NA)
exp<-c(NA)
for(j in 1:i){
p<-c(p,pSub[[j]][,p])
exp<-c(exp,expSub[[j]])
}

p<-p[-1]
exp<-exp[-1]
num.snp<-length(p)
qqplot(exp, -log10(p), xlim=c(0, max), ylim=c(0, max),main=t, ylab='Observed', xlab='Expected')

abline(0,1)

}
}

### Function 4: Manhattan Plot  ##

manhattan <- function(data, sig, min.p, scale,PlotAll,step2){
     if(PlotAll==1){
      cutoff=0
     } else {
      cutoff=1
     }
	#colnames(data)	<- c('p', 'Chr', 'MapInfo')
	setnames(data,c(1:3),c('p', 'Chr', 'MapInfo'))
	data$y <- -log10(data$p)
      Chr<-list()
      for(i in 1:22){
	  t <- data[data$Chr==i,]
	  Chr[[i]] <- t[t$y>cutoff,]
       }    


## Scale mapinfo for each Chromosome to range between 0-1 and add a unit increase for successive Chr
        x<-list()
        for(i in 1:22){

	    x[[i]] <- (Chr[[i]]$MapInfo-min(Chr[[i]]$MapInfo))/((max(Chr[[i]]$MapInfo)+0.0001)-min(Chr[[i]]$MapInfo))*scale+i-1

        }    

colset<-c(rep(c("blue","olivedrab4","red","darkorange3","purple","orange"),3),"blue","olivedrab4","red","darkorange3")

	plot(x[[1]], Chr[[1]]$y, col="blue", xlab="Chromosome", ylab="-log10(p)", xlim=c(0,22), ylim=c(cutoff,min.p), axes=F, pch=20)
      for(i in 2:22){
	   points(x[[i]], Chr[[i]]$y, col=colset[i], pch=20)
      }

if(step2==1){		
abline(-1*log10(sig), 0, lwd=1)
}

axis(1,at=seq(scale/2-1,21+scale/2,2), label=seq(0,22,2), cex.axis=0.6)

axis(2,at=c(cutoff:min.p),label=c(cutoff:min.p),cex.axis=0.65)

}



### Function 6: Weighted testing Plot ## needs to be sorted by step 2 P-value and keyed by grp

wtplot <- function(results, min.p,title,scale,last.sig,num,PlotAll){

     if(PlotAll==1){
      cutoff=0
     } else {
      cutoff=1
     }

	setnames(results,c(1:4), c('p', 'grp','wt','MapInfo'))
	results[,y:= -log10(results[,p])]

      glist<-list()
      for(i in 1:num){
	t <- results[J(i)]
      t[,ref:=-1*log10(min(t[,wt]))]
      t[,x:=(t[,MapInfo]-min(t[,MapInfo]))/((max(t[,MapInfo])+0.0001)-min(t[,MapInfo]))*scale+i-1]
      glist[[i]]<-t
      rm(t)
      }


## Scale mapinfo for each Bin to range between 0-1 and add a unit increase for successive Bin


      color<-rep(c("blue","olivedrab4"),100)
      
      plot(glist[[1]][,x], glist[[1]][,y], col="blue", xlab="Bin # for step1 p-values", ylab="-log10(step2 p-values)", xlim=c(0,num), ylim=c(0,min.p), axes=F,pch=19)
      lines(glist[[1]][,x], glist[[1]][,ref],col="black",lwd=2)

      if(PlotAll==0){
      for(i in 2:num-1){
      if(nrow(glist[[i]])>=5115){
	points(glist[[i]][y>cutoff][,x], glist[[i]][y>cutoff][,y], col=color[i],pch=19)
      rect(-.1+i-1,-0.07,0.8+i-1,min(glist[[i]][y>cutoff][,y])+.1,col=color[i],border=color[i])
      } else {
	points(glist[[i]][,x], glist[[i]][,y], col=color[i],pch=19)      
      }
      lines(glist[[i]][,x], glist[[i]][,ref],col="black",lwd=2)
      }
      
       
	points(glist[[num]][y>cutoff][,x], glist[[num]][y>cutoff][,y], col=color[num],pch=19)
      rect(-.1+num-1,-0.07,0.8+num-1,min(glist[[num]][y>cutoff][,y])+.1,col=color[num],border=color[num])
      lines(glist[[num]][,x], last.sig,col="black",lwd=2)
} else {
      for(i in 2:num-1){
	points(glist[[i]][,x], glist[[i]][,y], col=color[i],pch=19)      
      lines(glist[[i]][,x], glist[[i]][,ref],col="black",lwd=2)
      }
      
       
	points(glist[[num]][,x], glist[[num]][,y], col=color[num],pch=19)
      lines(glist[[num]][,x], last.sig,col="black",lwd=2)

}	

decix<-(num/2)-floor(num/2)
if(decix>0){
axis(1,at=c(-1.5,seq(.5,num-0.5,2)), label=c(0,seq(1,num,2)),cex.axis=0.6)
} else {

axis(1,at=seq(-.5,num-0.5,2), label=seq(0,num,2),cex.axis=0.6)
}

deciy<-(min.p/2)-floor(min.p/2)

axis(2,at=c(0:floor(min.p)),label=c(0:min.p),cex.axis=0.65)

title (main=title,sub="Bin Size=5")


}


## Function 7: Set Title ##

savetitle <- function(wb,sheetname,ref,srow,scol) {

mergeCells(wb,sheet=sheetname,reference=ref)
cs <- createCellStyle(wb)
setFillForegroundColor(cs, color = XLC$"COLOR.GREY_25_PERCENT")
setFillPattern(cs, fill = XLC$"FILL.SOLID_FOREGROUND")
setBorder(cs, side = c("bottom"), type = XLC$"BORDER.THICK",color = c(XLC$"COLOR.BLACK"))
setCellStyle(wb, sheet = sheetname, row = srow, col = scol, cellstyle = cs)

}


## Function 8: save objects to Excel 

save2excel <- function(wb,x, sheetname,srow,scol,h) {

writeWorksheet(wb, x, sheet = sheetname, startRow = srow, startCol = scol,header=h)

}

## Function 9: wrap

wrap <- function(wb,x, sheetname,srow,scol,h,wrap) {

writeWorksheet(wb, x, sheet = sheetname, startRow = srow, startCol = scol,header=h)
cs <- createCellStyle(wb)
setWrapText(cs, wrap =wrap)
setCellStyle(wb, sheet = sheetname, row = srow, col = scol,cellstyle = cs)

}


## Function 10: Set column width
setwidth <- function(wb,sheetname,col,w) {

setColumnWidth(wb, sheet = sheetname, column = col, width = w)

} 

## Function 11: Create Bin and significance level for weighted testing approach (needs to be sorted by p-value)
wt<-function(pv,alpha,num){

rk.pv<-c(1:nrow(pv))
grp=ceiling(log(rk.pv/num+1,base=2))
pv[,Bin:=grp]
setkey(pv,Bin)
for(i in 1:max(grp))
{
pv[J(i),Threshold:=alpha*2^(-i)/nrow(pv[J(i)])]
}

}

## Function 12: Fit Cocktail method

FitCocktail <- function(snpinfo,marg,EG,EB,CC,alpha1)
   {
       Step1_P <- ifelse(marg<alpha1,marg,EG)
       Step2_P <- ifelse(marg<alpha1,EB,CC)
       Method <- ifelse(marg<alpha1,"DG|EB","EG|GxE")
	 data.table(snpinfo,Step1_P,Method,Step2_P)
   }



############################## GWIS Test ##################################
## Note: All top hits (subset testing) table are sorted by step-2 P-value
##       All top hits (weighted tesing) table are sorted by step-1 P-values
##       All Ranks are sorted by corresponding P-values   
###########################################################################

## Read in the SNP info ##

info <-fread(paste(mytest,".snpinfo",sep=""), header=T, sep='\t',stringsAsFactors=F)
#info <-fread("testfam.bim", header=F, sep='\t',stringsAsFactors=F)
#info[,c(3,6):=NULL]
#setnames(info,c(1:4),c("CHR","SNP","BP","A1"))
M.info<-nrow(info)
setkey(info,SNPID)

########## If PlotAll==1 then draw all poins in QQ and Manhattan plot ###########
if(M.info<1e5){
PlotAll=1
} else {
PlotAll=0
}

########### Method 1: Control-Only #############

if( RunCNTL==1){
print("Control Only")

## Read in data

#if(cov.type=="binary"){
#cntl<-ReadMplink("CntlOnly.assoc.logistic")
#} else {
#cntl<-ReadMplink("CntlOnly.assoc.linear")
#}

cntl<-ReadMplink("Cntl_GE.gxeout")

## Check whether or not we have significant hits
M.cntl<-nrow(cntl)
sig_cntl <- alpha/M.cntl

if(M.cntl==0) {
stop("All SNPs are not converged from Control Only analyses")
}

#if(M.cntl<M.info) print("NA found in control-only gxescan output file")

########## If PlotAll==1 then draw all poins in QQ and Manhattan plot ###########
if(nrow(cntl)<1e5){
PlotAll.cntl=1
} else {
PlotAll.cntl=0
}

### find the axis limit for plots
cntlp.min <- min(cntl[,P])
maxlimqq.cntl <- ceiling(1-log10(cntlp.min))
cntl.minman<-min(cntlp.min,sig_cntl)
maxlimman.cntl <-ceiling(1-log10(cntl.minman))

## Create CntlOnly rank if "IncludeRanks==1"

if(IncludeRanks==1){
cntl.rank<-cntl[,list(SNPID)]
cntl.rank[,CntlOnly:=1:M.cntl]
}

## Subset to pick the top hits

cntl.top<- cntl[1:topnum]

## Specify significance indicator for the top hits

cntl.top[,SigTemp:=(cntl.top[,P]<sig_cntl)*1]
cntl.top[SigTemp==1,Sig:="***"]
cntl.top[SigTemp==0,Sig:=""]
cntl.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(cntl,SNPID)
setkey(cntl.top,SNPID)

cntl<-info[cntl]

cntl.top<-info[cntl.top]


## Make QQ Plot

png('CNTL_A.png')

QQ(cntl[,list(P)], 'Control-only: E-G correlation',maxlimqq.cntl,PlotAll.cntl)

dev.off()

## Make Manhattan Plot

png('CNTL_B.png')

manhattan(cntl[,list(P,CHR,BP)], sig_cntl, maxlimman.cntl,0.7,PlotAll.cntl,1)
title('Control-only: E-G correlation')

dev.off()

rm(cntl)
setnames(cntl.top,c("NMISS","BETA","Z","P"),c("N","Beta","tTest","Pvalue"))
} # end of "if(RunCNTL==1)"


######### Method 2: Case-Only ###########

if(RunCO==1){
print("Case Only")

## Read in data
#if(cov.type=="binary"){
#co<-ReadMplink("CO.assoc.logistic")
#} else {
#co<-ReadMplink("CO.assoc.linear")
#}

co<-ReadMplink("Case_GE.gxeout")

## Check whether or not we have significant hits
M.co<-nrow(co)
sig_co <- alpha/M.co

if(M.co==0) {
stop("All SNPs are not converged from Case Only analyses")
}

if(M.co<1e5){
PlotAll.co=1
} else {
PlotAll.co=0
}

#if(M.co<M.info) print("NA found in case-only gxescan output file")

### find the axis limit for plots
cop.min <- min(co[,P])
maxlimqq.co <- ceiling(1-log10(cop.min))
co.minman<-min(cop.min,sig_co)
maxlimman.co <-ceiling(1-log10(co.minman))


## Create Case-Only rank if "IncludeRanks==1"

if(IncludeRanks==1){
co.rank<-co[,list(SNPID)]
co.rank[,CO:=1:M.co]
}


## Subset to pick the top hits

co.top<- co[1:topnum]
co.top<-co.top[,-6]

## Specify significance indicator for the top hits

co.top[,SigTemp:=(co.top[,P]<sig_co)*1]
co.top[SigTemp==1,Sig:="***"]
co.top[SigTemp==0,Sig:=""]
co.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(co,SNPID)
setkey(co.top,SNPID)

co<-info[co]

co.top<-info[co.top]

## Make QQ Plot

png('CO_A.png')

QQ(co[,list(P)], 'Case-only: GxE',maxlimqq.co,PlotAll.co)

dev.off()

## Make Manhattan Plot

png('CO_B.png')

manhattan(co[,list(P,CHR,BP)], sig_co, maxlimman.co,0.7,PlotAll.co,1)
title('Case-only: GxE')

dev.off()

rm(co)

setnames(co.top,c("NMISS","BETA","Z","P"),c("N","Beta","tTest","Pvalue"))

} # end of "if(RunCO==1)"


######### Method 3: 3df (EG+2df) ###########

if(Run3df==1){
print("3 df test")
## Read in data
df3<-ReadMplink("CC_3df.gxeout")

## Check whether or not we have significant hits
M.3df<-nrow(df3)
sig_3df <- alpha/M.3df

if(M.3df==0) {
stop("All SNPs are not converged from 3df test")
}

if(M.3df<1e5){
PlotAll.df3=1
} else {
PlotAll.df3=0
}

#if(M.3df<M.info) print("NA found in 3df gxescan output file")

### find the axis limit for plots
df3p.min <- min(df3[,P])
maxlimqq.3df <- ceiling(1-log10(df3p.min))
df3.minman<-min(df3p.min,sig_3df)
maxlimman.3df <-ceiling(1-log10(df3.minman))

## Create 3df rank if "IncludeRanks==1"

if(IncludeRanks==1){
df3.rank<-df3[,list(SNPID)]
df3.rank[,df3:=1:M.3df]
}

## Subset to pick the top hits

df3.top<- df3[1:topnum]

## Specify significance indicator for the top hits

df3.top[,SigTemp:=(df3.top[,P]<sig_3df)*1]
df3.top[SigTemp==1,Sig:="***"]
df3.top[SigTemp==0,Sig:=""]
df3.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(df3,SNPID)
setkey(df3.top,SNPID)

df3<-info[df3]

df3.top<-info[df3.top]


## Make QQ Plot

png('df3_A.png')

QQ(df3[,list(P)], '3df: Joint G,GxE and E-G',maxlimqq.3df,PlotAll.df3)

dev.off()

## Make Manhattan Plot

png('df3_B.png')

manhattan(df3[,list(P,CHR,BP)], sig_3df, maxlimman.3df,0.7,PlotAll.df3,1)
title('3df: Joint G,GxE and E-G')

dev.off()

rm(df3)

setnames(df3.top,c("NMISS","CHISQ_3DF","P"),c("N","df3_Chisq","Pvalue"))

} # end of "if(Run3df==1)"



######### Method 4: 2df Test ###########

if(Run2df==1 | RunEG2df==1){
print("2 df test")
## Read in data

df2<-ReadMplink("CC_2df.gxeout")

if(Run2df==1){
## Check whether or not we have significant hits
M.2df<-nrow(df2)
sig_2df <- alpha/M.2df

if(M.2df==0) {
stop("All SNPs are not converged from 2df test")
}

if(M.2df<1e5){
PlotAll.df2=1
} else {
PlotAll.df2=0
}

#if(M.2df<M.info) print("NA found in 2df gxescan output file")

### find the axis limit for plots
df2p.min <- min(df2[,P])
maxlimqq.2df <- ceiling(1-log10(df2p.min))
df2.minman<-min(df2p.min,sig_2df)
maxlimman.2df <-ceiling(1-log10(df2.minman))

## Create 2df rank if "IncludeRanks==1"

if(IncludeRanks==1){
df2.rank<-df2[,list(SNPID)]
df2.rank[,df2:=1:M.2df]
}

## Subset to pick the top hits

df2.top<- df2[1:topnum]

## Specify significance indicator for the top hits

df2.top[,SigTemp:=(df2.top[,P]<sig_2df)*1]
df2.top[SigTemp==1,Sig:="***"]
df2.top[SigTemp==0,Sig:=""]
df2.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(df2,SNPID)
setkey(df2.top,SNPID)

df2<-info[df2]

df2.top<-info[df2.top]


## Make QQ Plot

png('df2_A.png')

QQ(df2[,list(P)], 'Case-control: 2df G,GxE',maxlimqq.2df,PlotAll.df2)

dev.off()

## Make Manhattan Plot

png('df2_B.png')

manhattan(df2[,list(P,CHR,BP)], sig_2df, maxlimman.2df,0.7,PlotAll.df2,1)
title('Case-control: 2df G,GxE')

dev.off()

setnames(df2.top,c("NMISS","CHISQ_2DF","P"),c("N","df2_Chisq","Pvalue"))

} # end of "if(Run2df==1)"


if(Run2df!=1 & RunEG2df==1){
setkey(df2,SNPID)

df2<-info[df2]

}

if(RunEG2df!=1){
rm(df2)
}

} # end of "if(Run2df==1 | RunEG2df==1)"




## E-G ##
if( RunEG2df==1 | RunEGGxE==1 | RunCocktail==1){
## Read in E-G
eg<-ReadMplink("CC_GE.gxeout")

setnames(eg,c("CHISQ_1DF","P"),c("St1_tTest","Step1_P"))

## EG|2df ##
if( RunEG2df==1){
print("EG|2df")

wt(eg,alpha,EG2df.bin)

setkey(eg,SNPID)


## Merge with 2df

setnames(df2,c("NMISS","CHISQ_2DF","P"),c("N","St2_Chisq","Step2_P"))

eg2df<-df2[eg,nomatch=0]
rm(df2)

M.eg2df<-nrow(eg2df)

if(M.eg2df==0) {
stop("All SNPs are not converged from EG|2df")
}

if(M.eg2df<1e5){
PlotAll.eg2df=1
} else {
PlotAll.eg2df=0
}

#if(M.eg2df<M.info) print("NA found in gxescan output file")

## Sort eg2df by step1 P-value
setkey(eg2df,Step1_P)

## Subset to pick the E-G top hits into step 2 (subset testing)
eg2df.s2<- eg2df[Step1_P<EG2dfscreen]
eg2df.s2[,c("Bin","Threshold"):=NULL]
eg2df[,SNPID:=NULL]
eg2df.s2.size<-nrow(eg2df.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
eg2df[,SigTemp:=(eg2df[,Step2_P]<eg2df[,Threshold])*1]
eg2df[SigTemp==1,Sig:="***"]
eg2df[SigTemp==0,Sig:=""]
eg2df[,SigTemp:=NULL]
eg2df[,Rank:=1:M.eg2df]

eg2df.wt.hits<-eg2df[Sig=="***"]
eg2df.wt.top<-eg2df[1:topnum]

eg2df.wt.top<-rbind(eg2df.wt.hits[!eg2df.wt.hits$Rank %in% eg2df.wt.top$Rank,],eg2df.wt.top)

rm(eg2df.wt.hits)


## Compute limit for step 1 plots
eg.min <- min(eg2df[,Step1_P])
maxlimqq.eg <- ceiling(1-log10(eg.min))

## Step1 QQ Plot ( will be used for both EG|GxE and EG|2df )

png('EG2df_A.png')
QQ(eg2df[,list(Step1_P)], 'EG | 2df: Step1 screen',maxlimqq.eg,PlotAll.eg2df)
dev.off()

## Step1 Manhattan Plot ( will be used for both EG|GxE and EG|2df )
png('EG2df_B.png')
manhattan(eg2df[,list(Step1_P, CHR, BP)],0,maxlimqq.eg,0.7,PlotAll.eg2df,0)
title('EG | 2df: Step1 screen')
dev.off()

eg2df[,c("CHR","SNP","BP","A1","N","St2_Chisq","St1_tTest","Step1_P","Sig"):=NULL]


if (eg2df.s2.size>0) {        ## If we have SNPs pass to step2
eg2df.min<-min(eg2df.s2[,Step2_P])
sig_eg2df_st2 <- alpha/eg2df.s2.size
eg2df.min2<-min(eg2df.min,sig_eg2df_st2)
maxlimqq.eg2df <- ceiling(1-log10(eg2df.min))
maxlimman.eg2df <- ceiling(1-log10(eg2df.min2))

if(eg2df.s2.size>1e5){
PlotAll.eg2df.stp2=0
} else {
PlotAll.eg2df.stp2=1
}

## Step2 QQ Plot
png('EG2df_C.png')
QQ(eg2df.s2[,list(Step2_P)], 'EG | 2df: Step 2 test',maxlimqq.eg2df,PlotAll.eg2df.stp2)
dev.off()

## Step2 Manhattan Plot
png('EG2df_D.png')
manhattan(eg2df.s2[,list(Step2_P, CHR, BP)],sig_eg2df_st2, maxlimman.eg2df,0.7,PlotAll.eg2df.stp2,1)
title('EG | 2df: Step 2 test')
dev.off()


## Compute significance indicator for subset testing
setkey(eg2df.s2,Step2_P)
eg2df.s2[,SigTemp:=(eg2df.s2[,Step2_P]<sig_eg2df_st2)*1]
eg2df.s2[SigTemp==1,Sig:="***"]
eg2df.s2[SigTemp==0,Sig:=""]
eg2df.s2[,SigTemp:=NULL]

## Create EG2df rank if "IncludeRanks==1"
if(IncludeRanks==1){
eg2df.rank<-eg2df.s2[,list(SNPID)]
eg2df.rank[,EG2df:=1:eg2df.s2.size]
eg2df.s2<-eg2df.s2[1:topnum]
}
} else {

sig_eg2df_st2<-NA

}



## Weighted testing
setkey(eg2df,Bin)
# Compute the limit for plot
eg2df.last<-max(eg2df[,Bin])
eg2df.k<-EG2df.bin*(2^(eg2df.last-1))
eg2df.wt.binsig<-alpha*((1/2)^eg2df.last)
eg2df.wt.sig<-eg2df.wt.binsig/eg2df.k
eg2df.wt.lnum<-nrow(eg2df[J(eg2df.last)])
eg2df.wt.ref<- rep(-1*log10(eg2df.wt.sig),eg2df.wt.lnum)

# Make plot
png('EG2df_E.png')

wteg2df.ref<-max(ceiling(1-log10(min(eg2df[,Step2_P]))),eg2df.wt.ref[1])
wtplot(eg2df,wteg2df.ref+1,"EG | 2df: weighted",0.7,eg2df.wt.ref,eg2df.last,PlotAll.eg2df)

dev.off()

rm(eg2df)
rm(eg2df.wt.ref)

eg[,c("Bin","Threshold"):=NULL]

} # End of RunEG2df 

if(RunEGGxE!=1 & RunCocktail!=1){
rm(eg)
} 


} # End of E-G


## case control (cc) ##

if( RunCC==1 | RunEDGxE==1 | RunDGGxE==1 |RunEGGxE==1 |RunCocktail==1){

## Read in data
cc<-ReadMplink("CC_GxE.gxeout")

if(RunCC==1){

print("Exhaustive CC scan")

## Check whether or not we have significant hits
M.cc<-nrow(cc)
sig_cc <- alpha/M.cc

if(M.cc==0) {
stop("All SNPs are not converged from case-control analyses")
}

if(M.cc<1e5){
PlotAll.cc=1
} else {
PlotAll.cc=0
}

#if(M.cc<M.info) print("NA found in case-control gxescan output file")

### find the axis limit for plots
ccp.min <- min(cc[,P])
maxlimqq.cc <- ceiling(1-log10(ccp.min))
cc.minman<-min(ccp.min,sig_cc)
maxlimman.cc <-ceiling(1-log10(cc.minman))

## Create CC rank if "IncludeRanks==1"

if(IncludeRanks==1){
cc.rank<-cc[,list(SNPID)]
cc.rank[,CC:=1:M.cc] ## Since SNPs are sorted by P-value in the original gxescan file 

}


## Subset to pick the top hits

cc.top<- cc[1:topnum]


## Specify significance indicator for the top hits

cc.top[,SigTemp:=(cc.top[,P]<sig_cc)*1]
cc.top[SigTemp==1,Sig:="***"]
cc.top[SigTemp==0,Sig:=""]
cc.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(cc,SNPID)
setkey(cc.top,SNPID)

cc<-info[cc]

cc.top<-info[cc.top]


## Make QQ Plot

png('CC_A.png')

QQ(cc[,list(P)], 'Case-control: GxE',maxlimqq.cc,PlotAll.cc)

dev.off()

## Make Manhattan Plot

png('CC_B.png')

manhattan(cc[,list(P,CHR,BP)], sig_cc, maxlimman.cc,0.7,PlotAll.cc,1)
title('Case-control: GxE')

dev.off()

setnames(cc.top,c("NMISS","BETA","Z","P"),c("N","Beta_GxE","tTest","Pvalue"))

} # end of "if(RunCC==1)"

## If CC is not implemented but any of the 2-step needs to be implemented

if(RunCC!=1 & (RunEDGxE==1 | RunDGGxE==1 | RunEGGxE==1 | RunCocktail==1  )){
setkey(cc,SNPID)
cc<-info[cc]

}


## Change the column names in CC to merge with step 1 later
setnames(cc,c("NMISS","BETA","Z","P"),c("N","CC_Beta","CC_STAT","CC_P"))

## EG+DG|GxE is implemented

## EG+DG ##
if( RunEDGxE==1 ){

print("EDGxE")

## Read in EG+DG
edge<-ReadMplink("CC_DGGE.gxeout")

wt(edge,alpha,EDGxE.bin)

setkey(edge,SNPID)


## Merge with GxE


setnames(edge,c("CHISQ_2DF","P"),c("St1_Chisq","Step1_P"))

edge<-cc[edge,nomatch=0]

M.edge<-nrow(edge)

if(M.edge==0) {
stop("All SNPs are not converged from EDGxE")
}

if(M.edge<1e5){
PlotAll.edge=1
} else {
PlotAll.edge=0
}

#if(M.edge<M.info) print("NA found in EDGxE gxescan output file")

## Sort EDGxE by step1 P-value
setkey(edge,Step1_P)

## Subset to pick the E-G top hits into step 2 (subset testing)
edge.s2<- edge[Step1_P<EDGxEscreen]
edge.s2[,c("Bin","Threshold"):=NULL]
edge[,SNPID:=NULL]
edge.s2.size<-nrow(edge.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
edge[,SigTemp:=(edge[,CC_P]<edge[,Threshold])*1]
edge[SigTemp==1,Sig:="***"]
edge[SigTemp==0,Sig:=""]
edge[,SigTemp:=NULL]
edge[,Rank:=1:M.edge]

edge.wt.hits<-edge[Sig=="***"]
edge.wt.top<-edge[1:topnum]

edge.wt.top<-rbind(edge.wt.hits[!edge.wt.hits$Rank %in% edge.wt.top$Rank,],edge.wt.top)

rm(edge.wt.hits)


## Compute limit for step 1 plots
edge.min <- min(edge[,Step1_P])
maxlimqq.edge <- ceiling(1-log10(edge.min))

## Step1 QQ Plot

png('EDGE_A.png')
QQ(edge[,list(Step1_P)], 'EDGxE: Step1 screen',maxlimqq.edge,PlotAll.edge)
dev.off()

## Step1 Manhattan Plot 
png('EDGE_B.png')
manhattan(edge[,list(Step1_P, CHR, BP)],0,maxlimqq.edge,0.7,PlotAll.edge,0)
title('EDGxE: Step1 screen')
dev.off()

edge[,c("CHR","SNP","BP","A1","N","CC_Beta","CC_STAT","St1_Chisq","Step1_P","Sig"):=NULL]


if (edge.s2.size>0) {        ## If we have SNPs pass to step2
edge.min<-min(edge.s2[,CC_P])
sig_edge_st2 <- alpha/edge.s2.size
edge.min2<-min(edge.min,sig_edge_st2)
maxlimqq.edge <- ceiling(1-log10(edge.min))
maxlimman.edge <- ceiling(1-log10(edge.min2))

if(edge.s2.size>1e5){
PlotAll.edge.stp2=0
} else {
PlotAll.edge.stp2=1
}


## Step2 QQ Plot
png('EDGE_C.png')
QQ(edge.s2[,list(CC_P)], 'EDGxE: Step 2 test',maxlimqq.edge,PlotAll.edge.stp2)
dev.off()

## Step2 Manhattan Plot
png('EDGE_D.png')
manhattan(edge.s2[,list(CC_P, CHR, BP)], sig_edge_st2, maxlimman.edge,0.7,PlotAll.edge.stp2,1)
title('EDGxE: Step 2 test')
dev.off()


## Compute significance indicator for subset testing
setkey(edge.s2,CC_P)
edge.s2[,SigTemp:=(edge.s2[,CC_P]<sig_edge_st2)*1]
edge.s2[SigTemp==1,Sig:="***"]
edge.s2[SigTemp==0,Sig:=""]
edge.s2[,SigTemp:=NULL]

## Create EDGE rank if "IncludeRanks==1"
if(IncludeRanks==1){
edge.rank<-edge.s2[,list(SNPID)]
edge.rank[,EDGxE:=1:edge.s2.size]
edge.s2<-edge.s2[1:topnum]
}
} else {

sig_edge_st2<-NA

}

## Weighted testing
setkey(edge,Bin)
# Compute the limit for plot
edge.last<-max(edge[,Bin])
edge.k<-EDGxE.bin*(2^(edge.last-1))
edge.wt.binsig<-alpha*((1/2)^edge.last)
edge.wt.sig<-edge.wt.binsig/edge.k
edge.wt.lnum<-nrow(edge[J(edge.last)])
edge.wt.ref<- rep(-1*log10(edge.wt.sig),edge.wt.lnum)

# Make plot
png('EDGE_E.png')

wtedge.ref<-max(ceiling(1-log10(min(edge[,CC_P]))),edge.wt.ref[1])
wtplot(edge,wtedge.ref+1,"EDGxE: weighted",0.7,edge.wt.ref,edge.last,PlotAll.edge)

dev.off()

rm(edge)
rm(edge.wt.ref)

setnames(edge.s2,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))
setnames(edge.wt.top,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))

} # End of RunEDGxE


## EG|GxE is implemented

if( RunEGGxE==1 ){
print("EG|GxE")
setkey(eg,Step1_P)
wt(eg,alpha,EGGxE.bin)

setkey(eg,SNPID)

## Merge with GxE

egGxE<-cc[eg,nomatch=0]

M.egGxE<-nrow(egGxE)

if(M.egGxE==0) {
stop("All SNPs are not converged from EG|GxE")
}

if(M.egGxE<1e5){
PlotAll.egGxE=1
} else {
PlotAll.egGxE=0
}

#if(M.egGxE<M.info) print("NA found in gxescan output file")

## Sort EDGxE by step1 P-value
setkey(egGxE,Step1_P)

## Subset to pick the E-G top hits into step 2 (subset testing)
egGxE.s2<- egGxE[Step1_P<EGscreen]
egGxE.s2[,c("Bin","Threshold"):=NULL]
egGxE[,SNPID:=NULL]
egGxE.s2.size<-nrow(egGxE.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
egGxE[,SigTemp:=(egGxE[,CC_P]<egGxE[,Threshold])*1]
egGxE[SigTemp==1,Sig:="***"]
egGxE[SigTemp==0,Sig:=""]
egGxE[,SigTemp:=NULL]
egGxE[,Rank:=1:M.egGxE]

egGxE.wt.hits<-egGxE[Sig=="***"]
egGxE.wt.top<-egGxE[1:topnum]

egGxE.wt.top<-rbind(egGxE.wt.hits[!egGxE.wt.hits$Rank %in% egGxE.wt.top$Rank,],egGxE.wt.top)

rm(egGxE.wt.hits)


if(RunEG2df!=1){
## Compute limit for step 1 plots
eg.min <- min(egGxE[,Step1_P])
maxlimqq.eg <- ceiling(1-log10(eg.min))
}

## Step1 QQ Plot ( will be used for both EG|GxE and EG|2df )

png('EGGxE_A.png')
QQ(egGxE[,list(Step1_P)], 'EG | GxE: Step1 screen',maxlimqq.eg,PlotAll.egGxE)
dev.off()

## Step1 Manhattan Plot ( will be used for both EG|GxE and EG|2df )
png('EGGxE_B.png')
manhattan(egGxE[,list(Step1_P, CHR, BP)],0,maxlimqq.eg,0.7,PlotAll.egGxE,0)
title('EG | GxE: Step1 screen')
dev.off()


egGxE[,c("CHR","SNP","BP","A1","N","CC_Beta","CC_STAT","St1_tTest","Step1_P","Sig"):=NULL]


if (egGxE.s2.size>0) {        ## If we have SNPs pass to step2
egGxE.min<-min(egGxE.s2[,CC_P])
sig_egGxE_st2 <- alpha/egGxE.s2.size
egGxE.min2<-min(egGxE.min,sig_egGxE_st2)
maxlimqq.egGxE <- ceiling(1-log10(egGxE.min))
maxlimman.egGxE<- ceiling(1-log10(egGxE.min2))

if(egGxE.s2.size>1e5){
PlotAll.egGxE.stp2=0
} else {
PlotAll.egGxE.stp2=1
}

## Step2 QQ Plot
png('EGGxE_C.png')
QQ(egGxE.s2[,list(CC_P)], 'EG | GxE: Step 2 test',maxlimqq.egGxE,PlotAll.egGxE.stp2)
dev.off()

## Step2 Manhattan Plot
png('EGGxE_D.png')
manhattan(egGxE.s2[,list(CC_P, CHR, BP)], sig_egGxE_st2, maxlimman.egGxE,0.7,PlotAll.egGxE.stp2,1)
title('EG | GxE: Step 2 test')
dev.off()


## Compute significance indicator for subset testing
setkey(egGxE.s2,CC_P)
egGxE.s2[,SigTemp:=(egGxE.s2[,CC_P]<sig_egGxE_st2)*1]
egGxE.s2[SigTemp==1,Sig:="***"]
egGxE.s2[SigTemp==0,Sig:=""]
egGxE.s2[,SigTemp:=NULL]

## Create EG|GxE rank if "IncludeRanks==1"
if(IncludeRanks==1){
egGxE.rank<-egGxE.s2[,list(SNPID)]
egGxE.rank[,EGGxE:=1:egGxE.s2.size]
egGxE.s2<-egGxE.s2[1:topnum]
}
} else {

sig_egGxE_st2<-NA

}



## Weighted testing
setkey(egGxE,Bin)
# Compute the limit for plot
egGxE.last<-max(egGxE[,Bin])
egGxE.k<-EGGxE.bin*(2^(egGxE.last-1))
egGxE.wt.binsig<-alpha*((1/2)^egGxE.last)
egGxE.wt.sig<-egGxE.wt.binsig/egGxE.k
egGxE.wt.lnum<-nrow(egGxE[J(egGxE.last)])
egGxE.wt.ref<- rep(-1*log10(egGxE.wt.sig),egGxE.wt.lnum)

# Make plot
png('EGGxE_E.png')

wtegGxE.ref<-max(ceiling(1-log10(min(egGxE[,CC_P]))),egGxE.wt.ref[1])
wtplot(egGxE,wtegGxE.ref+1,"EG | GxE: weighted",0.7,egGxE.wt.ref,egGxE.last,PlotAll.egGxE)

dev.off()

rm(egGxE)
rm(egGxE.wt.ref)

setnames(egGxE.s2,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))
setnames(egGxE.wt.top,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))
egGxE.s2<-egGxE.s2[,-12]
} # End of RunEGGxE

if(RunCocktail!=1 & exists("eg")){
rm(eg)
} else if (RunCocktail==1 & exists("eg")){

eg[,c("St1_tTest","Bin","Threshold"):=NULL]

setnames(eg,"Step1_P","EG.P")

}

if(RunDGGxE!=1 & RunCocktail!=1){
rm(cc)
} 

} # End of "if( RunCC==1 | RunEDGxE==1 | RunDGGxE==1 |RunEGGxE==1 |RunCocktail==1)"



## D-G ##

if( RunDG==1 | RunDGGxE==1 | RunDGEB==1 |RunCocktail==1){

## Read in data
dg<-ReadMplink("CC_DG.gxeout")

if(RunDG==1){
print("Marginal G scan")
## Check whether or not we have significant hits
M.dg<-nrow(dg)
sig_dg <- alpha/M.dg

if(M.dg==0) {
stop("All SNPs are not converged from marginal G scan")
}

if(M.dg<1e5){
PlotAll.dg=1
} else {
PlotAll.dg=0
}

#if(M.dg<M.info) print("NA found in marginal G (DG) gxescan output file")

### find the axis limit for plots
dgp.min <- min(dg[,P])
maxlimqq.dg <- ceiling(1-log10(dgp.min))
dg.minman<-min(dgp.min,sig_dg)
maxlimman.dg <-ceiling(1-log10(dg.minman))

## Create DG rank if "IncludeRanks==1"

if(IncludeRanks==1){
dg.rank<-dg[,list(SNPID)]
dg.rank[,DG:=1:M.dg] ## Since SNPs are sorted by P-value in the original gxescan file 

}


## Subset to pick the top hits

dg.top<- dg[1:topnum]

## Specify significance indicator for the top hits

dg.top[,SigTemp:=(dg.top[,P]<sig_dg)*1]
dg.top[SigTemp==1,Sig:="***"]
dg.top[SigTemp==0,Sig:=""]
dg.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(dg,SNPID)
setkey(dg.top,SNPID)

dg<-info[dg]

dg.top<-info[dg.top]


## Make QQ Plot   (This will also be used for DG|GxE and DG|EB

png('DG_A.png')

QQ(dg[,list(P)], 'MarG: Marginal G Scan',maxlimqq.dg,PlotAll.dg)

dev.off()

## Make Manhattan Plot

png('DG_B.png')

manhattan(dg[,list(P,CHR,BP)], sig_dg, maxlimman.dg,0.7,PlotAll.dg,1)
title('MarG: Marginal G Scan')

dev.off()

setnames(dg.top,c("NMISS","BETA","Z","P"),c("N","Beta","tTest","Pvalue"))

setkey(dg,P)

dg[,c("CHR","SNP","BP","A1","NMISS","BETA"):=NULL]

} # end of "if(RunDG==1)"
 
if(RunDG!=1 & (RunDGGxE==1 | RunDGEB==1 |RunCocktail==1)){
dg[,c("NMISS","BETA"):=NULL]
}

setnames(dg,c("Z","P"),c("St1_tTest","Step1_P"))

## DG|GxE is implemented

if( RunDGGxE==1 ){
print("DG|GxE")
wt(dg,alpha,DGGxE.bin)

setkey(dg,SNPID)

## Merge with GxE

dgGxE<-cc[dg,nomatch=0]

M.dgGxE<-nrow(dgGxE)

if(M.dgGxE==0) {
stop("All SNPs are not converged from DG|GxE")
}

if(M.dgGxE<1e5){
PlotAll.dgGxE=1
} else {
PlotAll.dgGxE=0
}

#if(M.dgGxE<M.info) print("NA found in gxescan output file")

## Sort DG|GxE by step1 P-value
setkey(dgGxE,Step1_P)

## Subset to pick the D-G top hits into step 2 (subset testing)
dgGxE.s2<- dgGxE[Step1_P<DGscreen]
dgGxE.s2[,c("Bin","Threshold"):=NULL]
dgGxE[,SNPID:=NULL]
dgGxE.s2.size<-nrow(dgGxE.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
dgGxE[,SigTemp:=(dgGxE[,CC_P]<dgGxE[,Threshold])*1]
dgGxE[SigTemp==1,Sig:="***"]
dgGxE[SigTemp==0,Sig:=""]
dgGxE[,SigTemp:=NULL]
dgGxE[,Rank:=1:M.dgGxE]

dgGxE.wt.hits<-dgGxE[Sig=="***"]
dgGxE.wt.top<-dgGxE[1:topnum]

dgGxE.wt.top<-rbind(dgGxE.wt.hits[!dgGxE.wt.hits$Rank %in% dgGxE.wt.top$Rank,],dgGxE.wt.top)

rm(dgGxE.wt.hits)


if(RunDG!=1){
## Compute limit for step 1 plots
dgp.min <- min(dgGxE[,Step1_P])
maxlimqq.dg <- ceiling(1-log10(dgp.min))
}

## Step1 QQ Plot ( will be used for both EG|GxE and EG|2df )

png('DGGxE_A.png')
QQ(dgGxE[,list(Step1_P)], 'DG | GxE: Step1 screen',maxlimqq.dg,PlotAll.dgGxE)
dev.off()

## Step1 Manhattan Plot ( will be used for both EG|GxE and EG|2df )
png('DGGxE_B.png')
manhattan(dgGxE[,list(Step1_P, CHR, BP)],0,maxlimqq.dg,0.7,PlotAll.dgGxE,0)
title('DG | GxE: Step1 screen')
dev.off()


dgGxE[,c("CHR","SNP","BP","A1","N","CC_Beta","CC_STAT","St1_tTest","Step1_P","Sig"):=NULL]


if (dgGxE.s2.size>0) {        ## If we have SNPs pass to step2
dgGxE.min<-min(dgGxE.s2[,CC_P])
sig_dgGxE_st2 <- alpha/dgGxE.s2.size
dgGxE.min2<-min(dgGxE.min,sig_dgGxE_st2)
maxlimqq.dgGxE <- ceiling(1-log10(dgGxE.min))
maxlimman.dgGxE <- ceiling(1-log10(dgGxE.min2))

if(dgGxE.s2.size>1e5){
PlotAll.dgGxE.stp2=0
} else {
PlotAll.dgGxE.stp2=1
}

## Step2 QQ Plot
png('DGGxE_C.png')
QQ(dgGxE.s2[,list(CC_P)], 'DG | GxE: Step 2 test',maxlimqq.dgGxE,PlotAll.dgGxE.stp2)
dev.off()

## Step2 Manhattan Plot
png('DGGxE_D.png')
manhattan(dgGxE.s2[,list(CC_P, CHR, BP)], sig_dgGxE_st2, maxlimman.dgGxE,0.7,PlotAll.dgGxE.stp2,1)
title('DG | GxE: Step 2 test')
dev.off()


## Compute significance indicator for subset testing
setkey(dgGxE.s2,CC_P)
dgGxE.s2[,SigTemp:=(dgGxE.s2[,CC_P]<sig_dgGxE_st2)*1]
dgGxE.s2[SigTemp==1,Sig:="***"]
dgGxE.s2[SigTemp==0,Sig:=""]
dgGxE.s2[,SigTemp:=NULL]

## Create DG|GxE rank if "IncludeRanks==1"
if(IncludeRanks==1){
dgGxE.rank<-dgGxE.s2[,list(SNPID)]
dgGxE.rank[,DGGxE:=1:dgGxE.s2.size]
dgGxE.s2<-dgGxE.s2[1:topnum]
}
} else {

sig_dgGxE_st2<-NA

}


## Weighted testing
setkey(dgGxE,Bin)
# Compute the limit for plot
dgGxE.last<-max(dgGxE[,Bin])
dgGxE.k<-DGGxE.bin*(2^(dgGxE.last-1))
dgGxE.wt.binsig<-alpha*((1/2)^dgGxE.last)
dgGxE.wt.sig<-dgGxE.wt.binsig/dgGxE.k
dgGxE.wt.lnum<-nrow(dgGxE[J(dgGxE.last)])
dgGxE.wt.ref<- rep(-1*log10(dgGxE.wt.sig),dgGxE.wt.lnum)

# Make plot
png('DGGxE_E.png')

wtdgGxE.ref<-max(ceiling(1-log10(min(dgGxE[,CC_P]))),dgGxE.wt.ref[1])
wtplot(dgGxE,wtdgGxE.ref+1,"DG | GxE: weighted",0.7,dgGxE.wt.ref,dgGxE.last,PlotAll.dgGxE)

dev.off()

rm(dgGxE)
rm(dgGxE.wt.ref)

setnames(dgGxE.s2,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))
setnames(dgGxE.wt.top,c("CC_Beta","CC_STAT","CC_P"),c("Beta_GxE","St2_tTest","Step2_P"))

dg[,c("Bin","Threshold"):=NULL]
} # End of RunDGGxE

if(RunDGEB!=1 & RunCocktail!=1){
rm(dg)
}

if(RunCocktail!=1 & exists("cc")){
rm(cc)
} else if(RunCocktail==1 & exists("cc")){
cc[,c("CHR","SNP","BP","A1","N","CC_Beta","CC_STAT"):=NULL]
}


} # End of if( RunDG==1 | RunDGGxE==1 | RunDGEB==1 |RunCocktail==1)




## Empirical Bayes (EB) ##

if( RunEB==1 | RunDGEB==1 |RunCocktail==1){

## Read in data
eb<-ReadMplink("EB.gxeout")

setnames(eb,"NMISS","N")

if(RunEB==1){
print("Exhaustive EB")
## Check whether or not we have significant hits
M.eb<-nrow(eb)
sig_eb <- alpha/M.eb

if(M.eb==0) {
stop("All SNPs are not converged from Empirical Bayes analyses")
}

if(M.eb<1e5){
PlotAll.eb=1
} else {
PlotAll.eb=0
}

#if(M.eb<M.info) print("NA found in Empirical Bayes (EB) gxescan output file")

### find the axis limit for plots
ebp.min <- min(eb[,P])
maxlimqq.eb <- ceiling(1-log10(ebp.min))
eb.minman<-min(ebp.min,sig_eb)
maxlimman.eb <-ceiling(1-log10(eb.minman))

## Create EB rank if "IncludeRanks==1"

if(IncludeRanks==1){
eb.rank<-eb[,list(SNPID)]
eb.rank[,EB:=1:M.eb] ## Since SNPs are sorted by P-value in the original gxescan file 

}


## Subset to pick the top hits

eb.top<- eb[1:topnum]

## Specify significance indicator for the top hits

eb.top[,SigTemp:=(eb.top[,P]<sig_eb)*1]
eb.top[SigTemp==1,Sig:="***"]
eb.top[SigTemp==0,Sig:=""]
eb.top[,SigTemp:=NULL]


## Merge with SNP info (need to sort by line first)

setkey(eb,SNPID)
setkey(eb.top,SNPID)

eb<-info[eb]

eb.top<-info[eb.top]

## Make QQ Plot

png('EB_A.png')

QQ(eb[,list(P)], 'EB: Empirical Bayes',maxlimqq.eb,PlotAll.eb)

dev.off()

## Make Manhattan Plot

png('EB_B.png')

manhattan(eb[,list(P,CHR,BP)], sig_eb, maxlimman.eb,0.7,PlotAll.eb,1)
title('EB: Empirical Bayes')

dev.off()

setnames(eb.top,c("BETA","STAT","P"),c("Beta_GxE","tTest","Pvalue"))

} # end of "if(RunEB==1)"

## If EB is not implemented but any of the 2-step needs to be implemented

if(RunEB!=1 & (RunDGEB==1 | RunCocktail==1  )){
setkey(eb,SNPID)
eb<-info[eb]
}


## DG|EB is implemented

if( RunDGEB==1 ){
print("DG|EB")
setkey(dg,Step1_P)
wt(dg,alpha,DGEB.bin)

setkey(dg,SNPID)

## Merge with EB

dgEB<-dg[eb,nomatch=0]
rm(dg)
rm(eb)

setnames(dgEB,c("BETA","STAT","P"),c("Beta_GxE","St2_tTest","Step2_P"))

M.dgEB<-nrow(dgEB)

if(M.dgEB==0) {
stop("All SNPs are not converged from DG|EB")
}

if(M.dgEB<1e5){
PlotAll.dgEB=1
} else {
PlotAll.dgEB=0
}

#if(M.dgEB<M.info) print("NA found in gxescan output file")

## Sort DG|GxE by step1 P-value
setkey(dgEB,Step1_P)

## Subset to pick the D-G top hits into step 2 (subset testing)
dgEB.s2<- dgEB[Step1_P<DGEBscreen]
dgEB.s2[,c("Bin","Threshold"):=NULL]
dgEB.s2.size<-nrow(dgEB.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
dgEB[,SigTemp:=(dgEB[,Step2_P]<dgEB[,Threshold])*1]
dgEB[SigTemp==1,Sig:="***"]
dgEB[SigTemp==0,Sig:=""]
dgEB[,SigTemp:=NULL]
dgEB[,Rank:=1:M.dgEB]

dgEB.wt.hits<-dgEB[Sig=="***"]
dgEB.wt.top<-dgEB[1:topnum]
dgEB.wt.top<-rbind(dgEB.wt.hits[!dgEB.wt.hits$Rank %in% dgEB.wt.top$Rank,],dgEB.wt.top)

rm(dgEB.wt.hits)
dgEB.wt.top[,SNPID:=NULL]

if(RunDG!=1){
## Compute limit for step 1 plots
dgp.min <- min(dgEB[,Step1_P])
maxlimqq.dg <- ceiling(1-log10(dgp.min))
}

## Step1 QQ Plot ( will be used for both EG|GxE and EG|2df )

png('DGEB_A.png')
QQ(dgEB[,list(Step1_P)], 'DG | EB: Step1 screen',maxlimqq.dg,PlotAll.dgEB)
dev.off()

## Step1 Manhattan Plot ( will be used for both EG|GxE and EG|2df )
png('DGEB_B.png')
manhattan(dgEB[,list(Step1_P, CHR, BP)],0,maxlimqq.dg,0.7,PlotAll.dgEB,0)
title('DG | EB: Step1 screen')
dev.off()



if (dgEB.s2.size>0) {        ## If we have SNPs pass to step2
dgEB.min<-min(dgEB.s2[,Step2_P])
sig_dgEB_st2 <- alpha/dgEB.s2.size
dgEB.min2<-min(dgEB.min,sig_dgEB_st2)
maxlimqq.dgEB <- ceiling(1-log10(dgEB.min))
maxlimman.dgEB <- ceiling(1-log10(dgEB.min2))

if(dgEB.s2.size>1e5){
PlotAll.dgEB.stp2=0
} else {
PlotAll.dgEB.stp2=1
}

## Step2 QQ Plot
png('DGEB_C.png')
QQ(dgEB.s2[,list(Step2_P)], 'DG | EB: Step 2 test',maxlimqq.dgEB,PlotAll.dgEB.stp2)
dev.off()

## Step2 Manhattan Plot
png('DGEB_D.png')
manhattan(dgEB.s2[,list(Step2_P, CHR, BP)], sig_dgEB_st2, maxlimman.dgEB,0.7,PlotAll.dgEB.stp2,1)
title('DG | EB: Step 2 test')
dev.off()


## Compute significance indicator for subset testing
setkey(dgEB.s2,Step2_P)
dgEB.s2[,SigTemp:=(dgEB.s2[,Step2_P]<sig_dgEB_st2)*1]
dgEB.s2[SigTemp==1,Sig:="***"]
dgEB.s2[SigTemp==0,Sig:=""]
dgEB.s2[,SigTemp:=NULL]

## Create DGEB rank if "IncludeRanks==1"
if(IncludeRanks==1){
dgEB.rank<-dgEB.s2[,list(SNPID)]
dgEB.rank[,DGEB:=1:dgEB.s2.size]
dgEB.s2<-dgEB.s2[1:topnum]
}
} else {

sig_dgEB_st2<-NA

}

## Weighted testing
setkey(dgEB,Bin)
# Compute the limit for plot
dgEB.last<-max(dgEB[,Bin])
dgEB.k<-DGEB.bin*(2^(dgEB.last-1))
dgEB.wt.binsig<-alpha*((1/2)^dgEB.last)
dgEB.wt.sig<-dgEB.wt.binsig/dgEB.k
dgEB.wt.lnum<-nrow(dgEB[J(dgEB.last)])
dgEB.wt.ref<- rep(-1*log10(dgEB.wt.sig),dgEB.wt.lnum)

# Make plot
png('DGEB_E.png')

wtdgEB.ref<-max(ceiling(1-log10(min(dgEB[,Step2_P]))),dgEB.wt.ref[1])
wtplot(dgEB[,list(Step2_P,Bin,Threshold,Rank)],wtdgEB.ref+1,"DG | EB: weighted",0.7,dgEB.wt.ref,dgEB.last,PlotAll.dgEB)

dev.off()

rm(dgEB.wt.ref)

} # End of RunDGEB

} # End of if( RunEB==1 | RunDGEB==1 |RunCocktail==1)

if(IncludeRanks==1){
## Rank contains ranks from all methods implemented
Rank<-data.table(info[,SNPID])
setnames(Rank,1,"SNPID")
setkey(Rank,SNPID)
}

#rm(info)


if( RunDGEB==1 & RunCocktail!=1 ){
rm(dgEB)
} else if ( RunDGEB==1 & RunCocktail==1) {
dgEB[,c("St1_tTest","Bin","Threshold","Beta_GxE","St2_tTest","Sig","Rank"):=NULL]
setnames(dgEB,c("Step1_P","Step2_P"),c("DG.P","EB.P"))
} else if(RunDGEB!=1 & RunCocktail==1){
setkey(dg,SNPID)
dgEB<-dg[eb,nomatch=0]
rm(dg)
rm(eb)
dgEB[,c("St1_tTest","BETA","STAT"):=NULL]
setnames(dgEB,c("Step1_P","P"),c("DG.P","EB.P"))
} 


## Cocktail

if( RunCocktail==1 ){
print("Cocktail")
## Merge with cc and eg
setkey(dgEB,SNPID)
cocktail<-dgEB[cc,nomatch=0]
rm(cc)
rm(dgEB)
cocktail<-cocktail[eg,nomatch=0]
rm(eg)

## Fit Cocktail

cocktail<-FitCocktail(cocktail[,list(SNPID,CHR,SNP,BP,A1,N)],cocktail[,DG.P],cocktail[,EG.P],cocktail[,EB.P],cocktail[,CC_P],cocktail.c)
setkey(cocktail,Step1_P)

wt(cocktail,alpha,Cocktail.bin)

M.cocktail<-nrow(cocktail)

if(M.cocktail==0) {
stop("All SNPs are not converged from Cocktail")
}

if(M.cocktail<1e5){
PlotAll.cocktail=1
} else {
PlotAll.cocktail=0
}

#if(M.cocktail<M.info) print("NA found in gxescan output file")
## Subset to pick the cocktail top hits into step 2 (subset testing)
cocktail.s2<- cocktail[Step1_P<Cocktailscreen]
cocktail.s2[,c("Bin","Threshold"):=NULL]
cocktail[,SNPID:=NULL]
cocktail.s2.size<-nrow(cocktail.s2)

## Subset to pick the top hits for weighted testing(weighted testing)
cocktail[,SigTemp:=(cocktail[,Step2_P]<cocktail[,Threshold])*1]
cocktail[SigTemp==1,Sig:="***"]
cocktail[SigTemp==0,Sig:=""]
cocktail[,SigTemp:=NULL]
cocktail[,Rank:=1:M.cocktail]

cocktail.wt.hits<-cocktail[Sig=="***"]
cocktail.wt.top<-cocktail[1:topnum]

cocktail.wt.top<-rbind(cocktail.wt.hits[!cocktail.wt.hits$Rank %in% cocktail.wt.top$Rank,],cocktail.wt.top)

rm(cocktail.wt.hits)


## Compute limit for step 1 plots
cocktail.min <- min(cocktail[,Step1_P])
maxlimqq.cocktail <- ceiling(1-log10(cocktail.min))

## Step1 QQ Plot

png('Cocktail_A.png')
QQ(cocktail[,list(Step1_P)], 'Cocktail I: Step1 screen',maxlimqq.cocktail,PlotAll.cocktail)
dev.off()

## Step1 Manhattan Plot 
png('Cocktail_B.png')
manhattan(cocktail[,list(Step1_P, CHR, BP)],0,maxlimqq.cocktail,0.7,PlotAll.cocktail,0)
title('Cocktail I: Step1 screen')
dev.off()

cocktail[,c("CHR","SNP","BP","A1","N","Step1_P","Method","Sig"):=NULL]


if (cocktail.s2.size>0) {        ## If we have SNPs pass to step2
cocktail.min<-min(cocktail.s2[,Step2_P])
sig_cocktail_st2 <- alpha/cocktail.s2.size
cocktail.min2<-min(cocktail.min,sig_cocktail_st2)
maxlimqq.cocktail <- ceiling(1-log10(cocktail.min))
maxlimman.cocktail <- ceiling(1-log10(cocktail.min2))

if(cocktail.s2.size>1e5){
PlotAll.cocktail.stp2=0
} else {
PlotAll.cocktail.stp2=1
}

## Step2 QQ Plot
png('Cocktail_C.png')
QQ(cocktail.s2[,list(Step2_P)], 'Cocktail I: Step1 screen',maxlimqq.cocktail,PlotAll.cocktail.stp2)
dev.off()

## Step2 Manhattan Plot
png('Cocktail_D.png')
manhattan(cocktail.s2[,list(Step2_P, CHR, BP)], sig_cocktail_st2, maxlimman.cocktail,0.7,PlotAll.cocktail.stp2,1)
title('Cocktail I: Step1 screen')
dev.off()


## Compute significance indicator for subset testing
setkey(cocktail.s2,Step2_P)
cocktail.s2[,SigTemp:=(cocktail.s2[,Step2_P]<sig_cocktail_st2)*1]
cocktail.s2[SigTemp==1,Sig:="***"]
cocktail.s2[SigTemp==0,Sig:=""]
cocktail.s2[,SigTemp:=NULL]

## Create cocktail rank if "IncludeRanks==1"
if(IncludeRanks==1){
cocktail.rank<-cocktail.s2[,list(SNPID)]
cocktail.rank[,Cocktail:=1:cocktail.s2.size]
cocktail.s2<-cocktail.s2[1:topnum]
}
} else {

sig_cocktail_st2<-NA

}

## Weighted testing
setkey(cocktail,Bin)
# Compute the limit for plot
cocktail.last<-max(cocktail[,Bin])
cocktail.k<-Cocktail.bin*(2^(cocktail.last-1))
cocktail.wt.binsig<-alpha*((1/2)^cocktail.last)
cocktail.wt.sig<-cocktail.wt.binsig/cocktail.k
cocktail.wt.lnum<-nrow(cocktail[J(cocktail.last)])
cocktail.wt.ref<- rep(-1*log10(cocktail.wt.sig),cocktail.wt.lnum)

# Make plot
png('Cocktail_E.png')

wtcocktail.ref<-max(ceiling(1-log10(min(cocktail[,Step2_P]))),cocktail.wt.ref[1])
wtplot(cocktail[,list(Step2_P,Bin,Threshold,Rank)],wtcocktail.ref+1,"Cocktail I: weighted",0.7,cocktail.wt.ref,cocktail.last,PlotAll.cocktail)

dev.off()

rm(cocktail)
rm(cocktail.wt.ref)

} # End of RunCocktail

rm(manhattan,QQ,wtplot)
##########################################
#
#               End of GWIS
#
##########################################
print("Done testing...formatting output")
## Format the output columns ##

## Cocktail

if (RunCocktail==1){
if(IncludeRanks==1){
setkey(cocktail.rank,SNPID)
Rank<-cocktail.rank[Rank]
rm(cocktail.rank)
}
cocktail.s2[,Step1_P:=signif(cocktail.s2[,Step1_P],2)]
cocktail.s2[,Step2_P:=signif(cocktail.s2[,Step2_P],2)]

cocktail.wt.top[,Step1_P:=signif(cocktail.wt.top[,Step1_P],2)]
cocktail.wt.top[,Step2_P:=signif(cocktail.wt.top[,Step2_P],2)]

}

## EG|2df

if (RunEG2df==1){
if(IncludeRanks==1){
setkey(eg2df.rank,SNPID)
Rank<-eg2df.rank[Rank]
rm(eg2df.rank)
}
eg2df.s2[,St1_tTest:=as.numeric(formatC(eg2df.s2[,St1_tTest],digits=2,format="f"))]
eg2df.s2[,St2_Chisq:=as.numeric(formatC(eg2df.s2[,St2_Chisq],digits=2,format="f"))]

eg2df.s2[,Step1_P:=signif(eg2df.s2[,Step1_P],2)]
eg2df.s2[,Step2_P:=signif(eg2df.s2[,Step2_P],2)]

eg2df.wt.top[,St1_tTest:=as.numeric(formatC(eg2df.wt.top[,St1_tTest],digits=2,format="f"))]
eg2df.wt.top[,St2_Chisq:=as.numeric(formatC(eg2df.wt.top[,St2_Chisq],digits=2,format="f"))]
eg2df.wt.top[,Step1_P:=signif(eg2df.wt.top[,Step1_P],2)]
eg2df.wt.top[,Step2_P:=signif(eg2df.wt.top[,Step2_P],2)]

}


## EG|GxE

if (RunEGGxE==1){
if(IncludeRanks==1){
setkey(egGxE.rank,SNPID)
Rank<-egGxE.rank[Rank]
rm(egGxE.rank)
}
egGxE.s2[,Beta_GxE:=as.numeric(formatC(egGxE.s2[,Beta_GxE],digits=2,format="f"))]
egGxE.s2[,St1_tTest:=as.numeric(formatC(egGxE.s2[,St1_tTest],digits=2,format="f"))]
egGxE.s2[,St2_tTest:=as.numeric(formatC(egGxE.s2[,St2_tTest],digits=2,format="f"))]
egGxE.s2[,Step1_P:=signif(egGxE.s2[,Step1_P],2)]
egGxE.s2[,Step2_P:=signif(egGxE.s2[,Step2_P],2)]

egGxE.wt.top[,Beta_GxE:=as.numeric(formatC(egGxE.wt.top[,Beta_GxE],digits=2,format="f"))]
egGxE.wt.top[,St1_tTest:=as.numeric(formatC(egGxE.wt.top[,St1_tTest],digits=2,format="f"))]
egGxE.wt.top[,St2_tTest:=as.numeric(formatC(egGxE.wt.top[,St2_tTest],digits=2,format="f"))]
egGxE.wt.top[,Step1_P:=signif(egGxE.wt.top[,Step1_P],2)]
egGxE.wt.top[,Step2_P:=signif(egGxE.wt.top[,Step2_P],2)]

}


## DG|EB

if (RunDGEB==1){
if(IncludeRanks==1){
setkey(dgEB.rank,SNPID)
Rank<-dgEB.rank[Rank]
rm(dgEB.rank)
}
dgEB.s2[,Beta_GxE:=as.numeric(formatC(dgEB.s2[,Beta_GxE],digits=2,format="f"))]
dgEB.s2[,St1_tTest:=as.numeric(formatC(dgEB.s2[,St1_tTest],digits=2,format="f"))]
dgEB.s2[,St2_tTest:=as.numeric(formatC(dgEB.s2[,St2_tTest],digits=2,format="f"))]
dgEB.s2[,Step1_P:=signif(dgEB.s2[,Step1_P],2)]
dgEB.s2[,Step2_P:=signif(dgEB.s2[,Step2_P],2)]

dgEB.wt.top[,Beta_GxE:=as.numeric(formatC(dgEB.wt.top[,Beta_GxE],digits=2,format="f"))]
dgEB.wt.top[,St1_tTest:=as.numeric(formatC(dgEB.wt.top[,St1_tTest],digits=2,format="f"))]
dgEB.wt.top[,St2_tTest:=as.numeric(formatC(dgEB.wt.top[,St2_tTest],digits=2,format="f"))]
dgEB.wt.top[,Step1_P:=signif(dgEB.wt.top[,Step1_P],2)]
dgEB.wt.top[,Step2_P:=signif(dgEB.wt.top[,Step2_P],2)]

}


## DG|GxE

if (RunDGGxE==1){
if(IncludeRanks==1){
setkey(dgGxE.rank,SNPID)
Rank<-dgGxE.rank[Rank]
rm(dgGxE.rank)
}
dgGxE.s2[,Beta_GxE:=as.numeric(formatC(dgGxE.s2[,Beta_GxE],digits=2,format="f"))]
dgGxE.s2[,St1_tTest:=as.numeric(formatC(dgGxE.s2[,St1_tTest],digits=2,format="f"))]
dgGxE.s2[,St2_tTest:=as.numeric(formatC(dgGxE.s2[,St2_tTest],digits=2,format="f"))]
dgGxE.s2[,Step1_P:=signif(dgGxE.s2[,Step1_P],2)]
dgGxE.s2[,Step2_P:=signif(dgGxE.s2[,Step2_P],2)]

dgGxE.wt.top[,Beta_GxE:=as.numeric(formatC(dgGxE.wt.top[,Beta_GxE],digits=2,format="f"))]
dgGxE.wt.top[,St1_tTest:=as.numeric(formatC(dgGxE.wt.top[,St1_tTest],digits=2,format="f"))]
dgGxE.wt.top[,St2_tTest:=as.numeric(formatC(dgGxE.wt.top[,St2_tTest],digits=2,format="f"))]
dgGxE.wt.top[,Step1_P:=signif(dgGxE.wt.top[,Step1_P],2)]
dgGxE.wt.top[,Step2_P:=signif(dgGxE.wt.top[,Step2_P],2)]

}


## EG+DG|GxE

if (RunEDGxE==1){
if(IncludeRanks==1){
setkey(edge.rank,SNPID)
Rank<-edge.rank[Rank]
rm(edge.rank)
}
edge.s2[,Beta_GxE:=as.numeric(formatC(edge.s2[,Beta_GxE],digits=2,format="f"))]
edge.s2[,St1_Chisq:=as.numeric(formatC(edge.s2[,St1_Chisq],digits=2,format="f"))]
edge.s2[,St2_tTest:=as.numeric(formatC(edge.s2[,St2_tTest],digits=2,format="f"))]
edge.s2[,Step1_P:=signif(edge.s2[,Step1_P],2)]
edge.s2[,Step2_P:=signif(edge.s2[,Step2_P],2)]

edge.wt.top[,Beta_GxE:=as.numeric(formatC(edge.wt.top[,Beta_GxE],digits=2,format="f"))]
edge.wt.top[,St1_Chisq:=as.numeric(formatC(edge.wt.top[,St1_Chisq],digits=2,format="f"))]
edge.wt.top[,St2_tTest:=as.numeric(formatC(edge.wt.top[,St2_tTest],digits=2,format="f"))]
edge.wt.top[,Step1_P:=signif(edge.wt.top[,Step1_P],2)]
edge.wt.top[,Step2_P:=signif(edge.wt.top[,Step2_P],2)]

}

## 3df
if (Run3df==1){
if(IncludeRanks==1){
setkey(df3.rank,SNPID)
Rank<-df3.rank[Rank]
rm(df3.rank)
}
df3.top[,df3_Chisq:= as.numeric(formatC(df3.top[,df3_Chisq],digits=2,format="f"))]
df3.top[,Pvalue:=signif(df3.top[,Pvalue],2)]
setnames(df3.top,"df3_Chisq","3df_Chisq")
}

## 2df
if (Run2df==1){
if(IncludeRanks==1){
setkey(df2.rank,SNPID)
Rank<-df2.rank[Rank]
rm(df2.rank)
}
df2.top[,df2_Chisq:= as.numeric(formatC(df2.top[,df2_Chisq],digits=2,format="f"))]
df2.top[,Pvalue:=signif(df2.top[,Pvalue],2)]
setnames(df2.top,"df2_Chisq","2df_Chisq")
}

## EB
if (RunEB==1){
if(IncludeRanks==1){
setkey(eb.rank,SNPID)
Rank<-eb.rank[Rank]
rm(eb.rank)
}
eb.top[,Beta_GxE:=as.numeric(formatC(eb.top[,Beta_GxE],digits=2,format="f"))]
eb.top[,tTest:= as.numeric(formatC(eb.top[,tTest],digits=2,format="f"))]
eb.top[,Pvalue:=signif(eb.top[,Pvalue],2)]

}

## CntlOnly
if (RunCNTL==1){
if(IncludeRanks==1){
setkey(cntl.rank,SNPID)
Rank<-cntl.rank[Rank]
rm(cntl.rank)
}
cntl.top[,Beta:=as.numeric(formatC(cntl.top[,Beta],digits=2,format="f"))]
cntl.top[,tTest:= as.numeric(formatC(cntl.top[,tTest],digits=2,format="f"))]
cntl.top[,Pvalue:=signif(cntl.top[,Pvalue],2)]
}


## CO
if (RunCO==1){
if(IncludeRanks==1){
setkey(co.rank,SNPID)
Rank<-co.rank[Rank]
rm(co.rank)
}
co.top[,Beta:=as.numeric(formatC(co.top[,Beta],digits=2,format="f"))]
co.top[,tTest:= as.numeric(formatC(co.top[,tTest],digits=2,format="f"))]
co.top[,Pvalue:=signif(co.top[,Pvalue],2)]

}


## CC
if (RunCC==1){
if(IncludeRanks==1){
setkey(cc.rank,SNPID)
Rank<-cc.rank[Rank]
rm(cc.rank)
}
cc.top[,Beta_GxE:=as.numeric(formatC(cc.top[,Beta_GxE],digits=2,format="f"))]
cc.top[,tTest:= as.numeric(formatC(cc.top[,tTest],digits=2,format="f"))]
cc.top[,Pvalue:=signif(cc.top[,Pvalue],2)]

}

## Mar G
if (RunDG==1){
if(IncludeRanks==1){
setkey(dg.rank,SNPID)
setnames(dg.rank,"DG","MarG")
Rank<-dg.rank[Rank]
rm(dg.rank)
}
dg.top[,Beta:=as.numeric(formatC(dg.top[,Beta],digits=2,format="f"))]
dg.top[,tTest:= as.numeric(formatC(dg.top[,tTest],digits=2,format="f"))]
dg.top[,Pvalue:=signif(dg.top[,Pvalue],2)]

}


## Add Ranks if checked
ImpMethod<- c("MarG","CC","CO","CntlOnly","EB","df2","df3","EDGxE","DGGxE","DGEB","EGGxE","EG2df","Cocktail")
ImpMethod<-ImpMethod[c(RunDG==1,RunCC==1,RunCO==1,RunCNTL==1,RunEB==1,Run2df==1,Run3df==1,RunEDGxE==1,RunDGGxE==1,RunDGEB==1,RunEGGxE==1,RunEG2df==1,RunCocktail==1)]

##MarG
if(RunDG==1){

if(IncludeRanks==1){
ImpMethod.marG<-ImpMethod[-which(ImpMethod=="MarG")]
dg.top<-Rank[dg.top]
dg.top[,SNPID:=NULL]
setcolorder(dg.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig","MarG",ImpMethod.marG))
setkey(dg.top,MarG)
} else {
dg.top[,SNPID:=NULL]
setcolorder(dg.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig"))
setkey(dg.top,Pvalue)
}


}


##CC
if(RunCC==1){
if(IncludeRanks==1){
ImpMethod.cc<-ImpMethod[-which(ImpMethod=="CC")]
cc.top<-Rank[cc.top]
cc.top[,SNPID:=NULL]
setcolorder(cc.top,c("SNP","CHR","BP","A1","N","Beta_GxE","tTest","Pvalue","Sig","CC",ImpMethod.cc))
setkey(cc.top,CC)
} else {
cc.top[,SNPID:=NULL]
setcolorder(cc.top,c("SNP","CHR","BP","A1","N","Beta_GxE","tTest","Pvalue","Sig"))
setkey(cc.top,Pvalue)
}

}


##CO
if(RunCO==1){
if(IncludeRanks==1){
ImpMethod.co<-ImpMethod[-which(ImpMethod=="CO")]
co.top<-Rank[co.top]
co.top[,SNPID:=NULL]
setcolorder(co.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig","CO",ImpMethod.co))
setkey(co.top,CO)
} else {
co.top[,SNPID:=NULL]
setcolorder(co.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig"))
setkey(co.top,Pvalue)
}
}

##CntlOnly
if(RunCNTL==1){
if(IncludeRanks==1){
ImpMethod.cntl<-ImpMethod[-which(ImpMethod=="CntlOnly")]
cntl.top<-Rank[cntl.top]
cntl.top[,SNPID:=NULL]
setcolorder(cntl.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig","CntlOnly",ImpMethod.cntl))
setkey(cntl.top,CntlOnly)
} else {
cntl.top[,SNPID:=NULL]
setcolorder(cntl.top,c("SNP","CHR","BP","A1","N","Beta","tTest","Pvalue","Sig"))
setkey(cntl.top,Pvalue)
}
}

## 2df
if(Run2df==1){
if(IncludeRanks==1){
ImpMethod.2df<-ImpMethod[-which(ImpMethod=="df2")]
df2.top<-Rank[df2.top]
df2.top[,SNPID:=NULL]
setcolorder(df2.top,c("SNP","CHR","BP","A1","N","2df_Chisq","Pvalue","Sig","df2",ImpMethod.2df))
setkey(df2.top,df2)
} else {
df2.top[,SNPID:=NULL]
setcolorder(df2.top,c("SNP","CHR","BP","A1","N","2df_Chisq","Pvalue","Sig"))
setkey(df2.top,Pvalue)
}
}


## 3df
if(Run3df==1){
if(IncludeRanks==1){
ImpMethod.3df<-ImpMethod[-which(ImpMethod=="df3")]
df3.top<-Rank[df3.top]
df3.top[,SNPID:=NULL]
setcolorder(df3.top,c("SNP","CHR","BP","A1","N","3df_Chisq","Pvalue","Sig","df3",ImpMethod.3df))
setkey(df3.top,df3)
} else {
df3.top[,SNPID:=NULL]
setcolorder(df3.top,c("SNP","CHR","BP","A1","N","3df_Chisq","Pvalue","Sig"))
setkey(df3.top,Pvalue)
}
}


## EB
if(RunEB==1){
if(IncludeRanks==1){
ImpMethod.eb<-ImpMethod[-which(ImpMethod=="EB")]
eb.top<-Rank[eb.top]
eb.top[,SNPID:=NULL]
setcolorder(eb.top,c("SNP","CHR","BP","A1","N","Beta_GxE","tTest","Pvalue","Sig","EB",ImpMethod.eb))
setkey(eb.top,EB)
} else {
eb.top[,SNPID:=NULL]
setcolorder(eb.top,c("SNP","CHR","BP","A1","N","Beta_GxE","tTest","Pvalue","Sig"))
setkey(eb.top,Pvalue)
}
}


## EDGxE

if(RunEDGxE==1){
if(IncludeRanks==1){
ImpMethod.edge<-ImpMethod[-which(ImpMethod=="EDGxE")]
setkey(edge.s2,SNPID)
edge.s2<-Rank[edge.s2]
edge.s2[,SNPID:=NULL]
setcolorder(edge.s2,c("SNP","CHR","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig","EDGxE",ImpMethod.edge))
setkey(edge.s2,EDGxE)
} else {
edge.s2[,SNPID:=NULL]
setcolorder(edge.s2,c("SNP","CHR","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig"))
setkey(edge.s2,Step2_P)
}


setcolorder(edge.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_tTest","Step2_P","Threshold","Bin","Sig"))

setkey(edge.wt.top,Step1_P)

}

## DG|GxE

if(RunDGGxE==1){
if(IncludeRanks==1){
ImpMethod.dgGxE<-ImpMethod[-which(ImpMethod=="DGGxE")]
setkey(dgGxE.s2,SNPID)
dgGxE.s2<-Rank[dgGxE.s2]
dgGxE.s2[,SNPID:=NULL]
setcolorder(dgGxE.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig","DGGxE",ImpMethod.dgGxE))
setkey(dgGxE.s2,DGGxE)
} else {
dgGxE.s2[,SNPID:=NULL]
setcolorder(dgGxE.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig"))
setkey(dgGxE.s2,Step2_P)
}

setcolorder(dgGxE.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Threshold","Bin","Sig"))

setkey(dgGxE.wt.top,Step1_P)

}


## DG|EB

if(RunDGEB==1){
if(IncludeRanks==1){
ImpMethod.dgEB<-ImpMethod[-which(ImpMethod=="DGEB")]
setkey(dgEB.s2,SNPID)
dgEB.s2<-Rank[dgEB.s2]
dgEB.s2[,SNPID:=NULL]
setcolorder(dgEB.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig","DGEB",ImpMethod.dgEB))
setkey(dgEB.s2,DGEB)
} else {
dgEB.s2[,SNPID:=NULL]
setcolorder(dgEB.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig"))
setkey(dgEB.s2,Step2_P)
}

setcolorder(dgEB.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Threshold","Bin","Sig"))

setkey(dgEB.wt.top,Step1_P)

}


egGxE.wt.top<-egGxE.wt.top[,-11]
## EG|GxE

if(RunEGGxE==1){
if(IncludeRanks==1){
ImpMethod.egGxE<-ImpMethod[-which(ImpMethod=="EGGxE")]
setkey(egGxE.s2,SNPID)
egGxE.s2<-Rank[egGxE.s2]
egGxE.s2[,SNPID:=NULL]
setcolorder(egGxE.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig","EGGxE",ImpMethod.egGxE))
setkey(egGxE.s2,EGGxE)
} else {
egGxE.s2[,SNPID:=NULL]
setcolorder(egGxE.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Sig"))
setkey(egGxE.s2,Step2_P)
}

setcolorder(egGxE.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_tTest","Step1_P","Beta_GxE","St2_tTest","Step2_P","Threshold","Bin","Sig"))

setkey(egGxE.wt.top,Step1_P)

}

eg2df.s2<-eg2df.s2[,-11]
eg2df.wt.top<-eg2df.wt.top[,-10]
## EG|2df

if(RunEG2df==1){
if(IncludeRanks==1){
ImpMethod.eg2df<-ImpMethod[-which(ImpMethod=="EG2df")]
setkey(eg2df.s2,SNPID)
eg2df.s2<-Rank[eg2df.s2]
eg2df.s2[,SNPID:=NULL]
setcolorder(eg2df.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","St2_Chisq","Step2_P","Sig","EG2df",ImpMethod.eg2df))
setkey(eg2df.s2,EG2df)
} else {
eg2df.s2[,SNPID:=NULL]
setcolorder(eg2df.s2,c("SNP","CHR","BP","A1","N","St1_tTest","Step1_P","St2_Chisq","Step2_P","Sig"))
setkey(eg2df.s2,Step2_P)
}

setcolorder(eg2df.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_tTest","Step1_P","St2_Chisq","Step2_P","Threshold","Bin","Sig"))

setkey(eg2df.wt.top,Step1_P)

}


## Cocktail

if(RunCocktail==1){
if(IncludeRanks==1){
ImpMethod.cocktail<-ImpMethod[-which(ImpMethod=="Cocktail")]
setkey(cocktail.s2,SNPID)
cocktail.s2<-Rank[cocktail.s2]
cocktail.s2[,SNPID:=NULL]
setcolorder(cocktail.s2,c("SNP","CHR","BP","A1","N","Step1_P","Method","Step2_P","Sig","Cocktail",ImpMethod.cocktail))
setkey(cocktail.s2,Cocktail)
} else {
cocktail.s2[,SNPID:=NULL]
setcolorder(cocktail.s2,c("SNP","CHR","BP","A1","N","Step1_P","Method","Step2_P","Sig"))
setkey(cocktail.s2,Step2_P)
}

setcolorder(cocktail.wt.top,c("Rank","CHR","SNP","BP","A1","N","Step1_P","Method","Step2_P","Threshold","Bin","Sig"))

setkey(cocktail.wt.top,Step1_P)

}

if(IncludeRanks==1){
rm(Rank)
}

## Write results to excel

FileName<-paste(mytest,"_output.xlsx",sep="")

# Delete file with the same name
id <- grep(FileName, dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

################################################## Save Test Summary ########################
print("Writing output to Excel file...")
num.method<- sum(RunDG,RunCC,RunCO,RunCNTL,RunEB,Run2df,Run3df,RunEDGxE,RunDGGxE,RunDGEB,RunEGGxE,RunEG2df,RunCocktail)
num.exaus<-sum(RunDG,RunCC,RunCO,RunCNTL,RunEB,Run2df,Run3df)
num.2step<-sum(RunEDGxE,RunDGGxE,RunDGEB,RunEGGxE,RunEG2df,RunCocktail)

diff.st1<- 7-num.exaus
diff.st2<- 13-num.method

## Create the excel file 
wb <- loadWorkbook(FileName, create = TRUE)
createSheet(wb, name = "summary")

setwidth( wb,"summary",1,2600)

##Exhaustive tests

## Total number of SNPs (M)
M <- c(NA)

if(RunDG==1){
M<-cbind(M,M.dg)
}

if(RunCC==1){
M<-cbind(M,M.cc)
}

if(RunCO==1){
M<-cbind(M,M.co)
}
if(RunCNTL==1){
M<-cbind(M,M.cntl)
}

if(RunEB==1){
M<-cbind(M,M.eb)
}

if(Run2df==1){
M<-cbind(M,M.2df)
}

if(Run3df==1){
M<-cbind(M,M.3df)
}

M<-M[-1]

M<-data.frame(M)

## Genomewide significance level

sig1 <- c(NA)

if(RunDG==1){
sig1<-cbind(sig1,sig_dg)
}

if(RunCC==1){
sig1<-cbind(sig1,sig_cc)
}

if(RunCO==1){
sig1<-cbind(sig1,sig_co)
}

if(RunCNTL==1){
sig1<-cbind(sig1,sig_cntl)
}

if(RunEB==1){
sig1<-cbind(sig1,sig_eb)
}

if(Run2df==1){
sig1<-cbind(sig1,sig_2df)
}

if(Run3df==1){
sig1<-cbind(sig1,sig_3df)
}

sig1<-sig1[-1]

sig1<-data.frame(sig1)


## alpha
alpha<-signif(alpha,3)
galpha <- rep(alpha,num.exaus)
galpha <- data.frame(galpha)
colnames(galpha) <- "alpha"

## Method names
test.exhaus <- c(NA)

if(RunDG==1){
test.exhaus.dg<- "MarG"
test.exhaus<-cbind(test.exhaus,test.exhaus.dg)
}


if(RunCC==1){
test.exhaus.cc<- "CC"
test.exhaus<-cbind(test.exhaus,test.exhaus.cc)
}

if(RunCO==1){
test.exhaus.ca<- "CO"
test.exhaus<-cbind(test.exhaus,test.exhaus.ca)
}

if(RunCNTL==1){
test.exhaus.con<- "CntlOnly"
test.exhaus<-cbind(test.exhaus,test.exhaus.con)
}

if(RunEB==1){
test.exhaus.eb<- "EB"
test.exhaus<-cbind(test.exhaus,test.exhaus.eb)
}

if(Run2df==1){
test.exhaus.2df<- "2df"
test.exhaus<-cbind(test.exhaus,test.exhaus.2df)
}

if(Run3df==1){
test.exhaus.3df<- "3df"
test.exhaus<-cbind(test.exhaus,test.exhaus.3df)
}

test.exhaus<-test.exhaus[-1]
test.exhaus<-data.frame(test.exhaus)
colnames(test.exhaus)<-"Method"


exhaustive<-data.frame(test.exhaus,M,galpha,sig1)
colnames(exhaustive)<-c("Method","M","alpha","alpha/M")

te0<-"Exhaustive tests"

wrap(wb,te0, "summary",1,1,F,F)
save2excel(wb,exhaustive, "summary",2,1,T)
rm(exhaustive)

######### 2-step tests
## Method names
te1<-"2-step tests"
if(num.2step>0){
test.2stp <- c(NA)

if(RunEDGxE==1){
test.2stp.egdg<- "EDGxE"
test.2stp<-cbind(test.2stp,test.2stp.egdg)
}

if(RunDGGxE==1){
test.2stp.dgGxE<- "DG|GxE"
test.2stp<-cbind(test.2stp,test.2stp.dgGxE)
}

if(RunDGEB==1){
test.2stp.dgeb<- "DG|EB"
test.2stp<-cbind(test.2stp,test.2stp.dgeb)
}

if(RunEGGxE==1){
test.2stp.egGxE<- "EG|GxE"
test.2stp<-cbind(test.2stp,test.2stp.egGxE)
}

if(RunEG2df==1){
test.2stp.eg2df<- "EG|2df"
test.2stp<-cbind(test.2stp,test.2stp.eg2df)
}

if(RunCocktail==1){
test.2stp.cocktail<- "Cocktail"
test.2stp<-cbind(test.2stp,test.2stp.cocktail)
}

test.2stp<-test.2stp[-1]
test.2stp<-data.frame(test.2stp)
colnames(test.2stp)<-"Method"

## Total number of SNPs (M)
M <- c(NA)

if(RunEDGxE==1){
M<-cbind(M,M.edge)
}

if(RunDGGxE==1){
M<-cbind(M,M.dgGxE)
}

if(RunDGEB==1){
M<-cbind(M,M.dgEB)
}

if(RunEGGxE==1){
M<-cbind(M,M.egGxE)
}

if(RunEG2df==1){
M<-cbind(M,M.eg2df)
}

if(RunCocktail==1){
M<-cbind(M,M.cocktail)
}

M<-M[-1]

M<-data.frame(M)

## Step1 significance level
s1_alpha1<-c(NA)

if(RunEDGxE==1){
s1_alpha1.egdg<-signif(EDGxEscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.egdg)
}

if(RunDGGxE==1){
s1_alpha1.dg<-signif(DGscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.dg)
}

if(RunDGEB==1){
s1_alpha1.dg<-signif(DGEBscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.dg)
}

if(RunEGGxE==1){
s1_alpha1.eg<-signif(EGscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.eg)
}

if(RunEG2df==1){
s1_alpha1.eg<-signif(EG2dfscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.eg)
}

if(RunCocktail==1){
s1_alpha1.cocktail<-signif(Cocktailscreen,3)
s1_alpha1<-cbind(s1_alpha1,s1_alpha1.cocktail)
}

s1_alpha1<-s1_alpha1[-1]

s1_alpha1<-data.frame(s1_alpha1)

colnames(s1_alpha1) <- "alpha1"

## Number of SNPs pass to step 2

ns2<-c(NA)

if(RunEDGxE==1){
ns2<-cbind(ns2,edge.s2.size)
}

if(RunDGGxE==1){
ns2<-cbind(ns2,dgGxE.s2.size)
}

if(RunDGEB==1){
ns2<-cbind(ns2,dgEB.s2.size)
}

if(RunEGGxE==1){
ns2<-cbind(ns2,egGxE.s2.size)
}

if(RunEG2df==1){
ns2<-cbind(ns2,eg2df.s2.size)
}

if(RunCocktail==1){
ns2<-cbind(ns2,cocktail.s2.size)
}

ns2<- ns2[-1]
ns2 <- data.frame(ns2)
colnames(ns2) <- "m"

## significance level at step 2

sig2<-c(NA)

if(RunEDGxE==1){
sig2<-cbind(sig2,signif(sig_edge_st2,3))
}

if(RunDGGxE==1){
sig2<-cbind(sig2,signif(sig_dgGxE_st2,3))
}

if(RunDGEB==1){
sig2<-cbind(sig2,signif(sig_dgEB_st2,3))
}

if(RunEGGxE==1){
sig2<-cbind(sig2,signif(sig_egGxE_st2,3))
}

if(RunEG2df==1){
sig2<-cbind(sig2,signif(sig_eg2df_st2,3))
}

if(RunCocktail==1){
sig2<-cbind(sig2,signif(sig_cocktail_st2,3))
}

sig2<-sig2[-1]
sig2 <- data.frame(sig2)
colnames(sig2) <- "sig2"

## Bin size for each method
binsize<-c(NA)

if(RunEDGxE==1){
binsize<-cbind(binsize,EDGxE.bin)
}

if(RunDGGxE==1){
binsize<-cbind(binsize,DGGxE.bin)
}

if(RunDGEB==1){
binsize<-cbind(binsize,DGEB.bin)
}

if(RunEGGxE==1){
binsize<-cbind(binsize,EGGxE.bin)
}

if(RunEG2df==1){
binsize<-cbind(binsize,EG2df.bin)
}

if(RunCocktail==1){
binsize<-cbind(binsize,Cocktail.bin)
}

binsize<-binsize[-1]
binsize<-data.frame(binsize)
colnames(binsize)<-"Bin Size"

step2 <- data.frame(test.2stp,M,s1_alpha1,ns2,sig2,binsize)
colnames(step2)<-c("Method","M","alpha1","m","alpha/m","Bin Size")

te1.row<- 11-diff.st1
te1.row2<- te1.row+1
wrap(wb,te1, "summary",te1.row,1,F,F)
save2excel(wb,step2, "summary",te1.row2,1,T)
rm(step2)
}


#####Descriptions

if(num.2step>0){
te2.row<- 20-diff.st2
} else {
te2.row<- 11-diff.st1
}


te2<-"M:  total number of SNPs tested ( After removing unconverged SNPs )"
wrap(wb,te2, "summary",te2.row,1,F,F)

te3.row<- te2.row+1
te3<-"alpha/M:  Bonferroni corrected significance level for exhaustive testing"
wrap(wb,te3, "summary",te3.row,1,F,F)

te4.row<- te3.row+1
te4<-"alpha1:  Step1 significance threshold for 2-step methods"
wrap(wb,te4, "summary",te4.row,1,F,F)

te5.row<- te4.row+1
te5<-"m: Number of SNPs that passed Step1"
wrap(wb,te5, "summary",te5.row,1,F,F)

te6.row<- te5.row+1
te6<-"alpha1/m:  Significance level for SNPs tested in Step 2; subset testing"
wrap(wb,te6, "summary",te6.row,1,F,F)

te7.row<- te6.row+1
te7<-"Bin Size:  Initial Bin size for weighted hypothesis testing in Step 2"
wrap(wb,te7, "summary",te7.row,1,F,F)

te0.row<- te7.row+2
wrap(wb,te0, "summary",te0.row,1,F,F)


te00.row<- te0.row+1
te00<-"MarG:  Standard GWAS analysis of marginal G effect for all M SNPs"
wrap(wb,te00, "summary",te00.row,1,F,F)


te8.row<- te00.row+1
te8<-"CC:  Standard case-control analysis of GxE interaction for all M SNPs"
wrap(wb,te8, "summary",te8.row,1,F,F)

te9.row<- te8.row+1
te9<-"CO: Standard case-only analysis of GxE interaction for all M SNPs"
wrap(wb,te9, "summary",te9.row,1,F,F)

te10.row<- te9.row+1
te10<-"CntlOnly: Analysis of E-G association in controls only"
wrap(wb,te10, "summary",te10.row,1,F,F)

te11.row<- te10.row+1
te11<-"EB: Empirical Bayes Analysis of GxE interaction for all M SNPs (Mukherjee &Chatterjee, 2008)"
wrap(wb,te11, "summary",te11.row,1,F,F)

te12.row<- te11.row+1
te12<-"2df: Two degree of freedom joint test of G,GxE for all M SNPs( Kraft et al.2007)"
wrap(wb,te12, "summary",te12.row,1,F,F)

te13.row<- te12.row+1
te13<-"3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs (Gauderman et al., 2013)"
wrap(wb,te13, "summary",te13.row,1,F,F)

te1.2.row<- te13.row+2
wrap(wb,te1, "summary",te1.2.row,1,F,F)

te14.row<- te1.2.row+1
te14<-"EDGxE: 2-step with Step-1 screen based on joint test of E-G and D-G association, Step-2 test of GxE using CC analysis (Gauderman et al., 2013)"
wrap(wb,te14, "summary",te14.row,1,F,F)

te15.row<- te14.row+1
te15<-"DG|GxE: 2-step with Step-1 screen based on D-G association, Step-2 test of GxE using CC analysis (Kooperberg and LeBlanc, 2009)"
wrap(wb,te15, "summary",te15.row,1,F,F)

te16.row<- te15.row+1
te16<-"DG|EB: 2-step with Step-1 screen based on D-G association, Step-2 test of GxE using EB analysis (Hsu et al., 2012)"
wrap(wb,te16, "summary",te16.row,1,F,F)

te17.row<- te16.row+1
te17<-"EG|GxE: 2-step with Step-1 screen based on E-G association, Step-2 test of GxE using CC analysis (Murcray et al., 2009)"
wrap(wb,te17, "summary",te17.row,1,F,F)

te18.row<- te17.row+1
te18<-"EG|2df: 2-step with Step-1 screen based on E-G association, Step-2 is 2df joint test of G,GxE (Gauderman et al., 2013)"
wrap(wb,te18, "summary",te18.row,1,F,F)

te19.row<- te18.row+1
te19<-"Cocktail: Cocktail I Analysis of GxE interaction (Hsu et al., 2012)"
wrap(wb,te19, "summary",te19.row,1,F,F)

rm(te0,te00,te2,te3,te4,te5,te6,te7,te8,te9,te10,te11,te12,te13,te14,te15,te16,te17,te18,te19)
##################### Write out top hits and plots #################################3

alphabet<- c("I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")


if(IncludeRanks==1){

## MarG
if(RunDG==1){
createSheet(wb, name = "MarG")
########## 1st panel

# Marginal G QQ

saveimage(wb,"graph","MarG","DG_A.png",1,1)

# Marginal G Manhattan

saveimage(wb,"graph","MarG","DG_B.png",1,11)

setwidth( wb,"MarG",1,3000)
setwidth( wb,"MarG",6,2600)

tMarG<-"MarG:  Standard GWAS analysis of marginal G effect for all M SNPs"
wrap(wb,tMarG, "MarG",26,1,F,F)

sigMarG<-paste("Number of SNPs tested = ",M.dg,"; significance threshold = ",sig_dg,sep="")
wrap(wb,sigMarG, "MarG",27,1,F,F)

RankDescrip<-"A blank under Ranks means this SNP does not pass through the screening step for the corresponding 2-step approach"
wrap(wb,RankDescrip, "MarG",28,1,F,F)

rMarG<-"Ranks"
save2excel(wb,rMarG, "MarG",29,10,F)

MarG.colend<- alphabet[1+num.method]
MarG.letter<- paste("J29:",MarG.colend,"29",sep="")
MarG.col<- 9+num.method
savetitle( wb,"MarG",MarG.letter,29,c(10:MarG.col))
save2excel(wb,dg.top, "MarG",30,1,T)
rm(dg.top)
}

##CC
if(RunCC==1){
createSheet(wb, name = "CC")
########## 1st panel

# Case-Control QQ

saveimage(wb,"graph","CC","CC_A.png",1,1)

# Case-Control Manhattan

saveimage(wb,"graph","CC","CC_B.png",1,11)

setwidth( wb,"CC",1,3000)
setwidth( wb,"CC",6,2600)

tcc<-"CC:  Standard case-control analysis of GxE interaction for all M SNPs"
wrap(wb,tcc, "CC",26,1,F,F)

sigcc<-paste("Number of SNPs tested = ",M.cc,"; significance threshold = ",sig_cc,sep="")
wrap(wb,sigcc, "CC",27,1,F,F)

wrap(wb,RankDescrip, "CC",28,1,F,F)

rcc<-"Ranks"
save2excel(wb,rcc, "CC",29,10,F)

cc.colend<- alphabet[1+num.method]
cc.letter<- paste("J29:",cc.colend,"29",sep="")
cc.col<- 9+num.method
savetitle( wb,"CC",cc.letter,29,c(10:cc.col))
save2excel(wb,cc.top, "CC",30,1,T)
rm(cc.top)
}

##CO
if(RunCO==1){
createSheet(wb, name = "CO")
########## 2nd panel
# Case-Only QQ
saveimage(wb,"graph","CO","CO_A.png",1,1)
# Case-Only Manhattan
saveimage(wb,"graph","CO","CO_B.png",1,11)

setwidth( wb,"CO",1,3000)

tco<-"CO: Standard case-only analysis of GxE interaction for all M SNPs"
wrap(wb,tco, "CO",26,1,F,F)

sigco<-paste("Number of SNPs tested = ",M.co,"; significance threshold = ",sig_co,sep="")
wrap(wb,sigco, "CO",27,1,F,F)

wrap(wb,RankDescrip, "CO",28,1,F,F)

rco<-"Ranks"
save2excel(wb,rco, "CO",29,10,F)

co.colend<- alphabet[1+num.method]
co.letter<- paste("J29:",co.colend,"29",sep="")
co.col<- 9+num.method

savetitle( wb,"CO",co.letter,29,c(10:co.col))

save2excel(wb,co.top, "CO",30,1,T)
rm(co.top)
}

## CntlyOnly
if(RunCNTL==1){
createSheet(wb, name = "CntlOnly")
########## 3rd panel
# Control-Only QQ
saveimage(wb,"graph","CntlOnly","CNTL_A.png",1,1)

# Control-Only Manhattan
saveimage(wb,"graph","CntlOnly","CNTL_B.png",1,11)

setwidth( wb,"CntlOnly",1,3000)

tcntl<-"CntlOnly: Analysis of G-E association in controls only"
wrap(wb,tcntl, "CntlOnly",26,1,F,F)

wrap(wb,RankDescrip, "CntlOnly",27,1,F,F)

rcntl<-"Ranks"
save2excel(wb,rcntl, "CntlOnly",28,10,F)
cntl.colend<- alphabet[1+num.method]
cntl.letter<- paste("J28:",cntl.colend,"28",sep="")
cntl.col<- 9+num.method
savetitle( wb,"CntlOnly",cntl.letter,28,c(10:cntl.col))
save2excel(wb,cntl.top, "CntlOnly",29,1,T)
rm(cntl.top)
}

##EB
if(RunEB==1){
createSheet(wb, name = "EB")
####### 4th panel
# EB QQ
saveimage(wb,"graph","EB","EB_A.png",1,1)
# EB Manhattan
saveimage(wb,"graph","EB","EB_B.png",1,11)

setwidth( wb,"EB",1,3000)
setwidth( wb,"EB",6,2600)

teb<-"EB: Empirical Bayes Analysis of GxE interaction for all M SNPs"
wrap(wb,teb, "EB",26,1,F,F)
sigeb<-paste("Number of SNPs tested = ",M.eb,"; significance threshold = ",sig_eb,sep="")
wrap(wb,sigeb, "EB",27,1,F,F)

wrap(wb,RankDescrip, "EB",28,1,F,F)

reb<-"Ranks"
save2excel(wb,reb, "EB",29,10,F)
eb.colend<- alphabet[1+num.method]
eb.letter<- paste("J29:",eb.colend,"29",sep="")
eb.col<- 9+num.method
savetitle( wb,"EB",eb.letter,29,c(10:eb.col))
save2excel(wb,eb.top, "EB",30,1,T)
rm(eb.top)
}

##2df
if(Run2df==1){
createSheet(wb, name = "df2")
####### 5th panel
# 2df QQ
saveimage(wb,"graph","df2","df2_A.png",1,1)
 
# 2df Manhattan
saveimage(wb,"graph","df2","df2_B.png",1,11)


setwidth( wb,"df2",6,2600)
setwidth( wb,"df2",1,3000)

t2df<-"2df: Two degree of freedom joint test of G,GxE for all M SNPs"
wrap(wb,t2df, "df2",26,1,F,F)

sig2df<-paste("Number of SNPs tested = ",M.2df,"; significance threshold = ",sig_2df,sep="")
wrap(wb,sig2df, "df2",27,1,F,F)

wrap(wb,RankDescrip, "df2",28,1,F,F)

r2df<-"Ranks"
save2excel(wb,r2df, "df2",29,9,F)
twodf.colend<- alphabet[num.method]
twodf.letter<- paste("I29:",twodf.colend,"29",sep="")
twodf.col<- 8+num.method
savetitle( wb,"df2",twodf.letter,29,c(9:twodf.col))
save2excel(wb,df2.top, "df2",30,1,T)
rm(df2.top)
}

##3df
if(Run3df==1){
createSheet(wb, name = "df3")
######## 6th panel 3df QQ
saveimage(wb,"graph","df3","df3_A.png",1,1)

# Control-Only Manhattan
saveimage(wb,"graph","df3","df3_B.png",1,11)

setwidth( wb,"df3",6,2600)
setwidth( wb,"df3",1,3000)

t3df<-"3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs"
wrap(wb,t3df, "df3",26,1,F,F)

sig3df<-paste("Number of SNPs tested = ",M.3df,"; significance threshold = ",sig_3df,sep="")
wrap(wb,sig3df, "df3",27,1,F,F)

wrap(wb,RankDescrip, "df3",28,1,F,F)

r3df<-"Ranks"
threedf.colend<- alphabet[num.method]
threedf.letter<- paste("I29:",threedf.colend,"29",sep="")
threedf.col<- 8+num.method
save2excel(wb,r3df, "df3",29,9,F)
savetitle( wb,"df3",threedf.letter,29,c(9:threedf.col))

save2excel(wb,df3.top, "df3",30,1,T)
rm(df3.top)
}

##EDGxE
if(RunEDGxE==1){
createSheet(wb, name = "EDGxE_subset")
########## 7th panel
saveimage(wb,"graph","EDGxE_subset","EDGE_A.png",1,1)
# CS/T step1 Manhattan
saveimage(wb,"graph","EDGxE_subset","EDGE_B.png",1,11)

# CS/T step2 QQ
if (edge.s2.size>0) {
saveimage(wb,"graph","EDGxE_subset","EDGE_C.png",26,1)

# CS/T step2 Manhattan
saveimage(wb,"graph","EDGxE_subset","EDGE_D.png",26,11)


setwidth( wb,"EDGxE_subset",6,2600)
setwidth( wb,"EDGxE_subset",8,2600)
setwidth( wb,"EDGxE_subset",9,2600)
setwidth( wb,"EDGxE_subset",1,3000)

tdgeg<-"EDGxE: 2-step with screening based on joint test of E-G and D-G association; subset testing of GxE in Step 2"
wrap(wb,tdgeg, "EDGxE_subset",51,1,F,F)

sigdgeg<-paste("Number of SNPs tested in Step 2= ",edge.s2.size,"; Step2 significance threshold = ",sig_edge_st2,sep="")
wrap(wb,sigdgeg, "EDGxE_subset",52,1,F,F)

wrap(wb,RankDescrip, "EDGxE_subset",53,1,F,F)

rdgeg<-"Ranks"
save2excel(wb,rdgeg, "EDGxE_subset",54,12,F)
egdg.colend<- alphabet[3+num.method]
egdg.letter<- paste("L54:",egdg.colend,"54",sep="")
egdg.col<- 11+num.method
savetitle( wb,"EDGxE_subset",egdg.letter,54,c(12:egdg.col))

save2excel(wb,edge.s2, "EDGxE_subset",55,1,T)
rm(edge.s2)
}

#EDGxE wgted

########## 8th panel
createSheet(wb, name = "EDGxE_wgted")
saveimage(wb,"graph","EDGxE_wgted","EDGE_E.png",1,1)

##dgeg wt
setwidth( wb,"EDGxE_wgted",7,2600)
setwidth( wb,"EDGxE_wgted",9,2600)
setwidth( wb,"EDGxE_wgted",10,2600)
setwidth( wb,"EDGxE_wgted",12,2600)
setwidth( wb,"EDGxE_wgted",3,3000)

tdgegwt<-"EDGxE: 2-step with screening based on joint test of E-G and D-G association; weighted testing of GxE in Step 2"
wrap(wb,tdgegwt, "EDGxE_wgted",26,1,F,F)
save2excel(wb,edge.wt.top, "EDGxE_wgted",27,1,T)
rm(edge.wt.top)
}

##DG|GxE
if(RunDGGxE==1){
createSheet(wb, name = "DGGxE_subset")
#######9th panel
# D-G step1 QQ
saveimage(wb,"graph","DGGxE_subset","DGGxE_A.png",1,1)

# D-G step1 Manhattan
saveimage(wb,"graph","DGGxE_subset","DGGxE_B.png",1,11)

##### D-G step2 QQ
if (dgGxE.s2.size>0) {
saveimage(wb,"graph","DGGxE_subset","DGGxE_C.png",26,1)

######D-G step2 Manhattan

saveimage(wb,"graph","DGGxE_subset","DGGxE_D.png",26,11)


setwidth( wb,"DGGxE_subset",6,2600)
setwidth( wb,"DGGxE_subset",8,2600)
setwidth( wb,"DGGxE_subset",9,2600)
setwidth( wb,"DGGxE_subset",1,3000)

tdg<-"DG|GxE: 2-step with screening based on D-G association; subset testing of GxE in Step 2"
wrap(wb,tdg, "DGGxE_subset",51,1,F,F)

sigdg<-paste("Number of SNPs tested in Step 2= ",dgGxE.s2.size,"; Step2 significance threshold = ",sig_dgGxE_st2,sep="")
wrap(wb,sigdg, "DGGxE_subset",52,1,F,F)

wrap(wb,RankDescrip, "DGGxE_subset",53,1,F,F)


rdg<-"Ranks"
save2excel(wb,rdg, "DGGxE_subset",54,12,F)
dg.colend<- alphabet[3+num.method]
dg.letter<- paste("L54:",dg.colend,"54",sep="")
dg.col<- 11+num.method
savetitle( wb,"DGGxE_subset",dg.letter,54,c(12:dg.col))

save2excel(wb,dgGxE.s2, "DGGxE_subset",55,1,T)
rm(dgGxE.s2)
}


##dg|GxE wt

#######10th panel

##### DG wgted 
createSheet(wb, name = "DGGxE_wgted")
saveimage(wb,"graph","DGGxE_wgted","DGGxE_E.png",1,1)

setwidth( wb,"DGGxE_wgted",12,2800)
setwidth( wb,"DGGxE_wgted",7,2600)
setwidth( wb,"DGGxE_wgted",9,2600)
setwidth( wb,"DGGxE_wgted",10,2600)
setwidth( wb,"DGGxE_wgted",3,3000)

tdgwt<-"DG|GxE: 2-step with screening based on D-G association; weighted testing of GxE in Step 2"
wrap(wb,tdgwt, "DGGxE_wgted",26,1,F,F)

save2excel(wb,dgGxE.wt.top, "DGGxE_wgted",27,1,T)
rm(dgGxE.wt.top)
}

##dg|EB
if(RunDGEB==1){
#######11th panel
createSheet(wb, name = "DGEB_subset")
# D-G |EB step1 QQ

saveimage(wb,"graph","DGEB_subset","DGEB_A.png",1,1)

# D-G step1 Manhattan
saveimage(wb,"graph","DGEB_subset","DGEB_B.png",1,11)

##### D-G step2 QQ
if (dgEB.s2.size>0) {
saveimage(wb,"graph","DGEB_subset","DGEB_C.png",26,1)


######D-G|EB step2 Manhattan

saveimage(wb,"graph","DGEB_subset","DGEB_D.png",26,11)

setwidth( wb,"DGEB_subset",6,2600)
setwidth( wb,"DGEB_subset",8,2600)
setwidth( wb,"DGEB_subset",9,2600)
setwidth( wb,"DGEB_subset",1,3000)

tdgeb<-"DGEB: 2-step with screening based on D-G association; subset testing of GxE using EB in Step 2"
wrap(wb,tdgeb, "DGEB_subset",51,1,F,F)

sigdgeb<-paste("Number of SNPs tested in Step 2= ",dgEB.s2.size,"; Step2 significance threshold = ",sig_dgEB_st2,sep="")
wrap(wb,sigdgeb, "DGEB_subset",52,1,F,F)

wrap(wb,RankDescrip, "DGEB_subset",53,1,F,F)

rdgeb<-"Ranks"
save2excel(wb,rdgeb, "DGEB_subset",54,12,F)
dgeb.colend<- alphabet[3+num.method]
dgeb.letter<- paste("L54:",dgeb.colend,"54",sep="")
dgeb.col<- 11+num.method
savetitle( wb,"DGEB_subset",dgeb.letter,54,c(12:dgeb.col))

save2excel(wb,dgEB.s2, "DGEB_subset",55,1,T)
rm(dgEB.s2)
}

##dg|EB wt

####### 12th panel
createSheet(wb, name = "DGEB_wgted")
##### DG wgted 
saveimage(wb,"graph","DGEB_wgted","DGEB_E.png",1,1)

setwidth( wb,"DGEB_wgted",12,2800)
setwidth( wb,"DGEB_wgted",7,2600)
setwidth( wb,"DGEB_wgted",9,2600)
setwidth( wb,"DGEB_wgted",10,2600)
setwidth( wb,"DGEB_wgted",3,3000)

tdgebwt<-"DGEB: 2-step with screening based on D-G association; weighted testing of GxE using EB in Step 2"
wrap(wb,tdgebwt, "DGEB_wgted",26,1,F,F)

save2excel(wb,dgEB.wt.top, "DGEB_wgted",27,1,T)
rm(dgEB.wt.top)
}


##eg
if(RunEGGxE==1){
createSheet(wb, name ="EGGxE_subset")
# E-G step1 QQ
saveimage(wb,"graph","EGGxE_subset","EGGxE_A.png",1,1)

# E-G step1 Manhattan
saveimage(wb,"graph","EGGxE_subset","EGGxE_B.png",1,11)


#E-G step2 QQ
if (egGxE.s2.size>0) {
saveimage(wb,"graph","EGGxE_subset","EGGxE_C.png",26,1)


#E-G step2 Manhattan

saveimage(wb,"graph","EGGxE_subset","EGGxE_D.png",26,11)


setwidth(wb, "EGGxE_subset",6,2600)
setwidth(wb, "EGGxE_subset",8,2600)
setwidth(wb, "EGGxE_subset",9,2600)
setwidth(wb, "EGGxE_subset",1,3000)

teg<-"EG|GxE: 2-step with screening based on E-G association; subset testing of GxE in Step 2"
wrap(wb,teg, "EGGxE_subset",51,1,F,F)

sigeg<-paste("Number of SNPs tested in Step 2= ",egGxE.s2.size,"; Step2 significance threshold = ",sig_egGxE_st2,sep="")
wrap(wb,sigeg, "EGGxE_subset",52,1,F,F)

wrap(wb,RankDescrip, "EGGxE_subset",53,1,F,F)

reg<-"Ranks"
eg.colend<- alphabet[3+num.method]
eg.letter<- paste("L54:",eg.colend,"54",sep="")
eg.col<- 11+num.method
save2excel(wb,reg, "EGGxE_subset",54,12,F)
savetitle( wb,"EGGxE_subset",eg.letter,54,c(12:eg.col))

save2excel(wb,egGxE.s2, "EGGxE_subset",55,1,T)
rm(egGxE.s2)
}

##eg wt
createSheet(wb, name = "EGGxE_wgted")
saveimage(wb,"graph","EGGxE_wgted","EGGxE_E.png",1,1)

setwidth( wb,"EGGxE_wgted",7,2600)
setwidth( wb,"EGGxE_wgted",9,2600)
setwidth( wb,"EGGxE_wgted",10,2600)
setwidth( wb,"EGGxE_wgted",12,2600)
setwidth( wb,"EGGxE_wgted",3,3000)

tegwt<-"EG|GxE: 2-step with screening based on E-G association; weighted testing of GxE in Step 2"
wrap(wb,tegwt, "EGGxE_wgted",26,1,F,F)

save2excel(wb,egGxE.wt.top, "EGGxE_wgted",27,1,T)
rm(egGxE.wt.top)
}

##eg2df
if(RunEG2df==1){
createSheet(wb, name = "EG2df_subset")
# EG-2df step1 QQ
saveimage(wb,"graph","EG2df_subset","EG2df_A.png",1,1)

saveimage(wb,"graph","EG2df_subset","EG2df_B.png",1,11)

#EG-2df Step2 QQ
if (eg2df.s2.size>0) {
saveimage(wb,"graph","EG2df_subset","EG2df_C.png",26,1)

#EG-2df step2 Manhattan
saveimage(wb,"graph","EG2df_subset","EG2df_D.png",26,11)

setwidth( wb,"EG2df_subset",6,2600)
setwidth( wb,"EG2df_subset",8,2600)
setwidth( wb,"EG2df_subset",1,3000)

teg2df<-"EG|2df: 2-step with screening based on E-G association; subset testing of G, GxE in Step 2"
wrap(wb,teg2df, "EG2df_subset",51,1,F,F)

sigeg2df<-paste("Number of SNPs tested in Step 2= ",eg2df.s2.size,"; Step2 significance threshold = ",sig_eg2df_st2,sep="")
wrap(wb,sigeg2df, "EG2df_subset",52,1,F,F)

wrap(wb,RankDescrip, "EG2df_subset",53,1,F,F)

reg2df<-"Ranks"
save2excel(wb,reg2df, "EG2df_subset",54,11,F)
eg2df.colend<- alphabet[2+num.method]
eg2df.letter<- paste("K54:",eg2df.colend,"54",sep="")
eg2df.col<- 10+num.method
savetitle( wb,"EG2df_subset",eg2df.letter,54,c(11:eg2df.col))

save2excel(wb,eg2df.s2, "EG2df_subset",55,1,T)
rm(eg2df.s2)
}

##eg2df wt
createSheet(wb, name = "EG2df_wgted")
saveimage(wb,"graph","EG2df_wgted","EG2df_E.png",1,1)

setwidth( wb,"EG2df_wgted",7,2600)
setwidth( wb,"EG2df_wgted",9,2600)
setwidth( wb,"EG2df_wgted",11,2600)
setwidth( wb,"EG2df_wgted",3,3000)


teg2dfwt<-"EG|2df: 2-step with screening based on E-G association; weighted testing of G, GxE in Step 2"
wrap(wb,teg2dfwt, "EG2df_wgted",26,1,F,F)
save2excel(wb,eg2df.wt.top, "EG2df_wgted",27,1,T)
rm(eg2df.wt.top)
}

##cocktail subset
if(RunCocktail==1){
createSheet(wb, name = "Cocktail_subset")
# Cocktail subset step1 QQ
saveimage(wb,"graph","Cocktail_subset","Cocktail_A.png",1,1)

# Cocktail subset step1 Manhattan

saveimage(wb,"graph","Cocktail_subset","Cocktail_B.png",1,11)

#Cocktail subset Step2 QQ
if (cocktail.s2.size>0) {

saveimage(wb,"graph","Cocktail_subset","Cocktail_C.png",26,1)

#Cocktail subset step2 Manhattan

saveimage(wb,"graph","Cocktail_subset","Cocktail_D.png",26,11)


setwidth( wb,"Cocktail_subset",1,3000)

tcocktail<-"Cocktail subset: 2-step with Cocktail I screening in Step 1 and subset testing in Step 2"
wrap(wb,tcocktail, "Cocktail_subset",51,1,F,F)

sigcocktail<-paste("Number of SNPs tested in Step 2= ",cocktail.s2.size,"; Step2 significance threshold = ",sig_cocktail_st2,sep="")
wrap(wb,sigcocktail, "Cocktail_subset",52,1,F,F)

wrap(wb,RankDescrip, "Cocktail_subset",53,1,F,F)

rcocktail<-"Ranks"
save2excel(wb,rcocktail, "Cocktail_subset",54,10,F)
cocktail.colend<- alphabet[1+num.method]
cocktail.letter<- paste("J54:",cocktail.colend,"54",sep="")
cocktail.col<- 9+num.method
savetitle( wb,"Cocktail_subset",cocktail.letter,54,c(10:cocktail.col))

save2excel(wb,cocktail.s2, "Cocktail_subset",55,1,T)
rm(cocktail.s2)
}

##cocktail wt
createSheet(wb, name = "Cocktail_wgted")
saveimage(wb,"graph","Cocktail_wgted","Cocktail_E.png",1,1)

setwidth( wb,"Cocktail_wgted",10,2600)
setwidth( wb,"Cocktail_wgted",3,3000)

tcocktailwt<-"Cocktail weighted: 2-step with Cocktail I screening in Step 1 and weighted testing in Step 2"
wrap(wb,tcocktailwt, "Cocktail_wgted",26,1,F,F)

save2excel(wb,cocktail.wt.top, "Cocktail_wgted",27,1,T)
rm(cocktail.wt.top)
}
saveWorkbook(wb)
rm(wb)
} else {


## MarG
if(RunDG==1){
createSheet(wb, name = "MarG")
########## 1st panel

# Marginal G QQ

saveimage(wb,"graph","MarG","DG_A.png",1,1)

# Marginal G Manhattan

saveimage(wb,"graph","MarG","DG_B.png",1,11)

setwidth( wb,"MarG",1,3000)
setwidth( wb,"MarG",6,2600)

tMarG<-"MarG:  Standard GWAS analysis of marginal G effect for all M SNPs"
wrap(wb,tMarG, "MarG",26,1,F,F)

sigMarG<-paste("Number of SNPs tested = ",M.dg,"; significance threshold = ",sig_dg,sep="")
wrap(wb,sigMarG, "MarG",27,1,F,F)

save2excel(wb,dg.top, "MarG",29,1,T)
rm(dg.top)
}

##CC
if(RunCC==1){
createSheet(wb, name = "CC")
########## 1st panel

# Case-Control QQ

saveimage(wb,"graph","CC","CC_A.png",1,1)

# Case-Control Manhattan

saveimage(wb,"graph","CC","CC_B.png",1,11)

setwidth( wb,"CC",1,3000)
setwidth( wb,"CC",6,2600)

tcc<-"CC:  Standard case-control analysis of GxE interaction for all M SNPs"
wrap(wb,tcc, "CC",26,1,F,F)

sigcc<-paste("Number of SNPs tested = ",M.cc,"; significance threshold = ",sig_cc,sep="")
wrap(wb,sigcc, "CC",27,1,F,F)

save2excel(wb,cc.top, "CC",29,1,T)
rm(cc.top)
}

##CO
if(RunCO==1){
createSheet(wb, name = "CO")
########## 2nd panel
# Case-Only QQ
saveimage(wb,"graph","CO","CO_A.png",1,1)
# Case-Only Manhattan
saveimage(wb,"graph","CO","CO_B.png",1,11)

setwidth( wb,"CO",1,3000)

tco<-"CO: Standard case-only analysis of GxE interaction for all M SNPs"
wrap(wb,tco, "CO",26,1,F,F)

sigco<-paste("Number of SNPs tested = ",M.co,"; significance threshold = ",sig_co,sep="")
wrap(wb,sigco, "CO",27,1,F,F)

save2excel(wb,co.top, "CO",29,1,T)
rm(co.top)
}

## CntlyOnly
if(RunCNTL==1){
createSheet(wb, name = "CntlOnly")
########## 3rd panel
# Control-Only QQ
saveimage(wb,"graph","CntlOnly","CNTL_A.png",1,1)

# Control-Only Manhattan
saveimage(wb,"graph","CntlOnly","CNTL_B.png",1,11)

setwidth( wb,"CntlOnly",1,3000)

tcntl<-"CntlOnly: Analysis of G-E association in controls only"
wrap(wb,tcntl, "CntlOnly",26,1,F,F)

save2excel(wb,cntl.top, "CntlOnly",28,1,T)
rm(cntl.top)
}

##EB
if(RunEB==1){
createSheet(wb, name = "EB")
####### 4th panel
# EB QQ
saveimage(wb,"graph","EB","EB_A.png",1,1)
# EB Manhattan
saveimage(wb,"graph","EB","EB_B.png",1,11)

setwidth( wb,"EB",1,3000)
setwidth( wb,"EB",6,2600)

teb<-"EB: Empirical Bayes Analysis of GxE interaction for all M SNPs"
wrap(wb,teb, "EB",26,1,F,F)
sigeb<-paste("Number of SNPs tested = ",M.eb,"; significance threshold = ",sig_eb,sep="")
wrap(wb,sigeb, "EB",27,1,F,F)

save2excel(wb,eb.top, "EB",29,1,T)
rm(eb.top)
}

##2df
if(Run2df==1){
createSheet(wb, name = "df2")
####### 5th panel
# 2df QQ
saveimage(wb,"graph","df2","df2_A.png",1,1)
 
# 2df Manhattan
saveimage(wb,"graph","df2","df2_B.png",1,11)


setwidth( wb,"df2",6,2600)
setwidth( wb,"df2",1,3000)

t2df<-"2df: Two degree of freedom joint test of G,GxE for all M SNPs"
wrap(wb,t2df, "df2",26,1,F,F)

sig2df<-paste("Number of SNPs tested = ",M.2df,"; significance threshold = ",sig_2df,sep="")
wrap(wb,sig2df, "df2",27,1,F,F)

save2excel(wb,df2.top, "df2",29,1,T)
rm(df2.top)
}

##3df
if(Run3df==1){
createSheet(wb, name = "df3")
######## 6th panel 3df QQ
saveimage(wb,"graph","df3","df3_A.png",1,1)

# Control-Only Manhattan
saveimage(wb,"graph","df3","df3_B.png",1,11)

setwidth( wb,"df3",6,2600)
setwidth( wb,"df3",1,3000)

t3df<-"3df: Three degree of freedom joint test of G,GxE and E-G association for all M SNPs"
wrap(wb,t3df, "df3",26,1,F,F)

sig3df<-paste("Number of SNPs tested = ",M.3df,"; significance threshold = ",sig_3df,sep="")
wrap(wb,sig3df, "df3",27,1,F,F)


save2excel(wb,df3.top, "df3",29,1,T)
rm(df3.top)
}

##EDGxE
if(RunEDGxE==1){
createSheet(wb, name = "EDGxE_subset")
########## 7th panel
saveimage(wb,"graph","EDGxE_subset","EDGE_A.png",1,1)
# CS/T step1 Manhattan
saveimage(wb,"graph","EDGxE_subset","EDGE_B.png",1,11)

# CS/T step2 QQ
if (edge.s2.size>0) {
saveimage(wb,"graph","EDGxE_subset","EDGE_C.png",26,1)

# CS/T step2 Manhattan
saveimage(wb,"graph","EDGxE_subset","EDGE_D.png",26,11)


setwidth( wb,"EDGxE_subset",6,2600)
setwidth( wb,"EDGxE_subset",8,2600)
setwidth( wb,"EDGxE_subset",9,2600)
setwidth( wb,"EDGxE_subset",1,3000)

tdgeg<-"EDGxE: 2-step with screening based on joint test of E-G and D-G association; subset testing of GxE in Step 2"
wrap(wb,tdgeg, "EDGxE_subset",51,1,F,F)

sigdgeg<-paste("Number of SNPs tested in Step 2= ",edge.s2.size,"; Step2 significance threshold = ",sig_edge_st2,sep="")
wrap(wb,sigdgeg, "EDGxE_subset",52,1,F,F)

save2excel(wb,edge.s2, "EDGxE_subset",54,1,T)
rm(edge.s2)
}

#EDGxE wgted

########## 8th panel
createSheet(wb, name = "EDGxE_wgted")
saveimage(wb,"graph","EDGxE_wgted","EDGE_E.png",1,1)

##dgeg wt
setwidth( wb,"EDGxE_wgted",7,2600)
setwidth( wb,"EDGxE_wgted",9,2600)
setwidth( wb,"EDGxE_wgted",10,2600)
setwidth( wb,"EDGxE_wgted",12,2600)
setwidth( wb,"EDGxE_wgted",3,3000)

tdgegwt<-"EDGxE: 2-step with screening based on joint test of E-G and D-G association; weighted testing of GxE in Step 2"
wrap(wb,tdgegwt, "EDGxE_wgted",26,1,F,F)
save2excel(wb,edge.wt.top, "EDGxE_wgted",27,1,T)
rm(edge.wt.top)
}

##DG|GxE
if(RunDGGxE==1){
createSheet(wb, name = "DGGxE_subset")
#######9th panel
# D-G step1 QQ
saveimage(wb,"graph","DGGxE_subset","DGGxE_A.png",1,1)

# D-G step1 Manhattan
saveimage(wb,"graph","DGGxE_subset","DGGxE_B.png",1,11)

##### D-G step2 QQ
if (dgGxE.s2.size>0) {
saveimage(wb,"graph","DGGxE_subset","DGGxE_C.png",26,1)

######D-G step2 Manhattan

saveimage(wb,"graph","DGGxE_subset","DGGxE_D.png",26,11)


setwidth( wb,"DGGxE_subset",6,2600)
setwidth( wb,"DGGxE_subset",8,2600)
setwidth( wb,"DGGxE_subset",9,2600)
setwidth( wb,"DGGxE_subset",1,3000)

tdg<-"DG|GxE: 2-step with screening based on D-G association; subset testing of GxE in Step 2"
wrap(wb,tdg, "DGGxE_subset",51,1,F,F)

sigdg<-paste("Number of SNPs tested in Step 2= ",dgGxE.s2.size,"; Step2 significance threshold = ",sig_dgGxE_st2,sep="")
wrap(wb,sigdg, "DGGxE_subset",52,1,F,F)

save2excel(wb,dgGxE.s2, "DGGxE_subset",54,1,T)
rm(dgGxE.s2)
}


##dg|GxE wt

#######10th panel

##### DG wgted 
createSheet(wb, name = "DGGxE_wgted")
saveimage(wb,"graph","DGGxE_wgted","DGGxE_E.png",1,1)

setwidth( wb,"DGGxE_wgted",12,2800)
setwidth( wb,"DGGxE_wgted",7,2600)
setwidth( wb,"DGGxE_wgted",9,2600)
setwidth( wb,"DGGxE_wgted",10,2600)
setwidth( wb,"DGGxE_wgted",3,3000)

tdgwt<-"DG|GxE: 2-step with screening based on D-G association; weighted testing of GxE in Step 2"
wrap(wb,tdgwt, "DGGxE_wgted",26,1,F,F)

save2excel(wb,dgGxE.wt.top, "DGGxE_wgted",27,1,T)
rm(dgGxE.wt.top)
}

##dg|EB
if(RunDGEB==1){
#######11th panel
createSheet(wb, name = "DGEB_subset")
# D-G |EB step1 QQ

saveimage(wb,"graph","DGEB_subset","DGEB_A.png",1,1)

# D-G step1 Manhattan
saveimage(wb,"graph","DGEB_subset","DGEB_B.png",1,11)

##### D-G step2 QQ
if (dgEB.s2.size>0) {
saveimage(wb,"graph","DGEB_subset","DGEB_C.png",26,1)


######D-G|EB step2 Manhattan

saveimage(wb,"graph","DGEB_subset","DGEB_D.png",26,11)

setwidth( wb,"DGEB_subset",6,2600)
setwidth( wb,"DGEB_subset",8,2600)
setwidth( wb,"DGEB_subset",9,2600)
setwidth( wb,"DGEB_subset",1,3000)

tdgeb<-"DGEB: 2-step with screening based on D-G association; subset testing of GxE using EB in Step 2"
wrap(wb,tdgeb, "DGEB_subset",51,1,F,F)

sigdgeb<-paste("Number of SNPs tested in Step 2= ",dgEB.s2.size,"; Step2 significance threshold = ",sig_dgEB_st2,sep="")
wrap(wb,sigdgeb, "DGEB_subset",52,1,F,F)


save2excel(wb,dgEB.s2, "DGEB_subset",54,1,T)
rm(dgEB.s2)
}

##dg|EB wt

####### 12th panel
createSheet(wb, name = "DGEB_wgted")
##### DG wgted 
saveimage(wb,"graph","DGEB_wgted","DGEB_E.png",1,1)

setwidth( wb,"DGEB_wgted",12,2800)
setwidth( wb,"DGEB_wgted",7,2600)
setwidth( wb,"DGEB_wgted",9,2600)
setwidth( wb,"DGEB_wgted",10,2600)
setwidth( wb,"DGEB_wgted",3,3000)

tdgebwt<-"DGEB: 2-step with screening based on D-G association; weighted testing of GxE using EB in Step 2"
wrap(wb,tdgebwt, "DGEB_wgted",26,1,F,F)

save2excel(wb,dgEB.wt.top, "DGEB_wgted",27,1,T)
rm(dgEB.wt.top)
}


##eg
if(RunEGGxE==1){
createSheet(wb, name ="EGGxE_subset")
# E-G step1 QQ
saveimage(wb,"graph","EGGxE_subset","EGGxE_A.png",1,1)

# E-G step1 Manhattan
saveimage(wb,"graph","EGGxE_subset","EGGxE_B.png",1,11)


#E-G step2 QQ
if (egGxE.s2.size>0) {
saveimage(wb,"graph","EGGxE_subset","EGGxE_C.png",26,1)


#E-G step2 Manhattan

saveimage(wb,"graph","EGGxE_subset","EGGxE_D.png",26,11)


setwidth(wb, "EGGxE_subset",6,2600)
setwidth(wb, "EGGxE_subset",8,2600)
setwidth(wb, "EGGxE_subset",9,2600)
setwidth(wb, "EGGxE_subset",1,3000)

teg<-"EG|GxE: 2-step with screening based on E-G association; subset testing of GxE in Step 2"
wrap(wb,teg, "EGGxE_subset",51,1,F,F)

sigeg<-paste("Number of SNPs tested in Step 2= ",egGxE.s2.size,"; Step2 significance threshold = ",sig_egGxE_st2,sep="")
wrap(wb,sigeg, "EGGxE_subset",52,1,F,F)


save2excel(wb,egGxE.s2, "EGGxE_subset",54,1,T)
rm(egGxE.s2)
}

##eg wt
createSheet(wb, name = "EGGxE_wgted")
saveimage(wb,"graph","EGGxE_wgted","EGGxE_E.png",1,1)

setwidth( wb,"EGGxE_wgted",7,2600)
setwidth( wb,"EGGxE_wgted",9,2600)
setwidth( wb,"EGGxE_wgted",10,2600)
setwidth( wb,"EGGxE_wgted",12,2600)
setwidth( wb,"EGGxE_wgted",3,3000)

tegwt<-"EG|GxE: 2-step with screening based on E-G association; weighted testing of GxE in Step 2"
wrap(wb,tegwt, "EGGxE_wgted",26,1,F,F)

save2excel(wb,egGxE.wt.top, "EGGxE_wgted",27,1,T)
rm(egGxE.wt.top)
}

##eg2df
if(RunEG2df==1){
createSheet(wb, name = "EG2df_subset")
# EG-2df step1 QQ
saveimage(wb,"graph","EG2df_subset","EG2df_A.png",1,1)

saveimage(wb,"graph","EG2df_subset","EG2df_B.png",1,11)

#EG-2df Step2 QQ
if (eg2df.s2.size>0) {
saveimage(wb,"graph","EG2df_subset","EG2df_C.png",26,1)

#EG-2df step2 Manhattan
saveimage(wb,"graph","EG2df_subset","EG2df_D.png",26,11)

setwidth( wb,"EG2df_subset",6,2600)
setwidth( wb,"EG2df_subset",8,2600)
setwidth( wb,"EG2df_subset",1,3000)

teg2df<-"EG|2df: 2-step with screening based on E-G association; subset testing of G, GxE in Step 2"
wrap(wb,teg2df, "EG2df_subset",51,1,F,F)

sigeg2df<-paste("Number of SNPs tested in Step 2= ",eg2df.s2.size,"; Step2 significance threshold = ",sig_eg2df_st2,sep="")
wrap(wb,sigeg2df, "EG2df_subset",52,1,F,F)


save2excel(wb,eg2df.s2, "EG2df_subset",54,1,T)
rm(eg2df.s2)
}

##eg2df wt
createSheet(wb, name = "EG2df_wgted")
saveimage(wb,"graph","EG2df_wgted","EG2df_E.png",1,1)

setwidth( wb,"EG2df_wgted",7,2600)
setwidth( wb,"EG2df_wgted",9,2600)
setwidth( wb,"EG2df_wgted",11,2600)
setwidth( wb,"EG2df_wgted",3,3000)


teg2dfwt<-"EG|2df: 2-step with screening based on E-G association; weighted testing of G, GxE in Step 2"
wrap(wb,teg2dfwt, "EG2df_wgted",26,1,F,F)
save2excel(wb,eg2df.wt.top, "EG2df_wgted",27,1,T)
rm(eg2df.wt.top)
}

##cocktail subset
if(RunCocktail==1){
createSheet(wb, name = "Cocktail_subset")
# Cocktail subset step1 QQ
saveimage(wb,"graph","Cocktail_subset","Cocktail_A.png",1,1)

# Cocktail subset step1 Manhattan

saveimage(wb,"graph","Cocktail_subset","Cocktail_B.png",1,11)

#Cocktail subset Step2 QQ
if (cocktail.s2.size>0) {

saveimage(wb,"graph","Cocktail_subset","Cocktail_C.png",26,1)

#Cocktail subset step2 Manhattan

saveimage(wb,"graph","Cocktail_subset","Cocktail_D.png",26,11)


setwidth( wb,"Cocktail_subset",1,3000)

tcocktail<-"Cocktail subset: 2-step with Cocktail I screening in Step 1 and subset testing in Step 2"
wrap(wb,tcocktail, "Cocktail_subset",51,1,F,F)

sigcocktail<-paste("Number of SNPs tested in Step 2= ",cocktail.s2.size,"; Step2 significance threshold = ",sig_cocktail_st2,sep="")
wrap(wb,sigcocktail, "Cocktail_subset",52,1,F,F)

save2excel(wb,cocktail.s2, "Cocktail_subset",54,1,T)
rm(cocktail.s2)
}

##cocktail wt
createSheet(wb, name = "Cocktail_wgted")
saveimage(wb,"graph","Cocktail_wgted","Cocktail_E.png",1,1)

setwidth( wb,"Cocktail_wgted",10,2600)
setwidth( wb,"Cocktail_wgted",3,3000)

tcocktailwt<-"Cocktail weighted: 2-step with Cocktail I screening in Step 1 and weighted testing in Step 2"
wrap(wb,tcocktailwt, "Cocktail_wgted",26,1,F,F)

save2excel(wb,cocktail.wt.top, "Cocktail_wgted",27,1,T)
rm(cocktail.wt.top)
}
saveWorkbook(wb)
rm(wb)
}



###### Delete pics
### pics

###### Delete pics
### pics

id <- grep("CC_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("CC_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("CO_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("CO_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("CNTL_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("CNTL_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EB_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EB_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("df2_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("df2_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("df3_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("df3_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EDGE_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EDGE_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EDGE_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EDGE_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("EDGE_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGGxE_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGGxE_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGGxE_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGGxE_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("DGGxE_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGEB_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGEB_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGEB_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DGEB_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("DGEB_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EGGxE_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EGGxE_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EGGxE_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EGGxE_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("EGGxE_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EG2df_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EG2df_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EG2df_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("EG2df_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("EG2df_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("Cocktail_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("Cocktail_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("Cocktail_C.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("Cocktail_D.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("Cocktail_E.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

id <- grep("DG_A.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)


id <- grep("DG_B.png", dir(homedir))
todelete <- dir(homedir, full.names = TRUE)[id]
unlink(todelete)

}


GWIS()
