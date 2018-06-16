## Update ##
## 08/26/2015: creat "plotall"for each individual method
#  09/18/2015: fix the problem by sorting the snps first by raw p-values then by rounded p-values ##
#              Add one line to excel describing blanks under Ranks
#              Change the reading routine to read in .snpinfo file (SNPID now is created by GxEScan)  
#              Change .gxescan to gxeout for mplink output files 







rm(list=ls(all=TRUE))
library("data.table")
options(java.parameters = "-Xmx1024m")
suppressMessages(library(XLConnect))
library(rJava)

GWIS<-function(){
  
  #####################################################################
  #args <- commandArgs(trailingOnly = TRUE)
  args<-c( "--home-dir","C:/Users/Alfred/Documents/UKB/gxescan/EduYearsHigh/",  	   # Path to GxEScan output files 
           "--alpha","0.05",                 # Family-wise error rate
           "--gxescan","EduYearsHigh_avMSE",            # Base filename of GxEScan/GxEMerge output
           "--marginal-g",                   # Output Marginal G scan, Section B.3.1
           "--GxE",                          # Output CC GxE test, Section B.3.2
           "--2df",                          # Output Case only analysis, Section 3.3.3
           ## Each of the following four lines specifies the available 2-step methods (see Section 3.2.2.4 of the documentation).  The first value (e.g. 0.05) is the step-1 significance threshold required for subset testing, and the second value (e.g. 5) is the initial bin size for weighted hypothesis testing (Section 3.2.2.5)
           "--JointGxE","0.05","5",          # Joint|GxE analysis, Section B.3.4.3
           "--YGGxE","0.05","5",             # YG|GxE, Section B.3.4.1
           "--rVarGxE","0.05","5",           # rVar|GxE, Section B.3.4.2
           ## Output options
           "--include-ranks",           # Output ranks for subset testing
           "--top","500")                # Output top 50 SNPs from each method}
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
        RunYG<- 1
      } else if(CurrentArgName=="--GxE"){
        RunGxE<- 1
      } else if(CurrentArgName=="--2df"){
        Run2df<- 1
      } else if(CurrentArgName=="--YGGxE"){
        RunYGGxE<-1
        YGGxEscreen<- as.numeric(args[j])
        YGGxE.bin<- as.integer(args[k])
      } else if(CurrentArgName=="--rVarGxE"){
        RunVarGxE<-1
        Varscreen<- as.numeric(args[j])
        VarGxE.bin<- as.integer(args[k])
      } else if(CurrentArgName=="--JointGxE"){
        RunJointGxE<-1
        JointGxEscreen<- as.numeric(args[j])
        JointGxE.bin<- as.integer(args[k])
      } else if(CurrentArgName=="--include-ranks"){
        IncludeRanks<- 1
      } else if(CurrentArgName=="--top"){
        topnum<- as.integer(args[j])
      } 
      
    }
  }
  
  if(!exists("IncludeRanks")){IncludeRanks=0}
  if(!exists("RunYG")){RunYG=0}
  if(!exists("RunGxE")){RunGxE=0}
  if(!exists("Run2df")){Run2df=0}
  if(!exists("RunVarGxE")){RunVarGxE=0}
  if(!exists("RunYGGxE")){RunYGGxE=0}
  if(!exists("RunJointGxE")){RunJointGxE=0}
  if(!exists("topnum")){topnum=50}
  
  
  
  
  ##########################################################
  
  
  ### Set Working Directory to Location of Mplink output Files
  setwd(homedir);
  
  
  ########################### Functions ###############################
  
  ### Function 1: Read in MPLINK output files and remove NAs ##
  
  ReadMplink<- function(filename){
    
    if(!exists("mytest")){
      res <-fread(paste("QT_",filename,sep=""), header=T, sep='\t',stringsAsFactors=F)
    } else {
      res <-fread(paste(mytest,"_QT_",filename,sep=""), header=T, sep='\t',stringsAsFactors=F)
    }
    
    
    # Remove missing values
    res.sub <- res[!is.na(P)]
    
    # Remove value 0
    res.sub <- res.sub[P!=0]
    
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
      abline(-1*log10(0.00000005), 0, lwd=1)
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
  
  
  
  
  ############################## GWIS Test ##################################
  ## Note: All top hits (subset testing) table are sorted by step-2 P-value
  ##       All top hits (weighted tesing) table are sorted by step-1 P-values
  ##       All Ranks are sorted by corresponding P-values   
  ###########################################################################
  
  ## Read in the SNP info ##
  
  info <-fread(paste(mytest,".snpinfo",sep=""), header=T, sep='\t',stringsAsFactors=F)
  M.info<-nrow(info)
  setkey(info,SNPID)
  
  ########## If PlotAll==1 then draw all poins in QQ and Manhattan plot ###########
  #if(M.info<1e5){
  #PlotAll=1
  #} else {
  #PlotAll=0
  #}
  
  ######### Method 1: 2df Test ###########
  
  if(Run2df==1){
    print("2 df test")
    ## Read in data
    
    df2<-ReadMplink("2df.gxeout")
    
    ## Check whether or not we have significant hits
    M.2df<-nrow(df2)
    sig_2df <- alpha/M.2df
    
    if(M.2df==0) {
      stop("All SNPs are not converged from 2df test")
    }
    
    #if(M.2df<M.info) print("NA found in 2df Mplink output file")
    
    ########## If PlotAll==1 then draw all poins in QQ and Manhattan plot ###########
    if(M.2df<1e5){
      PlotAll.2df=1
    } else {
      PlotAll.2df=0
    }
    
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
    
    df2<-merge(df2,info)
    
    df2.top<-merge(df2.top,info)
    
    
    ## Make QQ Plot
    
    png('df2_A.png')
    
    QQ(df2[,list(P)], '2df: Joint G,GxE',maxlimqq.2df,PlotAll.2df)
    
    dev.off()
    
    ## Make Manhattan Plot
    
    png('df2_B.png')
    
    manhattan(df2[,list(P,CHR,BP)], sig_2df, maxlimman.2df,0.7,PlotAll.2df,1)
    title('2df: Joint G,GxE')
    
    dev.off()
    
    setnames(df2.top,c("NMISS","CHISQ_2DF","P"),c("N","df2_Chisq","Pvalue"))
    
    rm(df2)
    
  } # end of "if(Run2df==1)"
  
  
  
  ## Standard GxE ##
  
  if( RunGxE==1 | RunYGGxE==1 |RunVarGxE==1 |RunJointGxE==1){
    
    ## Read in data
    gxe<-ReadMplink("GxE.gxeout")
    
    if(RunGxE==1){
      
      print("Exhaustive GxE scan")
      
      ## Check whether or not we have significant hits
      M.gxe<-nrow(gxe)
      sig_gxe <- alpha/M.gxe
      
      if(M.gxe==0) {
        stop("All SNPs are not converged from GxE analysis")
      }
      
      #if(M.gxe<M.info) print("NA found in GxE Mplink output file")
      
      ########## If PlotAll==1 then draw all points in QQ and Manhattan plot ###########
      if(M.gxe<1e5){
        PlotAll.gxe=1
      } else {
        PlotAll.gxe=0
      }
      
      ### find the axis limit for plots
      gxep.min <- min(gxe[,P])
      maxlimqq.gxe <- ceiling(1-log10(gxep.min))
      gxe.minman<-min(gxep.min,sig_gxe)
      maxlimman.gxe <-ceiling(1-log10(gxe.minman))
      
      ## Create GxE rank if "IncludeRanks==1"
      
      if(IncludeRanks==1){
        gxe.rank<-gxe[,list(SNPID)]
        gxe.rank[,GxE:=1:M.gxe] ## Since SNPs are sorted by P-value in the original Mplink file 
        
      }
      
      
      ## Subset to pick the top hits
      
      gxe.top<- gxe[1:topnum]
      
      
      ## Specify significance indicator for the top hits
      
      gxe.top[,SigTemp:=(gxe.top[,P]<sig_gxe)*1]
      gxe.top[SigTemp==1,Sig:="***"]
      gxe.top[SigTemp==0,Sig:=""]
      gxe.top[,SigTemp:=NULL]
      
      
      ## Merge with SNP info (need to sort by line first)
      
      setkey(gxe,SNPID)
      setkey(gxe.top,SNPID)
      
      gxe<-merge(info,gxe)
      
      gxe.top<-merge(info,gxe.top)
      
      
      ## Make QQ Plot
      
      png('GxE_A.png')
      
      QQ(gxe[,list(P)], 'GxE',maxlimqq.gxe,PlotAll.gxe)
      
      dev.off()
      
      ## Make Manhattan Plot
      
      png('GxE_B.png')
      
      manhattan(gxe[,list(P,CHR,BP)], sig_gxe, maxlimman.gxe,0.7,PlotAll.gxe,1)
      title('GxE')
      
      dev.off()
      
      setnames(gxe.top,c("NMISS","BETA","Z","P"),c("N","Beta_GxE","zTest","Pvalue"))
      
    } # end of "if(RunGxE==1)"
    
    ## If GxE is not implemented but any of the 2-step needs to be implemented
    
    if(RunGxE!=1 & (RunYGGxE==1 | RunVarGxE==1 | RunJointGxE==1  )){
      setkey(gxe,SNPID)
      gxe<-merge(info,gxe)
      
    }
    
    ## Change the column names in CC to merge with step 1 later
    setnames(gxe,c("NMISS","BETA","Z","P"),c("N","GxE_Beta","GxE_STAT","GxE_P"))
    
    ## YG+Var|GxE is implemented
    
    ## YG+Var ##
    if( RunJointGxE==1 ){
      
      print("Joint|GxE")
      
      ## Read in YG+Var
      joint<-ReadMplink("YGVH.gxeout")
      
      wt(joint,alpha,JointGxE.bin)
      
      setkey(joint,SNPID)
      
      
      ## Merge with GxE
      
      
      setnames(joint,c("CHISQ_4DF","P"),c("St1_Chisq","Step1_P"))
      
      joint<-gxe[joint,nomatch=0]
      
      M.joint<-nrow(joint)
      
      if(M.joint==0) {
        stop("All SNPs are not converged from joint test of Y-G association and residual variance heterogeneity")
      }
      
      #if(M.joint<M.info) print("NA found in joint test Mplink output file")
      
      ########## If PlotAll==1 then draw all points in QQ and Manhattan plot ###########
      if(M.joint<1e5){
        PlotAll.j1=1
      } else {
        PlotAll.j1=0
      }
      
      ## Sort Joint|GxE by step1 P-value
      setkey(joint,Step1_P)
      
      ## Subset to pick the Joint: YG+Var top hits into step 2 (subset testing)
      joint.s2<- joint[Step1_P<JointGxEscreen]
      joint.s2[,c("Bin","Threshold"):=NULL]
      joint[,SNPID:=NULL]
      joint.s2.size<-nrow(joint.s2)
      
      ## Subset to pick the top hits for weighted testing(weighted testing)
      joint[,SigTemp:=(joint[,GxE_P]<joint[,Threshold])*1]
      joint[SigTemp==1,Sig:="***"]
      joint[SigTemp==0,Sig:=""]
      joint[,SigTemp:=NULL]
      joint[,Rank:=1:M.joint]
      
      joint.wt.hits<-joint[Sig=="***"]
      joint.wt.top<-joint[1:topnum]
      
      joint.wt.top<-rbind(joint.wt.hits[!joint.wt.hits$Rank %in% joint.wt.top$Rank,],joint.wt.top)
      
      rm(joint.wt.hits)
      
      
      ## Compute limit for step 1 plots
      joint.min <- min(joint[,Step1_P])
      maxlimqq.joint <- ceiling(1-log10(joint.min))
      
      ## Step1 QQ Plot
      
      png('Joint_A.png')
      QQ(joint[,list(Step1_P)], 'Joint|GxE: Step1 screen',maxlimqq.joint,PlotAll.j1)
      dev.off()
      
      ## Step1 Manhattan Plot 
      png('Joint_B.png')
      manhattan(joint[,list(Step1_P, CHR, BP)],0,maxlimqq.joint,0.7,PlotAll.j1,0)
      title('Joint: Step1 screen')
      dev.off()
      
      joint[,c("CHR","SNP","BP","A1","N","GxE_Beta","GxE_STAT","St1_Chisq","Step1_P","Sig"):=NULL]
      
      
      if (joint.s2.size>0) {        ## If we have SNPs pass to step2
        joint.min<-min(joint.s2[,GxE_P])
        sig_joint_st2 <- alpha/joint.s2.size
        joint.min2<-min(joint.min,sig_joint_st2)
        maxlimqq.joint <- ceiling(1-log10(joint.min))
        maxlimman.joint <- ceiling(1-log10(joint.min2))
        
        if(joint.s2.size>1e5){
          PlotAll.j2=0
        } else {
          PlotAll.j2=1
        }
        
        
        ## Step2 QQ Plot
        png('Joint_C.png')
        QQ(joint.s2[,list(GxE_P)], 'Joint|GxE: Step 2 test',maxlimqq.joint,PlotAll.j2)
        dev.off()
        
        ## Step2 Manhattan Plot
        png('Joint_D.png')
        manhattan(joint.s2[,list(GxE_P, CHR, BP)], sig_joint_st2, maxlimman.joint,0.7,PlotAll.j2,1)
        title('Joint|GxE: Step 2 test')
        dev.off()
        
        
        ## Compute significance indicator for subset testing
        setkey(joint.s2,GxE_P)
        joint.s2[,SigTemp:=(joint.s2[,GxE_P]<sig_joint_st2)*1]
        joint.s2[SigTemp==1,Sig:="***"]
        joint.s2[SigTemp==0,Sig:=""]
        joint.s2[,SigTemp:=NULL]
        
        ## Create Joint|GxE rank if "IncludeRanks==1"
        if(IncludeRanks==1){
          joint.rank<-joint.s2[,list(SNPID)]
          joint.rank[,JointGxE:=1:joint.s2.size]
          joint.s2<-joint.s2[1:topnum]
        }
      } else {
        
        sig_joint_st2<-NA
        
      }
      
      ## Weighted testing
      setkey(joint,Bin)
      # Compute the limit for plot
      joint.last<-max(joint[,Bin])
      joint.k<-JointGxE.bin*(2^(joint.last-1))
      joint.wt.binsig<-alpha*((1/2)^joint.last)
      joint.wt.sig<-joint.wt.binsig/joint.k
      joint.wt.lnum<-nrow(joint[J(joint.last)])
      joint.wt.ref<- rep(-1*log10(joint.wt.sig),joint.wt.lnum)
      
      # Make plot
      png('Joint_E.png')
      
      wtedge.ref<-max(ceiling(1-log10(min(joint[,GxE_P]))),joint.wt.ref[1])
      wtplot(joint,wtedge.ref+1,"Joint|GxE: weighted",0.7,joint.wt.ref,joint.last,PlotAll.j1)
      
      dev.off()
      
      rm(joint)
      rm(joint.wt.ref)
      
      setnames(joint.s2,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      setnames(joint.wt.top,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      
      
    } # End of RunJointGxE
    
    
    ## Var|GxE is implemented
    
    if( RunVarGxE==1 ){
      
      print("rVar|GxE")
      
      ## Read in Var
      var<-ReadMplink("VH.gxeout")
      
      wt(var,alpha,VarGxE.bin)
      
      setkey(var,SNPID)
      
      
      ## Merge with GxE
      
      
      setnames(var,c("F","P"),c("St1_Lev","Step1_P"))
      
      var<-gxe[var,nomatch=0]
      
      M.var<-nrow(var)
      
      if(M.var==0) {
        stop("All SNPs are not converged from test of residual variance heterogeneity")
      }
      
      #if(M.var<M.info) print("NA found in test of variance heterogeneity Mplink output file")
      
      ########## If PlotAll==1 then draw all points in QQ and Manhattan plot ###########
      if(M.var<1e5){
        PlotAll.var1=1
      } else {
        PlotAll.var1=0
      }
      
      
      ## Sort Var|GxE by step1 P-value
      setkey(var,Step1_P)
      
      ## Subset to pick the levene's test  top hits into step 2 (subset testing)
      var.s2<- var[Step1_P<Varscreen]
      var.s2[,c("Bin","Threshold"):=NULL]
      var[,SNPID:=NULL]
      var.s2.size<-nrow(var.s2)
      
      ## Subset to pick the top hits for weighted testing(weighted testing)
      var[,SigTemp:=(var[,GxE_P]<var[,Threshold])*1]
      var[SigTemp==1,Sig:="***"]
      var[SigTemp==0,Sig:=""]
      var[,SigTemp:=NULL]
      var[,Rank:=1:M.var]
      
      var.wt.hits<-var[Sig=="***"]
      var.wt.top<-var[1:topnum]
      
      var.wt.top<-rbind(var.wt.hits[!var.wt.hits$Rank %in% var.wt.top$Rank,],var.wt.top)
      
      rm(var.wt.hits)
      
      
      ## Compute limit for step 1 plots
      var.min <- min(var[,Step1_P])
      maxlimqq.var <- ceiling(1-log10(var.min))
      
      ## Step1 QQ Plot
      
      png('levene_A.png')
      QQ(var[,list(Step1_P)], 'rVar|GxE: Step1 screen',maxlimqq.var,PlotAll.var1)
      dev.off()
      
      ## Step1 Manhattan Plot 
      png('levene_B.png')
      manhattan(var[,list(Step1_P, CHR, BP)],0,maxlimqq.var,0.7,PlotAll.var1,0)
      title('rVar|GxE: Step1 screen')
      dev.off()
      
      var[,c("CHR","SNP","BP","A1","N","GxE_Beta","GxE_STAT","St1_Lev","Step1_P","Sig"):=NULL]
      
      
      if (var.s2.size>0) {        ## If we have SNPs pass to step2
        var.min<-min(var.s2[,GxE_P])
        sig_var_st2 <- alpha/var.s2.size
        var.min2<-min(var.min,sig_var_st2)
        maxlimqq.var <- ceiling(1-log10(var.min))
        maxlimman.var <- ceiling(1-log10(var.min2))
        
        if(var.s2.size>1e5){
          PlotAll.var2=0
        } else {
          PlotAll.var2=1
        }
        
        
        ## Step2 QQ Plot
        png('levene_C.png')
        QQ(var.s2[,list(GxE_P)], 'rVar|GxE: Step 2 test',maxlimqq.var,PlotAll.var2)
        dev.off()
        
        ## Step2 Manhattan Plot
        png('levene_D.png')
        manhattan(var.s2[,list(GxE_P, CHR, BP)], sig_var_st2, maxlimman.var,0.7,PlotAll.var2,1)
        title('rVar|GxE: Step 2 test')
        dev.off()
        
        
        ## Compute significance indicator for subset testing
        setkey(var.s2,GxE_P)
        var.s2[,SigTemp:=(var.s2[,GxE_P]<sig_var_st2)*1]
        var.s2[SigTemp==1,Sig:="***"]
        var.s2[SigTemp==0,Sig:=""]
        var.s2[,SigTemp:=NULL]
        
        ## Create Var|GxE rank if "IncludeRanks==1"
        if(IncludeRanks==1){
          var.rank<-var.s2[,list(SNPID)]
          var.rank[,rVarGxE:=1:var.s2.size]
          var.s2<-var.s2[1:topnum]
        }
      } else {
        
        sig_var_st2<-NA
        
      }
      
      ## Weighted testing
      setkey(var,Bin)
      # Compute the limit for plot
      var.last<-max(var[,Bin])
      var.k<-JointGxE.bin*(2^(var.last-1))
      var.wt.binsig<-alpha*((1/2)^var.last)
      var.wt.sig<-var.wt.binsig/var.k
      var.wt.lnum<-nrow(var[J(var.last)])
      var.wt.ref<- rep(-1*log10(var.wt.sig),var.wt.lnum)
      
      # Make plot
      png('levene_E.png')
      
      wtedge.ref<-max(ceiling(1-log10(min(var[,GxE_P]))),var.wt.ref[1])
      wtplot(var,wtedge.ref+1,"rVar|GxE: weighted",0.7,var.wt.ref,var.last,PlotAll.var1)
      
      dev.off()
      
      rm(var)
      rm(var.wt.ref)
      
      setnames(var.s2,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      setnames(var.wt.top,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      
    } # End of RunVarGxE==1
    
    if(RunYGGxE!=1){
      rm(gxe)
    } 
    
  } # End of if( RunGxE==1 | RunYGGxE==1 |RunVarGxE==1 |RunJointGxE==1)
  
  ## Y-G ##
  
  if( RunYG==1 | RunYGGxE==1){
    
    ## Read in data
    yg<-ReadMplink("YG.gxeout")
    
    if(RunYG==1){
      print("Marginal G scan")
      ## Check whether or not we have significant hits
      M.yg<-nrow(yg)
      sig_yg <- alpha/M.yg
      
      if(M.yg==0) {
        stop("All SNPs are not converged from test of Y-G association")
      }
      
      #if(M.yg<M.info) print("NA found in marginal G Mplink output file")
      
      ########## If PlotAll==1 then draw all points in QQ and Manhattan plot ###########
      if(M.yg<1e5){
        PlotAll.yg=1
      } else {
        PlotAll.yg=0
      }
      
      
      
      ### find the axis limit for plots
      ygp.min <- min(yg[,P])
      maxlimqq.yg <- ceiling(1-log10(ygp.min))
      yg.minman<-min(ygp.min,sig_yg)
      maxlimman.yg <-ceiling(1-log10(yg.minman))
      
      ## Create YG rank if "IncludeRanks==1"
      
      if(IncludeRanks==1){
        yg.rank<-yg[,list(SNPID)]
        yg.rank[,YG:=1:M.yg] ## Since SNPs are sorted by P-value in the original Mplink file 
        
      }
      
      
      ## Subset to pick the top hits
      
      yg.top<- yg[1:topnum]
      
      ## Specify significance indicator for the top hits
      
      yg.top[,SigTemp:=(yg.top[,P]<sig_yg)*1]
      yg.top[SigTemp==1,Sig:="***"]
      yg.top[SigTemp==0,Sig:=""]
      yg.top[,SigTemp:=NULL]
      
      
      ## Merge with SNP info (need to sort by line first)
      
      setkey(yg,SNPID)
      setkey(yg.top,SNPID)
      
      yg<-merge(info,yg)
      
      yg.top<-merge(info,yg.top)
      
      
      ## Make QQ Plot 
      
      png('YG_A.png')
      
      QQ(yg[,list(P)], 'YG: Marginal G Scan',maxlimqq.yg,PlotAll.yg)
      
      dev.off()
      
      ## Make Manhattan Plot
      
      png('YG_B.png')
      
      manhattan(yg[,list(P,CHR,BP)], sig_yg, maxlimman.yg,0.7,PlotAll.yg,1)
      title('YG: Marginal G Scan')
      
      dev.off()
      
      setnames(yg.top,c("NMISS","BETA","Z","P"),c("N","Beta","zTest","Pvalue"))
      
      setkey(yg,P)
      
      yg[,c("CHR","SNP","BP","A1","NMISS","BETA"):=NULL]
      
    } # end of "if(RunYG==1)"
    
    if(RunYG!=1 & RunYGGxE==1){
      yg[,c("NMISS","BETA"):=NULL]
    }
    
    setnames(yg,c("Z","P"),c("St1_zTest","Step1_P"))
    
    ## YG|GxE is implemented
    
    if( RunYGGxE==1 ){
      print("YG|GxE")
      wt(yg,alpha,YGGxE.bin)
      
      setkey(yg,SNPID)
      
      ## Merge with GxE
      
      ygGxE<-gxe[yg,nomatch=0]
      
      M.ygGxE<-nrow(ygGxE)
      
      ## Sort YG|GxE by step1 P-value
      setkey(ygGxE,Step1_P)
      
      ## Subset to pick the D-G top hits into step 2 (subset testing)
      ygGxE.s2<- ygGxE[Step1_P<YGGxEscreen]
      ygGxE.s2[,c("Bin","Threshold"):=NULL]
      ygGxE[,SNPID:=NULL]
      ygGxE.s2.size<-nrow(ygGxE.s2)
      
      ## Subset to pick the top hits for weighted testing(weighted testing)
      ygGxE[,SigTemp:=(ygGxE[,GxE_P]<ygGxE[,Threshold])*1]
      ygGxE[SigTemp==1,Sig:="***"]
      ygGxE[SigTemp==0,Sig:=""]
      ygGxE[,SigTemp:=NULL]
      ygGxE[,Rank:=1:M.ygGxE]
      
      ygGxE.wt.hits<-ygGxE[Sig=="***"]
      ygGxE.wt.top<-ygGxE[1:topnum]
      
      ygGxE.wt.top<-rbind(ygGxE.wt.hits[!ygGxE.wt.hits$Rank %in% ygGxE.wt.top$Rank,],ygGxE.wt.top)
      
      rm(ygGxE.wt.hits)
      
      
      if(RunYG!=1){
        ## Compute limit for step 1 plots
        ygp.min <- min(ygGxE[,Step1_P])
        maxlimqq.yg <- ceiling(1-log10(ygp.min))
      }
      
      ## Step1 QQ Plot 
      
      png('YGGxE_A.png')
      QQ(ygGxE[,list(Step1_P)], 'YG | GxE: Step1 screen',maxlimqq.yg,PlotAll.yg)
      dev.off()
      
      ## Step1 Manhattan Plot ( will be used for both EG|GxE and EG|2df )
      png('YGGxE_B.png')
      manhattan(ygGxE[,list(Step1_P, CHR, BP)],0,maxlimqq.yg,0.7,PlotAll.yg,0)
      title('YG | GxE: Step1 screen')
      dev.off()
      
      
      ygGxE[,c("CHR","SNP","BP","A1","N","GxE_Beta","GxE_STAT","St1_zTest","Step1_P","Sig"):=NULL]
      
      
      if (ygGxE.s2.size>0) {        ## If we have SNPs pass to step2
        ygGxE.min<-min(ygGxE.s2[,GxE_P])
        sig_ygGxE_st2 <- alpha/ygGxE.s2.size
        ygGxE.min2<-min(ygGxE.min,sig_ygGxE_st2)
        maxlimqq.ygGxE <- ceiling(1-log10(ygGxE.min))
        maxlimman.ygGxE <- ceiling(1-log10(ygGxE.min2))
        
        if(ygGxE.s2.size>1e5){
          PlotAll.ygGxE=0
        } else {
          PlotAll.ygGxE=1
        }
        
        ## Step2 QQ Plot
        png('YGGxE_C.png')
        QQ(ygGxE.s2[,list(GxE_P)], 'YG | GxE: Step 2 test',maxlimqq.ygGxE,PlotAll.ygGxE)
        dev.off()
        
        ## Step2 Manhattan Plot
        png('YGGxE_D.png')
        manhattan(ygGxE.s2[,list(GxE_P, CHR, BP)], sig_ygGxE_st2, maxlimman.ygGxE,0.7,PlotAll.ygGxE,1)
        title('YG | GxE: Step 2 test')
        dev.off()
        
        
        ## Compute significance indicator for subset testing
        setkey(ygGxE.s2,GxE_P)
        ygGxE.s2[,SigTemp:=(ygGxE.s2[,GxE_P]<sig_ygGxE_st2)*1]
        ygGxE.s2[SigTemp==1,Sig:="***"]
        ygGxE.s2[SigTemp==0,Sig:=""]
        ygGxE.s2[,SigTemp:=NULL]
        
        ## Create YG|GxE rank if "IncludeRanks==1"
        if(IncludeRanks==1){
          ygGxE.rank<-ygGxE.s2[,list(SNPID)]
          ygGxE.rank[,YGGxE:=1:ygGxE.s2.size]
          ygGxE.s2<-ygGxE.s2[1:topnum]
        }
      } else {
        
        sig_ygGxE_st2<-NA
        
      }
      
      
      ## Weighted testing
      setkey(ygGxE,Bin)
      # Compute the limit for plot
      ygGxE.last<-max(ygGxE[,Bin])
      ygGxE.k<-YGGxE.bin*(2^(ygGxE.last-1))
      ygGxE.wt.binsig<-alpha*((1/2)^ygGxE.last)
      ygGxE.wt.sig<-ygGxE.wt.binsig/ygGxE.k
      ygGxE.wt.lnum<-nrow(ygGxE[J(ygGxE.last)])
      ygGxE.wt.ref<- rep(-1*log10(ygGxE.wt.sig),ygGxE.wt.lnum)
      
      # Make plot
      png('YGGxE_E.png')
      
      wtygGxE.ref<-max(ceiling(1-log10(min(ygGxE[,GxE_P]))),ygGxE.wt.ref[1])
      wtplot(ygGxE,wtygGxE.ref+1,"YG | GxE: weighted",0.7,ygGxE.wt.ref,ygGxE.last,PlotAll.yg)
      
      dev.off()
      
      rm(ygGxE)
      rm(ygGxE.wt.ref)
      rm(gxe)
      
      setnames(ygGxE.s2,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      setnames(ygGxE.wt.top,c("GxE_Beta","GxE_STAT","GxE_P"),c("Beta_GxE","St2_zTest","Step2_P"))
      
      yg[,c("Bin","Threshold"):=NULL]
      
    } # End of RunYGGxE
    
    rm(yg)
    
    
  } # End of if( RunYG==1 | RunYGGxE==1 )
  
  if(IncludeRanks==1){
    ## Rank contains ranks from all methods implemented
    Rank<-data.table(info[,SNPID])
    setnames(Rank,1,"SNPID")
    setkey(Rank,SNPID)
  }
  
  rm(info)
  
  
  rm(manhattan,QQ,wtplot)
  ##########################################
  #
  #               End of GWIS
  #
  ##########################################
  print("Done testing...formatting output")
  
  ## Format the output columns (Create combined Rank)##
  
  ## Var|GxE
  
  if (RunVarGxE==1){
    if(IncludeRanks==1){
      setkey(var.rank,SNPID)
      Rank<-var.rank[Rank]
      rm(var.rank)
    }
    var.s2[,Beta_GxE:=as.numeric(formatC(var.s2[,Beta_GxE],digits=2,format="f"))]
    var.s2[,St1_Lev:=as.numeric(formatC(var.s2[,St1_Lev],digits=2,format="f"))]
    var.s2[,St2_zTest:=as.numeric(formatC(var.s2[,St2_zTest],digits=2,format="f"))]
    var.s2[,Step1_P:=signif(var.s2[,Step1_P],2)]
    var.s2[,Step2_P:=signif(var.s2[,Step2_P],2)]
    
    var.wt.top[,Beta_GxE:=as.numeric(formatC(var.wt.top[,Beta_GxE],digits=2,format="f"))]
    var.wt.top[,St1_Lev:=as.numeric(formatC(var.wt.top[,St1_Lev],digits=2,format="f"))]
    var.wt.top[,St2_zTest:=as.numeric(formatC(var.wt.top[,St2_zTest],digits=2,format="f"))]
    var.wt.top[,Step1_P:=signif(var.wt.top[,Step1_P],2)]
    var.wt.top[,Step2_P:=signif(var.wt.top[,Step2_P],2)]
    
  }
  
  ## YG|GxE
  
  if (RunYGGxE==1){
    if(IncludeRanks==1){
      setkey(ygGxE.rank,SNPID)
      Rank<-ygGxE.rank[Rank]
      rm(ygGxE.rank)
    }
    ygGxE.s2[,Beta_GxE:=as.numeric(formatC(ygGxE.s2[,Beta_GxE],digits=2,format="f"))]
    ygGxE.s2[,St1_zTest:=as.numeric(formatC(ygGxE.s2[,St1_zTest],digits=2,format="f"))]
    ygGxE.s2[,St2_zTest:=as.numeric(formatC(ygGxE.s2[,St2_zTest],digits=2,format="f"))]
    ygGxE.s2[,Step1_P:=signif(ygGxE.s2[,Step1_P],2)]
    ygGxE.s2[,Step2_P:=signif(ygGxE.s2[,Step2_P],2)]
    
    ygGxE.wt.top[,Beta_GxE:=as.numeric(formatC(ygGxE.wt.top[,Beta_GxE],digits=2,format="f"))]
    ygGxE.wt.top[,St1_zTest:=as.numeric(formatC(ygGxE.wt.top[,St1_zTest],digits=2,format="f"))]
    ygGxE.wt.top[,St2_zTest:=as.numeric(formatC(ygGxE.wt.top[,St2_zTest],digits=2,format="f"))]
    ygGxE.wt.top[,Step1_P:=signif(ygGxE.wt.top[,Step1_P],2)]
    ygGxE.wt.top[,Step2_P:=signif(ygGxE.wt.top[,Step2_P],2)]
    
  }
  
  
  ## Joint|GxE
  
  if (RunJointGxE==1){
    if(IncludeRanks==1){
      setkey(joint.rank,SNPID)
      Rank<-joint.rank[Rank]
      rm(joint.rank)
    }
    joint.s2[,Beta_GxE:=as.numeric(formatC(joint.s2[,Beta_GxE],digits=2,format="f"))]
    joint.s2[,St1_Chisq:=as.numeric(formatC(joint.s2[,St1_Chisq],digits=2,format="f"))]
    joint.s2[,St2_zTest:=as.numeric(formatC(joint.s2[,St2_zTest],digits=2,format="f"))]
    joint.s2[,Step1_P:=signif(joint.s2[,Step1_P],2)]
    joint.s2[,Step2_P:=signif(joint.s2[,Step2_P],2)]
    
    joint.wt.top[,Beta_GxE:=as.numeric(formatC(joint.wt.top[,Beta_GxE],digits=2,format="f"))]
    joint.wt.top[,St1_Chisq:=as.numeric(formatC(joint.wt.top[,St1_Chisq],digits=2,format="f"))]
    joint.wt.top[,St2_zTest:=as.numeric(formatC(joint.wt.top[,St2_zTest],digits=2,format="f"))]
    joint.wt.top[,Step1_P:=signif(joint.wt.top[,Step1_P],2)]
    joint.wt.top[,Step2_P:=signif(joint.wt.top[,Step2_P],2)]
    
  }
  
  
  ## GxE
  if (RunGxE==1){
    if(IncludeRanks==1){
      setkey(gxe.rank,SNPID)
      Rank<-gxe.rank[Rank]
      rm(gxe.rank)
    }
    gxe.top[,Beta_GxE:=as.numeric(formatC(gxe.top[,Beta_GxE],digits=2,format="f"))]
    gxe.top[,zTest:= as.numeric(formatC(gxe.top[,zTest],digits=2,format="f"))]
    gxe.top[,Pvalue:=signif(gxe.top[,Pvalue],2)]
    
  }
  
  
  ## Mar G
  if (RunYG==1){
    if(IncludeRanks==1){
      setkey(yg.rank,SNPID)
      #setnames(yg.rank,"YG","MarG")
      Rank<-yg.rank[Rank]
      rm(yg.rank)
    }
    yg.top[,Beta:=as.numeric(formatC(yg.top[,Beta],digits=2,format="f"))]
    yg.top[,zTest:= as.numeric(formatC(yg.top[,zTest],digits=2,format="f"))]
    yg.top[,Pvalue:=signif(yg.top[,Pvalue],2)]
    
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
  
  
  
  ## Add Ranks if checked (Make sure the order here matches the layout in the excel file)
  ImpMethod<- c("YG","GxE","df2","JointGxE","YGGxE","rVarGxE")
  ImpMethod<-ImpMethod[c(RunYG==1,RunGxE==1,Run2df==1,RunJointGxE==1,RunYGGxE==1,RunVarGxE==1)]
  
  
  ##MarG
  if(RunYG==1){
    
    if(IncludeRanks==1){
      ImpMethod.marG<-ImpMethod[-which(ImpMethod=="YG")]
      yg.top<-Rank[yg.top]
      yg.top[,SNPID:=NULL]
      setcolorder(yg.top,c("SNP","CHR","BP","A1","N","Beta","zTest","Pvalue","Sig","YG",ImpMethod.marG))
      setkey(yg.top,YG)
    } else {
      yg.top[,SNPID:=NULL]
      setcolorder(yg.top,c("SNP","CHR","BP","A1","N","Beta","zTest","Pvalue","Sig"))
      setkey(yg.top,Pvalue)
    }
  }
  
  
  ##GxE
  
  if(RunGxE==1){
    if(IncludeRanks==1){
      ImpMethod.gxe<-ImpMethod[-which(ImpMethod=="GxE")]
      gxe.top<-Rank[gxe.top]
      gxe.top[,SNPID:=NULL]
      setcolorder(gxe.top,c("SNP","CHR","BP","A1","N","Beta_GxE","zTest","Pvalue","Sig","GxE",ImpMethod.gxe))
      setkey(gxe.top,GxE)
    } else {
      gxe.top[,SNPID:=NULL]
      setcolorder(gxe.top,c("SNP","CHR","BP","A1","N","Beta_GxE","zTest","Pvalue","Sig"))
      setkey(gxe.top,Pvalue)
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
  
  
  
  ## YG|GxE
  
  if(RunYGGxE==1){
    if(IncludeRanks==1){
      ImpMethod.ygGxE<-ImpMethod[-which(ImpMethod=="YGGxE")]
      setkey(ygGxE.s2,SNPID)
      ygGxE.s2<-Rank[ygGxE.s2]
      ygGxE.s2[,SNPID:=NULL]
      setcolorder(ygGxE.s2,c("SNP","CHR","BP","A1","N","St1_zTest","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig","YGGxE",ImpMethod.ygGxE))
      setkey(ygGxE.s2,YGGxE)
    } else {
      ygGxE.s2[,SNPID:=NULL]
      setcolorder(ygGxE.s2,c("SNP","CHR","BP","A1","N","St1_zTest","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig"))
      setkey(ygGxE.s2,Step2_P)
    }
    
    
    setcolorder(ygGxE.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_zTest","Step1_P","Beta_GxE","St2_zTest","Step2_P","Threshold","Bin","Sig"))
    
    setkey(ygGxE.wt.top,Step1_P)
    
  }
  
  
  
  
  ## Joint|GxE
  
  if(RunJointGxE==1){
    if(IncludeRanks==1){
      ImpMethod.joint<-ImpMethod[-which(ImpMethod=="JointGxE")]
      setkey(joint.s2,SNPID)
      joint.s2<-Rank[joint.s2]
      joint.s2[,SNPID:=NULL]
      setcolorder(joint.s2,c("SNP","CHR","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig","JointGxE",ImpMethod.joint))
      setkey(joint.s2,JointGxE)
    } else {
      joint.s2[,SNPID:=NULL]
      setcolorder(joint.s2,c("SNP","CHR","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig"))
      setkey(joint.s2,Step2_P)
    }
    
    setcolorder(joint.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_Chisq","Step1_P","Beta_GxE","St2_zTest","Step2_P","Threshold","Bin","Sig"))
    
    setkey(joint.wt.top,Step1_P)
    
  }
  
  
  
  ## Var|GxE
  
  if(RunVarGxE==1){
    if(IncludeRanks==1){
      ImpMethod.var<-ImpMethod[-which(ImpMethod=="rVarGxE")]
      setkey(var.s2,SNPID)
      var.s2<-Rank[var.s2]
      var.s2[,SNPID:=NULL]
      setcolorder(var.s2,c("SNP","CHR","BP","A1","N","St1_Lev","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig","rVarGxE",ImpMethod.var))
      setkey(var.s2,rVarGxE)
    } else {
      var.s2[,SNPID:=NULL]
      setcolorder(var.s2,c("SNP","CHR","BP","A1","N","St1_Lev","Step1_P","Beta_GxE","St2_zTest","Step2_P","Sig"))
      setkey(var.s2,Step2_P)
    }
    
    setcolorder(var.wt.top,c("Rank","CHR","SNP","BP","A1","N","St1_Lev","Step1_P","Beta_GxE","St2_zTest","Step2_P","Threshold","Bin","Sig"))
    
    setkey(var.wt.top,Step1_P)
    
  }
  
  if(IncludeRanks==1){
    rm(Rank)
  }
  
  
  ############################################
  ## 
  ##         Write results to excel
  ##
  ############################################
  
  
  FileName<-paste(mytest,"_output.xlsx",sep="")
  
  # Delete file with the same name
  id <- grep(FileName, dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  ####################### Save Test Summary ########################
  
  print("Writing output to Excel file...")
  
  num.method<- sum(RunYG,RunGxE,Run2df,RunJointGxE,RunYGGxE,RunVarGxE)
  num.exaus<-sum(RunYG,RunGxE,Run2df)
  num.2step<-sum(RunJointGxE,RunYGGxE,RunVarGxE)
  
  diff.st1<- 3-num.exaus
  diff.st2<- 6-num.method
  
  ## Create the excel file 
  wb <- loadWorkbook(FileName, create = TRUE)
  createSheet(wb, name = "summary")
  
  setwidth( wb,"summary",1,2600)
  
  ##Exhaustive tests
  
  ## Total number of SNPs (M)
  M <- c(NA)
  
  if(RunYG==1){
    M<-cbind(M,M.yg)
  }
  
  if(RunGxE==1){
    M<-cbind(M,M.gxe)
  }
  
  if(Run2df==1){
    M<-cbind(M,M.2df)
  }
  
  M<-M[-1]
  
  M<-data.frame(M)
  
  
  ## Genomewide significance level
  
  sig1 <- c(NA)
  
  if(RunYG==1){
    sig1<-cbind(sig1,sig_yg)
  }
  
  if(RunGxE==1){
    sig1<-cbind(sig1,sig_gxe)
  }
  
  if(Run2df==1){
    sig1<-cbind(sig1,sig_2df)
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
  
  if(RunYG==1){
    test.exhaus.yg<- "YG"
    test.exhaus<-cbind(test.exhaus,test.exhaus.yg)
  }
  
  if(RunGxE==1){
    test.exhaus.gxe<- "GxE"
    test.exhaus<-cbind(test.exhaus,test.exhaus.gxe)
  }
  
  if(Run2df==1){
    test.exhaus.2df<- "2df"
    test.exhaus<-cbind(test.exhaus,test.exhaus.2df)
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
    
    if(RunJointGxE==1){
      test.2stp.joint<- "Joint|GxE"
      test.2stp<-cbind(test.2stp,test.2stp.joint)
    }
    
    if(RunYGGxE==1){
      test.2stp.ygGxE<- "YG|GxE"
      test.2stp<-cbind(test.2stp,test.2stp.ygGxE)
    }
    
    
    if(RunVarGxE==1){
      test.2stp.var<- "rVar|GxE"
      test.2stp<-cbind(test.2stp,test.2stp.var)
    }
    
    
    test.2stp<-test.2stp[-1]
    test.2stp<-data.frame(test.2stp)
    colnames(test.2stp)<-"Method"
    
    ## Total number of SNPs (M)
    M <- c(NA)
    
    if(RunJointGxE==1){
      M<-cbind(M,M.joint)
    }
    
    if(RunYGGxE==1){
      M<-cbind(M,M.ygGxE)
    }
    
    if(RunVarGxE==1){
      M<-cbind(M,M.var)
    }
    
    M<-M[-1]
    
    M<-data.frame(M)
    
    ## Step1 significance level
    s1_alpha1<-c(NA)
    
    if(RunJointGxE==1){
      s1_alpha1.joint<-signif(JointGxEscreen,3)
      s1_alpha1<-cbind(s1_alpha1,s1_alpha1.joint)
    }
    
    if(RunYGGxE==1){
      s1_alpha1.yg<-signif(YGGxEscreen,3)
      s1_alpha1<-cbind(s1_alpha1,s1_alpha1.yg)
    }
    
    if(RunVarGxE==1){
      s1_alpha1.var<-signif(Varscreen,3)
      s1_alpha1<-cbind(s1_alpha1,s1_alpha1.var)
    }
    
    s1_alpha1<-s1_alpha1[-1]
    
    s1_alpha1<-data.frame(s1_alpha1)
    
    colnames(s1_alpha1) <- "alpha1"
    
    ## Number of SNPs pass to step 2
    
    ns2<-c(NA)
    
    if(RunJointGxE==1){
      ns2<-cbind(ns2,joint.s2.size)
    }
    
    if(RunYGGxE==1){
      ns2<-cbind(ns2,ygGxE.s2.size)
    }
    
    if(RunVarGxE==1){
      ns2<-cbind(ns2,var.s2.size)
    }
    
    ns2<- ns2[-1]
    ns2 <- data.frame(ns2)
    colnames(ns2) <- "m"
    
    ## significance level at step 2
    
    sig2<-c(NA)
    
    if(RunJointGxE==1){
      sig2<-cbind(sig2,signif(sig_joint_st2,3))
    }
    
    if(RunYGGxE==1){
      sig2<-cbind(sig2,signif(sig_ygGxE_st2,3))
    }
    
    if(RunVarGxE==1){
      sig2<-cbind(sig2,signif(sig_var_st2,3))
    }
    
    sig2<-sig2[-1]
    sig2 <- data.frame(sig2)
    colnames(sig2) <- "sig2"
    
    ## Bin size for each method
    binsize<-c(NA)
    
    if(RunJointGxE==1){
      binsize<-cbind(binsize,JointGxE.bin)
    }
    
    if(RunYGGxE==1){
      binsize<-cbind(binsize,YGGxE.bin)
    }
    
    if(RunVarGxE==1){
      binsize<-cbind(binsize,VarGxE.bin)
    }
    
    binsize<-binsize[-1]
    binsize<-data.frame(binsize)
    colnames(binsize)<-"Bin Size"
    
    step2 <- data.frame(test.2stp,M,s1_alpha1,ns2,sig2,binsize)
    colnames(step2)<-c("Method","M","alpha1","m","alpha/m","Bin Size")
    
    te1.row<- 7-diff.st1
    te1.row2<- te1.row+1
    wrap(wb,te1, "summary",te1.row,1,F,F)
    save2excel(wb,step2, "summary",te1.row2,1,T)
    rm(step2)
  }
  
  #####Descriptions
  
  if(num.2step>0){
    te2.row<- 13-diff.st2
  } else {
    te2.row<- 7-diff.st2
  }
  
  
  te2<-"M:  total number of SNPs tested"
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
  te00<-"YG:  Standard GWAS analysis of marginal G effect for all M SNPs"
  wrap(wb,te00, "summary",te00.row,1,F,F)
  
  
  te8.row<- te00.row+1
  te8<-"GxE:  Standard analysis of GxE interaction for all M SNPs"
  wrap(wb,te8, "summary",te8.row,1,F,F)
  
  te9.row<- te8.row+1
  te9<-"2df: Two degree of freedom joint test of G,GxE for all M SNPs( Kraft et al.2007)"
  wrap(wb,te9, "summary",te9.row,1,F,F)
  
  
  te1.2.row<- te9.row+2
  wrap(wb,te1, "summary",te1.2.row,1,F,F)
  
  te11.row<- te1.2.row+1
  te11<-"Joint|GxE: 2-step with Step-1 screen based on joint test of Y-G association and residual variance heterogeneity, Step-2 test of GxE using GxE analysis"
  wrap(wb,te11, "summary",te11.row,1,F,F)
  
  te12.row<- te11.row+1
  te12<-"YG|GxE: 2-step with Step-1 screen based on Y-G association, Step-2 test of GxE using GxE analysis (Kooperberg and LeBlanc, 2009)"
  wrap(wb,te12, "summary",te12.row,1,F,F)
  
  te13.row<- te12.row+1
  te13<-"rVar|GxE: 2-step with Step-1 screen based on residual variance heterogeneity, Step-2 test of GxE using GxE analysis"
  wrap(wb,te13, "summary",te13.row,1,F,F)
  
  rm(te0,te00,te2,te3,te4,te5,te6,te7,te8,te9,te11,te12,te13)
  
  
  ##################### Write out top hits and plots #################################3
  
  alphabet<- c("I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  
  
  if(IncludeRanks==1){
    
    ## MarG
    if(RunYG==1){
      createSheet(wb, name = "YG")
      ########## 1st panel
      
      # Marginal G QQ
      
      saveimage(wb,"graph","YG","YG_A.png",1,1)
      
      # Marginal G Manhattan
      
      saveimage(wb,"graph","YG","YG_B.png",1,11)
      
      setwidth( wb,"YG",1,3000)
      #setwidth( wb,"MarG",6,2600)
      
      tMarG<-"YG:  Standard GWAS analysis of marginal G effect for all M SNPs"
      wrap(wb,tMarG, "YG",26,1,F,F)
      
      sigMarG<-paste("Number of SNPs tested = ",M.yg,"; significance threshold = ",sig_yg,sep="")
      wrap(wb,sigMarG, "YG",27,1,F,F)
      
      RankDescrip<-"A blank under Ranks means this SNP does not pass through the screening step for the corresponding 2-step approach"
      wrap(wb,RankDescrip, "YG",28,1,F,F)
      
      rMarG<-"Ranks"
      save2excel(wb,rMarG, "YG",29,10,F)
      
      MarG.colend<- alphabet[1+num.method]
      MarG.letter<- paste("J29:",MarG.colend,"29",sep="")
      MarG.col<- 9+num.method
      savetitle( wb,"YG",MarG.letter,29,c(10:MarG.col))
      save2excel(wb,yg.top, "YG",30,1,T)
      rm(yg.top)
    }
    
    ## GxE
    if(RunGxE==1){
      createSheet(wb, name = "GxE")
      
      ########## 2nd panel
      
      # G x E QQ
      
      saveimage(wb,"graph","GxE","GxE_A.png",1,1)
      
      # GxE Manhattan
      
      saveimage(wb,"graph","GxE","GxE_B.png",1,11)
      
      setwidth( wb,"GxE",1,3000)
      setwidth( wb,"GxE",6,2600)
      
      tgxe<-"GxE:  Standard analysis of GxE interaction for all M SNPs"
      wrap(wb,tgxe, "GxE",26,1,F,F)
      
      siggxe<-paste("Number of SNPs tested = ",M.gxe,"; significance threshold = ",sig_gxe,sep="")
      wrap(wb,siggxe, "GxE",27,1,F,F)
      
      wrap(wb,RankDescrip, "GxE",28,1,F,F)
      
      rgxe<-"Ranks"
      save2excel(wb,rgxe, "GxE",29,10,F)
      
      gxe.colend<- alphabet[1+num.method]
      gxe.letter<- paste("J29:",gxe.colend,"29",sep="")
      gxe.col<- 9+num.method
      savetitle( wb,"GxE",gxe.letter,29,c(10:gxe.col))
      save2excel(wb,gxe.top, "GxE",30,1,T)
      rm(gxe.top)
    }
    
    
    ##2df
    
    if(Run2df==1){
      createSheet(wb, name = "df2")
      
      ####### 3rd panel
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
    
    
    
    ## Joint|GxE
    
    if(RunJointGxE==1){
      createSheet(wb, name = "JointGxE_subset")
      
      ########## 5th panel
      
      saveimage(wb,"graph","JointGxE_subset","Joint_A.png",1,1)
      
      # Joint step1 Manhattan
      
      saveimage(wb,"graph","JointGxE_subset","Joint_B.png",1,11)
      
      # Joint step2 QQ
      
      if (joint.s2.size>0) {
        saveimage(wb,"graph","JointGxE_subset","Joint_C.png",26,1)
        
        # Joint step2 Manhattan
        saveimage(wb,"graph","JointGxE_subset","Joint_D.png",26,11)
        
        
        setwidth( wb,"JointGxE_subset",6,2600)
        setwidth( wb,"JointGxE_subset",8,2600)
        setwidth( wb,"JointGxE_subset",9,2600)
        setwidth( wb,"JointGxE_subset",1,3000)
        
        tjoint<-"JointGxE: 2-step with screening based on joint test of Y-G association and variance heterogeneity; subset testing of GxE in Step 2"
        wrap(wb,tjoint, "JointGxE_subset",51,1,F,F)
        
        sigjoint<-paste("Number of SNPs tested in Step 2= ",joint.s2.size,"; Step2 significance threshold = ",sig_joint_st2,sep="")
        wrap(wb,sigjoint, "JointGxE_subset",52,1,F,F)
        
        wrap(wb,RankDescrip, "JointGxE_subset",53,1,F,F)
        
        rjoint<-"Ranks"
        save2excel(wb,rjoint, "JointGxE_subset",54,12,F)
        joint.colend<- alphabet[3+num.method]
        joint.letter<- paste("L54:",joint.colend,"54",sep="")
        joint.col<- 11+num.method
        savetitle( wb,"JointGxE_subset",joint.letter,54,c(12:joint.col))
        
        save2excel(wb,joint.s2, "JointGxE_subset",55,1,T)
        rm(joint.s2)
      }
      
      # JointGxE wgted
      
      ########## 6th panel
      createSheet(wb, name = "JointGxE_wgted")
      saveimage(wb,"graph","JointGxE_wgted","Joint_E.png",1,1)
      
      ##dgeg wt
      setwidth( wb,"JointGxE_wgted",7,2600)
      setwidth( wb,"JointGxE_wgted",9,2600)
      setwidth( wb,"JointGxE_wgted",10,2600)
      setwidth( wb,"JointGxE_wgted",12,2600)
      setwidth( wb,"JointGxE_wgted",3,3000)
      
      tjointwt<-"JointGxE: 2-step with screening based on joint test of Y-G association and variance heterogeneity; weighted testing of GxE in Step 2"
      wrap(wb,tjointwt, "JointGxE_wgted",26,1,F,F)
      save2excel(wb,joint.wt.top, "JointGxE_wgted",27,1,T)
      rm(joint.wt.top)
    }
    
    
    ## YG|GxE
    
    if(RunYGGxE==1){
      createSheet(wb, name = "YGGxE_subset")
      
      ####### 7th panel
      # Y-G step1 QQ
      saveimage(wb,"graph","YGGxE_subset","YGGxE_A.png",1,1)
      
      # Y-G step1 Manhattan
      saveimage(wb,"graph","YGGxE_subset","YGGxE_B.png",1,11)
      
      ##### Y-G step2 QQ
      if (ygGxE.s2.size>0) {
        saveimage(wb,"graph","YGGxE_subset","YGGxE_C.png",26,1)
        
        ###### Y-G step2 Manhattan
        
        saveimage(wb,"graph","YGGxE_subset","YGGxE_D.png",26,11)
        
        
        setwidth( wb,"YGGxE_subset",6,2600)
        setwidth( wb,"YGGxE_subset",8,2600)
        setwidth( wb,"YGGxE_subset",9,2600)
        setwidth( wb,"YGGxE_subset",1,3000)
        
        tyg<-"YG|GxE: 2-step with screening based on Y-G association; subset testing of GxE in Step 2"
        wrap(wb,tyg, "YGGxE_subset",51,1,F,F)
        
        sigyg<-paste("Number of SNPs tested in Step 2= ",ygGxE.s2.size,"; Step2 significance threshold = ",sig_ygGxE_st2,sep="")
        wrap(wb,sigyg, "YGGxE_subset",52,1,F,F)
        
        wrap(wb,RankDescrip, "YGGxE_subset",53,1,F,F)
        
        ryg<-"Ranks"
        save2excel(wb,ryg, "YGGxE_subset",54,12,F)
        yg.colend<- alphabet[3+num.method]
        yg.letter<- paste("L54:",yg.colend,"54",sep="")
        yg.col<- 11+num.method
        savetitle( wb,"YGGxE_subset",yg.letter,54,c(12:yg.col))
        
        save2excel(wb,ygGxE.s2, "YGGxE_subset",55,1,T)
        rm(ygGxE.s2)
      }
      
      
      ## YG|GxE wt
      
      ####### 8th panel
      
      ##### YG wgted 
      createSheet(wb, name = "YGGxE_wgted")
      saveimage(wb,"graph","YGGxE_wgted","YGGxE_E.png",1,1)
      
      setwidth( wb,"YGGxE_wgted",12,2800)
      setwidth( wb,"YGGxE_wgted",7,2600)
      setwidth( wb,"YGGxE_wgted",9,2600)
      setwidth( wb,"YGGxE_wgted",10,2600)
      setwidth( wb,"YGGxE_wgted",3,3000)
      
      tygwt<-"YG|GxE: 2-step with screening based on Y-G association; weighted testing of GxE in Step 2"
      wrap(wb,tygwt, "YGGxE_wgted",26,1,F,F)
      
      save2excel(wb,ygGxE.wt.top, "YGGxE_wgted",27,1,T)
      rm(ygGxE.wt.top)
    }
    
    
    ## Var|GxE
    
    if(RunVarGxE==1){
      createSheet(wb, name = "rVarGxE_subset")
      
      ########## 9th panel
      
      saveimage(wb,"graph","rVarGxE_subset","levene_A.png",1,1)
      
      # Var step1 Manhattan
      
      saveimage(wb,"graph","rVarGxE_subset","levene_B.png",1,11)
      
      # Var step2 QQ
      
      if (var.s2.size>0) {
        saveimage(wb,"graph","rVarGxE_subset","levene_C.png",26,1)
        
        # Var step2 Manhattan
        saveimage(wb,"graph","rVarGxE_subset","levene_D.png",26,11)
        
        
        setwidth( wb,"rVarGxE_subset",6,2600)
        setwidth( wb,"rVarGxE_subset",8,2600)
        setwidth( wb,"rVarGxE_subset",9,2600)
        setwidth( wb,"rVarGxE_subset",1,3000)
        
        tvar<-"rVarGxE: 2-step with screening based on residual variance heterogeneity; subset testing of GxE in Step 2"
        wrap(wb,tvar, "rVarGxE_subset",51,1,F,F)
        
        sigvar<-paste("Number of SNPs tested in Step 2= ",var.s2.size,"; Step2 significance threshold = ",sig_var_st2,sep="")
        wrap(wb,sigvar, "rVarGxE_subset",52,1,F,F)
        
        wrap(wb,RankDescrip, "rVarGxE_subset",53,1,F,F)
        
        rvar<-"Ranks"
        save2excel(wb,rvar, "rVarGxE_subset",54,12,F)
        var.colend<- alphabet[3+num.method]
        var.letter<- paste("L54:",var.colend,"54",sep="")
        var.col<- 11+num.method
        savetitle( wb,"rVarGxE_subset",var.letter,54,c(12:var.col))
        
        save2excel(wb,var.s2, "rVarGxE_subset",55,1,T)
        rm(var.s2)
      }
      
      # VarGxE wgted
      
      ########## 10th panel
      createSheet(wb, name = "rVarGxE_wgted")
      saveimage(wb,"graph","rVarGxE_wgted","levene_E.png",1,1)
      
      ## VarGxE wt
      setwidth( wb,"rVarGxE_wgted",7,2600)
      setwidth( wb,"rVarGxE_wgted",9,2600)
      setwidth( wb,"rVarGxE_wgted",10,2600)
      setwidth( wb,"rVarGxE_wgted",12,2600)
      setwidth( wb,"rVarGxE_wgted",3,3000)
      
      tvarwt<-"rVarGxE: 2-step with screening based on residual variance heterogeneity; weighted testing of GxE in Step 2"
      wrap(wb,tvarwt, "rVarGxE_wgted",26,1,F,F)
      save2excel(wb,var.wt.top, "rVarGxE_wgted",27,1,T)
      rm(var.wt.top)
    }
    
    saveWorkbook(wb)
    rm(wb)
  } else {
    
    ## MarG
    if(RunYG==1){
      createSheet(wb, name = "YG")
      ########## 1st panel
      
      # Marginal G QQ
      
      saveimage(wb,"graph","YG","YG_A.png",1,1)
      
      # Marginal G Manhattan
      
      saveimage(wb,"graph","YG","YG_B.png",1,11)
      
      setwidth( wb,"YG",1,3000)
      #setwidth( wb,"YG",6,2600)
      
      tMarG<-"YG:  Standard GWAS analysis of marginal G effect for all M SNPs"
      wrap(wb,tMarG, "YG",26,1,F,F)
      
      sigMarG<-paste("Number of SNPs tested = ",M.yg,"; significance threshold = ",sig_yg,sep="")
      wrap(wb,sigMarG, "YG",27,1,F,F)
      
      save2excel(wb,yg.top, "YG",29,1,T)
      
      rm(yg.top)
    }
    
    ## GxE
    if(RunGxE==1){
      createSheet(wb, name = "GxE")
      
      ########## 2nd panel
      
      # G x E QQ
      
      saveimage(wb,"graph","GxE","GxE_A.png",1,1)
      
      # GxE Manhattan
      
      saveimage(wb,"graph","GxE","GxE_B.png",1,11)
      
      setwidth( wb,"GxE",1,3000)
      setwidth( wb,"GxE",6,2600)
      
      tgxe<-"GxE:  Standard analysis of GxE interaction for all M SNPs"
      wrap(wb,tgxe, "GxE",26,1,F,F)
      
      siggxe<-paste("Number of SNPs tested = ",M.gxe,"; significance threshold = ",sig_gxe,sep="")
      wrap(wb,siggxe, "GxE",27,1,F,F)
      
      save2excel(wb,gxe.top, "GxE",29,1,T)
      rm(gxe.top)
    }
    
    
    ##2df
    
    if(Run2df==1){
      createSheet(wb, name = "df2")
      
      ####### 3rd panel
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
    
    
    
    
    
    ## Joint|GxE
    
    if(RunJointGxE==1){
      createSheet(wb, name = "JointGxE_subset")
      
      ########## 5th panel
      
      saveimage(wb,"graph","JointGxE_subset","Joint_A.png",1,1)
      
      # Joint step1 Manhattan
      
      saveimage(wb,"graph","JointGxE_subset","Joint_B.png",1,11)
      
      # Joint step2 QQ
      
      if (joint.s2.size>0) {
        saveimage(wb,"graph","JointGxE_subset","Joint_C.png",26,1)
        
        # Joint step2 Manhattan
        saveimage(wb,"graph","JointGxE_subset","Joint_D.png",26,11)
        
        
        setwidth( wb,"JointGxE_subset",6,2600)
        setwidth( wb,"JointGxE_subset",8,2600)
        setwidth( wb,"JointGxE_subset",9,2600)
        setwidth( wb,"JointGxE_subset",1,3000)
        
        tjoint<-"Joint|GxE: 2-step with screening based on joint test of Y-G association and residual variance heterogeneity; subset testing of GxE in Step 2"
        wrap(wb,tjoint, "JointGxE_subset",51,1,F,F)
        
        sigjoint<-paste("Number of SNPs tested in Step 2= ",joint.s2.size,"; Step2 significance threshold = ",sig_joint_st2,sep="")
        wrap(wb,sigjoint, "JointGxE_subset",52,1,F,F)
        
        save2excel(wb,joint.s2, "JointGxE_subset",54,1,T)
        rm(joint.s2)
      }
      
      # JointGxE wgted
      
      ########## 6th panel
      createSheet(wb, name = "JointGxE_wgted")
      saveimage(wb,"graph","JointGxE_wgted","Joint_E.png",1,1)
      
      ##dgeg wt
      setwidth( wb,"JointGxE_wgted",7,2600)
      setwidth( wb,"JointGxE_wgted",9,2600)
      setwidth( wb,"JointGxE_wgted",10,2600)
      setwidth( wb,"JointGxE_wgted",12,2600)
      setwidth( wb,"JointGxE_wgted",3,3000)
      
      tjointwt<-"Joint|GxE: 2-step with screening based on joint test of Y-G association and residual variance heterogeneity; weighted testing of GxE in Step 2"
      wrap(wb,tjointwt, "JointGxE_wgted",26,1,F,F)
      save2excel(wb,joint.wt.top, "JointGxE_wgted",27,1,T)
      rm(joint.wt.top)
    }
    
    
    ## YG|GxE
    
    if(RunYGGxE==1){
      createSheet(wb, name = "YGGxE_subset")
      
      ####### 7th panel
      # Y-G step1 QQ
      saveimage(wb,"graph","YGGxE_subset","YGGxE_A.png",1,1)
      
      # Y-G step1 Manhattan
      saveimage(wb,"graph","YGGxE_subset","YGGxE_B.png",1,11)
      
      ##### Y-G step2 QQ
      if (ygGxE.s2.size>0) {
        saveimage(wb,"graph","YGGxE_subset","YGGxE_C.png",26,1)
        
        ###### Y-G step2 Manhattan
        
        saveimage(wb,"graph","YGGxE_subset","YGGxE_D.png",26,11)
        
        
        setwidth( wb,"YGGxE_subset",6,2600)
        setwidth( wb,"YGGxE_subset",8,2600)
        setwidth( wb,"YGGxE_subset",9,2600)
        setwidth( wb,"YGGxE_subset",1,3000)
        
        tyg<-"YG|GxE: 2-step with screening based on D-G association; subset testing of GxE in Step 2"
        wrap(wb,tyg, "YGGxE_subset",51,1,F,F)
        
        sigyg<-paste("Number of SNPs tested in Step 2= ",ygGxE.s2.size,"; Step2 significance threshold = ",sig_ygGxE_st2,sep="")
        wrap(wb,sigyg, "YGGxE_subset",52,1,F,F)
        
        save2excel(wb,ygGxE.s2, "YGGxE_subset",54,1,T)
        rm(ygGxE.s2)
      }
      
      
      ## YG|GxE wt
      
      ####### 8th panel
      
      ##### YG wgted 
      createSheet(wb, name = "YGGxE_wgted")
      saveimage(wb,"graph","YGGxE_wgted","YGGxE_E.png",1,1)
      
      setwidth( wb,"YGGxE_wgted",12,2800)
      setwidth( wb,"YGGxE_wgted",7,2600)
      setwidth( wb,"YGGxE_wgted",9,2600)
      setwidth( wb,"YGGxE_wgted",10,2600)
      setwidth( wb,"YGGxE_wgted",3,3000)
      
      tygwt<-"YG|GxE: 2-step with screening based on Y-G association; weighted testing of GxE in Step 2"
      wrap(wb,tygwt, "YGGxE_wgted",26,1,F,F)
      
      save2excel(wb,ygGxE.wt.top, "YGGxE_wgted",27,1,T)
      rm(ygGxE.wt.top)
    }
    
    
    ## Var|GxE
    
    if(RunVarGxE==1){
      createSheet(wb, name = "rVarGxE_subset")
      
      ########## 9th panel
      
      saveimage(wb,"graph","rVarGxE_subset","levene_A.png",1,1)
      
      # Var step1 Manhattan
      
      saveimage(wb,"graph","rVarGxE_subset","levene_B.png",1,11)
      
      # Var step2 QQ
      
      if (var.s2.size>0) {
        saveimage(wb,"graph","rVarGxE_subset","levene_C.png",26,1)
        
        # Var step2 Manhattan
        saveimage(wb,"graph","rVarGxE_subset","levene_D.png",26,11)
        
        
        setwidth( wb,"rVarGxE_subset",6,2600)
        setwidth( wb,"rVarGxE_subset",8,2600)
        setwidth( wb,"rVarGxE_subset",9,2600)
        setwidth( wb,"rVarGxE_subset",1,3000)
        
        tvar<-"rVar|GxE: 2-step with screening based on residual variance heterogeneity; subset testing of GxE in Step 2"
        wrap(wb,tvar, "rVarGxE_subset",51,1,F,F)
        
        sigvar<-paste("Number of SNPs tested in Step 2= ",var.s2.size,"; Step2 significance threshold = ",sig_var_st2,sep="")
        wrap(wb,sigvar, "rVarGxE_subset",52,1,F,F)
        
        save2excel(wb,var.s2, "rVarGxE_subset",54,1,T)
        rm(var.s2)
      }
      
      # VarGxE wgted
      
      ########## 10th panel
      createSheet(wb, name = "rVarGxE_wgted")
      saveimage(wb,"graph","rVarGxE_wgted","levene_E.png",1,1)
      
      ## VarGxE wt
      setwidth( wb,"rVarGxE_wgted",7,2600)
      setwidth( wb,"rVarGxE_wgted",9,2600)
      setwidth( wb,"rVarGxE_wgted",10,2600)
      setwidth( wb,"rVarGxE_wgted",12,2600)
      setwidth( wb,"rVarGxE_wgted",3,3000)
      
      tvarwt<-"rVar|GxE: 2-step with screening based on residual variance heterogeneity; weighted testing of GxE in Step 2"
      wrap(wb,tvarwt, "rVarGxE_wgted",26,1,F,F)
      save2excel(wb,var.wt.top, "rVarGxE_wgted",27,1,T)
      rm(var.wt.top)
    }
    
    saveWorkbook(wb)
    rm(wb)
    
  }
  
  ###################### Delete pics ##########################
  
  
  id <- grep("GxE_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("GxE_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  id <- grep("YG_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("YG_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("df2_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("df2_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("mdf_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("mdf_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("YGGxE_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("YGGxE_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  id <- grep("YGGxE_C.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("YGGxE_D.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  id <- grep("YGGxE_E.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("Joint_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("Joint_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("Joint_C.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("Joint_D.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  id <- grep("Joint_E.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("levene_A.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("levene_B.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("levene_C.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  id <- grep("levene_D.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
  id <- grep("levene_E.png", dir(homedir))
  todelete <- dir(homedir, full.names = TRUE)[id]
  unlink(todelete)
  
  
}


GWIS()




