####################
##### Nearwork #####
####################

# Model1_threeway_interaction_Age_Read_SNP
#Significant P_TIME:SNP:NEARWORK hits - 5(rs7829127(ZMAT4), rs10880855(ARID2), rs10500355(RBFOX1), 
#                                                    rs6420484(TSPAN10), rs7253703(AP1M1_KLF2))


{var_names<-names(all[,c(12:141, 143:158, 160:161)])
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
results_matrix=matrix(nrow=(num_snps),ncol=25)
colnames(results_matrix)=c("SNP","Chr","Pos","Gene","RiskAllele","OtherAllele","SNPcheck",
                           "beta_SNP","se_SNP","P_SNP",
                           "beta_Nearwork","se_Nearwork","P_Nearwork",
                           "beta_Nearwork:SNP","se_Nearwork:SNP","P_Nearwork:SNP",
                           "beta_Time:SNP","se_Time:SNP","P_Time:SNP",
                           "beta_Time:Nearwork","se_Time:Nearwork","P_Time:Nearwork",
                           "beta_Time:Nearwork:SNP","se_Time:Nearwork:SNP","P_Time:Nearwork:SNP" 
)
results_matrix<-as.data.frame(results_matrix)
results_matrix[,1:6]<-cream2017_hits[c(2:131, 133:148, 150:151),c(1,3:7)]
# 
for (snp in 1:num_snps){
  #
  # write formula
  # 
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],  
                            " + ReadingBinary",
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + ReadingBinary:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):ReadingBinary",
                            " + I(age_at_visit - 7.5):ReadingBinary:",vars[snp_pos[snp]],
                            sep=""))
  # 
  # run model, store summary to "a"
  # 
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  # 
  # save fixed effects results
  # 
  results_matrix[snp,7] <-noquote(vars[snp_pos[snp]]) #SNP name check
  results_matrix[snp,8] <-a$tTable[2,1]               #beta_SNP
  results_matrix[snp,9] <-a$tTable[2,2]               #se_SNP
  results_matrix[snp,10]<-a$tTable[2,5]               #P_SNP
  results_matrix[snp,11]<-a$tTable[3,1]               #beta_Nearwork
  results_matrix[snp,12]<-a$tTable[3,2]               #se_Nearwork
  results_matrix[snp,13]<-a$tTable[3,5]               #P_Nearwork
  results_matrix[snp,14]<-a$tTable[8,1]               #beta_Nearwork:SNP
  results_matrix[snp,15]<-a$tTable[8,2]               #se_Nearwork:SNP
  results_matrix[snp,16]<-a$tTable[8,5]               #P_Nearwork:SNP
  results_matrix[snp,17]<-a$tTable[9,1]               #beta_Time:SNP
  results_matrix[snp,18]<-a$tTable[9,2]               #se_Time:SNP
  results_matrix[snp,19]<-a$tTable[9,5]               #P_Time:SNP
  results_matrix[snp,20]<-a$tTable[10,1]              #beta_Time:Nearwork
  results_matrix[snp,21]<-a$tTable[10,2]              #se_Time:Nearwork
  results_matrix[snp,22]<-a$tTable[10,5]              #P_Time:Nearwork
  results_matrix[snp,23]<-a$tTable[11,1]              #beta_Time:Nearwork:SNP
  results_matrix[snp,24]<-a$tTable[11,2]              #se_Time:Nearwork:SNP
  results_matrix[snp,25]<-a$tTable[11,5]              #P_Time:Nearwork:SNP
  # 
}} # next SNP
write.csv(results_matrix, file="E:/CSV/MyModel_tripleB.csv", row.names=FALSE)
#
# Plot SNPs with significant TIME:SNP:NEARWORK effects
graph_data=matrix(nrow=54,ncol=5)
pdata<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), ReadingBinary=seq(0,1,by=1), mse_at_visit="1")
snplist  <- c(59,98,122,135,138)
var_names<-names(all[,11:161])
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
curr_col <- 1
for (n in 1:5){
  snp<-snplist[n]
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + ReadingBinary:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):ReadingBinary",
                            " + I(age_at_visit - 7.5):ReadingBinary:",vars[snp_pos[snp]],
                            sep=""))
  
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  
  colnames(pdata)<-c("age_at_visit", vars[snp_pos[snp]], "ReadingBinary", "mse_at_visit" )
  pdata[,4]<-predict(a,pdata,level=0)
  graph_data[1:54,curr_col]<-pdata$mse_at_visit[1:54]
  curr_col<-curr_col+1
}

{pdata1<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), ReadingBinary=seq(0,1,by=1), mse_at_visit="1")
colnames(pdata1)<-c("age_at_visit", "SNP", "Nearwork", "mse_at_visit" )
pdata1$Nearwork <- as.factor(pdata1$Nearwork)
levels(pdata1$Nearwork) <- c("Low","High")
pdata1$SNP <- as.factor(pdata1$SNP)
levels(pdata1$SNP) <- c("0","1", "2")}

{pdata1$mse_at_visit<-graph_data[1:54,1]
pdata2<-pdata1;pdata2$mse_at_visit<-graph_data[1:54,2]
pdata3<-pdata1;pdata3$mse_at_visit<-graph_data[1:54,3]
pdata4<-pdata1;pdata4$mse_at_visit<-graph_data[1:54,4]
pdata5<-pdata1;pdata5$mse_at_visit<-graph_data[1:54,5]}

# For rs7253703(AP1M1_KLF2) 0 and 1 SNP copy number
{graph_data=matrix(nrow=36,ncol=1)
pdata<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,1,by=1), ReadingBinary=seq(0,1,by=1), mse_at_visit="1")
snplist  <- 138
formula<-as.formula(paste("mse_at_visit ~ ", 
                          vars[snp_pos[snplist]],
                          " + poly(I(age_at_visit - 7.5),4)",
                          " + ReadingBinary:",vars[snp_pos[snplist]],
                          " + I(age_at_visit - 7.5):",vars[snp_pos[snplist]],
                          " + I(age_at_visit - 7.5):ReadingBinary",
                          " + I(age_at_visit - 7.5):ReadingBinary:",vars[snp_pos[snplist]],
                          sep=""))
a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                             correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                             na.action = na.omit, 
                                             method="ML",
                                             data=all)))
colnames(pdata)<-c("age_at_visit", vars[snp_pos[snplist]], "ReadingBinary", "mse_at_visit" )
pdata[,4]<-predict(a,pdata,level=0)
graph_data[1:36,]<-pdata$mse_at_visit[1:36]

pdata6<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,1,by=1), ReadingBinary=seq(0,1,by=1), mse_at_visit="1")
  colnames(pdata6)<-c("age_at_visit", "SNP", "Nearwork", "mse_at_visit" )
  pdata6$Nearwork <- as.factor(pdata6$Nearwork)
  levels(pdata6$Nearwork) <- c("Low","High")
  pdata6$SNP <- as.factor(pdata6$SNP)
  levels(pdata6$SNP) <- c("0","1")
pdata6$mse_at_visit<-graph_data[1:36,1]}

# For rs7253703(AP1M1_KLF2) 0 SNP copy number
  pdata7<-expand.grid(age_at_visit=seq(7,15,by=1), SNP="0", ReadingBinary=seq(0,1,by=1), mse_at_visit="1")
  colnames(pdata7)<-c("age_at_visit", "SNP", "Nearwork", "mse_at_visit" )
  pdata7$Nearwork <- as.factor(pdata7$Nearwork)
  levels(pdata7$Nearwork) <- c("Low","High")
  pdata7$SNP <- as.factor(pdata7$SNP)
  levels(pdata7$SNP) <- c("0")
pdata7$mse_at_visit<-graph_data[c(1:9,19:27),1]


theme_fred <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      plot.margin = unit(c(0.25,0.25,0.5,0.75), "lines"),
      plot.background = element_rect(fill="white",colour="white"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1.3), vjust=-0.25),
      axis.title.y = element_text(colour = "black", size=rel(1.3), angle=90, vjust=1.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.key = element_rect(fill = "white",colour="white"),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size=rel(0.7), face="bold"),
      legend.text  = element_text(size=rel(0.7), face="plain"),
      legend.position=c(0.01,0), legend.justification=c(0.01,0),
      legend.box.just = "left"
    )   
}
theme_set(theme_fred())

{way1<-ggplot(pdata1, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs7829127
(ZMAT4)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
    geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    theme(panel.background = element_rect(fill="white",colour="black") ) + 
    annotate("text", x = 14.5, y=0.25,label = "a",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7"))
  
  way2<-ggplot(pdata2, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs10880855
(ARID2)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
    geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    annotate("text", x = 14.5, y=0.25,label = "b",cex=10, hjust = 0.5) + 
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way3<-ggplot(pdata3, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs10500355
(RBFOX1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
    geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) + 
    annotate("text", x = 14.5, y=0.25,label = "c",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way4<-ggplot(pdata4, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs6420484
(TSPAN10)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
    geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    annotate("text", x = 14.5, y=0.25,label = "d",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    guides(fill = guide_legend(order=1), linetype = guide_legend(order=1), shape = guide_legend(order=1)) +
    theme(panel.background = element_rect(fill="white",colour="black") )}

png("E:/CSV/Age_Nearwork_SNP.png", units=c("mm"), res=260, width=200, height=246)
grid.arrange(way1,way2,way3,way4,
             left=textGrob("Refractive error (D)", rot=90, vjust=1, gp=gpar(cex=1.7)),
             bottom=textGrob ("Age (years)", vjust=-0.5, gp=gpar(cex=1.7)),
             nrow=2, ncol=2, clip=FALSE, top = textGrob("Age:Nearwork:SNP", gp=gpar(cex=1.7)))
dev.off()
  
{way5<-ggplot(pdata5, aes(age_at_visit, mse_at_visit)) + 
  labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs7253703
(AP1M1_KLF2)") + 
  scale_linetype_manual(values=c(2,1)) +
  geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
  geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
  scale_shape_manual(values = c(21,21))  +
  scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
  scale_y_continuous(limits=c(-7.5, 0.8), breaks=seq(-7.5,0.8,0.5)) +
  scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  annotate("text", x = 14.5, y=0.8,label = "a",cex=10, hjust = 0.5) +
  theme(panel.background = element_rect(fill="white",colour="black") )
way6<-ggplot(pdata6, aes(age_at_visit, mse_at_visit)) + 
  labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs7253703
(AP1M1_KLF2)") + 
  scale_linetype_manual(values=c(2,1)) +
  geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
  geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
  scale_shape_manual(values = c(21,21))  +
  scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
  scale_y_continuous(limits=c(-4, 0.8), breaks=seq(-4,0.8,0.5)) +
  scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  annotate("text", x = 14.5, y=0.8,label = "b",cex=10, hjust = 0.5) +
  theme(panel.background = element_rect(fill="white",colour="black") )
way7<-ggplot(pdata7, aes(age_at_visit, mse_at_visit)) + 
  labs (x="", y="", shape="Nearwork", linetype="Nearwork", colour="rs7253703
(AP1M1_KLF2)") + 
  scale_linetype_manual(values=c(2,1)) +
  geom_line  (size=1.1, aes(colour=SNP, linetype=Nearwork)) + 
  geom_point (size=3, colour="black",aes(fill=Nearwork, shape=Nearwork)) +
  scale_shape_manual(values = c(21,21))  +
  scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
  scale_y_continuous(limits=c(-0.75, 0.25), breaks=seq(-0.75,0.25,0.25)) +
  scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
  annotate("text", x = 14.5, y=0.25,label = "c",cex=10, hjust = 0.5) +
  theme(panel.background = element_rect(fill="white",colour="black") )}

png("E:/CSV/Age_Nearwork_SNP_extra.png", units=c("mm"), res=150, width=200, height=140)
grid.arrange(way5,way6,way7,
             left=textGrob("Refractive error (D)", rot=90, vjust=1, gp=gpar(cex=1.7)),
             bottom=textGrob ("Age (years)", vjust=-0.5, gp=gpar(cex=1.7)),
             nrow=1, ncol=3, clip=FALSE, top = textGrob("Age:Nearwork:SNP", gp=gpar(cex=1.7)))
dev.off()

#######################
##### Outdoor1151 #####
#######################

# Model1_threeway_interaction_Age_OUTDOORS1151_SNP
#Significant P_TIME:SNP:Outdoors hits - 5(rs1556867(5S_rRNA_PBX1), rs17428076(HAT1_METAP1D), rs1064583(COL10A1), 
#                                                    rs28891973(NCOA2_TRAM1), rs12451582(NOG_C17orf67))


{var_names<-names(all[,c(160:161)])#12:141, 143:147, 149:158, 160:161
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
results_matrix=matrix(nrow=(num_snps),ncol=25)
colnames(results_matrix)=c("SNP","Chr","Pos","Gene","RiskAllele","OtherAllele","SNPcheck",
                           "beta_SNP","se_SNP","P_SNP",
                           "beta_Outdoors","se_Outdoors","P_Outdoors",
                           "beta_Outdoors:SNP","se_Outdoors:SNP","P_Outdoors:SNP",
                           "beta_Time:SNP","se_Time:SNP","P_Time:SNP",
                           "beta_Time:Outdoors","se_Time:Outdoors","P_Time:Outdoors",
                           "beta_Time:Outdoors:SNP","se_Time:Outdoors:SNP","P_Time:Outdoors:SNP" 
)
results_matrix<-as.data.frame(results_matrix)
results_matrix[,1:6]<-cream2017_hits[c(150:151),c(1,3:7)]#2:131, 133:137, 139:148, 150:151
# 
for (snp in 1:num_snps){
  #
  # write formula
  # 
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],  
                            " + Outdoor1151",
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + Outdoor1151:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):Outdoor1151",
                            " + I(age_at_visit - 7.5):Outdoor1151:",vars[snp_pos[snp]],
                            sep=""))
  # 
  # run model, store summary to "a"
  # 
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  # 
  # save fixed effects results
  # 
  results_matrix[snp,7] <-noquote(vars[snp_pos[snp]]) #SNP name check
  results_matrix[snp,8] <-a$tTable[2,1]               #beta_SNP
  results_matrix[snp,9] <-a$tTable[2,2]               #se_SNP
  results_matrix[snp,10]<-a$tTable[2,5]               #P_SNP
  results_matrix[snp,11]<-a$tTable[3,1]               #beta_Outdoors
  results_matrix[snp,12]<-a$tTable[3,2]               #se_Outdoors
  results_matrix[snp,13]<-a$tTable[3,5]               #P_Outdoors
  results_matrix[snp,14]<-a$tTable[8,1]               #beta_Outdoors:SNP
  results_matrix[snp,15]<-a$tTable[8,2]               #se_Outdoors:SNP
  results_matrix[snp,16]<-a$tTable[8,5]               #P_Outdoors:SNP
  results_matrix[snp,17]<-a$tTable[9,1]               #beta_Time:SNP
  results_matrix[snp,18]<-a$tTable[9,2]               #se_Time:SNP
  results_matrix[snp,19]<-a$tTable[9,5]               #P_Time:SNP
  results_matrix[snp,20]<-a$tTable[10,1]              #beta_Time:Outdoors
  results_matrix[snp,21]<-a$tTable[10,2]              #se_Time:Outdoors
  results_matrix[snp,22]<-a$tTable[10,5]              #P_Time:Outdoors
  results_matrix[snp,23]<-a$tTable[11,1]              #beta_Time:Outdoors:SNP
  results_matrix[snp,24]<-a$tTable[11,2]              #se_Time:Outdoors:SNP
  results_matrix[snp,25]<-a$tTable[11,5]              #P_Time:Outdoors:SNP
  # 
}} # next SNP
write.csv(results_matrix, file="E:/CSV/MyModel_tripleBB.csv", row.names=FALSE)
#
# Plot SNPs with significant TIME:SNP:OUTDOORS1151 effects
graph_data=matrix(nrow=54,ncol=5)
pdata<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), Outdoor1151=seq(0,1,by=1), mse_at_visit="1")
snplist  <- c(6,19,47,64,131)
var_names<-names(all[,11:161])
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
curr_col <- 1
for (n in 1:5){
  snp<-snplist[n]
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + Outdoor1151:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):Outdoor1151",
                            " + I(age_at_visit - 7.5):Outdoor1151:",vars[snp_pos[snp]],
                            sep=""))
  
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  
  colnames(pdata)<-c("age_at_visit", vars[snp_pos[snp]], "Outdoor1151", "mse_at_visit" )
  pdata[,4]<-predict(a,pdata,level=0)
  graph_data[1:54,curr_col]<-pdata$mse_at_visit[1:54]
  curr_col<-curr_col+1
}

{pdata1<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), Outdoor1151=seq(0,1,by=1), mse_at_visit="1")
  colnames(pdata1)<-c("age_at_visit", "SNP", "Outdoors", "mse_at_visit" )
  pdata1$Outdoors <- as.factor(pdata1$Outdoors)
  levels(pdata1$Outdoors) <- c("Low","High")
  pdata1$SNP <- as.factor(pdata1$SNP)
  levels(pdata1$SNP) <- c("0","1", "2")}

{pdata1$mse_at_visit<-graph_data[1:54,1]
  pdata2<-pdata1;pdata2$mse_at_visit<-graph_data[1:54,2]
  pdata3<-pdata1;pdata3$mse_at_visit<-graph_data[1:54,3]
  pdata4<-pdata1;pdata4$mse_at_visit<-graph_data[1:54,4]
  pdata5<-pdata1;pdata5$mse_at_visit<-graph_data[1:54,5]}

theme_fred <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      plot.margin = unit(c(0.25,0.25,0.5,0.75), "lines"),
      plot.background = element_rect(fill="white",colour="white"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1.3), vjust=-0.25),
      axis.title.y = element_text(colour = "black", size=rel(1.3), angle=90, vjust=1.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.key = element_rect(fill = "white",colour="white"),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size=rel(0.7), face="bold"),
      legend.text  = element_text(size=rel(0.7), face="plain"),
      legend.position=c(0.01,0), legend.justification=c(0.01,0),
      legend.box.just = "left"
    )   
}
theme_set(theme_fred())

{way1<-ggplot(pdata1, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs1556867
          (5S_rRNA_PBX1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    theme(panel.background = element_rect(fill="white",colour="black") ) + 
    annotate("text", x = 15, y=0.25,label = "a",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7"))
  
  way2<-ggplot(pdata2, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs17428076
          (HAT1_METAP1D)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    annotate("text", x = 15, y=0.25,label = "b",cex=10, hjust = 0.5) + 
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way3<-ggplot(pdata3, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs1064583
(COL10A1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) + 
    annotate("text", x = 15, y=0.25,label = "c",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way4<-ggplot(pdata4, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs1064583
(NCOA2_TRAM1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    annotate("text", x = 15, y=0.25,label = "d",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    guides(fill = guide_legend(order=1), linetype = guide_legend(order=1), shape = guide_legend(order=1)) +
    theme(panel.background = element_rect(fill="white",colour="black") )
  
  way5<-ggplot(pdata5, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs12451582
(NOG_C17orf67)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    annotate("text", x = 15, y=0.25,label = "e",cex=10, hjust = 0.5) +
    theme(panel.background = element_rect(fill="white",colour="black") )}

png("E:/CSV/Age_Outdoors_SNP.png", units=c("mm"), res=300, width=300, height=246)
grid.arrange(way1,way2,way3,way4,way5,
             left=textGrob("Refractive error (D)", rot=90, vjust=1, gp=gpar(cex=1.7)),
             bottom=textGrob ("Age (years)", vjust=-0.5, gp=gpar(cex=1.7)),
             nrow=2, ncol=3, clip=FALSE, top = textGrob("Age:Outdoors:SNP", gp=gpar(cex=1.7)))
dev.off()



#######################
##### Outdoor1171 #####
#######################

# Model1_threeway_interaction_Age_OUTDOORS1151_SNP
#Significant P_TIME:SNP:Outdoors hits - 8 rs11210537(HIVEP3), rs2745953(CD34), rs478304(RNASEH2C_AP5B1), 
#                                         rs1954761(GRIA4), rs811156(BMP4_CDKN3), rs11160044(NDUFB1),
#                                         rs62070229(MYO1D_TMEM98), rs9606967(AK123891_SYN3)   


{var_names<-names(all[,c(39:141, 143:147, 149:158, 160:161)])#12:37, 39:141, 143:147, 149:158, 160:161
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
results_matrix=matrix(nrow=(num_snps),ncol=25)
colnames(results_matrix)=c("SNP","Chr","Pos","Gene","RiskAllele","OtherAllele","SNPcheck",
                           "beta_SNP","se_SNP","P_SNP",
                           "beta_Outdoors","se_Outdoors","P_Outdoors",
                           "beta_Outdoors:SNP","se_Outdoors:SNP","P_Outdoors:SNP",
                           "beta_Time:SNP","se_Time:SNP","P_Time:SNP",
                           "beta_Time:Outdoors","se_Time:Outdoors","P_Time:Outdoors",
                           "beta_Time:Outdoors:SNP","se_Time:Outdoors:SNP","P_Time:Outdoors:SNP" 
)
results_matrix<-as.data.frame(results_matrix)
results_matrix[,1:6]<-cream2017_hits[c(29:131, 133:137, 139:148, 150:151),c(1,3:7)]#2:27, 29:131, 133:137, 139:148, 150:151
# 
for (snp in 1:num_snps){
  #
  # write formula
  # 
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],  
                            " + Outdoor1171",
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + Outdoor1171:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):Outdoor1171",
                            " + I(age_at_visit - 7.5):Outdoor1171:",vars[snp_pos[snp]],
                            sep=""))
  # 
  # run model, store summary to "a"
  # 
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  # 
  # save fixed effects results
  # 
  results_matrix[snp,7] <-noquote(vars[snp_pos[snp]]) #SNP name check
  results_matrix[snp,8] <-a$tTable[2,1]               #beta_SNP
  results_matrix[snp,9] <-a$tTable[2,2]               #se_SNP
  results_matrix[snp,10]<-a$tTable[2,5]               #P_SNP
  results_matrix[snp,11]<-a$tTable[3,1]               #beta_Outdoors
  results_matrix[snp,12]<-a$tTable[3,2]               #se_Outdoors
  results_matrix[snp,13]<-a$tTable[3,5]               #P_Outdoors
  results_matrix[snp,14]<-a$tTable[8,1]               #beta_Outdoors:SNP
  results_matrix[snp,15]<-a$tTable[8,2]               #se_Outdoors:SNP
  results_matrix[snp,16]<-a$tTable[8,5]               #P_Outdoors:SNP
  results_matrix[snp,17]<-a$tTable[9,1]               #beta_Time:SNP
  results_matrix[snp,18]<-a$tTable[9,2]               #se_Time:SNP
  results_matrix[snp,19]<-a$tTable[9,5]               #P_Time:SNP
  results_matrix[snp,20]<-a$tTable[10,1]              #beta_Time:Outdoors
  results_matrix[snp,21]<-a$tTable[10,2]              #se_Time:Outdoors
  results_matrix[snp,22]<-a$tTable[10,5]              #P_Time:Outdoors
  results_matrix[snp,23]<-a$tTable[11,1]              #beta_Time:Outdoors:SNP
  results_matrix[snp,24]<-a$tTable[11,2]              #se_Time:Outdoors:SNP
  results_matrix[snp,25]<-a$tTable[11,5]              #P_Time:Outdoors:SNP
  # 
}} # next SNP
write.csv(results_matrix, file="E:/CSV/MyModel_tripleBB.csv", row.names=FALSE)
#
# Plot SNPs with significant TIME:SNP:OUTDOORS1151 effects
graph_data=matrix(nrow=54,ncol=8)
pdata<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), Outdoor1171=seq(0,1,by=1), mse_at_visit="1")
snplist  <- c(2,9,87,89,110,113,127,150)
var_names<-names(all[,11:161])
vars<-as.matrix(var_names)
num_vars<-nrow(vars)
num_snps<-0
snp_pos<-matrix(nrow=num_vars,ncol=1)
for (n in 1:num_vars){
  if (substr(vars[n],1,2)=="rs"){
    num_snps<-num_snps+1
    snp_pos[num_snps]<-n
  } 
}
curr_col <- 1
for (n in 1:8){
  snp<-snplist[n]
  formula<-as.formula(paste("mse_at_visit ~ ", 
                            vars[snp_pos[snp]],
                            " + poly(I(age_at_visit - 7.5),4)",
                            " + Outdoor1171:",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):",vars[snp_pos[snp]],
                            " + I(age_at_visit - 7.5):Outdoor1171",
                            " + I(age_at_visit - 7.5):Outdoor1171:",vars[snp_pos[snp]],
                            sep=""))
  
  a<-summary(model<-do.call("lme", args = list(formula, random=~I(age_at_visit - 7.5) | alfred_ID1, 
                                               correlation = corCAR1(form = ~ visit| alfred_ID1), 
                                               na.action = na.omit, 
                                               method="ML",
                                               data=all)))
  
  colnames(pdata)<-c("age_at_visit", vars[snp_pos[snp]], "Outdoor1171", "mse_at_visit" )
  pdata[,4]<-predict(a,pdata,level=0)
  graph_data[1:54,curr_col]<-pdata$mse_at_visit[1:54]
  curr_col<-curr_col+1
}

{pdata1<-expand.grid(age_at_visit=seq(7,15,by=1), SNP=seq(0,2,by=1), Outdoor1171=seq(0,1,by=1), mse_at_visit="1")
  colnames(pdata1)<-c("age_at_visit", "SNP", "Outdoors", "mse_at_visit" )
  pdata1$Outdoors <- as.factor(pdata1$Outdoors)
  levels(pdata1$Outdoors) <- c("Low","High")
  pdata1$SNP <- as.factor(pdata1$SNP)
  levels(pdata1$SNP) <- c("0","1", "2")}

{pdata1$mse_at_visit<-graph_data[1:54,1]
  pdata2<-pdata1;pdata2$mse_at_visit<-graph_data[1:54,2]
  pdata3<-pdata1;pdata3$mse_at_visit<-graph_data[1:54,3]
  pdata4<-pdata1;pdata4$mse_at_visit<-graph_data[1:54,4]
  pdata5<-pdata1;pdata5$mse_at_visit<-graph_data[1:54,5]
  pdata6<-pdata1;pdata6$mse_at_visit<-graph_data[1:54,6]
  pdata7<-pdata1;pdata7$mse_at_visit<-graph_data[1:54,7]
  pdata8<-pdata1;pdata8$mse_at_visit<-graph_data[1:54,8]}

theme_fred <- function (base_size = 12, base_family = "") {
  theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      plot.margin = unit(c(0.25,0.25,0.5,0.75), "lines"),
      plot.background = element_rect(fill="white",colour="white"),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(colour = "black", size=rel(1.3), vjust=-0.25),
      axis.title.y = element_text(colour = "black", size=rel(1.3), angle=90, vjust=1.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.key = element_rect(fill = "white",colour="white"),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.2, "cm"),
      legend.title = element_text(size=rel(0.7), face="bold"),
      legend.text  = element_text(size=rel(0.7), face="plain"),
      legend.position=c(0.01,0), legend.justification=c(0.01,0),
      legend.box.just = "left"
    )   
}
theme_set(theme_fred())

{way1<-ggplot(pdata1, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs11210537
          (HIVEP3)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    theme(panel.background = element_rect(fill="white",colour="black") ) + 
    annotate("text", x = 15, y=0.25,label = "a",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7"))
  
  way2<-ggplot(pdata2, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs2745953
          (CD34)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.25,0.25)) +
    annotate("text", x = 15, y=0.25,label = "b",cex=10, hjust = 0.5) + 
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way3<-ggplot(pdata3, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs478304
          (RNASEH2C_AP5B1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) + 
    annotate("text", x = 15, y=0.25,label = "c",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    theme(panel.background = element_rect(fill="white",colour="black") ) 
  
  way4<-ggplot(pdata4, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs1954761
          (GRIA4)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    annotate("text", x = 15, y=0.25,label = "d",cex=10, hjust = 0.5) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    guides(fill = guide_legend(order=1), linetype = guide_legend(order=1), shape = guide_legend(order=1)) +
    theme(panel.background = element_rect(fill="white",colour="black") )
  
  way5<-ggplot(pdata5, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs811156
          (BMP4_CDKN3)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    annotate("text", x = 15, y=0.25,label = "e",cex=10, hjust = 0.5) +
    theme(panel.background = element_rect(fill="white",colour="black") )
  
  way6<-ggplot(pdata6, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs11160044
          (NDUFB1)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    annotate("text", x = 15, y=0.25,label = "f",cex=10, hjust = 0.5) +
    theme(panel.background = element_rect(fill="white",colour="black") )
  
  way7<-ggplot(pdata7, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs62070229
          (MYO1D_TMEM98)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    annotate("text", x = 15, y=0.25,label = "g",cex=10, hjust = 0.5) +
    theme(panel.background = element_rect(fill="white",colour="black") )
  
  way8<-ggplot(pdata8, aes(age_at_visit, mse_at_visit)) + 
    labs (x="", y="", shape="Outdoors", linetype="Outdoors", colour="rs9606967
          AK123891_SYN3)") + 
    scale_linetype_manual(values=c(2,1)) +
    geom_line  (size=1.1, aes(colour=SNP, linetype=Outdoors)) + 
    geom_point (size=3, colour="black",aes(fill=Outdoors, shape=Outdoors)) +
    scale_shape_manual(values = c(21,21))  +
    scale_x_continuous(limits=c(7.5, 15.5), breaks=seq(8, 15, 1)) +
    scale_y_continuous(limits=c(-0.8, 0.3), breaks=seq(-0.75,0.5,0.25)) +
    scale_colour_manual(values=c("#0000FF", "#00FF00", "#990000", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    scale_fill_manual(values=c("#000000", "#FFFFFF", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00","#CC79A7")) +
    annotate("text", x = 15, y=0.25,label = "h",cex=10, hjust = 0.5) +
    theme(panel.background = element_rect(fill="white",colour="black") )}

png("E:/CSV/Age_Outdoor1171_SNP.png", units=c("mm"), res=300, width=300, height=246)
grid.arrange(way1,way2,way3,way4,way5, way6, way7, way8,
             left=textGrob("Refractive error (D)", rot=90, vjust=1, gp=gpar(cex=1.7)),
             bottom=textGrob ("Age (years)", vjust=-0.5, gp=gpar(cex=1.7)),
             nrow=2, ncol=4, clip=FALSE, top = textGrob("Age:Outdoors:SNP", gp=gpar(cex=1.7)))
dev.off()
