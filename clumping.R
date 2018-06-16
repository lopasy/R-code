# Creating clumps
ukb <- read_delim("~/Documents/UKB/gxescan/ukb_full_0.05.snpinfo", 
                  +     "\t", escape_double = FALSE, trim_ws = TRUE)
YG <- read_delim("~/Documents/UKB/gxescan/ukb_full_0.05_QT_YG.gxeout", 
                 +     "\t", escape_double = FALSE, trim_ws = TRUE)
VH <- read_delim("~/Documents/UKB/gxescan/ukb_full_0.05_QT_VH.gxeout", 
                 +     "\t", escape_double = FALSE, trim_ws = TRUE)
YGVH <- read_delim("~/Documents/UKB/gxescan/ukb_full_0.05_QT_YGVH.gxeout", 
                   +     "\t", escape_double = FALSE, trim_ws = TRUE)

zvh=merge(ukb, VH,by="SNPID")
zygvh=merge(ukb, YGVH,by="SNPID")
zyg=merge(ukb, YG,by="SNPID")


write.table(zvh, file = "clump_ukb_full_vh_5e8.txt", row.names = F, quote = F)
write.table(zygvh, file = "clump_ukb_full_ygvh_5e8.txt", row.names = F, quote = F)
write.table(zyg, file = "clump_ukb_full_yg_5e8.txt", row.names = F, quote = F)