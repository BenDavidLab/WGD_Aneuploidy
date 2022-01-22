ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/Final checks"

setwd(ROOT_DIR)

data <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_data.txt", header=T, sep='\t', as.is=T)

head(data$Genome_doublings)
TYPES<-unique(data$Type)

library("stringr")

#subsetting data based on WGD status
CA <- data
CA_WGD0<-CA[CA$Genome_doublings==0,]
CA_WGD1<-CA[CA$Genome_doublings==1,]
CA_WGD2<-CA[CA$Genome_doublings==2,]
CA_WGD<-CA[CA$Genome_doublings>0,]

#all tumor types with >=20 events in each group
g_20_types <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC")



WGD_fractions <- c()
WGD0_fractions <- c()

for (t in 1:length(g_20_types)){
  
  library("MASS")
  
  g_20_types = c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC")
  
  
  df <- data.frame(frac_top3 = c(), wgd_status = as.character(c()), type = as.character(c()))
  
  
  ty <- g_20_types[t]
  
  CA<-data[data$Type==ty,]
  CA_WGD0<-CA[CA$Genome_doublings==0,]
  CA_WGD1<-CA[CA$Genome_doublings==1,]
  CA_WGD2<-CA[CA$Genome_doublings==2,]
  CA_WGD<-CA[CA$Genome_doublings>0,]
  
  data_plus_minus <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_plus_minus_data.txt", header=T, sep='\t', as.is=T)
  
  
  CA_PM<-data_plus_minus[data_plus_minus$Type==ty,]
  CA_WGD0_plus_minus<-CA_PM[CA_PM$Genome_doublings==0,]
  CA_WGD1_plus_minus<-CA_PM[CA_PM$Genome_doublings==1,]
  CA_WGD2_plus_minus<-CA_PM[CA_PM$Genome_doublings==2,]
  CA_WGD_plus_minus<-CA_PM[CA_PM$Genome_doublings>0,]
  
  
  CA_WGD_mat <- CA_WGD[,14:(ncol(CA_WGD)-17)]
  CA_WGD0_mat <- CA_WGD0[,14:(ncol(CA_WGD0)-17)]
  
  CA_PM_WGD0<-c()
  CA_PM_WGD<-c()
  
  
  for (i in 14:91){
    WGD0_s <- sum(CA_WGD0_plus_minus[,i]==1, na.rm=T)
    WGD_s <- sum(CA_WGD_plus_minus[,i]==1, na.rm=T)
    CA_PM_WGD0[i]<- WGD0_s
    CA_PM_WGD[i]<- WGD_s
  }
  
  
  WGD_sum_desc <- sort(CA_PM_WGD, decreasing = T)
  WGD0_sum_desc <- sort(CA_PM_WGD0, decreasing = T)
  
  WGD_sum_top_3 <- WGD_sum_desc[1:3]
  WGD0_sum_top_3 <- WGD0_sum_desc[1:3]
  WGD_sum_total <- sum(WGD_sum_top_3, na.rm=T)
  WGD0_sum_total <- sum(WGD0_sum_top_3, na.rm=T)
  WGD_total <- sum(CA_PM_WGD, na.rm=T)
  WGD0_total <- sum(CA_PM_WGD0, na.rm=T)
  WGD_frac <- WGD_sum_total/WGD_total
  WGD0_frac <- WGD0_sum_total/WGD0_total
  WGD_fractions[t] <- WGD_frac
  WGD0_fractions[t] <- WGD0_frac
  
  
}

frac_top3_t <- c(WGD0_fractions, WGD_fractions)
wgd_status_t <- c(rep(as.character("non-WGD"), length(WGD0_fractions)),rep("WGD", length(WGD_fractions)))
type_t <- c(rep(as.character(g_20_types), 2))

new_data <- cbind(frac_top3_t, wgd_status_t, type_t)


df <- rbind(df, new_data)


plot.new()

#par(mar = c(4.4,2.5,1.5,1.5))
color_vec <- c("chartreuse3", "mediumpurple3", "violetred1", "turquoise3", "blue4", "slateblue2", "gray75", "firebrick4", "chocolate1", "gray49","springgreen4","bisque" , "lavenderblush2", "cyan4", "purple3", "gray32", "deepskyblue1", "blue", "gold", "gray14", "lightslateblue", "salmon")

#pdf(paste0(ROOT_DIR, "/", "pan_4.3_boxplot_11.pdf"), width = 12, height = 5.5)
pdf(paste0(ROOT_DIR, "/", "Pan-cancer top 3_ax_sz2.pdf"), width = 12, height = 5.5)
par(mar = c(4.4,5,1.5,12))
par(cex.lab=1.05)
par(cex.axis=1.5)

boxplot(WGD0_fractions, WGD_fractions, col=c("white"),names = c("WGD-", "WGD+"), main = "Pan-cancer fractions of the top three aberrations", ylab = "Cumulative fraction of\n3 most recurrent aneuploidies", outline = FALSE, ylim = c(0, 0.42), boxwex=0.8)
for (i in c(1:22)){
  stripchart(WGD0_fractions[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35)
  stripchart(WGD_fractions[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35, at = 2)
}

legend("topright", legend=unique(sort(g_20_types)),
       pt.bg=color_vec, col = c("black"),cex=1.08, pch = 21, xpd = TRUE, inset=c(-0.12,0), pt.cex = 1)

par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos<-yrange[2]-diff(yrange)/10
#segments(1,ypos,2,ypos)
polygon(x = c(1,2,2,1), y = c(ypos,ypos,ypos+0.00125,ypos+0.00125), col = "black", border = NA)
text((1+2)/2,ypos+diff(yrange)/40,"*",cex=3)
par(xpd=FALSE)

dev.off()

t.test(WGD0_fractions,WGD_fractions, paired = T)