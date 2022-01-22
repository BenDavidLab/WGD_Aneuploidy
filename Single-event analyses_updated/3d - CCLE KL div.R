ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Paper outline/Part 2.5/Fig 3"
setwd(ROOT_DIR)

data <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_data.txt", header=T, sep='\t', as.is=T)

tumor_type <- c()
wgd_loss <- c()
non_wgd_loss <- c()
wgd_gain <- c()
non_wgd_gain <- c()

TYPES <- unique(data$Type)

TYPES <- c("BONE", "BREAST", "CENTRAL_NERVOUS_SYSTEM", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
           "KIDNEY", "LARGE_INTESTINE", "LIVER", "LUNG", "OVARY", "PANCREAS", "SKIN", "STOMACH",
           "UPPER_AERODIGESTIVE_TRACT", 'URINARY_TRACT')

for (t in 1:length(TYPES)){
  #   
  TYPE <- TYPES[t]
  #   
  cat(TYPE, " ", t, "\n")
  input_df <- read.table("C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/CCLE analysis/Prelim/Type_wise_df_no-ke97.tsv", header=T, sep='\t', as.is=T)
  
  wgd_minus <- input_df[(input_df$ploidy < 2.5) & (grepl(TYPE, input_df$CCLE_ID, fixed = TRUE)),]
  wgd_plus <- input_df[(input_df$ploidy > 3) & (grepl(TYPE, input_df$CCLE_ID, fixed = TRUE)),]
  
  wgd_minus_loss <- c()
  wgd_plus_loss <- c()
  
  wgd_minus_gain <- c()
  wgd_plus_gain <- c()
  
  for(i in c(5:43))
  {
    wgd_minus_loss <- c(wgd_minus_loss, sum(wgd_minus[, i]==-1))
    wgd_plus_loss <- c(wgd_plus_loss, sum(wgd_plus[, i]==-1))
    
    wgd_minus_gain <- c(wgd_minus_gain, sum(wgd_minus[, i]==1))
    wgd_plus_gain <- c(wgd_plus_gain, sum(wgd_plus[, i]==1))
  }
  
  
  KL_wgd_loss <- (wgd_plus_loss)[1:39]
  uniform_dist_wgd_loss <- rep((sum(wgd_plus_loss)/39), 39)
  
  KL_non_wgd_loss <- (wgd_minus_loss)[1:39]
  uniform_dist_non_wgd_loss <- rep((sum(wgd_minus_loss)/39), 39)
  
  KL_wgd_gain <- (wgd_plus_gain)[1:39]
  uniform_dist_wgd_gain <- rep((sum(wgd_plus_gain)/39), 39)
  
  KL_non_wgd_gain <- (wgd_minus_gain)[1:39]
  uniform_dist_non_wgd_gain <- rep((sum(wgd_minus_gain)/39), 39)
  

  library(entropy)

  
  KL_calc_wgd_loss <- KL.plugin(KL_wgd_loss, uniform_dist_wgd_loss)
  KL_calc_non_wgd_loss <- KL.plugin(KL_non_wgd_loss, uniform_dist_non_wgd_loss)
  KL_calc_wgd_gain <- KL.plugin(KL_wgd_gain, uniform_dist_wgd_gain)
  KL_calc_non_wgd_gain <- KL.plugin(KL_non_wgd_gain, uniform_dist_non_wgd_gain)
  
  wgd_loss <- c(wgd_loss, KL_calc_wgd_loss)
  wgd_gain <- c(wgd_gain, KL_calc_wgd_gain)
  non_wgd_loss <- c(non_wgd_loss, KL_calc_non_wgd_loss)
  non_wgd_gain <- c(non_wgd_gain, KL_calc_non_wgd_gain)
  tumor_type <- c(tumor_type, TYPE)
  
}

df_1 <- data.frame(WGD_loss = wgd_loss,
                   Non_WGD_loss = non_wgd_loss,
                   WGD_gain = wgd_gain,
                   Non_WGD_gain = non_wgd_gain,
                   Tumor_type = tumor_type)


#The next two lines are used to create spaces between boxes of different cancer types
# n <- prod(dim(with(df_1, table(wgd_status_fil, type_fil))))
# VEC <- seq(0.5, n/2, length.out=n)*2 - c(0, .18)

#sig_types = c("CESC", "BRCA", "COAD", "ESCA", "HNSC", "KIRC", "LIHC", "LUAD", "LUSC","OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "UCEC") #this vector stores the names of those types whose t-tests indicate a significant difference between WGD and non-WGD groups - I individually looked at the t-test for each type to compile this (was not done programmatically) 


#This function is used to introduce transparency to an input colour
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
  
  return(t.col)
}
## END



color_vec <- c("gold", "violetred1", "gray75", "darkorange", "#ff4302", "blue4", "springgreen4", "lavenderblush2", "cyan4", "purple3", "gray14", "lightslateblue", "firebrick4", "yellow")
color_vec_trans <- c()

for (i in c(1:14)){
  color_vec_trans <- c(color_vec_trans, color_vec[i],t_col(color_vec[i] , 40))
  
}


draw_sig_bars <- function(x1, x2)
{par(xpd=TRUE)
  yrange<-par("usr")[3:4]
  ypos<-yrange[2]-diff(yrange)/10
  segments(x1,ypos,x2,ypos)
  text((x1+x2)/2,ypos+diff(yrange)/40,"**",cex=1)
  par(xpd=FALSE)}



plot.new()

#width 12, inset = -0.07, boxwex = 0.8, par(mar = c(4.4,4.5,1.5,12))

pdf(paste0(ROOT_DIR, "/", "KL_arms_loss_new_ccle_1_no-ke97_names_sz_check.pdf"), width = 13.4, height = 5.5) #Renamed to "Pan cancer KL div of relative losses WGD vs nonWGD.pdf"

par(mar = c(4.4,4.5,1.5,13))
par(cex.lab=1.5)
par(cex.axis=1.5)

boxplot(df_1$Non_WGD_loss, df_1$WGD_loss, col = c("white"),names = c("WGD-", "WGD+"), main = "Pan-cancer comparison of KL divergence values", ylab = "KL divergence score (losses)", outline = FALSE, ylim = c(0, 1.4),boxwex = 0.8)
for (i in c(1:14)){
  stripchart(df_1$Non_WGD_loss[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35)
  stripchart(df_1$WGD_loss[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35, at = 2)
}

TYPES_names <- c("BONE", "BREAST", "CNS", "HEMATO",
                 "KIDNEY", "COLON", "LIVER", "LUNG", "OVARY", "PANCREAS", "SKIN", "STOMACH",
                 "UPPER AIRWAY", 'URINARY TRACT')

legend("topright", legend=unique(sort(TYPES_names)),
       pt.bg=color_vec, col = c("black"),cex=1.08, pch = 21, xpd = TRUE, inset=c(-0.22,0), pt.cex = 1)

par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos<-yrange[2]-diff(yrange)/10
#segments(1,ypos,2,ypos)
polygon(x = c(1,2,2,1), y = c(ypos,ypos,ypos+0.005,ypos+0.005), col = "black", border = NA)
text((1+2)/2,ypos+diff(yrange)/40,"*",cex=3)
par(xpd=FALSE)


dev.off()



#t.test(x = df_1$WGD_gain,df_1$Non_WGD_gain, paired = T ) ##p val is p-value = 4.163e-07 so significance bar is drawn in the plot
