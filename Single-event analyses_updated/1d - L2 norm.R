ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/Rani_suggestions"
setwd(ROOT_DIR)

data <- read.table("Taylor_et_al._Arm-Level_WGD_TCGA_data.txt", header=T, sep='\t', as.is=T)

tumor_type <- c()
wgd_loss <- c()
non_wgd_loss <- c()
wgd_gain <- c()
non_wgd_gain <- c()

TYPES <- unique(data$Type)
TYPES <- c("ACC", "BLCA", "BRCA", "CESC", "COAD", "ESCA", "GBM", "HNSC", "KIRC", "LGG", "LIHC", "LUAD", "LUSC", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "UCEC")
for (t in 1:length(TYPES)){
  TYPE <- TYPES[t]
  #   
  cat(TYPE, " ", t, "\n")
  #frac_table <- read.table(paste0(TYPE,"combined_table.txt"), header=T, sep='\t', as.is=T)
  type_data <- data[data$Type == TYPE,]
  data_WGD_minus <- type_data[type_data$Genome_doublings == 0, c(14:52)]
  data_WGD_plus <- type_data[type_data$Genome_doublings > 0, c(14:52)]
  
  chr_name <- c()
  WGD_minus_gains <- c()
  WGD_minus_losses <- c()
  WGD_plus_gains <- c()
  WGD_plus_losses <- c()
  
  for(j in c(1:39))
  {
    chr_name <- c(chr_name, colnames(data_WGD_plus)[j])
    WGD_minus_gains <- c(WGD_minus_gains, sum(data_WGD_minus[j] == 1, na.rm = T))
    WGD_minus_losses <- c(WGD_minus_losses, sum(data_WGD_minus[j] == -1, na.rm = T))
    
    WGD_plus_gains <- c(WGD_plus_gains, sum(data_WGD_plus[j] == 1, na.rm = T))
    WGD_plus_losses <- c(WGD_plus_losses, sum(data_WGD_plus[j] == -1, na.rm = T))
  }
  
  frac_table <- cbind.data.frame(chr_name, WGD_minus_gains, WGD_minus_losses, WGD_plus_gains, WGD_plus_losses)
  frac_table$WGD_plus_frac_gains <- frac_table$WGD_plus_gains/sum(frac_table$WGD_plus_gains)
  frac_table$WGD_plus_frac_losses <- frac_table$WGD_plus_losses/sum(frac_table$WGD_plus_losses)
  frac_table$WGD_minus_frac_gains <- frac_table$WGD_minus_gains/sum(frac_table$WGD_minus_gains)
  frac_table$WGD_minus_frac_losses <- frac_table$WGD_minus_losses/sum(frac_table$WGD_minus_losses)
  
  
  L2_wgd_loss <- (frac_table$WGD_plus_frac_losses)[1:39]
  uniform_dist_wgd_loss <- rep((sum(frac_table$WGD_plus_frac_losses)/39), 39)
  
  L2_non_wgd_loss <- (frac_table$WGD_minus_frac_losses)[1:39]
  uniform_dist_non_wgd_loss <- rep((sum(frac_table$WGD_minus_frac_losses)/39), 39)
  
  L2_wgd_gain <- (frac_table$WGD_plus_frac_gains)[1:39]
  uniform_dist_wgd_gain <- rep((sum(frac_table$WGD_plus_frac_gains)/39), 39)
  
  L2_non_wgd_gain <- (frac_table$WGD_minus_frac_gains)[1:39]
  uniform_dist_non_wgd_gain <- rep((sum(frac_table$WGD_minus_frac_gains)/39), 39)
  

  L2_calc_wgd_loss <- l2norm(L2_wgd_loss, uniform_dist_wgd_loss)
  L2_calc_non_wgd_loss <- l2norm(L2_non_wgd_loss, uniform_dist_non_wgd_loss)
  L2_calc_wgd_gain <- l2norm(L2_wgd_gain, uniform_dist_wgd_gain)
  L2_calc_non_wgd_gain <- l2norm(L2_non_wgd_gain, uniform_dist_non_wgd_gain)
  
  wgd_loss <- c(wgd_loss, L2_calc_wgd_loss)
  wgd_gain <- c(wgd_gain, L2_calc_wgd_gain)
  non_wgd_loss <- c(non_wgd_loss, L2_calc_non_wgd_loss)
  non_wgd_gain <- c(non_wgd_gain, L2_calc_non_wgd_gain)
  tumor_type <- c(tumor_type, TYPE)
  
}

df_1l2 <- data.frame(WGD_loss = wgd_loss,
                        Non_WGD_loss = non_wgd_loss,
                        WGD_gain = wgd_gain,
                        Non_WGD_gain = non_wgd_gain,
                        Tumor_type = tumor_type)
#write.table(df_1, file = "Pan-cancer L2 norm gain loss.tsv", sep = "\t")





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

#The vector of colours for each of the 22 cancer types
#color_vec <- c("darkorchid4", "blueviolet","darkblue","lightslateblue","cornflowerblue", "cyan2", "aquamarine2", "chartreuse4","darkolivegreen4","chartreuse1", "khaki1","gold1", "sienna1", "orangered1", "red3" ,"indianred3","tan4","peachpuff3","lightcoral","palevioletred3", "maroon3", "lavenderblush3")
color_vec <- c("chartreuse3", "mediumpurple3", "violetred1", "turquoise3", "blue4", "slateblue2", "gray75", "firebrick4", "chocolate1", "gray49","springgreen4","bisque" , "lavenderblush2", "cyan4", "purple3", "gray32", "deepskyblue1", "blue", "gold", "gray14", "lightslateblue", "salmon")
color_vec_trans <- c() #will include a colour and its transparent version


for (i in c(1:22)){
  color_vec_trans <- c(color_vec_trans, color_vec[i],t_col(color_vec[i] , 40))
  
}

draw_sig_bars <- function(x1, x2)
{par(xpd=TRUE)
  yrange<-par("usr")[3:4]
  ypos<-yrange[2]-diff(yrange)/10
  segments(x1,ypos,x2,ypos)
  text((x1+x2)/2,ypos+diff(yrange)/40,"*",cex=1)
  par(xpd=FALSE)}







plot.new()

par(mar = c(4.4,2.5,1.5,1.5))

pdf(paste0(ROOT_DIR, "/", "L2norm_arms_loss_new_corrected.pdf"), width = 12, height = 5.5) #Renamed to "Pan cancer L2 norm of relative losses WGD vs nonWGD.pdf"

par(mar = c(4.4,4.5,1.5,12))

boxplot(df_1l2$Non_WGD_loss, df_1l2$WGD_loss,col=c("white") ,names = c("WGD-", "WGD+"), main = "Pan-cancer comparison of L2 norm values for gains", ylab = "L2 norm", outline = FALSE, ylim = c(0, 0.4))
for (i in c(1:22)){
  stripchart(df_1l2$Non_WGD_loss[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35)
  stripchart(df_1l2$WGD_loss[i],vertical=TRUE, method="jitter",add=TRUE, pch=21,bg=color_vec[i], col = c("black"), cex = 1, jitter = 0.35, at = 2)
}

legend("topright", legend=unique(sort(TYPES)),
       pt.bg=color_vec, col = c("black"),cex=0.65, pch = 21, xpd = TRUE, inset=c(-0.07,0), pt.cex = 1)

par(xpd=TRUE)
yrange<-par("usr")[3:4]
ypos<-yrange[2]-diff(yrange)/10
segments(1,ypos,2,ypos)
text((1+2)/2,ypos+diff(yrange)/40,"*",cex=1)
par(xpd=FALSE)

dev.off()



t.test(x = df_1l2$WGD_loss,df_1l2$Non_WGD_loss, paired = T ) 
t.test(x = df_1l2$WGD_gain,df_1l2$Non_WGD_gain, paired = T )
