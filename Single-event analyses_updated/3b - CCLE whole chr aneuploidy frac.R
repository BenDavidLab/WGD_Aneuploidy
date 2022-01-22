ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Paper outline/Part 2.5/Fig 3"

setwd(ROOT_DIR)


library("stringr")

data <- read.csv("C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/CCLE analysis/Prelim/wca_samplewise_no-ke97.tsv", sep = "\t", header = T)


input_box_1_fil_2 <- data.frame(count_anp_fil = data$Fraction, type_fil = data$Type, wgd_status_fil = data$`WGD.status`)
write.table(input_box_1_fil_2, file = "Table-CCLE-Pan-cancer whole chr anp frac_no-ke97.tsv", sep = "\t")

#The next two lines are used to create spaces between boxes of different cancer types
n <- prod(dim(with(input_box_1_fil_2, table(wgd_status_fil, type_fil))))
VEC <- seq(0.5, n/2, length.out=n)*2 - c(0, .18)


sig_types <- c("SKIN", "URINARY_TRACT", "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "UPPER_AERODIGESTIVE_TRACT") ##upon doing t-test (end of this code), only these 4 types were significant
all_types <- c("BONE", "BREAST", "CENTRAL_NERVOUS_SYSTEM","HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "KIDNEY",  "LARGE_INTESTINE", "LIVER", "LUNG", "OVARY", "PANCREAS" ,"SKIN","STOMACH", "UPPER_AERODIGESTIVE_TRACT", "URINARY_TRACT")


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
color_vec_trans <- c() #will include a colour and its transparent version


for (i in c(1:22)){
  color_vec_trans <- c(color_vec_trans, t_col(color_vec[i] , 40), color_vec[i])
  
}


p_vals <- c()
for(t in sort(unique(input_box_1_fil_2$type_fil)))
{
  cat(t, "\n")
  p_vals <- c(p_vals, t.test(input_box_1_fil_2[input_box_1_fil_2$type_fil==t & input_box_1_fil_2$wgd_status_fil=="WGD-","count_anp_fil"], input_box_1_fil_2[input_box_1_fil_2$type_fil==t & input_box_1_fil_2$wgd_status_fil=="WGD+","count_anp_fil"])$p.value)
  q_vals <- p.adjust(p_vals, method = "BH")
}

####Storing asterisks
asterisks <- c()
for(qval in q_vals)
{
  if( qval <0.0001){asterisks <- c(asterisks, "****")}
  else if( qval <0.001){asterisks <- c(asterisks, "***")}
  else if( qval <0.01){asterisks <- c(asterisks, "**")}
  else if( qval <0.05){asterisks <- c(asterisks, "*")}
  else if( qval >= 0.05){asterisks <- c(asterisks, "n.s.")}
}


draw_sig_bars <- function(x1, x2,index)
{par(xpd=TRUE)
  yrange<-par("usr")[3:4]
  ypos<-yrange[2]-diff(yrange)/10
  segments(x1,ypos,x2,ypos)
  text((x1+x2)/2,ypos+diff(yrange)/40,asterisks[index],cex=1)
  par(xpd=FALSE)}

####setting x-labels for boxplot
tumor_abbrevs <- c("BONE", "BRST", "CENS", "HELT", "KDNY", "LINT",
                   "LIVR", "LUNG", "OVRY", "PNCR", "SKIN", "STMC",
                   "UAET", "URIN")
box_labels <- c()
textbox_labels <- c()
for(abbrev in tumor_abbrevs)
{
  #name <- paste0("+      -\n", abbrev)
  box_labels <- c(box_labels, "-", "+")
  textbox_labels <- c(abbrev, "        ")
}



#Sets appropriate margins for the plot chart


pdf(paste0(ROOT_DIR, "/", "CCLE-pan cancer wholechr anp frac_no-ke97_names1_ns_qvals.pdf"), width = 13, height = 7)
par(mar = c(14.2,4,1.5,1.5))
boxplot(count_anp_fil~wgd_status_fil*type_fil, col=color_vec_trans, data=input_box_1_fil_2, las = 1, ylab = "Fraction of whole-chromosome aneuploidies", cex.lab = 1, cex.axis = 0.8, xlab = "", cex.main = 0.9, main = "Pan-cancer cell line whole chromosome aneuploidy fraction", outline = FALSE, at = VEC, ylim = c(0,1.2), mar = c(8,4,1.5,1.5), names = box_labels)
mtext('BONE', side=1, line=2, at=1.5, cex = 0.85)
mtext('BREAST', side=1, line=2, at=3.5, cex = 0.85)
mtext('CNS', side=1, line=2, at=5.5, cex = 0.85)
mtext('HEMATO', side=1, line=2, at=7.5, cex = 0.85)
mtext('KIDNEY', side=1, line=2, at=9.5, cex = 0.85)
mtext('COLON', side=1, line=2, at=11.5, cex = 0.85)
mtext('LIVER', side=1, line=2, at=13.5, cex = 0.85)
mtext('LUNG', side=1, line=2, at=15.5, cex = 0.85)
mtext('OVARY', side=1, line=2, at=17.5, cex = 0.85)
mtext('PANCREAS', side=1, line=2, at=19.5, cex = 0.85)
mtext('SKIN', side=1, line=2, at=21.5, cex = 0.85)
mtext('STOMACH', side=1, line=2, at=23.5, cex = 0.85)
mtext('\n\nUPPER\nAIRWAY', side=1, line=2.5, at=25.5, cex = 0.85)
mtext('\n\nURINARY\nTRACT', side=1, line=2.5, at=27.5, cex = 0.85)


for (i in 1 : length(all_types))
{t <- all_types[i]
idx <- which(t == sort(as.character(unique(input_box_1_fil_2$type_fil))))
print(idx)
x1 <- 0.2 + idx + 0.972*(idx-1) #0.2 is the offset from the y axis
x2 <- x1 + 1
#cat(x1," ", x2)
draw_sig_bars(x1, x2,i)
}

dev.off()


