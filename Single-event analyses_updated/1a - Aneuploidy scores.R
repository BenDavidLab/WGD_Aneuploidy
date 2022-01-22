ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/Pan-cancer-graphs"

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

CA_WGD0_no_acro<-cbind(CA_WGD0[,1:37], CA_WGD0[,41:50])
CA_WGD_no_acro<-cbind(CA_WGD[,1:37], CA_WGD[,41:50])
CA_WGD0_no_acro[is.na(CA_WGD0_no_acro)]<-0
CA_WGD_no_acro[is.na(CA_WGD_no_acro)]<-0
CA_WGD0_chr_arm_ratio<-c()

CA_WGD0_rearr<-cbind(CA_WGD0[,14:37], CA_WGD0[,41:50], CA_WGD0[,38:40], CA_WGD0[, 51:52])
CA_WGD_rearr<-cbind(CA_WGD[,14:37], CA_WGD[,41:50], CA_WGD[,38:40], CA_WGD[, 51:52])
#############################

count_anp_WGD0 <- c()
for (j in 1:nrow(CA_WGD0_rearr)){
  count_anp_sample <- 0
  #print(CA_WGD0[j,])
  for (i in seq(1, 34, 2)) ###thr first 34 columns correspond to non-acrocentric chr in our rearrnged df
  {##counting whole chr anp
    #print(CA_WGD0_rearr[j,i])
    if(!is.na(CA_WGD0_rearr[j,i]) && !is.na(CA_WGD0_rearr[j,i+1]) && (CA_WGD0_rearr[j,i] == -1 | CA_WGD0_rearr[j,i] == 1) && (CA_WGD0_rearr[j,i] == CA_WGD0_rearr[j,i+1])) 
    {
      count_anp_sample <- count_anp_sample + 1
    }
    ##counting the arms as separate aneuploidies
    else
    {
      if(!is.na(CA_WGD0_rearr[j,i]) && (CA_WGD0_rearr[j,i] == -1 | CA_WGD0_rearr[j,i] == 1))
      {
        count_anp_sample <- count_anp_sample + 1 
      }
      if(!is.na(CA_WGD0_rearr[j,i+1]) && (CA_WGD0_rearr[j,i+1] == -1 | CA_WGD0_rearr[j,i+1] == 1))
      {
        count_anp_sample <- count_anp_sample + 1 
      }
    }
  }
  for (i in seq(35, 39)) ###thr first 34 columns correspond to non-acrocentric chr in our rearrnged df
  {##counting whole chr anp
    #print(CA_WGD0_rearr[j,i])

    ##counting the arms as separate aneuploidies

    if(!is.na(CA_WGD0_rearr[j,i]) && (CA_WGD0_rearr[j,i] == -1 | CA_WGD0_rearr[j,i] == 1))
    {
      count_anp_sample <- count_anp_sample + 1 
    }

    
  }
  count_anp_WGD0 <- c(count_anp_WGD0, count_anp_sample)
}


count_anp_WGD <- c()
for (j in 1:nrow(CA_WGD_rearr)){
  count_anp_sample <- 0
  #print(CA_WGD0[j,])
  for (i in seq(1, 34, 2)) ###thr first 47 columns correspond to non-acrocentric chr in our rearrnged df
  {##counting whole chr anp
    if(!is.na(CA_WGD_rearr[j,i]) && !is.na(CA_WGD_rearr[j,i+1]) && (CA_WGD_rearr[j,i] == -1 | CA_WGD_rearr[j,i] == 1) && (CA_WGD_rearr[j,i] == CA_WGD_rearr[j,i+1])) 
    {
      count_anp_sample <- count_anp_sample + 1
    }
    ##counting the arms as separate aneuploidies
    else
    {
      if(!is.na(CA_WGD_rearr[j,i]) && (CA_WGD_rearr[j,i] == -1 | CA_WGD_rearr[j,i] == 1))
      {
        count_anp_sample <- count_anp_sample + 1 
      }
      if(!is.na(CA_WGD_rearr[j,i+1]) && (CA_WGD_rearr[j,i+1] == -1 | CA_WGD_rearr[j,i+1] == 1))
      {
        count_anp_sample <- count_anp_sample + 1 
      }
    }
  }
  for (i in seq(35, 39)) ###thr first 34 columns correspond to non-acrocentric chr in our rearrnged df
  {##counting whole chr anp
    #print(CA_WGD0_rearr[j,i])
    
    ##counting the arms as separate aneuploidies
    
    if(!is.na(CA_WGD_rearr[j,i]) && (CA_WGD_rearr[j,i] == -1 | CA_WGD_rearr[j,i] == 1))
    {
      count_anp_sample <- count_anp_sample + 1 
    }
    
    
  }
  count_anp_WGD <- c(count_anp_WGD, count_anp_sample)
}


##################################################


count_anp <- as.numeric(c(count_anp_WGD0, count_anp_WGD))
wgd_status <- c(CA_WGD0$Genome_doublings, CA_WGD$Genome_doublings)
type <- c(CA_WGD0$Type, CA_WGD$Type)

#initial data frame which will serve as input to the boxplot fn
input_box_1 <- data.frame(count_anp = count_anp, wgd_status = wgd_status, type = type)
#input_box_1[1:20, ]

#replacing 0,1,2 with "non WGD", "WGD"
input_box_1$wgd_status <- with(input_box_1, replace(wgd_status, wgd_status==0, "WGD-"))
input_box_1$wgd_status <- with(input_box_1, replace(wgd_status, wgd_status==1, "WGD+"))
input_box_1$wgd_status <- with(input_box_1, replace(wgd_status, wgd_status==2, "WGD+"))



input_box_1_fil <- input_box_1

#filtering out those cancer types with less than 20  samples in WGD/non WGD category
for (t in TYPES){
  num_wgd = nrow(input_box_1_fil[input_box_1_fil$type == t & input_box_1_fil$wgd_status == "WGD+", ])
  num_non_wgd = nrow(input_box_1_fil[input_box_1_fil$type == t & input_box_1_fil$wgd_status == "WGD-", ])
  cat( t, " ", num_wgd, " ", num_non_wgd, "\n")
  if (num_wgd < 20 | num_non_wgd < 20){
    input_box_1_fil <- input_box_1_fil[input_box_1_fil$type != t, ]
    cat(t, " filtered out \n")
  }
}
count_anp_fil <- input_box_1_fil$count_anp
wgd_status_fil <- as.factor(input_box_1_fil$wgd_status)
type_fil <- input_box_1_fil$type


nrow(input_box_1_fil)
#input_box_1_fil <- input_box_1_fil[is.na(input_box_1_fil$count_anp) != TRUE, ]
#nrow(input_box_1_fil)


#plot.window(xlim = c(0,1000), ylim = c(0,1000)) 


type_fil <- as.factor(type_fil)
#Without dropping unused levels, the names of all cancer types, even those filtered out, appear along the x axis in the boxplot
type_fil <- droplevels(type_fil)
#wgd_status_fil <- as.factor(wgd_status_fil)


input_box_1_fil_2 <- data.frame(count_anp_fil = count_anp_fil, type_fil = type_fil, wgd_status_fil = wgd_status_fil )
write.table(input_box_1_fil_2, "Pan-cancer aneuploidy scores-final-wca.tsv", sep = "\t")
#The next two lines are used to create spaces between boxes of different cancer types
n <- prod(dim(with(input_box_1_fil_2, table(wgd_status_fil, type_fil))))
VEC <- seq(0.5, n/2, length.out=n)*2 - c(0, .18)

#sig_types = c("CESC", "BRCA", "COAD", "ESCA", "HNSC", "KIRC", "LIHC", "LUAD", "LUSC","OV", "PAAD", "PCPG", "PRAD", "SARC", "SKCM", "STAD", "UCEC") #this vector stores the names of those types whose t-tests indicate a significant difference between WGD and non-WGD groups - I individually looked at the t-test for each type to compile this (was not done programmatically) 
sig_types <- sort(unique(input_box_1_fil_2$type_fil)) ##upon doing t-test (end of this code), all 22 tumor types were significant

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
  color_vec_trans <- c(color_vec_trans, t_col(color_vec[i] , 40), color_vec[i])
  
}


#my_color_vec <- c("lightpink2","lightpink3", "plum2", "plum3", "palevioletred", "palevioletred4", "mediumorchid", "purple4", "thistle1", "thistle2", "slateblue1", "slateblue4", "skyblue3", "skyblue4", "cyan3", "cyan4","cornflowerblue", "blue4","slategray1", "slategray3", "paleturquoise1", "paleturquoise3", "turquoise", "turquoise4", "olivedrab2", "olivedrab4", "chartreuse3", "chartreuse4","lightgoldenrod1", "lightgoldenrod3", "darkgoldenrod1", "darkgoldenrod3", "yellow2", "yellow3", "tan1", "tan2", "darkorange2", "darkorange3", "tomato", "tomato3", "red3", "red4", "indianred3", "indianred4","wheat1", "wheat3", "tan3", "tan4", "peachpuff3", "peru", "bisque3", "bisque4", "salmon1", "salmon3","pink3", "pink4", "mistyrose1", "mistyrose3","lavenderblush3" ,"lavenderblush4","lightsteelblue3", "lightsteelblue4" )
#my_color_vec <- c("lightpink2","lightpink3", "plum2", "plum3", "palevioletred", "palevioletred4", "mediumorchid", "purple4", "thistle1", "thistle2", "slateblue1", "slateblue4", "skyblue3", "skyblue4", "cyan3", "cyan4","cornflowerblue", "blue4","slategray1", "slategray3", "paleturquoise1", "paleturquoise3", "turquoise", "turquoise4", "olivedrab2", "olivedrab4", "chartreuse3", "chartreuse4","lightgoldenrod1", "lightgoldenrod3", "darkgoldenrod1", "darkgoldenrod3", "yellow2", "yellow3", "tan1", "tan2", "darkorange2", "darkorange3", "tomato", "tomato3", "red3", "red4", "indianred3", "indianred4","wheat1", "wheat3", "tan3", "tan4", "bisque3", "bisque4", "salmon1", "salmon3","pink3", "pink4", "mistyrose1", "mistyrose3","lightsteelblue3", "lightsteelblue4" )

####conducting t-test for each tumor type
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
tumor_abbrevs <- sort(unique(input_box_1_fil_2$type_fil))
box_labels <- c()
textbox_labels <- c()
for(abbrev in tumor_abbrevs)
{
  #name <- paste0("+      -\n", abbrev)
  box_labels <- c(box_labels, "-", "+")
  textbox_labels <- c(abbrev, "        ")
}



#Sets appropriate margins for the plot chart
par(mar = c(4.4,4,1.5,1.5))

pdf(paste0(ROOT_DIR, "/", "pcas1_wca_qvals.pdf"), width = 13, height = 4.3)
boxplot(count_anp_fil~wgd_status_fil*type_fil, col=color_vec_trans, data=input_box_1_fil_2, las = 1, ylab = "Number of aneuploidies", cex.lab = 1, cex.axis = 0.8, xlab = "", cex.main = 0.9, main = "Pan-cancer aneuploidy scores", outline = FALSE, at = VEC, ylim = c(0,43), names = box_labels)

mtext('ACC', side=1, line=2, at=1.5, cex = 0.85)
mtext('BLCA', side=1, line=2, at=3.5, cex = 0.85)
mtext('BRCA', side=1, line=2, at=5.5, cex = 0.85)
mtext('CESC', side=1, line=2, at=7.5, cex = 0.85)
mtext('COAD', side=1, line=2, at=9.5, cex = 0.85)
mtext('ESCA', side=1, line=2, at=11.5, cex = 0.85)
mtext('GBM', side=1, line=2, at=13.5, cex = 0.85)
mtext('HNSC', side=1, line=2, at=15.5, cex = 0.85)
mtext('KIRC', side=1, line=2, at=17.5, cex = 0.85)
mtext('LGG', side=1, line=2, at=19.5, cex = 0.85)
mtext('LIHC', side=1, line=2, at=21.5, cex = 0.85)
mtext('LUAD', side=1, line=2, at=23.5, cex = 0.85)
mtext('LUSC', side=1, line=2, at=25.5, cex = 0.85)
mtext('OV', side=1, line=2, at=27.5, cex = 0.85)
mtext('PAAD', side=1, line=2, at=29.5, cex = 0.85)
mtext('PCPG', side=1, line=2, at=31.5, cex = 0.85)
mtext('PRAD', side=1, line=2, at=33.5, cex = 0.85)
mtext('READ', side=1, line=2, at=35.5, cex = 0.85)
mtext('SARC', side=1, line=2, at=37.5, cex = 0.85)
mtext('SKCM', side=1, line=2, at=39.5, cex = 0.85)
mtext('STAD', side=1, line=2, at=41.5, cex = 0.85)
mtext('UCEC', side=1, line=2, at=43.5, cex = 0.85)

for (i in 1 : length(sig_types))
{t <- sig_types[i]
idx <- which(t == sort(as.character(unique(input_box_1_fil_2$type_fil))))
print(idx)
x1 <- 0.2 + idx + 0.972*(idx-1) #0.2 is the offset from the y axis
x2 <- x1 + 1
#cat(x1," ", x2)
draw_sig_bars(x1, x2,i)
}


# legend("topleft", legend=unique(sort(type_fil)),
#        col=color_vec, cex=0.51, pch = 15, inset=c(1,0.1), xpd = TRUE, bty = 'o')
dev.off()

# ####conducting t-test for each tumor type
# p_vals <- c()
# for(t in unique(input_box_1_fil_2$type_fil))
# {
#   cat(t, "\n")
#   p_vals <- c(p_vals, t.test(input_box_1_fil_2[input_box_1_fil_2$type_fil==t & input_box_1_fil_2$wgd_status_fil=="WGD-","count_anp_fil"], input_box_1_fil_2[input_box_1_fil_2$type_fil==t & input_box_1_fil_2$wgd_status_fil=="WGD+","count_anp_fil"])$p.value)
# }
