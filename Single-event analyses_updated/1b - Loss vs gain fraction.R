ROOT_DIR = "C:/Users/kiwii/OneDrive/Documents/Kavya/Uri/Work/Tasks/Pan-cancer-graphs"
setwd(ROOT_DIR)

data <- read.table("TCGA_mastercalls.abs_tables_JSedit.fixed.txt", header=T, sep='\t', as.is=T)

library(dplyr)
library(ggplot2)

loss_more <- c()
gain_more <- c()
wgd_status <- c()
tumor_type <- c()
chi_square_val <- c()
chi_square_p <- c()
fisher_p <- c()

#TYPES <- unique(data$Type)
#TYPES <- c("ACC")
#g_20_types <- c()
g_20_types = c("ACC","BLCA","BRCA", "CESC", "COAD", "ESCA","GBM" ,"HNSC", "KIRC","LGG" ,"LIHC", "LUAD", "LUSC","OV", "PAAD", "PCPG", "PRAD", "READ","SARC", "SKCM", "STAD", "UCEC")


data_tumor <- data 

#CA_WGD0 <- data %>%  filter(Genome_doublings == 0) %>% select(c(53:69), c(38:40), c(51:52))#change BRCA later
#CA_WGD1 <- data %>%  filter(Genome_doublings > 0) %>% select(c(53:69), c(38:40), c(51:52))#change BRCA later

CA_WGD0 <- data_tumor %>%  filter(Genome.doublings == 0) %>% select(c(5:6))#change BRCA later
CA_WGD1 <- data_tumor %>%  filter(Genome.doublings == 1) %>% select(c(5:6))#change BRCA later
CA_WGD2 <- data_tumor %>%  filter(Genome.doublings == 2) %>% select(c(5:6))

CA_WGD0_loss <- CA_WGD0 %>% filter(ploidy < 2)
CA_WGD0_gain <- CA_WGD0 %>% filter(ploidy > 2)

CA_WGD1_loss <- CA_WGD1 %>% filter(ploidy < 4)
CA_WGD1_gain <- CA_WGD1 %>% filter(ploidy > 4)

CA_WGD2_loss <- CA_WGD2 %>% filter(ploidy < 8)
CA_WGD2_gain <- CA_WGD2 %>% filter(ploidy > 8)

CA_WGD0_total <- nrow(CA_WGD0_loss) + nrow(CA_WGD0_gain)
CA_WGD1_total <- nrow(CA_WGD1_loss) + nrow(CA_WGD1_gain)
CA_WGD2_total <- nrow(CA_WGD2_loss) + nrow(CA_WGD2_gain)

plot_data <- data.frame(WGD_status = c("WGD-", "WGD-", "WGD_1", "WGD_1", "WGD_2", "WGD_2"),
                        Loss_gain_status = rep(c('Loss', "Gain"), 3),
                        Count = c(nrow(CA_WGD0_loss)/CA_WGD0_total, nrow(CA_WGD0_gain)/CA_WGD0_total,
                                  nrow(CA_WGD1_loss)/CA_WGD1_total, nrow(CA_WGD1_gain)/CA_WGD1_total,
                                  nrow(CA_WGD2_loss)/CA_WGD2_total, nrow(CA_WGD2_gain)/CA_WGD2_total))

ggplot(plot_data, aes(fill=Loss_gain_status, y=Count, x=WGD_status)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#ed1e1e", "#1e50ed"))


# ggplot(data_tumor,aes(x=`ploidy`)) + 
#   geom_histogram(data=CA_WGD0,fill = "red", alpha = 0.2) +
#   geom_histogram(data=CA_WGD1,fill = "blue", alpha = 0.2) +
#   geom_histogram(data=CA_WGD2,fill = "green", alpha = 0.2)

#ggsave("WGD ploidy distribution bar.pdf", width = 7, height = 7)

chisq <- matrix(data = rep(NA, 6), nrow = 2)
chisq[1,1] <- nrow(CA_WGD0_loss)
chisq[2,1] <- nrow(CA_WGD0_gain)
chisq[1,2] <- nrow(CA_WGD1_loss)
chisq[2,2] <- nrow(CA_WGD1_gain)
chisq[1,3] <- nrow(CA_WGD2_loss)
chisq[2,3] <- nrow(CA_WGD2_gain)

chisq.test(chisq)
