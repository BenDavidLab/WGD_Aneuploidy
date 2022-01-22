library(dplyr)
library(SuperExactTest)

##### this is for filtering the data for each tumor type
raw_data<- read.csv2("C://Users//yelie//OneDrive//Desktop//Phd//kavya project//superexacttest//taylor_et_al._arm-level_WGD_TCGA_plus_minus.csv",sep = ",")
reccuring_events_table<- read.table("C://Users//yelie//OneDrive//Desktop//Phd//kavya project//superexacttest//reccuring_event_from_code.csv",sep ="," )
reccuring_events_table <- reccuring_events_table %>% mutate_all(na_if,"")
collumnames<- colnames(raw_data)
collumnames[15:93]<-gsub("[.][.]1","_loss",collumnames[15:93])
collumnames[15:93]<- gsub("[.]","_gain",collumnames[15:93])
collumnames[15:93]<- gsub("X","",collumnames[15:93])
colnames(raw_data)<-collumnames
tumor_type_list <- unique(reccuring_events_table[,1])
for (tumor in tumor_type_list){
  tumor_type = tumor
  raw_data_filtered <- filter(raw_data,raw_data$Type == tumor_type)
  filtered_reccuring_events_table<- filter(reccuring_events_table,reccuring_events_table[,1] == tumor_type )
  reccuring_events <- na.omit(as.character(filtered_reccuring_events_table[1,][2:length(filtered_reccuring_events_table[1,])]))
  filtered<- raw_data_filtered[,reccuring_events]
  reccuring_arm_calls<- cbind(raw_data_filtered[,1],filtered)
  print(tumor_type)
  dir = paste("C://Users//yelie//OneDrive//Desktop//Phd//kavya project//superexacttest//filtered_reccuring_arm_calls//recurring_events//",tumor_type,"_reccuring_events.csv")
  write.csv(reccuring_arm_calls,file = dir,row.names = FALSE)
  super_data<- reccuring_arm_calls
  super_data[is.na(super_data)] <- 0
  for (row in 1:nrow(super_data))
    {
    relevant_row <- super_data[row,]
    replacment<- c("1")
    relevant_row<- replace(relevant_row,relevant_row %in% replacment,relevant_row[1,1])
    super_data[row,] <-relevant_row
    }
  super_data<- na_if(super_data, 0)
  reccuring_list <- as.list(super_data[2:ncol(super_data)])
  super_results<- supertest(reccuring_list,n=nrow(filtered),degree = 2:5)
  dir2 <- paste("C://Users//yelie//OneDrive//Desktop//Phd//kavya project//superexacttest//filtered_reccuring_arm_calls//super_exact_result//",tumor_type,"_super_results.csv")
  write.csv(summary(super_results)$Table,file=dir2,row.names = FALSE)
}

