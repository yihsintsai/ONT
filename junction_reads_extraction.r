#/mnt/nas2/yh/R/R-4.4.0/bin/R

lib <- c("CNEr","ape","dplyr","tidyverse","rdist","bioseq","parallel","foreach")
  

for(f in lib){
    if(!require(f,character.only = TRUE)) install.packages(f)
    library(f,character.only = TRUE)
  }

args <- commandArgs(trailingOnly = TRUE)
pathway <- args[1]
sample.table_path <- args[2]
data_file <- args[3]
reference_file <- args[4]

print("prepare data")

table <- read.table(sample.table_path,header=T)%>%`colnames<-`(c("strand","chromosome","position"))
data <- read.table(data_file)%>%`colnames<-`(c("sample_name","chromosome","str_site","end_site","strand","length","qual","cigar"))
reference <- read.table(reference_file)%>%`colnames<-`(c("sample_name","sequence"))
file_name <- tail((str_split(sample.table_path,pattern="/")[[1]]),n=1) %>% gsub(".breakpoints.xls","",.)

print("seq_extraction")
numCores <- parallel::detectCores()-1
cl <- parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl)

dist.matrix <- foreach::foreach(n = (1:nrow(table)), .combine = 'rbind',.packages = c("CNEr","ape","dplyr","tidyverse","rdist","bioseq")) %dopar%  { 
  pos <- table[n,3]
  chromosome <- table[n,2]
  chrpos <- paste(chromosome,pos,sep="_")
  df <- data[(data$chromosome == chromosome) & (data$str_site <= pos ) & (data$end_site >= pos), ]
  
  ##calculate needed position
  df$specific_position_length <- ""
  for ( read in 1:(nrow(df))) {
    cigar <- df[read,8]
    strand <- df[read,5]
    df[read, 8] <- ifelse(strand == "-", reverseCigar(cigar), cigar)   
    df[read, 9] <- ifelse(strand == "-", df[read,4] - pos, pos - df[read,3])
  }
  df <- df[,c(1:3,5,7,8,9)]
  df$specific_len <- ""

 #calcluate position in cigar  
  for ( read in 1:(nrow(df))) {
    cigar <- df[read,6]
    specific_pos <- as.numeric(df[read,7])
    specific_len <- df[read,8]
    cigar <- gsub("(?<=\\D)", ";", cigar , perl=TRUE)
    cigar_split <- strsplit(cigar,";")
   
    sum_num=0
    for (f in 1:length(cigar_split[[1]])){
      if(grepl("I|S|H",cigar_split[[1]][f])){
        next
       }else{
         num=gsub("[A-Z]","",cigar_split[[1]][f])
         sum_num = as.numeric(sum_num) + as.numeric(num)
        }
      if (sum_num >=  specific_pos){
        extra <- sum_num - specific_pos
        colnum <- f
        break }
    } 

   ##calculate real position in cigar 
    sum_num = 0
    for (i in 1:colnum){
      if(grepl("D",cigar_split[[1]][i])){
        next
       }else{
         num=gsub("[A-Z]","",cigar_split[[1]][i])
         sum_num = as.numeric(sum_num) + as.numeric(num)
        }
      if (colnum == 1) {
         df[read,8] <- specific_pos
       }else {
         df[read,8] <- sum_num - extra
        }
      }
    }

   #extract sequence 
   df$upstream <- ""
   df$downstream <- ""
   df$length_seq<- ""
   for ( seq in 1:(nrow(df))) {
    name <- df[seq,1]
    len_str <- as.numeric(df[seq,8])
    reference_seq <- reference[grep(paste(name,"$",sep=""),reference$sample_name),]
    df[seq, 11] <- str_length(reference_seq$sequence)
    strand <-  df[seq, 4] 
    if (strand == "+") { 
      up_str_pos <- len_str - 1000 
      up_end_pos <- len_str - 1
      down_str_pos <- len_str 
      down_end_pos <- len_str + 999
      if ( up_str_pos < 0 || down_end_pos > df[read, 11]){
        if ( up_str_pos < 0 ) { 
          up_str_pos = 1
        }else {
          down_end_pos == df[seq, 11]
         }
       }
      df[seq, 9] <- substr(reference_seq$sequence, up_str_pos , up_end_pos)
      df[seq, 10] <- substr(reference_seq$sequence, down_str_pos , down_end_pos)
     }else{
      up_str_pos <- len_str + 2 
      up_end_pos <- len_str + 1001
      down_str_pos <- len_str - 998  
      down_end_pos <- len_str + 1
      if ( down_end_pos < 0 || up_str_pos > df[seq, 11]){
        if ( down_end_pos < 0 ) { 
          up_str_pos = 1
        }else {
          up_str_pos == df[seq, 11]
         }
       }
      df[seq, 9] <- substr(reference_seq$sequence, up_str_pos , up_end_pos)
      df[seq, 9] <- seq_reverse(dna(df[seq,9]))
      df[seq, 9] <- seq_complement(dna(df[seq,9]))
      
      df[seq, 10] <- substr(reference_seq$sequence, down_str_pos , down_end_pos)
      df[seq, 10] <- seq_reverse(dna(df[seq,10]))
      df[seq, 10] <- seq_complement(dna(df[seq,10]))
     } 
   }
  df$breakpoint <- chrpos
  df <- df
}
df_table<- data.frame(dist.matrix)
df_table <- df_table %>%
  mutate_all(~ replace(., . == "" | is.na(.), NA))

write.table(df_table,paste(pathway,file_name,".bp.txt",sep=""),row.names = T)
parallel::stopCluster(cl)
print("finish")
