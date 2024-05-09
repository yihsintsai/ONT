
 lib <- c("CNEr","bioseq","ape","dplyr","tidyverse","rdist","TraMineR","seqinr","msa")
  for(f in lib){
    if(!require(f,character.only = TRUE)) install.packages(f)
    library(f,character.only = TRUE)
  }

#### input 
table <- data.frame(pos_name=c("A40_chr7_54658974","A40_chr7_54753591","A40_chr7_55477848","A40_chr7_55636272"),
             breakpoint_strand = c("+","+","-","-"))

pathway <- "/mnt/nas2/yh/14.ONT/soft_clip/A40/"

for (n in (1:nrow(table))){
    pos_name <- table[n,1]
    pos <- as.numeric(gsub(".*_","",pos_name))
    data <- read.table(paste(pathway,pos_name,".bed",sep=""))%>%`colnames<-`(c("sample_name","chromosome","str_site","end_site","strand","length","qual","cigar"))
    reference <- read.table(paste(pathway,pos_name,"_sample.sequence",sep=""))%>%`colnames<-`(c("sample_name","sequence"))
    breakpoint_strand <- table[n,2]
    df <- data
      df$specific_position_length <- ""
  
##extract 1000bps upstream and downstream of breakpoint    
    for ( read in 1:(nrow(df))) {
      cigar <- df[read,8]
      strand <- df[read,5]
      len <- df[read,9]
      df[read, 8] <- ifelse(strand == "-", reverseCigar(cigar), cigar)
      df[read, 9] <- ifelse(strand == "-", df[read,4] - pos, pos - df[read,3]  )
    }
  
    df <- df[,c(1:3,5,7:9)]
      df$specific_len <- ""
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
          print(colnum)
          print(cigar_split[[1]][f])
          print(sum_num)
          break}
      }
      sum_num = 0
        for (i in 1:colnum){
          if(grepl("D",cigar_split[[1]][i])){
            next
          }else{
            num=gsub("[A-Z]","",cigar_split[[1]][i])
            sum_num = as.numeric(sum_num) + as.numeric(num)
          }
        }
          if (colnum == 1) {
            df[read,8] <- sum_num + specific_pos
          } else {
            df[read,8] <- sum_num - extra 
          }
    }

  df$upstream <- ""
  df$downstream <- ""
  df$length_seq<- ""
  for ( read in 1:(nrow(df))) {
    name <- df[read,1]
    len_str <- as.numeric(df[read,8])
    reference_seq <- reference[grep(name,reference$sample_name),]
    df[read, 11] <- str_length(reference_seq$sequence)
    strand <-  df[read, 4] 
    if (strand == "+") { 
      up_str_pos <- len_str - 1000 
      up_end_pos <- len_str - 1
      down_str_pos <- len_str 
      down_end_pos <- len_str + 999
      if ( up_str_pos < 0 || down_end_pos > df[read, 11]){
        if ( up_str_pos < 0 ) { 
          up_str_pos = 1
        } else {
          down_end_pos == df[read, 11]
        }
      }
      df[read, 9] <- substr(reference_seq$sequence, up_str_pos , up_end_pos)
      df[read, 10] <- substr(reference_seq$sequence, down_str_pos , down_end_pos)
    } else {
      up_str_pos <- len_str + 2 
      up_end_pos <- len_str + 1001
      down_str_pos <- len_str - 998  
      down_end_pos <- len_str + 1
      if ( down_end_pos < 0 || up_str_pos > df[read, 11]){
        if ( down_end_pos < 0 ) { 
          up_str_pos = 1
        } else {
          up_str_pos == df[read, 11]
        }
      }
      df[read, 9] <- substr(reference_seq$sequence, up_str_pos , up_end_pos)
      df[read, 9] <- seq_reverse(dna(df[read,9]))
      df[read, 9] <- seq_complement(dna(df[read,9]))
      df[read, 10] <- substr(reference_seq$sequence, down_str_pos , down_end_pos)
      df[read, 10] <- seq_reverse(dna(df[read,10]))
      df[read, 10] <- seq_complement(dna(df[read,10]))
    }
  } 
  df$breakpoint <- pos
  write.table(df,paste(pathway,pos_name,"_table.txt",sep=""),sep="\t",row.names = F)
  }
