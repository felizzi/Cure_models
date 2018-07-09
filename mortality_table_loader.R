install.packages("RCurl")
library("RCurl")

table_gen <- function(str_Country_Code, str_Country, sex, parent_dir){
  
  if (sex == 'M'){
    url_name <- paste("http://www.mortality.org/hmd",str_Country_Code,"STATS/mltper_1x1.txt",sep="/")
  }else{
    url_name <- paste("http://www.mortality.org/hmd",str_Country_Code,"STATS/fltper_1x1.txt",sep="/")
  }
  aut_data <- getURL(url_name,userpwd = "federico.felizzi@roche.com:duilio01#")
  fdata <- "libraries/app_aut.txt"
  write(aut_data, file = "libraries/app_aut.txt")
  
  con=file(fdata,open="r")
  line=readLines(con) 
  long=length(line)
  #line with names
  entries <- unlist(strsplit(line[3], ' '));
  entry_names_id <- which(entries != "")
  entry_names <- entries[entry_names_id]
  
  #initialize an empty data-frame 
  aut_frame <- data.frame(t(rep(NA,length(entry_names_id))))
  names(aut_frame) <- entry_names
  aut_frame <- aut_frame[-1,]
  
  for (i in 4:long){
    entries <- unlist(strsplit(line[i], ' '));
    ids <- which(entries != "")
    if (length(ids) > 0){
      entry_values <- entries[ids];
      for (j in 1:length(ids)){
        aut_frame[i-3,entry_names[j]] <- entry_values[j]
      }
    }
  }
  close(con)
  #write the dataframe  
  if (sex == 'M'){
    write.csv(aut_frame, file = paste(parent_dir,str_Country,"Male","Mortality.csv", sep ="/"),sep =",");
  }else{
    write.csv(aut_frame, file = paste(parent_dir,str_Country,"Female","Mortality.csv", sep ="/"),sep =",");
  }
}
