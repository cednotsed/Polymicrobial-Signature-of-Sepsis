setwd("~/git_repos/Polymicrobial-Signature-of-Sepsis/")
require(data.table); require(tidyverse)
dir_list <- list.dirs("datasets/kraken2_reports")
dir_list <- dir_list[grepl("maxi", dir_list)]
data_names <- str_split(dir_list, "/", simplify = T)[, 3]
data_names <- str_split(data_names, "_", simplify = T)[, 1]

perc_data <- list()

for (i in seq(length(dir_list))) {
  base_dir <- dir_list[i]
  perc_list <- c()
  file_list <- list.files(base_dir)
  
  for (j in seq(length(file_list))) {
    file <- file_list[j]
    report <- fread(paste0(base_dir, "/", file)) %>%
      filter(V6 == "unclassified")
    
    if (nrow(report) == 0) {
      perc_unclassified <- 0
    } else {
      perc_unclassified <- report$V1
    }
    
    perc_list <- c(perc_list, perc_unclassified)
  }
  
  perc_data[[data_names[i]]] <- perc_list 
}


final_df <- data.frame()

for (n in names(perc_data)) {
  percs <- perc_data[[n]]
  print(length(percs))
  avg <- mean(percs)
  range_low <- range(percs)[1]
  range_high <- range(percs)[2]
  morsel <- data.frame(avg = c(avg), range_low = c(range_low), range_high = c(range_high))
  final_df <- rbind(final_df, morsel)
}

rownames(final_df) <- names(perc_data)

fwrite(final_df, "results/proportion_unclassified_reads.csv", row.names = T)

