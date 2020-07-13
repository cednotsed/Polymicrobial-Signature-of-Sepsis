setwd("~/git_repos/Polymicrobial-Signature-of-Sepsis/datasets")

# Grumaz 2019
df <- read.csv(file = "PRJEB30958.txt", sep = "\t", stringsAsFactors = F)
df <- df[, c("run_accession", "library_name", "fastq_ftp")]
df <- df[grep("_T0", df$library_name), ]
df$study <- "Grumaz_2019"
df$y <- "septic"

df2 <- read.csv(file = "PRJEB21872.txt", sep = "\t", stringsAsFactors = F)
df2 <- df2[, c("run_accession", "library_name", "fastq_ftp")]
df2 <- df2[grep("_T0", df2$library_name), ]
df2 <- df2[-grep("POP", df2$library_name), ]
df2$study <- "Grumaz_2019"
df2$y <- "septic"

# Grumaz 2016 
df3 <- read.csv(file = "PRJEB13247.txt", sep = "\t", stringsAsFactors = F)
df3 <- df3[, c("run_accession", "sample_alias", "fastq_ftp")]
colnames(df3) <- c("run_accession", "library_name", "fastq_ftp")
df3 <- df3[-grep("T1|T2|T3|T4|T5|T6|POP", df3$library_name), ]
df3 <- df3[order(df3$library_name), ]
df3$study <- "Grumaz_2016"
df3$y <- "-1"
df3[grep("G", df3$library_name), "y"] <- "healthy"
df3[grep("S", df3$library_name), "y"] <- "septic"

# Combine
final <- rbind(df, df2, df3)
final <- final[order(final$library_name), ]

table(final$y)

write.csv(final, 'grumaz_pooled_metadata.csv', row.names = F)
# write.table(final$fastq_ftp, 'grumaz_accessions.csv', row.names = F, quote = F, col.names = F)

df <- read.csv("karius_parsed_metadata.csv", stringsAsFactors = F)
write.table(df$run, "karius_accessions.csv", quote = F, row.names = F, col.names = F)
