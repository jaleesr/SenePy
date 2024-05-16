#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)



library(org.Hs.eg.db)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")


#print('i working')

#fpath <- args[1]
#opath <- args[2]


f <- file("stdin")
open(f)



# print('reading...')
# print(fpath)
# print('saving to...')
# print(opath)

genes <- readLines(f)



mm <- org.Hs.eg.db

lookup <- select(mm, 
                 keys = genes,
                 columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")


#write.csv(lookup, stdout(), row.names = FALSE, quote=FALSE)



df <- as.data.frame(lookup)

#print(df)

entrez_ids <- as.character(df$ENTREZID)
entrez_sym<- as.character(df$SYMBOL)


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
cols <- c("TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
temp <- select(txdb, keys=entrez_ids, columns=cols, keytype="GENEID")
temp_pos <- temp[temp$TXSTRAND == '+',]
temp_neg <- temp[temp$TXSTRAND == '-',]


library(plyr)
library(dplyr)



temp_pos <- ddply(temp_pos, .(GENEID), function(df){
  df$TXSTART <- as.numeric(df$TXSTART)
  df.arrange<-df %>% arrange(TXSTART) # arrange(desc(TXstart))
  val<-df.arrange[1,]
  return(val)
})

temp_neg <- ddply(temp_neg, .(GENEID), function(df){
  df$TXSTART <- as.numeric(df$TXSTART)
  df.arrange<-df %>% arrange(desc(TXSTART)) # arrange(desc(TXstart))
  val<-df.arrange[1,]
  return(val)
})


temp_pos$slice_start <- temp_pos$TXSTART - 1000
temp_pos$slice_end <- temp_pos$TXSTART + 500

temp_neg$slice_start <- temp_neg$TXEND - 500
temp_neg$slice_end <- temp_neg$TXEND  + 1000

master_df <- rbind(temp_pos, temp_neg)

master_df$symbol <- mapvalues(master_df$GENEID, from=entrez_ids, to=entrez_sym)



write.csv(master_df, stdout(), row.names = FALSE, quote=FALSE)














