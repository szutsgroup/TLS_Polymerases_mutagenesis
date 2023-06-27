library(ASCAT)
library(dplyr)

centro <- read.delim("~/Downloads/GRCh38_centro.bed")
centro$X.bin <- centro$name <- NULL
centro$chrom <- gsub("chr", "", centro$chrom)
group_by(centro, chrom) %>% summarize(Start = min(chromStart)-1e6, End = max(chromEnd)+1e6) ->centro

blacklist <- read.delim("~/Downloads/hg38-blacklist.v2.bed", header = FALSE)
names(blacklist) <- c("chrom", "Start", "End", "reason")
blacklist$chrom <- gsub("chr", "", blacklist$chrom)

ff <- list.files(path = "./", pattern = "*_getaf.txt", full.names = TRUE)
ff <- ff[1:2]

w <- read.delim(ff[1], sep = " ")
w$Chrom <- factor(w$Chrom, levels = c(1:22, "X", "Y"))
arrange(w, Chrom, Pos) -> w
w$Centro <- FALSE
w$Blacklist <- FALSE
for (i in levels(w$Chrom)) {
  w[which(w$Chrom == i & w$Pos > filter(centro, chrom == i)$Start & w$Pos <= filter(centro, chrom == i)$End), "Centro"] <- TRUE
}
gr0 <- GRanges(seqnames = w$Chrom, IRanges(start = w$Pos, width = 1))
olaps <- findOverlaps(query = gr0, subject = GRanges(blacklist))
w[queryHits(olaps), "Blacklist"] <- TRUE

baf <- logr <- data.frame(chrs = w$Chrom, pos = w$Pos, Centro = w$Centro, Blacklist = w$Blacklist) %>% filter(!Centro, !Blacklist) %>% dplyr::select(-Centro, -Blacklist)

# collect data in individual files into one big dataframe
for (f in ff) {
  fn = strsplit(f, split = "/") %>% sapply("[", 9) %>% gsub("_afs.txt", "", .)
  w <- read.delim(f, sep = " ")
  w$BAF = w[,5] / w[,6]
  w$Cov = log2(w[,6] / mean(w[,6], na.rm = TRUE, trim = 0.001))
  w$Chrom <- factor(w$Chrom, levels =c(1:22, "X", "Y"))
  arrange(w, Chrom, Pos) -> w
  w$Centro <- FALSE
  for (i in levels(w$Chrom)) {
    w[which(w$Chrom == i & w$Pos > filter(centro, chrom == i)$Start & w$Pos <= filter(centro, chrom == i)$End), "Centro"] <- TRUE
  }
  filter(w, !Centro) %>%
    dplyr::select(Chrom, Pos, Cov) %>%
    dplyr::rename("chrs" = Chrom, "pos" = Pos, !!fn := Cov) %>%
    merge(logr, ., by = c("chrs", "pos")) -> logr
  filter(w, !Centro) %>%
    dplyr::select(Chrom, Pos, BAF) %>%
    dplyr::rename("chrs" = Chrom, "pos" = Pos, !!fn := BAF) %>%
    merge(baf, ., by = c("chrs", "pos")) -> baf
  print(fn)
}
logr <- arrange(logr, chrs, pos)
baf <- arrange(baf, chrs, pos)

# calculate logR values
logr2 <- logr



# set up input & output
ascat_list = vector(mode = "list", length = length(ff))
names(ascat_list) <- strsplit(ff, split = "/") %>% sapply("[", 9) %>% gsub("_afs.txt", "", .)
bc_list = vector(mode = "list", length = length(ff))
names(bc_list) <- strsplit(ff, split = "/") %>% sapply("[", 9) %>% gsub("_afs.txt", "", .)

done <- list.files("./", pattern = "*ASCATprofile.png")
done <- substr(done, 1, 5)

# generate files for ASCAT
# run ASCAT, starting is RS009 if cell was created from P53, RS010 if from P53 BRCA1
for (f in ff) {
  
  
  fn = strsplit(f, split = "/") %>% sapply("[", 9) %>% gsub("_afs.txt", "", .)
 
  
  
  dplyr::select(baf, chrs, pos, all_of(fn)) %>%
    write.table("sample_BAF.txt", sep = "\t", quote = FALSE)
  dplyr::select(logr2, chrs, pos, all_of(fn)) %>%
    write.table("sample_logR.txt", sep = "\t", quote = FALSE)
  
  
  
    dplyr::select(baf, chrs, pos) %>%
      mutate(imaginary_normal = ifelse(baf[,fn] > 0.9, 1, 0.5)) %>%
      write.table("germline_BAF.txt", sep = "\t", quote = FALSE)
    dplyr::select(logr2, chrs, pos) %>%
      mutate(imaginary_normal = 0) %>%
      write.table("germline_logR.txt", sep = "\t", quote = FALSE)

  
  
  ascat.bc = ascat.loadData(Tumor_LogR_file = "sample_logR.txt",
                            Tumor_BAF_file =  "sample_BAF.txt",
                            Germline_LogR_file = "germline_logR.txt", 
                            Germline_BAF_file = "germline_BAF.txt")
  ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_custom_loci.txt", replictimingfile = "RT_custom_loci.txt")
  ascat.plotRawData(ascat.bc, img.dir = "plots/")
  ascat.bc = ascat.aspcf(ascat.bc)
  bc_list[[fn]] <- ascat.bc
  ascat.plotSegmentedData(ascat.bc, img.dir = "plots/")
  ascat_list[[fn]] = ascat.runAscat(ascat.bc, img.dir = "plots/", gamma = 1)
  for (i in list.files("./", path = "*txt", full.names = TRUE)) {
    file.remove(i)
  }
  
}


  





