library(tidyverse)
library(rpart)
library(BSgenome.Hsapiens.UCSC.hg38)

# A

id <- readxl::read_excel("~/Public/Szuts_csoport_megosztott/ALL_WGS_SAMPLES.xlsx", sheet = "WGS samples")

centro <- read.delim("~/Downloads/GRCh38_centro.bed")
centro$X.bin <- centro$name <- NULL
centro$chrom <- gsub("chr", "", centro$chrom)
group_by(centro, chrom) %>% summarize(Start = min(chromStart)-1e6, End = max(chromEnd)+1e6) ->centro
blacklist <- read.delim("~/Downloads/hg38-blacklist.v2.bed", header = FALSE)
names(blacklist) <- c("chrom", "Start", "End", "reason")
blacklist$chrom <- sub("chr", "", blacklist$chrom)

seqNames <- seqnames(Hsapiens)[1:24]
seqLengths <-  as.numeric(seqlengths(Hsapiens)[1:24])
cumulativeLengths <- Reduce(sum, seqLengths, accumulate = TRUE)
genomeLength <- sum(seqLengths)
normalizeGenomicPositions <- function(chr, pos) {
  pos <- as.numeric(pos)
  chr <- paste0("chr", chr)
  whichChrom <- which(seqNames %in% chr)
  earlierChromsLength <- sum(seqLengths[0:(whichChrom - 1)])
  return(as.numeric(earlierChromsLength) + pos)
}

get_info <- function(sample) {
  filter(id, `Sample ID` == sample) -> tmp
  return(c(tmp$Genotype, tmp$Treatment))
}

ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

coverage_difference <- function(sample, control) {
  a <- read.table(gzfile(paste0("~/Documents/202205_cnv/simple_cov/beds/", sample, "_bedcov_10k.bed.gz")))
  b <- read.table(gzfile(paste0("~/Documents/202205_cnv/simple_cov/beds/", control, "_bedcov_10k.bed.gz")))
  names(a) <- names(b) <- c("Chr", "Start", "End", "Cov")
  filter(a, Chr %in% c(1:22, "X", "Y")) -> a
  filter(b, Chr %in% c(1:22, "X", "Y")) -> b
  merge(a, b, by = c("Chr", "Start", "End")) -> oo
  names(oo)[4:5] <- c(sample, control)
  oo$Chr <- factor(oo$Chr, levels = c(1:22, "X", "Y"))
  arrange(oo, Chr, Start) -> oo
  oo[,4] <- ma(oo[,4]/median(oo[,4]), 10)
  oo[,5] <- ma(oo[,5]/median(oo[,5]), 10)
  oo[,4] <- as.numeric(oo[,4])
  oo[,5] <- as.numeric(oo[,5])
  
  return(oo)

}

get_sc <- function(sample, cell_type) {
  filter(id, `Sample ID` == sample) -> tmp
  filter(id, Genotype == tmp$Genotype, Treatment == "starting clone", Round != "R1", Type == cell_type) %>% 
    head(1) %>% .[,"Sample ID"] %>% unlist() %>% unname() %>% return()
}

gc <- read.delim("~/Documents/202205_cnv/simple_cov/GRCh38pa_EBV_10k_windows.gc")
names(gc) <- c("Chr", "Start", "End", "AT", "GC", "A", "C", "G", "T", "N", "other", "seq_len")
filter(gc, Chr %in% c(1:28, "X", "Y")) -> gc
gc2 <- gc[-which(gc$seq_len != "10000"),]
gc2 <- gc2[-which(gc2$N / gc2$seq_len > 0.05),]
gc2$Chr <- factor(gc2$Chr, levels = gsub("chr", "", seqNames))

dont_include <- c("RS009", "RS057", "RS018", paste0("RS0", c(49:52, 70:73)))

for (g in c("P53", "P53 REV1", "P53 PRIMPOL", "P53 REV1 PRIMPOL",
            "P53 K164R", "P53 K164R PRIMPOL", "P53 REV3L", "P53 K164R REV1")) {
  png(paste0("_", gsub(" ", "_", g), ".png"), width = 480, height = 480)
 
  filter(id, Type == "RPE-1", Genotype == g, 
         Treatment == "mock", !(`Sample ID` %in% dont_include))$`Sample ID` -> m
  filter(id, Type == "RPE-1", Genotype == g, 
         Treatment == "starting clone", !(`Sample ID` %in% dont_include))$`Sample ID` -> sc
  if (g == "P53") {
    sc <- "RS012"
    m <- setdiff(m, sc)
  }
  x <- read.delim(paste0("~/Documents/202205_cnv/simple_cov/SNP/", sc, "_afs.txt"), sep = " ")
  x$Chrom <- factor(x$Chrom, levels = c(1:22, "X", "Y"))
  x <- arrange(x, Chrom)
  x <- filter(x, Pos != "Pos")
  x[,5] <- as.numeric(x[,5])
  x[,6] <- as.numeric(x[,6])
  x$Pos <- as.numeric(x$Pos)
  x$Centro <- x$Blacklist <- FALSE
  for (i in levels(x$Chrom)) {
    x[which(x$Chrom == i & x$Pos > filter(centro, chrom == i)$Start & x$Pos <= filter(centro, chrom == i)$End), "Centro"] <- TRUE
  }
  gr0 <- GRanges(seqnames = x$Chrom, IRanges(start = x$Pos, width = 1))
  blo <- findOverlaps(query = gr0, subject = GRanges(blacklist))
  x[queryHits(blo), "Blacklist"] <- TRUE
  
  x$normPos <- apply(x, 1, function(x) normalizeGenomicPositions(x[1], x[2]))
  
  x2 <- filter(x, !Centro, !Blacklist)
  
  nf <- layout(matrix(c(1, 2, 3, 4, 5), 5, 1, byrow = TRUE),
               heights = c(0.75, 0.75, 1, 1, 1))
  par(mar = c(0, 4, 1, 0))
  kolor = scales::alpha("grey60", .05)
  plot(x2[,6] + rnorm(nrow(x2), mean = 0, sd = 1), x = x2$normPos, cex = .1, 
       ylim = c(0, 60), xlab = "", col = kolor, bg = kolor,
       ylab = "" ,pch = 21, xaxt = "n", yaxt = "n")
  axis(2, at = c(0, 20, 40, 60), cex.axis = 1.5, las = 2)
  for (i in names(table(x2$Chrom))) {
    if (i == "Y") next
    filter(x2, Chrom == i) -> tmp
    tmp <- tmp[which(tmp[,6] < median(tmp[,6]*2, na.rm = TRUE)),]
    tmp$MA <- ma(tmp[,6], 100)
    lines(tmp$MA, x = tmp$normPos, col = "purple", lwd = 4)
  }
  abline(v = cumulativeLengths, col = "red", lty = 3, lwd = 2, las = 2)
  
  par(mar = c(0, 4, 0.5, 0))
  plot(x2[,5] / x2[,6], x = x2$normPos, cex = .1, xaxt = "n", ylab = "", 
       col = kolor, bg = kolor, pch = 21, yaxt = "n")
  axis(2, at = c(0, .25, .5, .75, 1), cex.axis = 1.5, las = 2)
  abline(v = cumulativeLengths, col = "red", lty = 3, lwd = 2, las = 2)
  
  for (s in m) {
    oo <- coverage_difference(s, sc)
    oo <- oo[-which(gc$seq_len != "10000"),]
    oo <- oo[-which(gc$N / gc$seq_len > 0.05),]
    oo$normPos <- apply(oo, 1, function(x) normalizeGenomicPositions(x[1], x[2]))
    ss <- merge(oo, gc2) %>% arrange(Chr, Start)
    ss$Centro <- ss$Blacklist <- FALSE
    for (i in 1:nrow(centro)) {
      ss[which(ss$Chr == centro$chrom[i] & ss$Start >= centro$Start[i] & ss$End <= centro$End[i]), "Centro"] <- TRUE
    }
    gr0 <- GRanges(seqnames = ss$Chr, IRanges(start = ss$Start, end = ss$End))
    blo <- findOverlaps(query = gr0, subject = GRanges(blacklist))
    ss[queryHits(blo), "Blacklist"] <- TRUE
    
    
    ss <- ss[-union(which(ss$Centro), which(ss$Blacklist)),]
    ool <- findOverlaps(GRanges(ss), GRanges(x2$Chrom, IRanges(x2$Pos, width = 1)))
    ss <- ss[queryHits(ool),]
    lm1 <- lm(get(s) ~ get(sc) + GC, data = ss)
    
    par(mar = c(0, 4, 0, 0))
    plot(lm1$residuals, x = ss$normPos,  las = 2,
         cex = .1, ylim = c(-.7, .7), ylab = "", xaxt = "n", col = "grey75", yaxt = "n")
    axis(2, at = c(-.5, 0, .5), cex.axis = 1.5, las = 2)
    abline(v = cumulativeLengths, col = "red", lty = 3, lwd = 2, las = 2)
    tree <- rpart(lm1$residuals ~ ss$normPos, control = rpart.control(minsplit = 100, minbucket = 100, cp = 0.001))
    predict(tree, newdata = data.frame(ss$normPos)) -> n
    lines(n, x = ss$normPos, col = "black", lwd = 4)
  }
  
  
  
  
  dev.off()
}

# B

aa <- readxl::read_excel("ascat/ASCAT_summary.xlsx", sheet = "summary")
bb <- readxl::read_excel("ascat/ASCAT_summary.xlsx", sheet = "events")
mutate(bb, Sample = factor(Sample, levels = aa$Sample)) -> bb
merge(bb, aa[,1:3]) -> bb
bb$Chrom_percent <- bb$Width / seqlengths(Hsapiens)[paste0("chr", bb$Chr)]

group_by(bb, Sample, Type, .drop = FALSE) %>% 
  summarize(N = n(), .drop = FALSE) %>% 
  ggplot(aes(x = factor(Sample, levels = setdiff(aa$Sample, aa[grepl("ancestral", aa$Name),]$Sample)), 
             fill = Type, y = N, .drop = FALSE)) + 
  geom_col(col = "black", position = position_dodge()) + theme_bw() + 
  xlab("") + ylab("Width (bp)") + 
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(color = "black")) + 
  scale_x_discrete(drop = FALSE) + 
  scale_fill_discrete(drop = FALSE)
ggsave("CNV_counts.pdf", width = 15, height = 7)

# C
ggplot(bb, aes(x = Genotype, fill = Type, y = Width, group = Type)) + 
  geom_point(position = position_jitterdodge(jitter.width = .1, dodge.width = .5), 
             size = 5, color = "black", pch = 21) + theme_bw() + 
  scale_y_log10(breaks = c(1e5, 1e6, 1e7, 1e8), limits = c(1e5, 5e8)) + 
  xlab("") + ylab("Width (bp)") + 
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(color = "black"))
ggsave("CNV_lengths.pdf", width = 12, height = 7)
