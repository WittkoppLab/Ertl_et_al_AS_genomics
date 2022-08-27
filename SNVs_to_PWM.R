library(dplyr)
library(plyr)

# read in df
df <- read.table("./z30_zhr_variants_longer.vcf_SNPsonly_30minQ.gz.recode_dm6_Grh_motif_p4_overlap.bed", header = F, sep = "\t")

# split up positive and negative strands
df_pos <- df[df[,4] == "+",]
df_neg <- df[df[,4] == "-",]

df_neg$complement_ref <- NA
for (i in 1:nrow(df_neg)) {
if (df_neg[i,9] == "G"){
  df_neg$complement_ref[i] <- "C"
} else if (df_neg[i,9] == "C"){
  df_neg$complement_ref[i] <- "G"
} else if (df_neg[i,9] == "T"){
  df_neg$complement_ref[i] <- "A"
} else if (df_neg[i,9] == "A"){
  df_neg$complement_ref[i] <- "T"
}
}

df_neg$complement_alt <- NA
for (i in 1:nrow(df_neg)) {
if (df_neg[i,10] == "G"){
  df_neg$complement_alt[i] <- "C"
} else if (df_neg[i,10] == "C"){
  df_neg$complement_alt[i] <- "G"
} else if (df_neg[i,10] == "T"){
  df_neg$complement_alt[i] <- "A"
} else if (df_neg[i,10] == "A"){
  df_neg$complement_alt[i] <- "T"
}
}

df_neg <- df_neg %>% na.omit()
df_neg$V9 <- NULL
df_neg$V10 <- NULL

df_pos <- df_pos[df_pos$V9 != "." & df_pos$V10 != ".",]

colnames(df_pos)[9:10] <- c("V9", "V10")
colnames(df_neg)[9:10] <- c("V9", "V10")
df <- rbind(df_pos, df_neg)


# get relative snp location w/in motif
df$rel_snp_location <- df[,8] - df[,2]

#change all up uppercase
df[,5] <- toupper(df[,5])

df[,9] <- as.character(df[,9])
df[,10] <- as.character(df[,10])


# keep only rows with dinucleotide SNVs
df <- df[nchar(df[,9]) == 1 & nchar(df[,10]) == 1 & df[,9] != "*" & df[,10] != "*",]

# for loop to contruct reference seq
df$ref_seq <- as.character(df[,5])
df[,9] <- as.character(df[,9])
df[,11] <- as.character(df[,11])

for (i in 1:nrow(df)) {
substring(df$ref_seq[i], df[i,11]) <- df[i,9]
}

# for loop to contruct alt seq
df$alt_seq <- as.character(df[,5])
df[,10] <- as.character(df[,10])
df[,11] <- as.character(df[,11])

for (i in 1:nrow(df)) {
substring(df$alt_seq[i], df[i,11]) <- df[i,10]
}


# input function for PWM calculation
# useful reference for PWM code: https://davetang.org/muse/2013/10/01/position-weight-matrix/
pwm <- function(freq, total, bg=0.25){
  #using the formulae above
  p <- (freq + (sqrt(total) * 1/4)) / (total + (4 * (sqrt(total) * 1/4)))
  log2(p/bg)
}

#define the frequencies of nucleotides for Grh
A <- c(1170.00,1915.00,3090.00,3211.00,29.00,187.00,1914.00,20.00,49.00,58.00,358.00,546.00,754.00)
C <- c(730.00,477.00,27.00,20.00,3265.00,2825.00,194.00,57.00,34.00,96.00,569.00,786.00,933.00)
G <- c(883.00,535.00,89.00,57.00,18.00,199.00,838.00,3221.00,23.00,23.00,407.00,818.00,783.00)
T <- c(549.00,405.00,126.00,44.00,20.00,121.00,386.00,34.00,3226.00,3155.00,1998.00,1182.00,862.00)
m <- matrix(data=c(A,C,G,T),nrow=4,byrow=T,dimnames=list(c('A','C','G','T')))
mm <- pwm(m,sum(m[,1]))

#compute PWM score
df$ref_score <- NA
for (i in 1:nrow(df)){
seq <- df[i,12]
x <- strsplit(x=seq,split='')
#initialise vector
seq_score <- vector()
#get the corresponding values
for (k in 1:13){
  seq_score[k] <- mm[x[[1]][k],k]
}
df$ref_score[i] <- sum(seq_score)
}

df$alt_score <- NA
for (i in 1:nrow(df)){
seq <- df[i,13]
x <- strsplit(x=seq,split='')
#initialise vector
seq_score <- vector()
#get the corresponding values
for (k in 1:13){
  seq_score[k] <- mm[x[[1]][k],k]
}
df$alt_score[i] <- sum(seq_score)
}

df$pwm_delta <- df$ref_score - df$alt_score
a <- ggplot(df, aes(x=pwm_delta))+
geom_histogram(bins=30)

colnames(df) <- c("chrom", "start", "end", "strand", "ref", "chr_SNV", "start_CNV", "end_CNV", "ZHR_SNV", "Z30_SNV", "rel_snp_location", "ZHR_seq", "Z30_seq", "ZHR_score", "Z30_score", "pwm_delta")
