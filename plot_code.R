#### MAIN FIGURES ####
# FIG 1C
snps <- read.delim("./z30_zhr_variants_longer.vcf_SNPsonly_30minQ.gz.recode_dm6_Grh_motif_p4_overlap_PWMs.txt", header = T)
snps$rel_snp_location <- factor(snps$rel_snp_location,levels = c("1","2","3","4","5","6","7","8","9","10","11","12"))
snv_fig <- snps[snps$ZHR_score > 8,] %>% ggplot(aes(x=rel_snp_location)) + geom_bar() +
ylab("ZHR - Z30 # of SNVs") +
xlab("Grh motif position") +
theme_main()
ggsave(snv_fig, file = "./snv_fig.pdf", width = 15, height = 15)

# FIG 1D
#########################################
### command line w/ deepTools in PATH ###
sed 1d Full_results_output_CnR.txt > Full_results_output_CnR.bed
sed 1d Full_results_output_ATAC.txt > Full_results_output_ATAC.bed
sed 1d Full_results_output_RNA.txt > Full_results_output_RNA.bed

computeMatrix scale-regions -S ZHR_3_Z30_1_CnR.bw -R Full_results_output_CnR.bed -a 500 -b 500 -o Full_results_output_CnR_1000_mat
computeMatrix scale-regions -S ZHR_3_Z30_1_ATAC.bw -R Full_results_output_ATAC.bed -a 2000 -b 2000 -o Full_results_output_ATAC_4000_mat
computeMatrix scale-regions -S ZHR_3_Z30_1_RNA.bw -R Full_results_output_RNA.bed -a 100 -b 1500 -o Full_results_output_RNA_TSS_mat

plotHeatmap -m Full_results_output_CnR_1000_mat --colorMap inferno --whatToShow heatmap -o Full_results_output_CnR_1000.png
plotHeatmap -m Full_results_output_ATAC_4000_mat --colorMap inferno --whatToShow heatmap -o Full_results_output_ATAC_4000.png
plotHeatmap -m Full_results_output_RNA_TSS_mat --colorMap inferno --whatToShow heatmap -o Full_results_output_RNA_TSS.png
########################################

# FIG 2A
cnr <- read.table("./Full_results_output_Cnr.txt", header = T) %>% unique()
hist <- function(df,xx,yy){
  ggplot(df, aes(x=xx, fill=yy < 0.05))+
  geom_histogram()+
  theme_main()
  }
cnr_hist <- hist(cnr,P_est.mean,P_qvalue)
ggsave(cnr_hist, file = "./hist_CnR_P.pdf", width = 15, height = 15)

# FIG 2B
cnr_hist <- hist(cnr,H_est.mean,H_qvalue)
ggsave(cnr_hist, file = "./hist_CnR_H.pdf", width = 15, height = 15)

# FIG 2C
##function for iterative plot
cis_trans <- function(df){
  ggplot(df,aes(x=xx, y=yy))+
  theme_main()+
  geom_point(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + xlim(-1.5, 1.5) + ylim(-1.5, 1.5) +
  }
###plot
cnr_cistrans <- cis_trans(cnr,P_est.mean,H_est.mean)
ggsave(cnr_cistrans, file = "./cistrans_CnR_all.pdf", width = 15, height = 15)

# FIG 2F
all_snps <- read.delim("./z30_zhr_variants_longer.vcf_SNPsonly_30minQ.gz.recode_dm6.bed", header = F)
full_res_all_snps <- bedtoolsr::bt.intersect(a=cnr, b=all_snps, wa = TRUE, wb = TRUE) %>% unique()
full_res_all_snps <- full_res_all_snps[,4] %>% as.data.frame()
full_res_all_snps <- ddply(full_res_all_snps,.(.),nrow)
colnames(full_res_all_snps) <- c("Paste_locus", "num_snvs")
Full_results_output <- left_join(cnr, full_res_all_snps, by = "Paste_locus") %>% unique()
Full_results_output$snv_freq <- Full_results_output$num_snvs/(Full_results_output$end - Full_results_output$start)
snv_dens <- Full_results_output %>%
ggplot(aes(x=H_qvalue < 0.05, y=snv_freq*500, fill=H_qvalue < 0.05)) +
geom_boxplot(notch=T)+
theme_main() +
xlab("Grh binding") +
ylab("SNV density (per 500bp)") +
guides(fill = FALSE, size = FALSE)
ggsave(c, file = "./snv_density.pdf", width = 5, height = 8)

# FIG 3A
atac <- read.table("./Full_results_output_ATAC.txt", header = T) %>% unique()
cnr_atac <- bedtoolsr::bt.intersect(a=cnr, b=atac, wa=TRUE, wb=TRUE) %>% as.data.frame()
atac_hist <- hist(atac,H_est.mean,H_qvalue)
ggsave(atac_hist, file = "./hist_atac_all.pdf", width = 15, height = 15)

# FIG 3C
atac_cis_trans_grh <- cis_trans(cnr_atac,V54,V56)
ggsave(atac_cis_trans_grh), file = "./cistrans_atac_grh.pdf", width = 15, height = 15)

# FIG 3D
corr_plot <- function(df,xx,yy){
  ggplot(df,aes(x=xx, y=yy))+
  theme_main()+
  geom_point(size = 4)+
  geom_smooth(method="lm")
}
cnr_vs_atac <- corr_plot(cnr_atac,V21,V56)
ggsave(cnr_vs_atac), file = "./cnr_vs_atac.pdf", width = 15, height = 15)

# FIG 3E
#########################################
### command line w/ deepTools in PATH ###


# FIG 4A
rna_hist <- hist(rna)
ggsave(rna_hist, file = "./hist_rna_all.pdf", width = 15, height = 15)

# FIG 4B
rna_cistrans <- cis_trans(rna,P_est.mean,H_est.mean)
ggsave(rna_cistrans, file = "./cistrans_rna_all.pdf", width = 15, height = 15)

# FIG 4C
cnr_atac_rna <- bedtoolsr::bt.closest(a=cnr_atac, b=rna, wa=TRUE, wb=TRUE) %>% as.data.frame()
rna_cis_trans_grh <- cis_trans(cnr_atac_rna,V88,V90)
ggsave(rna_cis_trans_grh), file = "./cistrans_rna_grh.pdf", width = 15, height = 15)

# FIG 4D
atac_vs_rna <- corr_plot(cnr_atac_rna,V56,V90)
ggsave(atac_vs_rna), file = "./atac_vs_rna.pdf", width = 15, height = 15)

# FIG 5A
cnr_atac_rna_Hvals <- cbind(cnr_atac_rna$V4, cnr_atac_rna$V26,  cnr_atac_rna$V61, cnr_atac_rna$V95, cnr_atac_rna$V21,  cnr_atac_rna$V56, cnr_atac_rna$V90) %>% as.data.frame()
colnames(cnr_atac_rna_Hvals) <- c("locus", "cnrQ", "atacQ", "rnaQ", "cnr","atac","rna")
#shuffle rows of each column of q-vals: 1000 permutations
critical_value <- 0.05
mat <- matrix(ncol = 8, nrow = 1000) %>% as.data.frame()
for (j in 1:1000) {
	random <- NULL
	random <- cbind(cnr_atac_rna_Hvals$locus, sample(cnr_atac_rna_Hvals$cnrQ),  sample(cnr_atac_rna_Hvals$atacQ), sample(cnr_atac_rna_Hvals$rnaQ)) %>% as.data.frame()
	mat[j,1] <- nrow(random[random$V2 < critical_value & random$V3 < critical_value & random$V4 < critical_value,])
	mat[j,2] <- nrow(random[random$V2 < critical_value & random$V3 < critical_value & random$V4 > critical_value,])
	mat[j,3] <- nrow(random[random$V2 < critical_value & random$V3 > critical_value & random$V4 > critical_value,])
	mat[j,4] <- nrow(random[random$V2 > critical_value & random$V3 > critical_value & random$V4 > critical_value,])
	mat[j,5] <- nrow(random[random$V2 > critical_value & random$V3 < critical_value & random$V4 < critical_value,])
	mat[j,6] <- nrow(random[random$V2 > critical_value & random$V3 > critical_value & random$V4 < critical_value,])
	mat[j,7] <- nrow(random[random$V2 > critical_value & random$V3 < critical_value & random$V4 > critical_value,])
	mat[j,8] <- nrow(random[random$V2 < critical_value & random$V3 > critical_value & random$V4 < critical_value,])
}

# create full df with real and permuted values
mat <- rbind(mat, c(33,73,121,251,22,72,72,37)) #real counts of each category
temp <- cbind(mat, c(rep("permut",times=1000),"real"))
colnames(temp)[9] <- "class"

temp_melt <- (melt(temp))
temp_melt_real <- temp_melt[temp_melt$class == "real",]
temp_melt_perm <- temp_melt[temp_melt$class == "permut",]

temp_melt_perm$variable <- factor(temp_melt_perm$variable, levels=c("V1", "V2", "V8", "V3", "V5", "V7", "V6", "V4"))
temp_melt_real$variable <- factor(temp_melt_real$variable, levels=c("V1", "V2", "V8", "V3", "V5", "V7", "V6", "V4"))

permut_real <- ggplot()+
geom_violin(data=temp_melt_perm, aes(x=variable,y=value),color="black", fill="grey")+
geom_point(data=temp_melt_real, aes(x=variable,y=value), color="red", size = 5)+
theme_main()
ggsave(permut_real, file = "./all_permut_violin.pdf", width = 16, height = 5)

# FIG 5B
cnr_atac_rna_Hvals <- cnr_atac_rna_Hvals[,c(1,5:8)]
ab <- melt(cnr_atac_rna_Hvals, id.vars = c("locus", "class")) %>%
ggplot(aes(x=variable, y=abs(as.numeric(value)), group=factor(locus))) +
geom_smooth(aes(group = 1), method="loess", formula=y~x) +
  geom_line(alpha=0.08) + geom_point(alpha=0.08) +
  theme_main() +
  facet_wrap(~class) +
  theme(legend.position="none")+
  ylim(-0,1.5)
ggsave(ab, file = "./spagetti_abs.pdf", width = 15, height = 15)

######### SUPP FIGS ##########
# FIG S1A
df1 <- read.delim("./z30_zhr_variants_longer.vcf_SNPsonly_30minQ.gz.recode_dm6.bed", header = F) %>% unique()
df2 <- read.delim("./Full_results_output_RNA.bed", header = F) %>% unique()
df2$length <- df2$V3 - df2$V2
df3 <- read.delim("./Full_results_output_ATAC.bed", header = F) %>% unique()
df3$length <- df3$V3 - df3$V2
df4 <- read.delim("./Full_results_output_CnR.bed", header = F) %>% unique()
df4$length <- df4$V3 - df4$V2

df12 <- bedtoolsr::bt.intersect(a=df1, b=df2, wa = TRUE,wb=TRUE) %>% unique()
df12 <- df12[,c(4:6,38)]
df12$locus <- paste(df12$V4, df12$V5, df12$V6, sep="_")
df12 <- ddply(df12,.(locus,V38),nrow)
df12$class <- "rna"

df13 <- bedtoolsr::bt.intersect(a=df1, b=df3, wa = TRUE,wb=TRUE) %>% unique()
df13 <- df13[,c(4:6,38)]
df13$locus <- paste(df13$V4, df13$V5, df13$V6, sep="_")
df13 <- ddply(df13,.(locus,V38),nrow)
df13$class <- "atac"

df14 <- bedtoolsr::bt.intersect(a=df1, b=df4, wa = TRUE,wb=TRUE) %>% unique()
df14 <- df14[,c(4:6,39)]
df14$locus <- paste(df14$V4, df14$V5, df14$V6, sep="_")
df14 <- ddply(df14,.(locus,V39),nrow)
colnames(df14)[2] <- "V38"
df14$class <- "cnr"
df5 <- rbind(df12, df13, df14)
i <- ggplot(df5, aes(x=class,y=(V1/V38)*100))+
geom_boxplot(notch=T)+
geom_hline(yintercept=1.2,color="red")+
theme_main()+
xlab("")+
ylab("SNV frequency in region/gene (percent)")
ggsave(i, file = "./snv_freq_overall.pdf", width = 8, height = 11)

# FIG S1B
df_long <- read.table("./z30_zhr_variants_longer.vcf_SNPsonly_30minQ.gz.recode_dm6_Grh_motif_p4_overlap_plusminus20.bed", header = F, sep = "\t")
df_long[,c(9:10)] <- NULL
df_long$pos_rel_center <- (df_long[,8] - (df_long[,4]+6))
snp_fig2 <- df_long%>% ggplot(aes(x=pos_rel_center)) +
geom_line(aes(fill=..count..),stat="bin",binwidth=1, size = 2) +
ylab("ZHR - Z30 # of SNVs") +
xlab("Position relative to center of Grh motif") +
theme_main()
ggsave(snp_fig2, file = "./figs_rough/snp_fig_long.pdf", width = 15, height = 15)

# FIG S3
#functions
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

# read in data and plot
df <- read.table("./Full_results_output_Cnr.txt", header = T) %>% unique()
cormat <- round(cor(df[,c(4:ncol(df))]),2)
upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
l <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
  theme_minimal()+ # minimal theme
   theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
ggsave(l, file = "./cnr_corr_heatmap.pdf", width = 10, height = 10)

df <- read.table("./Full_results_output_ATAC.txt", header = T) %>% unique()
cormat <- round(cor(df[,c(4:ncol(df))]),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
m <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
  theme_minimal()+ # minimal theme
   theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
ggsave(m, file = "./atac_corr_heatmap.pdf", width = 10, height = 10)

df <- read.table("./Full_results_output_RNA.txt", header = T) %>% unique()cormat <- round(cor(df[,c(4:ncol(df))]),2)
df <- df[,c(1:7,9,8,11,10,13,12)]
cormat <- round(cor(df[,c(2:ncol(df))]),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
n <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
  theme_minimal()+ # minimal theme
   theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
ggsave(n, file = "./rna_corr_heatmap.pdf", width = 10, height = 10)

# FIG S4
atac_cistrans <- cis_trans(atac,P_est.mean,H_est.mean)
ggsave(atac_cistrans, file = "./cistrans_atac_all.pdf", width = 15, height = 15)

# FIG S5
rna <- read.table("./Full_results_output_RNA.txt", header = T) %>% unique()
atac_rna <- bedtoolsr::bt.closest(a=atac, b=rna, wa=TRUE, wb=TRUE) %>% as.data.frame()
atac_vs_rna_all <- corr_plot(atac_rna,V21,V56)
ggsave(atac_vs_rna), file = "./atac_vs_rna.pdf", width = 15, height = 15)

# FIG S6
cnr_atac_rna_Hvals <- cnr_atac_rna_Hvals[,c(1,5:8)]
no_ab <- melt(cnr_atac_rna_Hvals, id.vars = c("locus", "class")) %>%
ggplot(aes(x=variable, y=as.numeric(value), group=factor(locus))) +
  geom_line(alpha=0.08) + geom_point(alpha=0.08) +
  theme_main() +
  facet_wrap(~class) +
  theme(legend.position="none")+
  ylim(-1.5,1.5)
ggsave(no_ab, file = "./spagetti_noabs.pdf", width = 15, height = 15)

ab <- melt(cnr_atac_rna_Hvals, id.vars = c("locus", "class")) %>%
ggplot(aes(x=variable, y=abs(as.numeric(value)), group=factor(locus))) +
geom_smooth(aes(group = 1), method="loess", formula=y~x) +
  geom_line(alpha=0.08) + geom_point(alpha=0.08) +
  theme_main() +
  facet_wrap(~class) +
  theme(legend.position="none")+
  ylim(-0,1.5)
ggsave(ab, file = "./spagetti_abs.pdf", width = 15, height = 15)
