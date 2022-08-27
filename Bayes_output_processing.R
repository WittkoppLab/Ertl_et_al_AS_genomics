#############
##Libraries##
#############
library(plyr)
library(cowplot)
library(magrittr)
library(ggplot2)


#########################
##Set master plot theme##
#########################
theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 35),
  axis.title = element_text(size = 40),
  strip.text = element_text(size = 40),
  legend.text= element_text(size = 35),
  legend.title = element_text(size = 40),
  plot.title = element_text(size = 40, face = "bold")
)
}

# set wd
setwd("")

#####################
##Read primary data##
#####################
full_dataset <- read.delim("./CnR/CnR_grh_great20_1mismap_CPM.txt", header = T)
colnames(full_dataset) <- c("chrom", "start", "end", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3", "HYB_1_P2", "HYB_1_P1", "HYB_2_P2", "HYB_2_P1", "HYB_3_P2", "HYB_3_P1")
full_dataset$Paste_locus <- paste(full_dataset$chrom, full_dataset$start, full_dataset$end, sep = "_")

Parental_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "P1_1", "P1_2", "P1_3", "P2_1", "P2_2", "P2_3")]

Hybrid_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "HYB_1_P1", "HYB_2_P1", "HYB_3_P1", "HYB_1_P2", "HYB_2_P2", "HYB_3_P2")]


####################
##Combine datasets##
####################
Parental_results <- read.table("./CnR/CnR_grh_great20_1mismap_CPM_bayes_parents.txt", header = T) %>% unique()
Hybrid_results <- read.table("./CnR/CnR_grh_great20_1mismap_CPM_bayes_hybrids.txt", header = T) %>% unique()

Full_results_output <- join_all(list(Parental_data, Hybrid_data, Parental_results, Hybrid_results), by = 'Paste_locus', type = 'full')

Full_results_output <- na.omit(Full_results_output) %>% unique


##Apply FDR correction
##First plot distribution of p-values
Parent_P_plot <- ggplot(Full_results_output, aes(x = P_p_value)) + geom_histogram(bins = 100) + ggtitle("Parents")
Hybrid_P_plot <- ggplot(Full_results_output, aes(x = H_p_value)) + geom_histogram(bins = 100) + ggtitle("Hybrids")
p_vals <- plot_grid(Parent_P_plot, Hybrid_P_plot, nrow = 2)
##run FDR correction
Full_results_output$P_qvalue <- p.adjust(Full_results_output$P_p_value, method = "BY")
Full_results_output$H_qvalue <- p.adjust(Full_results_output$H_p_value, method = "BY")

####################################
##Get consistent allele directions##
####################################
Full_results_output$Direction_parent <- NA

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$P1_1[i] > Full_results_output$P2_1[i] & Full_results_output$P1_2[i] > Full_results_output$P2_2[i] & Full_results_output$P1_3[i] > Full_results_output$P2_3[i]){

	Full_results_output$Direction_parent[i] <- "P1"

} else if (Full_results_output$P1_1[i] < Full_results_output$P2_1[i] & Full_results_output$P1_2[i] < Full_results_output$P2_2[i] & Full_results_output$P1_3[i] < Full_results_output$P2_3[i]){

	Full_results_output$Direction_parent[i] <- "P2"

} else {Full_results_output$Direction_parent[i] <- "Ambig"}
}


Full_results_output$Direction_hybrid <- NA

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$HYB_1_P1[i] > Full_results_output$HYB_1_P2[i] & Full_results_output$HYB_2_P1[i] > Full_results_output$HYB_2_P2[i] & Full_results_output$HYB_3_P1[i] > Full_results_output$HYB_3_P2[i]){

	Full_results_output$Direction_hybrid[i] <- "P1"

} else if (Full_results_output$HYB_1_P1[i] < Full_results_output$HYB_1_P2[i] & Full_results_output$HYB_2_P1[i] < Full_results_output$HYB_2_P2[i] & Full_results_output$HYB_3_P1[i] < Full_results_output$HYB_3_P2[i]){

	Full_results_output$Direction_hybrid[i] <- "P2"

} else {Full_results_output$Direction_hybrid[i] <- "Ambig"}
}

##########################################################################
### Get Parental / Hybrid ratio to classify opposing vs same cis+trans ###
##########################################################################

Full_results_output$P_H_ratio <- abs(Full_results_output$P_est.mean / Full_results_output$H_est.mean)

####################################################
##Run classifier to establish class of each region##
####################################################
##Set qvalue cut-off
critical_value <- 0.05

##Run classifier
Full_results_output$Regulatory_class <- "Ambiguous"

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$P_qvalue[i] > critical_value & Full_results_output$H_qvalue[i] > critical_value & Full_results_output$P_H_qvalue[i] > critical_value){

	Full_results_output$Regulatory_class[i] <- "Conserved"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] > critical_value){

	Full_results_output$Regulatory_class[i] <- "Cis"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] > critical_value & Full_results_output$P_H_qvalue[i] < critical_value){

	Full_results_output$Regulatory_class[i] <- "Trans"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] == Full_results_output$Direction_hybrid[i] & Full_results_output$P_H_ratio[i] > 1 & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

	Full_results_output$Regulatory_class[i] <- "Cis_+_Trans,opposing"


} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] == Full_results_output$Direction_hybrid[i] & Full_results_output$P_H_ratio[i] < 1 & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

  Full_results_output$Regulatory_class[i] <- "Cis_+_Trans,same"

} else if (Full_results_output$P_qvalue[i] < critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value & Full_results_output$Direction_parent[i] != Full_results_output$Direction_hybrid[i] & Full_results_output$Direction_hybrid[i] != "Ambig" & Full_results_output$Direction_parent[i] != "Ambig"){

	Full_results_output$Regulatory_class[i] <- "Cis_*_Trans"

} else if (Full_results_output$P_qvalue[i] > critical_value & Full_results_output$H_qvalue[i] < critical_value & Full_results_output$P_H_qvalue[i] < critical_value){

	Full_results_output$Regulatory_class[i] <- "Compensatory"
}
}

##Run classifier for opposing and same
Full_results_output$Direction <- "NA"

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,opposing"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Cis_*_Trans"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Compensatory"){

	Full_results_output$Direction[i] <- "Opposing"

} else if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,same"){

	Full_results_output$Direction[i] <- "Reinforcing"
}
}

##Run classifier for sampler categorization
Full_results_output$Simple_class <- "Ambiguous"

for (i in 1:nrow(Full_results_output)) {

if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,opposing"){

	Full_results_output$Simple_class[i] <- "Cis-trans"

} else if (Full_results_output$Regulatory_class[i] == "Cis_*_Trans"){

	Full_results_output$Simple_class[i] <- "Cis-trans"

} else if (Full_results_output$Regulatory_class[i] == "Compensatory"){

	Full_results_output$Simple_class[i] <- "Cis-trans"

} else if (Full_results_output$Regulatory_class[i] == "Cis_+_Trans,same"){

	Full_results_output$Simple_class[i] <- "Cis-trans"

} else if (Full_results_output$Regulatory_class[i] == "Cis"){

	Full_results_output$Simple_class[i] <- "Cis"

} else if (Full_results_output$Regulatory_class[i] == "Trans"){

	Full_results_output$Simple_class[i] <- "Trans"

} else if (Full_results_output$Regulatory_class[i] == "Conserved"){

  Full_results_output$Simple_class[i] <- "Conserved"
}
}


#### Compute % CIS and TRANS ####
Full_results_output$trans_reg_diff <- Full_results_output$P_est.mean - Full_results_output$H_est.mean
Full_results_output$perc_cis <- (abs(Full_results_output$H_est.mean)/(abs(Full_results_output$H_est.mean) + abs(Full_results_output$trans_reg_diff))) * 100


##################################
##Write out full results to file##
##################################

write.table(Full_results_output, file = "/Users/henryertl/repos/MISC/Full_results_output_Cnr.txt", sep = "\t", row.names = F, quote = F)
