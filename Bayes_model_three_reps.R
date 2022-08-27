##Libraries
library(data.table)
library(plyr)
library(dplyr)
library(INLA)
library(parallel)
library(qvalue)
library(magrittr)


##############################
##Get command line arguments##
##############################
Arguments = (commandArgs(TRUE))

#########################
##Define test functions##
#########################
##################
##Parental model##
##################
Parental_model <- function(locus_ID, Parental_data){

	chrom <- Parental_data[Parental_data$Paste_locus == locus_ID,]$chrom
	pos_start <- Parental_data[Parental_data$Paste_locus == locus_ID,]$start
	pos_end <- Parental_data[Parental_data$Paste_locus == locus_ID,]$end
	Paste_locus <- Parental_data$Paste_locus[Parental_data$Paste_locus == locus_ID]

	##Reformat data for modelling
	reformed_matrix <- matrix(ncol = 2, nrow = 3) %>% as.data.frame()
	reformed_matrix[1,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(5, 8)])
	reformed_matrix[2,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(6, 9)])
	reformed_matrix[3,] <- c(Parental_data[Parental_data$Paste_locus == locus_ID, c(7, 10)])

	names(reformed_matrix)[1] <- "p1_reads"
	names(reformed_matrix)[2] <- "p2_reads"
	reformed_matrix$Total_reads <- reformed_matrix$p1_reads + reformed_matrix$p2_reads

	P_mod <- inla(p1_reads ~ 1, data = reformed_matrix , family = "binomial", Ntrials = Total_reads)

	coef <- P_mod$summary.fixed
	fixed_effect_posterior <- P_mod$marginals.fixed[[1]]
	lower_p <- inla.pmarginal(0, fixed_effect_posterior)
	upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
	post_pred_p <- 2 * (min(lower_p, upper_p))


	P_mod_output <- data.table(chrom = chrom, pos_start = pos_start, pos_end = pos_end, Paste_locus = Paste_locus, P_est = coef[1], P_p_value = post_pred_p)

	return(P_mod_output)
}

################
##Hybrid model##
################
Hybrid_model <- function(locus_ID, Hybrid_data){

	##ADJUST WHEN WE HAVE FULL DATA
	chrom <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$chrom
	pos_start <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$start
	pos_end <- Hybrid_data[Hybrid_data$Paste_locus == locus_ID,]$end
	Paste_locus <- Hybrid_data$Paste_locus[Hybrid_data$Paste_locus == locus_ID]


	##Reformat data for modelling
	reformed_matrix <- matrix(ncol = 3, nrow = 3) %>% as.data.frame()
	reformed_matrix[1,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(5, 8)])
	reformed_matrix[2,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(6, 9)])
	reformed_matrix[3,] <- c(Hybrid_data[Hybrid_data$Paste_locus == locus_ID, c(7, 10)])

	names(reformed_matrix)[1] <- "p1_reads"
	names(reformed_matrix)[2] <- "p2_reads"
	reformed_matrix$Total_reads <-  reformed_matrix$p1_reads + reformed_matrix$p2_reads

	H_mod <- inla(p1_reads ~ 1, data = reformed_matrix , family = "binomial", Ntrials = Total_reads)

	coef <- H_mod$summary.fixed
	fixed_effect_posterior <- H_mod$marginals.fixed[[1]]
	lower_p <- inla.pmarginal(0, fixed_effect_posterior)
	upper_p <- 1 - inla.pmarginal(0, fixed_effect_posterior)
	post_pred_p <- 2 * (min(lower_p, upper_p))


	H_mod_output <- data.table(chrom = chrom, pos_start = pos_start, pos_end = pos_end, Paste_locus = Paste_locus, H_est = coef[1], H_p_value = post_pred_p)

	return(H_mod_output)
}

######################
##Read primary data ##
######################
full_dataset <- read.delim("./CnR/CnR_grh_great20_1mismap_CPM.txt", header = T)
full_dataset$Paste_locus <- paste(full_dataset$chrom, full_dataset$start, full_dataset$end, sep = "_")

##Get collumns for Parent, hybrid and hybrid parent data
Parental_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "zhr1_zhr", "zhr2_zhr", "zhr3_zhr", "z301_z30", "z302_z30", "z303_z30")]

Hybrid_data <- full_dataset[, c("chrom", "start", "end", "Paste_locus", "zhrz301_zhr", "zhrz302_zhr", "zhrz303_zhr", "zhrz301_z30", "zhrz302_z30", "zhrz303_z30")]

#############################################
##run approrpriate test based on argument 1##
#############################################
if (Arguments[1] == "Parents") {

  Parental_results <- do.call(rbind, mclapply(Parental_data$Paste_locus, function(x) Parental_model(x, Parental_data), mc.cores = 4))

  write.table(Parental_results, file = "./CnR/CnR_grh_great20_5mismap_CPM_bayes_parents.txt", row.names = F, quote = F)

} else if (Arguments[1] == "Hybrids"){

  Hybrid_results <- do.call(rbind, mclapply(Hybrid_data$Paste_locus, function(x) Hybrid_model(x, Hybrid_data), mc.cores = 4))

  write.table(Hybrid_results, file = "./CnR/CnR_grh_great20_5mismap_CPM_bayes_hybrids.txt", row.names = F, quote = F)
