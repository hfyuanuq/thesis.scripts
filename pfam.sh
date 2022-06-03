## This code identifies Pfam domains that are enriched among a selected set of genes.
## Created 6-24-08 by Gregory Hather.
##
## The inputs are:
## "Aq_All.txt" -- A list of all the genes considered (AGI numbers). In our case, it was all the chromosomal genes represented on the ATH1 Affymetrix microarray.
## "Pfam_annotation_final.txt" -- A table with two columns. The first column is the AGI number, and the second column is the Pfam domain.
## "pfam_description.txt" -- A table with two columns. The first column is the Pfam domain, and the second column is the description of that domain.
## "genes-with-no-conversion.txt" -- A list of AGI numbers which should be excluded because there was AGI2Uniprot conversion.
## 1.txt -- A list of AGI numbers for the subset of gene of interest.
##
## The output is:
## output_file -- a table listing the enriched domains
##
##
## Note:
## "converted-swisspfam-arath.txt", "swisspfam-descriptions.txt", and "genes-with-noconversion.txt" were generated using perl and python scripts.
## The inputs to the scripts were swisspfam (from Pfam 22.0) and AGI2Uniprot.20080418 (from TAIR).
## All the genes in "converted-swisspfam-arath.txt" are in "ATH1all.txt".
## No genes from "converted-swisspfam-arath.txt" are in "genes-with-noconversion.txt".




# ----------------------------------------------
# %% before starting:
# %% change depending on analysis either 
# %%%% adjp0.05/0.1 (3 occ) 
# %%%% posFC / posFC (3 occ)
# %%%% or trt combination : nat (3 occ)
# %% save script with analysis in name of file. 
# 2017.08.04
# ----------------------------------------------


## read data

################################################ to change / 1 up from  output folder
setwd("/opsin/u/huifang/project/ATAC.analysis/allinone/10.pfam.enrichment");
################################################ to change



# refer to the two new files below - now with Aqu2.1 model naming system
#ATH1all = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_rmtext_Aq_all.txt", what = list(""))[[1]]

ATH1all = scan("/opsin/u/huifang/project/ATAC.analysis/allinone/10.pfam.enrichment/input/Aqu2.1_rmtext_Aq_all.txt", what = list(""))[[1]]


################################################ to change
# conduct a pfam enrichment using only the "unique" pfam domains within each Aqu2.1 model
#converted_swisspfam = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_Pfam_annotation_final.txt", what = list("",""), sep = "\t")
#
# total abundance of pfam domains 
#converted_swisspfam = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_Pfam_annotation_final.txt", what = list("",""), sep = "\t")

converted_swisspfam = scan("/opsin/u/huifang/project/ATAC.analysis/allinone/10.pfam.enrichment/input/Aqu2.1_Pfam_annotation_final.txt", what = list("",""), sep = "\t") # this has been renamed and in the same file path as # above

################################################ to change



# no need to rename below
swisspfam_descriptions = scan("/opsin/u/huifang/project/ATAC.analysis/allinone/10.pfam.enrichment/input/pfam_description.txt", what = list("",""), sep ="\t", quote = "")

# 2016.10.25
# input files- listAq - copied from below location (using new v R and ~group analysis, and lists fo Aq models were extracted from VENNY diagrams (hence incl. shared etc). 
######################
#/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/DESeq2_TS1015_h03.02_~group_2016.10.24_newRv/TS1015_VENNY_lists/TS1015_t1.t9_group_adjp0.05_FC/t1.t9_group_adjp_0.05_posFC_nat-nat.txt
######################



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# to change: DEG1248
# to change: uniq
################################################ to change
top_genes = scan("/opsin/u/huifang/project/ATAC.analysis/allinone/10.pfam.enrichment/list/FoxK.FoxD.brown.all.genes.txt", what = list(""))[[1]]
output_file = "FoxK.FoxD.brown.out.pfam.txt"
################################################ always write as .txt file (file layout)




## split genes by domains
domain_grouping = split(converted_swisspfam[[1]], converted_swisspfam[[2]])

## Do a hypergeometric test for enrichment of each domain among top_genes
## The balls in the urn are all the genes that are in ATH1all and not in genes_with_no_conversion.
## q: Genes with given domain that are in the top set are the white balls drawn from the urn.
## m: Genes with given domain are white balls.
## n: Genes without given domain are black balls.
## k: The balls drawn are the genes in the top set.

## loop over each domain
p_value_vector = c()
enrichment_vector = c()
q_vector = c()
m_vector = c()
for (index in 1:length(domain_grouping)){
domain = names(domain_grouping)[index]
genes_with_domain = domain_grouping[[index]]
q = length(as.vector(na.omit(match(genes_with_domain, top_genes))))
m = length(genes_with_domain)
n = length(ATH1all) - m
k = length(top_genes)
p_value = phyper(q - 1, m, n, k, lower.tail = FALSE)
p_value_vector[index] = p_value
enrichment = q/(m*k/(n+m))
enrichment_vector[index] = enrichment
q_vector[index] = q
m_vector[index] = m
}

## adjust p-values
library(multtest)
adjusted_p = mt.rawp2adjp(p_value_vector, proc = c("BH"))

## write result
p_cut_off = 0.1
significant_indices = as.vector(na.omit(ifelse(adjusted_p$adjp[,2] < p_cut_off, adjusted_p$index, NA)))
significant_enrichments = enrichment_vector[significant_indices]
indices_by_enrichment = significant_indices[sort(significant_enrichments, decreasing = TRUE, index.return = TRUE)[[2]]]
sorted_descriptions = swisspfam_descriptions[[2]][match(names(domain_grouping)[indices_by_enrichment], swisspfam_descriptions[[1]])]
sorted_Pfam_ID = names(domain_grouping)[indices_by_enrichment]
sorted_full_descriptions = paste(sorted_descriptions, " (", sorted_Pfam_ID, ")", sep = "")
sorted_enrichment = enrichment_vector[indices_by_enrichment]
sorted_p_values = adjusted_p$adjp[match(indices_by_enrichment, adjusted_p$index),2]
sorted_q = q_vector[indices_by_enrichment]
sorted_m = m_vector[indices_by_enrichment]
length(top_genes)
length(ATH1all)
sorted_gene_groups = c()
for (index in 1:length(indices_by_enrichment)){
temp_genes = top_genes[as.vector(na.omit(match(domain_grouping[[indices_by_enrichment[index]]], top_genes)))]
sorted_gene_groups[index] = paste(temp_genes, sep = ", ", collapse = ", ")
}
output = cbind(sorted_full_descriptions, round(sorted_enrichment, digits = 1), as.character(format(signif(sorted_p_values, digits = 1), scientific = TRUE)), sorted_q, sorted_m, sorted_gene_groups)

write(t(output), file = output_file, sep = "\t", ncolumns = 6)


writeLines(capture.output(sessionInfo()), "FoxK.FoxD.brown.out.sessionInfo.txt")

