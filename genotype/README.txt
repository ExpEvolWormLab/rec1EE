

* Note: 
** the samples were sequenced at 10X for lines and 60X for population
** sequencing and libraries were done at BGI tech solution Honk-Kong (WGS; random library selection; BGISEQ-500; 150bp paired reads).
** alignement to ws245 reference genome was done with bwa mem (default) and mapped read filtering was done with samtools view; -q 20 -f 0x0002 -F 0x0004 -F 0x0008 options.
** variant calling was with freebayes (option -@) calling only know SNPs from CeMEEv2 (https://doi.org/10.1093/g3journal/jkaa061)


#################################################################################################
### Pool-sequencing data during adaptation of the sister population (MR: rec-1 wild-type; or LR:Mutant)
### (Experiment II in article; code name = Rhom)
### The samples were sequenced at 60X at BGI tech solution Honk-Kong (BGIseq)


-> chr#_Rhom_snps_AllelicDepth_CeMEEv2_ws245.csv: SNPs information
  => ID: SNPs unique identifier Chrom_pos
  => chrom: chromosome (I,II,III,IV,V,X)
  => pos = physical position (bp, ws245)
  => ref & alt: SNP genotype
  => cM.wt & cM.mut: rec-1 wild-type and mutant genetic distance (cM) from paree et al,2024 (https://doi.org/10.1093/genetics/iyad205), corrected to yield 50cM.
  => rrdiff: difference between mutant and wild-type recombination rate (cM/Mb) within 1Mb genetic distance around that SNP (referred to as Δr in the article)
  => domain: recombination domains (left or right tips, arms, or centers) assuming boundaries of Rockman & Kruglyak, 2009 (https://doi.org/10.1371/journal.pgen.1000419)
  => domain_type: tip, arm, center (no specification of right/left)

-> chr#_Rhom_poolsed_AllelicDepth_CeMEEv2_ws245.csv: 
  => is the allelic depth (count) in format ref:alt
  => column names are the samples ID is format popname_GX (e.g, "SLR5_G33")
     where the pop name is LR (leveled recombination, rec-1 mutant) or HR (heterogenous recombination, rec-1 wild-type),
     with a S in front if evolution was done in the High salt environment (if not, in the domestication environment),
     and followed by a first number being the replicate number, and GX being the generation X.
  

-> create3Darray.R: R script to execute to obtain the two following files (used in the analysis):

* Rhom_poolsed_unfiltered.Rdata: same data as above but different format. Contains:
  => snps: data.frame with same info as Rhom_snps_AllelicDepth_CeMEEv2_ws245.csv
  => alt & ref: two 3D array containing reference and alternative allele count (dimensions = SNP ID, popname, generation)
  
  
* Rhom_poolsed.Rdata: 
  => same as Rhom_poolsed_unfiltered.Rdata but SNP with depth in the lower 5% quantile or fixed are flitered out
  => few outlier samples were excluded
  => This dataset is the one used for the genomic analysis 
  




###################################################################################
### Sequencing of inbred lines from SMR populations (rec-1 polymorphic) ###########
### (Experiment III; code name = Rpoly) ###########################################
The samples were sequenced at 10X at BGI tech solution Honk-Kong (BGIseq)

-> chr#_SMR_inbredlines_snps_CeMEEv2_ws245.csv.gz
-> chr#_SMR_inbredlines_genotype_CeMEEv2_ws245.csv.gz

Files contains snps and genotype tables of the SMR lines and the EEV1401 line. The EEV1403 line have the same genotype as the EEV1401 line (the two lines only differ for the CRISPR-Cas9 deletion in the rec-1 gene, see: Paree et al, 2024:  https://doi.org/10.1093/genetics/iyad205) .










