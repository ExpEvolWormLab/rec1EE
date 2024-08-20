
(!) for some of these analyses, you will need first to execute the genotype/create3Darray.R script

(!) some relevant resources are found at the https://github.com/ExpEvolWormLab/PAREE/ and https://github.com/lukemn/


*** USEFUL NOMENCLATURE (AND SYNONYM):

-> Environments
domestication (synonym: standard; ngm; NGM): domestication environment
salt (synonym: high salt; nacl): high salt environment (9-fold NaCl increase in the medium)

-> strains:
EEV1401: one of the strains used (rec-1 wt)
EEV1402: one of the strains used, with a GFP construct (rec-1 wt)
EEV1403: rec-1 mutant in EEV1401 background
EEV1404: rec-1 mutant in EEV1402 background

-> populations: 
(S)HR#: rec-1 WT genetically population (evolved in salt), replicate # 
(S)LR#: rec-1 mut genetically population (evolved in salt), replicate # 
(S)MR#: rec-1 polymorphic genetically population (evolved in salt), replicate # 

-> Experiment ID:
Experiment I: indirect selection experiment with rec-1 isogenic population
Experiment II (synonym: Rhom): adaptation experiment with rec-1 wt or mut sister population.
Experiment III (synonym: Rpoly): indirect selection experiment with rec-1 polymorphic population.




*** CONTENT: 

***** FILES:
#rec1EE_utils.R: redundant functions used in the different analysis


***** FOLDERS:

# genotype: contain the pool-sequencing and inbred lines genotype tables.
# ancestral_genomics: contains code relative to Figure 2 & S2 (SNV density, delta r, population genetics predictions)
# heritability: contains the code used to calculate the heritability for fitness (Fig 2F)
# adaptation_fitness: analysis of fitness before and after adaptation of HR & LR populations (Experiment II)
# caldidateSelection: candidate SNPs for significant differentiation during adaptation of HR and LR populations (Experiment II)
# direct_selection: rec-1 mutant frequency and selection coefficients in isogenic populations (Fig 1, Experiment I)
# observed_r2: allele frequency correlation between evolution replicates (Fig 3b; Experiment II)
# predicted_r2: expected linkage, from population genetic equations (Fig 2E)
# outcrossing&outliers: PCA of afc & outcrossing data
# salt: effect of high NaCl concentration
# recombinationMaps: delta r & rec-1 recombination maps 



