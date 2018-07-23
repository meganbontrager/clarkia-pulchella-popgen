library(BEDASSLE)

# full dataset ------------------------------------------------------------

load("data/allele_freqs.Rdata")
load("data/pairwise_matrices.Rdata")

rownames(counts.pop.matrix)
rownames(tave_prism.mat)

counts.pop.matrix[1:20, 1:20]
ss.pop.matrix[1:20, 1:20]

# unlike regional versions, these haven't been standardized yet
sd(geodist.mat) # 164.7215
geodist.mat.std = geodist.mat/sd(geodist.mat)

sd(ppt_spring_prism.mat) # 7.817393
ppt_spring_prism.mat.std = ppt_spring_prism.mat/sd(ppt_spring_prism.mat)

hist(geodist.mat.std)
hist(tave_prism.mat)
hist(ppt_spring_prism.mat.std)

# first run to get in the ballpark before reducing step size
MCMC(counts = counts.pop.matrix, 
     sample_sizes = ss.pop.matrix, 
     D = geodist.mat.std, 
     E = list(tave_prism.mat, ppt_spring_prism.mat.std),
     k = 32, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.25, 
     aE_stp = 0.005, 
     a2_stp= 0.01, 
     thetas_stp = 0.1, 
     mu_stp = 0.17, 
     ngen = 1e6, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     prefix = file.path("bedassle_output","bed_tp_run9_main_"))

make.continuing.params("bedassle_output/bed_tp_run9_main_MCMC_output1.Robj", file = "bedassle_output/bed_tp_run9_main_Continuing.parameters")

# second run from first's parameters
MCMC(counts = counts.pop.matrix, 
     sample_sizes = ss.pop.matrix, 
     D = geodist.mat.std, 
     E = list(tave_prism.mat, ppt_spring_prism.mat.std),
     k = 32, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.01, 
     aE_stp = 0.00001, 
     a2_stp= 0.001, 
     thetas_stp = 0.05, 
     mu_stp = 0.17, 
     ngen = 1e7, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     continue = TRUE,
     continuing.params = "bedassle_output/bed_tp_run9_main_Continuing.parameters",
     prefix = file.path("bedassle_output","bed_tp_run10_main_"))


# central populations -----------------------------------------------------

load("data/bed_inputs_center.Rdata")

hist(tave_prism.mat.center)
hist(ppt_spring_prism.mat.center.std)
hist(geodist.mat.center.std)

MCMC(counts = counts.pop.matrix.center, 
     sample_sizes = ss.pop.matrix.center, 
     D = geodist.mat.center.std, 
     E = list(tave_prism.mat.center, ppt_spring_prism.mat.center.std),
     k = 9, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.25, 
     aE_stp = 0.005, 
     a2_stp= 0.01, 
     thetas_stp = 0.1, 
     mu_stp = 0.17, 
     ngen = 1e6, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     prefix = file.path("bedassle_output","bed_tp_run13_center_"))

make.continuing.params("bedassle_output/bed_tp_run13_center_MCMC_output1.Robj", file = "bedassle_output/bed_tp_run13_center_Continuing.parameters")

MCMC(counts = counts.pop.matrix.center, 
     sample_sizes = ss.pop.matrix.center, 
     D = geodist.mat.center.std, 
     E = list(tave_prism.mat.center, ppt_spring_prism.mat.center.std),
     k = 9, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.01, 
     aE_stp = 0.00001, 
     a2_stp= 0.001, 
     thetas_stp = 0.05, 
     mu_stp = 0.17, 
     ngen = 1e7, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     continue = TRUE,
     continuing.params = "bedassle_output/bed_tp_run13_center_Continuing.parameters",
     prefix = file.path("bedassle_output","bed_tp_run14_center_"))


# northern populations ----------------------------------------------------

load("data/bed_inputs_north.Rdata")

MCMC(counts = counts.pop.matrix.north, 
     sample_sizes = ss.pop.matrix.north, 
     D = geodist.mat.north.std, 
     E = list(tave_prism.mat.north, ppt_spring_prism.mat.north.std),
     k = 14, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.25, 
     aE_stp = 0.005, 
     a2_stp= 0.01, 
     thetas_stp = 0.1, 
     mu_stp = 0.17, 
     ngen = 1e6, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     prefix = file.path("bedassle_output","bed_tp_run11_north_"))

make.continuing.params("bedassle_output/bed_tp_run11_north_MCMC_output1.Robj", file = "bedassle_output/bed_tp_run11_north_Continuing.parameters")

MCMC(counts = counts.pop.matrix.north, 
     sample_sizes = ss.pop.matrix.north, 
     D = geodist.mat.north.std, 
     E = list(tave_prism.mat.north, ppt_spring_prism.mat.north.std),
     k = 14, 
     loci = 2982, 
     delta = 0.0001, 
     aD_stp = 0.01, 
     aE_stp = 0.00001, 
     a2_stp= 0.001, 
     thetas_stp = 0.05, 
     mu_stp = 0.17, 
     ngen = 1e7, 
     printfreq = 10000,
     samplefreq = 1000,     
     savefreq = 1000,
     continue = TRUE,
     continuing.params = "bedassle_output/bed_tp_run11_north_Continuing.parameters",
     prefix = file.path("bedassle_output","bed_tp_run12_north_"))


