# install and load trio for it's read.pedfile function
# source("https://bioconductor.org/biocLite.R")
# biocLite("trio")

# prep data for BEDASSLE and conStruct

library(trio) # to import plink file
library(tidyverse)

# modified from github repo: QuentinRougemont/abc_inferences/02-additional_analyses/03-spacemix/00.scripts
# input files: 1) ped file 
# output: genotypic matrix (AC, TG, AA, CC, etc. with inds in row, snp in cols)
# some steps might be a bit circuitous now as I've modified it to fit my purposes

orig = read.pedfile("data/batch_1.plink.ped", first.row = TRUE) # load ped file, output from stacks::populations

# turn into matrix
orig.mat = as.matrix(orig) 

# vis - alleles are in separate columns at each locus
orig.mat[1:40,1:20]

# get rid of id columns for now
orig.mat.noid = orig.mat[, -c(1:6)] 

# make an interval for concatenating fields
start = seq(1, by = 2, length = ncol(orig.mat.noid) / 2) 

# concatenate contents of paired columns from original plink file
cat.loci = sapply(start, function(i, orig.mat.noid) paste(as.character(orig.mat.noid[,i]),as.character(orig.mat.noid[,i+1]), sep=""), orig.mat.noid = orig.mat.noid)

# now alleles are in one column at each locus
cat.loci[1:40, 1:20]

# reattach identities
cat.loci.id = cbind(orig.mat[,c(1:2)], as.data.frame(cat.loci)) 

# some rows are missing a lot of data
apply(as.data.frame(cat.loci), MARGIN=1, table)

# clean up names
cat.loci.id.clean = cat.loci.id %>% 
  separate(pid, into = c("pid", "remove"), extra = "merge") %>% 
  dplyr::select(-remove)

# this table has 1 row per individual, 1 column per locus
# each cell contains the genotype (i.e. AA, or 00 if data is missing)
# pop and ind ids are appended

# population ids here are different from the ones I usually use - get keys to use for geographic data:
# for each individual
ind.ids = cat.loci.id.clean %>% 
  dplyr::select(famid, pid) %>% 
  separate(pid, 2, into = c("popcode", "ind"))

# ids for each population
pop.ids = ind.ids %>% 
  dplyr::select(-ind) %>% 
  distinct()

# remove identifier columns from genotypic matrix
gen.mat = cat.loci.id.clean[, -c(1:2)] 

# make missing values NA
gen.mat.na = gen.mat
gen.mat.na[gen.mat.na == "00"] = NA

# transpose so that each column is an individual
gen.mat.tall.na = t(gen.mat.na)
gen.mat.tall.na[1:20, 1:20]

# make dataframe with numbers of individuals successfully genotyped at each locus
loci = data.frame(N = rowSums(!is.na(gen.mat.tall.na)))  
hist(loci$N)

# make dataframe with numbers of loci successfully genotyped in each ind
inds = data.frame(N = colSums(!is.na(gen.mat.tall.na)))  
hist(inds$N)

# get list of the different alleles at each loci
alleles.na = apply(cbind(substr(gen.mat.tall.na, 1, 1), substr(gen.mat.tall.na, 2, 2)), 1, unique) 
gen.mat.tall.na[1:20, 1:20]
alleles.na[1:3, 1:20]
alleles = apply(alleles.na, 2, sort) # sort alleles, alphabetical order; ignore NA
alleles[1:2, 1:20]

loci$num.alleles = apply(alleles, 2, length) # how many alleles are in each locus?
if(any(loci$num.alleles>2)) { stop("genot contains more than two alleles.\n") } # check that no locus has more than 2 different alleles

# transform the matrix as a class list and as characters 
if(is.matrix(alleles)) {alleles = lapply(apply(alleles, 1, as.list), as.character)}
alleles[[1]][1:20]
alleles[[2]][1:20]

loci$A1 = NA # make a new column that will contain the first allele
inds = which(loci$num.alleles > 0) # which rows have more than 0 alleles? (should be all of them)
loci$A1 = alleles[[1]] # put the first allele in column A1

loci[1:20, 1:3]

loci$A2 = NA # make a new column that will contain the second allele
inds = which(loci$num.alleles>1) # which rows have more than 1 alleles? (should be all of them)
loci$A2 = alleles[[2]] # put the second allele in column A2

loci[1:20, 1:4]

ref = loci$A1 # make a vector of first alleles
alt = NA # make an empty vector for alternate alleles
# inds = which(ref == loci$A1) 
alt[inds] = loci$A2[inds] # for loci in which the ref column = the primary allele, the alterate should be the second allele
# inds = which(ref == m$A2); alt[inds] = m$A1[inds] # for loci in which the ref column = the secondary allele, the alterate should be the primary allele
# this second case never happens, the alleles in ref are always the primary alleles

if(any(ref != loci$A1 & ref != loci$A2)) { warning("ref allele not present in genot for some m. Conversions for these m cannot be performed and will be coerced to NA.") }
# this checks whether, for any locus, there are genotypes missing the reference allele, I'm not sure how this could happen because ref = m$A1...

# columns of possible genotypes
loci$G2 = paste(ref, ref, sep = "")	# 2 copies of ref 
loci$G1.1 = paste(ref, alt, sep = "")	# 1 copy of ref, ref allele coded first
loci$G1.2 = paste(alt, ref, sep = "")	# 1 copy of ref, reversed coding
loci$G0 = paste(alt, alt, sep = "")	# 0 copy of ref
# loci$G2[is.na(ref)] = NA # inputting missing, this set doesn't seem to have any
# loci$G1.1[is.na(alt)] = NA
# loci$G1.2[is.na(alt)] = NA
# loci$G0[is.na(alt)] = NA

head(loci, 20)

# create empty final genotypic matrix
count.mat.na = matrix(NA, ncol = ncol(gen.mat.tall.na), nrow = nrow(gen.mat.tall.na), dimnames = dimnames(gen.mat.tall.na) )
# if both alleles are the reference allele, input 2
count.mat.na[gen.mat.tall.na == loci$G2] = 2 
# if one allele is the reference allele, input 1
count.mat.na[gen.mat.tall.na == loci$G1.1 | gen.mat.tall.na == loci$G1.2] = 1 
# if neither allele is the reference allele, input 0
count.mat.na[gen.mat.tall.na == loci$G0] = 0 

count.mat.na[1:20, 1:20]

count.mat = count.mat.na

# put ids back on
count.mat.ids = cbind((cat.loci.id.clean[, c(1:2)]), t(count.mat))

# get sample size by individual
# create empty matrix
ss.mat = matrix(NA, ncol = ncol(gen.mat.tall.na), nrow = nrow(gen.mat.tall.na), dimnames = dimnames(gen.mat.tall.na))
ss.mat[count.mat.na == 2] = 2 # if genotypic data is present, assign 2 
ss.mat[count.mat.na == 1] = 2
ss.mat[count.mat.na == 0] = 2
ss.mat[is.na(ss.mat)] = 0
ss.mat.ids = cbind((cat.loci.id.clean[,1:2]), t(ss.mat)) #  add sample ids
ss.mat.ids[1:20, 1:20]

# get population allele counts
counts.pop = count.mat.ids %>% 
  dplyr::select(-pid) %>% 
  group_by(famid) %>% 
  summarize_all(funs(sum(., na.rm = TRUE)))
dim(counts.pop)
counts.pop[1:20, 1:20]
counts.pop.rn = left_join(pop.ids, counts.pop) %>% select(-famid) %>% column_to_rownames("popcode")
counts.pop.rn[1:20, 1:20]
dim(counts.pop.rn)
counts.pop.matrix = as.matrix(counts.pop.rn)
counts.pop.matrix[1:20, 1:20]
dim(count.mat.ids)
count.mat.ids[,3]
count.mat.ids[1:20, 1:20]

# get population sample sizes
ss.mat.ids[1:20, 1:20]
ss.pop = ss.mat.ids %>% 
  dplyr::select(-pid) %>% 
  group_by(famid) %>% 
  summarize_all(funs(sum(., na.rm = TRUE)))
dim(ss.pop)
ss.pop[1:20, 1:20]
ss.pop.rn = left_join(pop.ids, ss.pop) %>% dplyr::select(-famid) %>% column_to_rownames("popcode")
ss.pop.rn[1:20, 1:20]
dim(ss.pop.rn)
ss.pop.matrix = as.matrix(ss.pop.rn)

counts.pop.matrix[1:20, 1:20]
ss.pop.matrix[1:20, 1:20]

# save data for bedassle
save(counts.pop.matrix, ss.pop.matrix, file = "data/allele_freqs.Rdata")

# get average genotyping success rate by population
ss.pop.average = data.frame(t(ss.pop.rn)) %>% 
  summarize_all(funs(mean))
ss.pop.max = data.frame(t(ss.pop.rn)) %>% 
  summarize_all(funs(max))
ss.pop.prop = ss.pop.average/ss.pop.max
ss.pop.print = data.frame(prop = t(ss.pop.prop)) %>% rownames_to_column("pop")

write.csv(ss.pop.print, "data/representation.csv", row.names = FALSE)


# prep for construct ------------------------------------------------------

# for conStruct, when ss is 0, replace allele freq with NA

na.pop.matrix = counts.pop.matrix
na.pop.matrix[ss.pop.matrix==0] = NA

na.pop.matrix[1:20, 1:20]
ss.pop.matrix[1:20, 1:20]

check = ss.pop.matrix - counts.pop.matrix
length(check[check < 0])
length(check[check >= 0])

construct.pop.matrix = na.pop.matrix/ss.pop.matrix
construct.pop.matrix[1:20, 1:20]

save(construct.pop.matrix, file = "data/construct_matrix.Rdata")


# now make datasets for each region ----------------------------------

load("data/allele_freqs.Rdata")
load("data/pairwise_matrices.Rdata")

# check that all pops are there
rownames(counts.pop.matrix)
rownames(tave_prism.mat)

counts.pop.matrix[1:20, 1:20]
ss.pop.matrix[1:20, 1:20]

# standardize to put geog and ppt on similar scale to temperature
sd(geodist.mat) # 164.7215
geodist.mat.std = geodist.mat/sd(geodist.mat)

sd(ppt_spring_prism.mat) # 7.817393
ppt_spring_prism.mat.std = ppt_spring_prism.mat/sd(ppt_spring_prism.mat)

plot(geodist.mat, tave_prism.mat)
plot(geodist.mat, fst.mat)
plot(tave_prism.mat, fst.mat)
# expect high effect size of dist relative to env

# make list of "northern" populations  
north = read.csv("data/population_data.csv") %>% 
  filter(cluster == "north") %>% 
  select(popcode) %>% 
  arrange(popcode)
north = north$popcode

# filter matrices to northern pops
counts.pop.matrix.north = counts.pop.matrix[north,]
ss.pop.matrix.north = ss.pop.matrix[north,]
geodist.mat.north.std = geodist.mat.std[north, north]
tave_prism.mat.north = tave_prism.mat[north, north]
ppt_spring_prism.mat.north.std = ppt_spring_prism.mat.std[north, north]

# save these inputs
save(counts.pop.matrix.north, ss.pop.matrix.north, geodist.mat.north.std, tave_prism.mat.north, ppt_spring_prism.mat.north.std, file = "data/bed_inputs_north.Rdata")

# make list of "central" populations  
center = read.csv("data/population_data.csv") %>% 
  filter(cluster == "central") %>% 
  select(popcode) %>% 
  arrange(popcode)
center = center$popcode

# filter matrices to central pops
counts.pop.matrix.center = counts.pop.matrix[center,]
ss.pop.matrix.center = ss.pop.matrix[center,]
geodist.mat.center.std = geodist.mat.std[center, center]
tave_prism.mat.center = tave_prism.mat[center, center]
ppt_spring_prism.mat.center.std = ppt_spring_prism.mat.std[center, center]

# save these inputs
save(counts.pop.matrix.center, ss.pop.matrix.center, geodist.mat.center.std, tave_prism.mat.center, ppt_spring_prism.mat.center.std, file = "data/bed_inputs_center.Rdata")

