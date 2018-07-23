# source("https://bioconductor.org/biocLite.R")
# biocLite("trio")

library(trio) # to import plink file
library(hierfstat)
library(tidyverse)

# calculate pairwise Fst and other basic popgen stats

# prepare data ------------------------------------------------------------

orig = read.pedfile("data/batch_1.plink.ped", first.row = TRUE) # load ped file, output from stacks::populationss

pops = orig %>% select(pid) %>% tidyr::separate(pid, 2, into = c("pop", "omit")) %>% select(-omit)

orig.mat = as.matrix(orig) # turn into matrix

# need to replace ACGT with 1234
orig.mat[orig.mat == "A"] = 1
orig.mat[orig.mat == "C"] = 2
orig.mat[orig.mat == "G"] = 3
orig.mat[orig.mat == "T"] = 4

orig.mat.noid = orig.mat[, -c(1:6)] # get rid of id columns for now

orig.mat[1:10, 1:10]
orig.mat.noid[1:10, 1:10]

start = seq(1, by = 2, length = ncol(orig.mat.noid) / 2) # make an interval for concatenating fields

cat.loci = sapply(start, function(i, orig.mat.noid) paste(as.character(orig.mat.noid[,i]),as.character(orig.mat.noid[,i+1]), sep=""), orig.mat.noid = orig.mat.noid) 
# concatenate contents of paired columns from original plink file

cat.loci.id = cbind(pops, as.data.frame(cat.loci)) # reattach identitites
# cat.loci.id = cbind(pop = orig.mat[,1], as.data.frame(cat.loci))

cat.loci.id[1:10, 1:10]

cat.loci.id[cat.loci.id == "00"] = NA
# convert 00 to NA

write.csv(cat.loci.id, "data/batch_1_intermediate.csv", row.names = FALSE)
num.loci = read.csv("data/batch_1_intermediate.csv", stringsAsFactors = FALSE)


# calculate fst -----------------------------------------------------------

# calculate fst - wc (takes a while, 1-2 hr)
fst.wc = pairwise.WCfst(num.loci, diploid = TRUE)

# make tall dataframe
fst.wc.df = rownames_to_column(data.frame(fst.wc)) %>% 
  rename(pop1 = rowname) %>% 
  gather("pop2", "fst_wc_hfs", 2:33) %>% 
  filter(!is.na(fst_wc_hfs))

# write out
write.csv(fst.wc, "data/fst_wc.csv")


# other basic stats -------------------------------------------------------

# fis, gene diversity (exp. het.)
num.loci[1:10, 1:10]

basic_stats = basic.stats(num.loci, diploid = TRUE)

dim(basic_stats$Hs)

summary(basic_stats$Hs)
summary(basic_stats$Fis)

hist(basic_stats$perloc$Fst)
fst_df = data.frame(basic_stats$perloc)

ggplot(fst_df, aes(x = Fst)) +
  geom_histogram(fill = "coral", color = "black", binwidth = 0.02) +
  ylab("Frequency") +
  xlab(expression(F["ST"]))
# ggsave("figs_for_paper/fst_dist.pdf", height = 3, width = 5)

fis = rownames_to_column(data.frame(apply(basic_stats$Fis, 2, mean, na.rm = TRUE)))
names(fis) = c("pop", "fis")

expected_het = rownames_to_column(data.frame(apply(basic_stats$Hs, 2, mean, na.rm = TRUE)))
names(expected_het) = c("pop", "exp_het")

obs_het = rownames_to_column(data.frame(apply(basic_stats$Ho, 2, mean, na.rm = TRUE)))
names(obs_het) = c("pop", "obs_het")


# merge and write out new values

both = left_join(fis, expected_het)

both = left_join(both, obs_het)

plot(both$exp_het~ both$obs_het)

plot(both$fis~both$obs_het)
plot(both$fis~both$exp_het)

ggplot(both, aes(exp_het, obs_het, color = fis)) +
  geom_point()

write.csv(both, "data/pops_stats.csv", row.names = FALSE)
