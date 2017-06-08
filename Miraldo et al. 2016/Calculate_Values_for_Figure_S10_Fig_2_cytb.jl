########################################################################
# Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Michael K. Borregaard
#########################################################################

#Load the necessary library
using DataFrames
include("GD_summary_functions.jl")

# read in the mammals and subset to valid sequence comparisons
mamcyt = readtable("pairwise_equalarea_mammals_cytb.csv")
mamcyt = mamcyt[!isna(mamcyt[:seq1]),:]
mamcyt = mamcyt[mamcyt[:overlap] .>= 0.5,:]

# read in the amphibians and subset to valid sequence comparisons
amphs = readtable("pairwise_equalarea_amphibians_cytb.csv")
amphs = amphs[!isna(amphs[:seq1]),:]
amphs = amphs[amphs[:overlap] .>= 0.5,:]

# Calculate the values and write them to disk
amp = totalbasepairs(amphs)
amp[:genespecies] = totalspecies(amphs)[:x1]
writetable("amphfig2.csv", amp)

mam = totalbasepairs(mamcyt)
mam[:genespecies] = totalspecies(mamcyt)[:x1]
writetable("mamfig2.csv", mam)

total = vcat(mamcyt, amphs)
fig2 = totalbasepairs(total)
fig2[:genespecies] = totalspecies(total)[:x1]
writetable("totalfig2.csv", fig2)
