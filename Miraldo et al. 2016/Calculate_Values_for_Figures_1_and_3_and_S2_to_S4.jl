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


for typ in ["gridcells", "anthromes", "biomes", "latbands", "equalarea"]
    # Read in the matrices for mammals and amphibians and concatenate them
    am = ready_data("pairwise_$(typ)_amphibians_cytb.csv")
    ma = ready_data("pairwise_$(typ)_mammals_cytb.csv")
    macoi = ready_data("pairwise_$(typ)_mammals_coi.csv")
    dats = vcat(am, ma)
    # restrict to pairs with at least 0.5 overlap value

    amphibs, amphibsfull = calc_cellvalues(am)
    mamms, mammalfull  = calc_cellvalues(ma)
    total, full = calc_cellvalues(dats)
    macois, fullmacois = calc_cellvalues(macoi)

    fig_data = join(total, mamms, on = :site, kind = :left)
    fig_data = join(fig_data, amphibs, on = :site, kind = :left)
    fig_data = join(fig_data, macois, on = :site, kind = :left)

    names!(fig_data, [:site, :totalGDval, :totalrichness, :totalbasepairs, :mammalsGDval, :mammalsrichness, :mammalsbasepairs, :amphibiansGDval, :amphibiansrichness, :amphibiansbasepairs, :macoiGDval, :macoirichness, :macoibasepairs])
    writetable("$(typ)_numbers.csv", fig_data)
    writetable("$(typ)_full_numbers.csv", full)
    writetable("$(typ)_amphibsfull_numbers.csv", amphibsfull)
    writetable("$(typ)_mammalfull_numbers.csv", mammalfull)
    writetable("$(typ)_macoifull_numbers.csv", fullmacois)
end



ma = readtable("pairwise_equalarea_mammals_coi.csv")
ma = ma[!isna(ma[:seq1]), :]
ma = ma[ma[:overlap] .> 0.5, :]
# restrict to pairs with at least 0.5 overlap value

mamms, mammalfull  = calc_cellvalues(ma)

names!(mamms, [:site, :mammalsGDval, :mammalsrichness, :mammalsbasepairs])
writetable("equalarea_coi_numbers.csv", mamms)
writetable("equalarea_coi_full_numbers.csv", mammalfull)
