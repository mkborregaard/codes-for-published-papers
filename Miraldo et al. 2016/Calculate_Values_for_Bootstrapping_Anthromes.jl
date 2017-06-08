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
using DataFrames, DataFramesMeta
include("GD_summary_functions.jl")

function boot_species(dat::DataFrame)
    by(dat, :cell) do df
        nrows = size(df, 1)
        rows = rand(1:nrows, nrows)
        mean(df[:x1][rows])
    end
end

function bootstrap_species(dat::DataFrame)
    boot_full = hcat([boot_species(dat)[:x1] for i in 1:100]...) #creates a matrix of cell values
    ret = by(df -> mean(df[:x1]), dat, :cell)
    names!(ret, [:site, :empirical])
    ret[:boot_sd] = vec(mapslices(std, boot_full, 2))
    ret[:boot_low] = vec(mapslices(x -> quantile(x,  [:0.025]), boot_full, 2))
    ret[:boot_high] = vec(mapslices(x -> quantile(x,  [:0.975]), boot_full, 2))
    ret
end



for typ in ["anthromes", "biomes", "latbands"]
    # Read in the matrices for mammals and amphibians and concatenate them
    am = ready_data("pairwise_$(typ)_amphibians_cytb.csv")
    ma = ready_data("pairwise_$(typ)_mammals_cytb.csv")
    macoi = ready_data("pairwise_$(typ)_mammals_coi.csv")
    dats = vcat(am, ma)
    # restrict to pairs with at least 0.5 overlap value

    amphibs, amphibsfull = calc_cellvalues(am)
    mamms, mammalfull  = calc_cellvalues(ma)
    macoi, macoifull = calc_cellvalues(macoi)
    total, full = calc_cellvalues(dats)

    bootfull = bootstrap_species(full)
    bootmams = bootstrap_species(mammalfull)
    bootamphibs = bootstrap_species(amphibsfull)
    bootmacoi = bootstrap_species(macoifull)

    boot_data = join(bootfull, bootmams, on = :site, kind = :left)
    boot_data = join(boot_data, bootamphibs, on = :site, kind = :left)
    boot_data = join(boot_data, bootmacoi, on = :site, kind = :left)
    names!(boot_data, [:site, :total_empirical, :total_SD, :total_0_025, :total_0_975, :mamm_empirical, :mamm_SD, :mamm_0_025, :mamm_0_975,:amph_empirical, :amph_SD, :amph_0_025, :amph_0_975,:macoi_empirical, :macoi_SD, :macoi_0_025, :macoi_0_975])

    writetable("$(typ)_bootstrap_numbers.csv", boot_data)
end
