######################################################################## # Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Science, 2016
# Code in this file by Sen Li and Michael K. Borregaard
#########################################################################

# Load the necessary library and functions
include("Compare_Pairwise_Function.jl")
using DataFrames

"""
A function to calculate the summary statistic for all sites (e.g. grid cell or biome) for one species.

**Parameters**
* 'species':        A string with the name of the species
* 'species_seqs':   A matrix where the rows are aligned genetic sequences, and columns are loci. Basepairs must be coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is absent from the alignment.
* 'sitenames':      A vector of strings with the names of the sites
"""
function calcspecies(species::String, species_seqs::Matrix{Int}, sitenames::Vector{String})
    res = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])
    uniq_grid = unique(sitenames, 1)
    for iter_grid = 1:size(uniq_grid,1) # Go through each site in turn
        seqs_grid = species_seqs[findin(sitenames, uniq_grid[iter_grid, :]), :]

        tot_muts = compare_pairwise(seqs_grid)
        lns = size(tot_muts, 1)
        tmp = DataFrame(species = repeat([species], inner= lns), cell = repeat(vec(uniq_grid[iter_grid, :]), inner = lns)) #Expand species and cell names to the length of the resulting DataFrame
        tot_muts = hcat(tmp, tot_muts)

        res = vcat(res, tot_muts)
    end

    for sym in [:seq1, :seq2, :length_seq1, :length_seq2, :overlap, :commons, :num_per_bp]
        res[sym][res[sym] .< 0.] = NA   #Replace empty values with an NA identifier
    end
    res
end


"""
A function to create master data matrices that are used to compute genetic diversity, assess data quality and do sensitivity analyses.

**Parameters**
* 'foldername' : The name of of a folder containing the data files. There must be files of two types (file extension 'fasta' and file extension 'coords') with  the same filename, e.g. the species names (e.g. folder contents could be 'Bufo_bufo.fasta, Bufo_bufo.coords, Rana_arvalis.fasta, Rana_arvalis.coords', etc.). The .fasta files contain the sequences as an m x n integer matrix, where m is the number of sequences and n is the length of the alignments. The .coords files contain the geographic coordinates of the sequences, as an m x 2 floating point matrix with latitude in the first column and longitude in the second.  TODO: A folder with the anthrome and biome of each sequence
"""
function create_master_matrices(foldername::String, sitesfolder::String)
    species_list = [x[1:(end-6)] for x in filter(st->contains(st, ".fasta"), readdir(foldername))]      #identify unique file names ignoring the extension
    num_files = size(species_list,1)

    #equalarea = latbands = biomes = anthrome = gridcells = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])   # Pre-initialize the DataFrame to ensure correct element types
    ret = DataFrame(species = String[], cell = String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])   # Pre-initialize the DataFrame to ensure correct element types

    for iter_file = 1:num_files
        # read in the sequence data for one species
        file_name_seq = joinpath(foldername, species_list[iter_file] * ".fasta")
        species_seqs = readdlm(file_name_seq, Int)
        length_seq = size(species_seqs,2)

        # read the coordinates for one species
        file_name_sites = joinpath(foldername, species_list[iter_file] * ".coords")
        species_sites = readdlm(file_name_sites, String)[:,1]
        #if eltype(species_sites) <: AbstractFloat
        #    species_sites = floor(Int, species_sites)
        #end

        ret = vcat(ret, calcspecies(species_list[iter_file], species_seqs, species_sites))
        println(iter_file)

    end
    writetable("pairwise_$(sitesfolder)_$(foldername).csv", ret)
end

function grid_coordinates(foldername = ".", xsize = 1, ysize = xsize, lowerleftoffset = (0,0))
    in("gridcells", readdir(foldername)) && error("There is already a folder named gridcells!")
    mkdir(joinpath(foldername, "gridcells"))
    species_list = [x[1:(end-7)] for x in filter(st->contains(st, ".coords"), readdir(foldername))]
    for species in species_list
        species_coords = readdlm(joinpath(foldername,"$(species).coords"))
        coords_grid = hcat(ceil(Int, (species_coords[:,1] + lowerleftoffset[1]) ./ xsize) .* xsize .- lowerleftoffset[1] .- 0.5xsize, ceil(Int, (species_coords[:,2] + lowerleftoffset[2]) ./ ysize) .* ysize .- lowerleftoffset[2] .- 0.5ysize)
        sitenames = mapslices(x -> "$(x[1])_$(x[2])", coords_grid, 2)
        writedlm(joinpath(foldername, "gridcells", species), sitenames)
    end
end

function latband_coordinates(foldername = ".", bandwidth = 1)
    in("latbands", readdir(foldername)) && error("There is already a folder named latbands!")
    mkdir(joinpath(foldername, "latbands"))
    species_list = [x[1:(end-7)] for x in filter(st->contains(st, ".coords"), readdir(foldername))]
    for species in species_list
        species_coords = readdlm(joinpath(foldername,"$(species).coords"))
        latband = bandwidth * floor(species_coords[:,1]/bandwidth)
        writedlm(joinpath(foldername, "latbands", species), latband)
    end
end
