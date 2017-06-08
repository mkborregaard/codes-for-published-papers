Analysis codes for the paper: An Anthropocene Map of Genetic Diversity
=======================================================================

Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Flórez-Rodríguez, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo

Submitted to Science, 2016.

This folder contains the codes to reproduce all the analyses in the paper, and some base data files.

## Data files
The included data files are:

* `Anthro1by1_world.csv`    : A file listing the prevalent anthrome in each terrestrial 1x1 lat-long grid cell of the world. See the Supplementary Information for how this file was created.
* `BIOMES_1BY1.csv`         : A file listing the prevalent terrestrial biome in each 1x1 lat-long grid cell of the world.
* `plotcols.txt`            : sRGB values for the color scheme used in most maps (e.g. Fig. 1 of the main text).

## Additional data
To run the codes, additional data is needed:

### Georeferenced genetic data
The genetic analyses depend on georeferenced and aligned species sequences for cytochrome B. They should be in one folder with files of two types (file extension 'fasta' and file extension 'coords'). The two types should have the same filenames, (the species names; e.g. folder contents should be 'Bufo_bufo.fasta, Bufo_bufo.coords, Rana_arvalis.fasta, Rana_arvalis.coords', etc.). The .fasta files contain the sequences as an m x n (row x column) integer matrix, where m is the number of sequences and n is the length of the alignments. DNA bases are coded as `[1, 2, 3, 4]`, with a `0` indicating a base missing in the alignment. The .coords files contain the geographic coordinates of the sequences, as an m x 2 floating-point-number matrix with latitude in the first column and longitude in the second. See the Methods section for information on how these files were generated. A list of accession numbers can be downloaded from [iMapGenes](http://www.macroecology.ku.dk/resources/imapgenes).

### IUCN distribution shapefiles
The file `Grid_IUCN_polygons.R` creates `R` data objects (in the format of the `nodiv` package) from the distribution shapefiles downloaded from IUCN. These shapefiles are freely available for download from: [IUCN](http://www.iucnredlist.org/technical-documents/spatial-data).


### Extra files used for plotting Figure 2
Figure 2 depends on all georeferenced sequences (not just those used for cytochrome B.) To get these data, we followed the procedure described in the Supplementary Information. For the codes, these were summarized into these files (not included):

* `species_grid_amphs`                : The number of amphibian species with any georeferenced sequences in a Behrmann grid cell
* `species_grid_mammals`              : The number of mammal species with any georeferenced sequences in a Behrmann grid cell
* `tot_length_grid_amphs`             : The total number of basepairs of georeferenced sequences for amphibians in a Behrmann grid cell
* `tot_length_grid_mammals`           : The total number of basepairs of georeferenced sequences for mammals in a Behrmann grid cell
* `corrected_basepairs_gridcell.csv`  : The total number of basepairs of georeferenced sequences for mammals and amphibians in a Behrmann grid cell


## Code files
The analytical codes are written in two different languages: computationally demanding operations were written in the very fast programming language `julia`, whereas most plotting operations were performed in `R`, which has very versatile plotting capabilities.


### Julia source files (in expected order of use)
* `Compare_Pairwise_Function.jl`    : An implementation of the first step in calculating GD (our Genetic Diversity metric), generating all pairwise comparisons of aligned sequences.
* `Create_Master_Matrices.jl`       : Uses the pairwise function to create master matrices from the genetic data files, used in all following calculations.
* `Create_matrices_script.jl`       : Creates the matrices for various grid cells or other spatial units
* `GD_summary_functions.jl`         : The basic functions for calculating GD  and other basic summaries of the master matrices.
* `Calculate_Values_for_Figures_1_and_3_and_S2_to_S4.jl` : Calculates the basic values used by the plotting functions
* `Calculate_Values_for_Figure_S10_Fig_2_cytb_.jl`: Calculates the basic values to plot Fig S10, the Supplementary Information version of Figure 2 (based just on CytB)
* `Calculate_Values_for_S5_to_S9.jl` : Calculates the numbers used in the sensitivity analyses in the Supplementary Information.


### R source files (in expected order of use)
* `Overlay_grid_files.R`            : Assigns the genetic sequences to the Behrmann grid based on their geographic coordinates
* `New_Grid_IUCN_polygons.R`        : The codes used to extract species' distributions from the IUCN shapefiles, apply them on a Behrmann grid and generative `R` objects (in the `nodiv` package format) from them.
* `Figure_1.R`                      : Uses the numbers generated by the julia scripts to plot the three maps constituting Figure 1.
* `Figure_2.R`                      : Script for generating Figure 2 based on the data files for all genes
* `Figure_S10_Supplementary_Fig_2_cytb_script.R`   : Script for generating Figure 2 just based on the Cytochrome B data (included in the supplementary materials)
* `Figure_3_and_S2_to_S4_script.R`  : Generate Fig. 3 and supplemental Figs.
* `Figures_S6_to_S9`                : Mainly plotting functions used to generate the maps in the sensitivity analysis in the Supplementary Information.
* `Figures_S11_to_S15`                : Mainly plotting functions and some statistic analyses used to generate the respective figures in the Supplementary Information.




## Analytical workflow:
The analytical workflow followed this procedure:

The first steps are the most time consuming and occur before using the codes presented here.

1. Generate files with georeferenced sequences of all genes for mammals and amphibians.

2. Generate files with georeferenced and aligned sequences of CytB and CO1 for mammals and and of CytB for amphibians.

3. Use the coordinates to assign each sequence to a biome and an anthrome based on published shapefiles.

The remaining steps use the codes presented here:

4. Generate gridded distributions for amphibians and mammals.
  * Download the IUCN shape files
  * Run `New_Grid_IUCN_polygons.R` to generate R files

5. Create master matrices for the analysis of the genetic data
  * Run `Create_matrices_script.jl`

6. Plot Figures 1 and 3
  * Run `Calculate_Values_for_Figures_1_and_3_and_S2_to_S4.jl`
  * Run `Figure_1.R`
  * Run `Figure_3_and_S2_to_S4_script.R`

7. Plot Figure 2
  * Run `Figure_2_script.R`

8. Plot Figure S10 - Supplementary version of Fig 2
  * Run `Calculate_Values_for_Figure_S10_Fig_2_cytb_.jl`
  * Run `Figure_S10_Supplementary_Fig_2_cytb_script.R`

9. Supplementary sensitivity analysis
  * Run `Calculate_Values_for_S5_to_S9.jl`
  * Run `Figures_S6_to_S9.R`

10. Supplementary analyses regarding Fig S11 to FigS15
  * Run `Figures_S11_to_S15.R`
