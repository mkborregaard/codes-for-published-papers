
##################
# Script         #
##################

include("Create_Master_Matrices.jl")

# This code can be run if you have .coords files in your directory
grid_coordinates("testfiles", 4, 4, (2,0))  #apply 4x4 grid cells to the coordinates before calculating GD
create_master_matrices("testfiles", "gridcells")

#apply 10-degree wide latitudinal bands to the coordinates
latband_coordinates("testfiles", 10)    # apply 10-degree wide latitudinal bands before calculating GD
create_master_matrices("testfiles", "latbands")

# and so forth if you have data for biomes, equal-area-cells etc
#create_master_matrices("testfiles", "anthromes")
#create_master_matrices("testfiles", "biomes")
#create_master_matrices("testfiles", "equal_area")
