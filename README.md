# Toolbox IME detection #

This toolbox contains the programs used to detect Island Mass Effects (IMEs) from a chlorophyll map and an island database.   
The method was developed within the [SAPPHIRE](https://sapphire.mio.osupytheas.fr) project funded by Horizon 2020 (Marie Skłodowska-Curie grant agreement No. 746530).  
The programs are written for Matlab and require the Image Processing toolbox.


### Get started

For examples on how to run the IME detection, see script:  

	start_IME_toolbox

Description of functions:  
`ime_detect_IMEs`: detects IMEs for a list of islands, based on a Chl map.  
`ime_combine_islands` (called by ime_detect_IMEs): combines islands within NaN patches (mostly useful on non-curated island databases and/or Chl maps with missing data).  


### Description of data inputs, available here for demonstration purposes

`Chl_climatology.mat` contains 12 climatological maps of MODIS chlorophyll in the tropical Pacific, originally downloaded from https://coastwatch.pfeg.noaa.gov/erddap/info/index.html (erdMH1chlamday product) with data within the shallow mask (above 30 m extended by 1 pixel in all directions) set to NaN.   
Also contains the mask for area of study (.mask_cont), the shallow pixel mask (.mask_30m) and pixel areas (.area). 
		  
`island_database.mat` contains the 437 emerged islands and 227 shallow reefs obtained after combining islands (detected from GSHHS coastline and GEBCO bathymetry) within the MODIS land pixels.  
For each island, the database contains its position (.lon, .lat), island area (.Iarea, =0 for shallow reefs), reef area (.Rarea), reef ID (.reefID, unique identifier corresponding to a set of connected pixels within the shallow mask, useful because several islands can share the same reef), and .is_NunnDB that indicates whether the island is part of the Nunn et al. (2016) database (https://doi.org/10.1186/s40562-016-0041-8).  
The island database is also available as a csv file at https://doi.org/10.5281/zenodo.6416130.


### Description of outputs generated by `start_IME_toolbox`

IME_and_REF_regions.jpg: example of IME detection for the Fiji/Tonga region in August  
mean_IME.jpg: map of IMEs detected from a mean chlorophyll map (Fig. 1 in paper)  
Marquesas_insert.jpg: idem for the Marquesas Islands in May (Fig. 1F in paper)  
IME_climato.mat: IME detection outputs obtained from climatological chlorophyll maps. The full IME database (including primary production and PHYSAT outputs) is available as a NetCDF file at https://doi.org/10.5281/zenodo.6416130.  


* * *

### Reference

Please refer this paper when using the toolbox (available upon request):  

Messié, M., A. Petrenko, A. Doglioli, E. Martinez, and S. Alvain (2022). Basin-scale biogeochemical and ecological impacts of islands in the tropical Pacific Ocean. Nature Geoscience, submitted.

* * *

### Contact

monique@mbari.org

Do not hesitate to contact me if you cannot run the code in `start_IME_toolbox`, if you notice bugs, or if you need help implementing the code for a custom application.