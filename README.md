# Surfmelt_mapping_tools
Google Earth Engine and Matlab codes to automatically map surface meltwater across the entire Antarctic Ice Sheet.

This code relates to the following paper preprint, which describes the continent-wide dataset and methodology:

Tuckett, P., Ely, J., Sole, A., Livingstone, S., Jones, J., Lea, J., & Gilbert, E. (in review). 
Continent-scale mapping reveals increasing sensitivity of East Antarctica to meltwater ponding.

The continent-wide dataset described here was produced using a scaled up version of the method descibed in Tuckett et al., 2021.
See Tuckett et al. (2021) for an in-depth description of the mapping method, image visibility assessments, and post-processing steps:

Tuckett, P. A., Ely, J. C., Sole, A. J., Lea, J. M., Livingstone, S. J., Jones, J. M., & van Wessem, J. M. (2021). 
Automated mapping of the seasonal evolution of surface meltwater and its links to climate on the Amery Ice Shelf, Antarctica. The Cryosphere. 1-35.

Please contact Pete Tuckett (University of York) for more information - pete.tuckett@york.ac.uk, petetuckett47@gmail.com


SOFTWARE REQUIREMENTS:

1) Google Earth Engine (requires a user account) - copy the contents of each .txt file into a new GEE script (https://code.earthengine.google.com).
2) Matlab (Used for post-processing stages. Code could be adapted for use in other platforms)
3) Use of a High Performace Computing (HPC) platform is beneficial, but not essential, for the post-processing Matlab script

Our continent-wide dataset was produced in 2022. All scripts have been updated to reflect system changes within GEE, but please be aware that the image collections
within GEE are constantly updated, and minor adaptions to scripts lines relating to image file pathways and image band names may be required.


CODE DESCRIPTIONS & INSTRUCTIONS:

Three key scripts are required to run our Continent-wide mapping methodology. A shapefile specfying the area to be mapped over is also required. Our continent-wide
dataset was created by mapping over one quarter of Antarctica at a time, to avoid memory limits within GEE. See below a description of the steps taken to map surface meltwater
over any given pre-defined area of Antarctica. Meltwater is mapped over a user-specified time period and temporal resolution.

1) Generate a shapefile polygon to define the area you intend to map over. If the area is larger than 100 x 100 km2, split the polygon into multiple tiles. See Amery.shp as an example for mapping
over the Amery region of Antarctica (Tuckett et al., 2021). Ensure that each tile has a unique Tile_ID number.
2) Import this shapefile into GEE as an asset.
3) Generate an ice mask for your given study area (and selected time period) by running the GEE script: 'Create_ice_mask'. Ensure that the assetpath is set to match the location of the shapefile asset
imported in step 2. This script will create an ice mask and save it as a new asset within GEE.
4) Run the lake mapping script: 'Automated_surface_meltwater_mapping'. Ensure that 'icemaskPath' is set to match the location of the ice mask within assets. Specify the date range
and temporal resolution (time window length) under 'SPECIFY CONFIG VARIABLES'. This script maps surface water for the given region and date range, and attaches visibility metadata.
It exports a GEOJSON file to Google Drive - set the output location and file names at the end of the script.
5) Download the GEOJSON file from Google Drive, and save it with your documents.
6) Generate mapped lake shapefiles for each time window by running the Matlab script: 'Post_processing_script_Matlab.m'. The script requires some functions that may need to be downloaded 
from Mathworks. File pathways/directory structure will need editing by the user. The script is setup to create three versions of shapefile outputs:
(i) Individual_tiles (raw shapefiles for each individual tile, per time window)
(ii) Combined (Merges tiles, to create a single shapefile per time window)
(iii) Filtered_shapefiles (Lake polygons are unioned over tile boundaries, lake centroids are recalculated, an area threshold is applied and metadata is assigned, including visibility scores)


