# Hutson_GLRCM_2023

This reposity accompanies the manuscript "Testing the Sensitivity of a WRF-based Great Lakes Regional Climate Model to Cumulus Parameterization and Spectral Nudging" by Abby Hutson (CIGLR), Ayumi Fujisaki-Manome (CIGLR/CLASP), and Brent Lofgren ("What Makes Climate Tick?"), submitted to the American Meteorological Society Journal of Hydrometeorology. 

The model in Hutson et al. (submitted) is version 4.3 of the Weather and Research Forecasting model (https://www2.mmm.ucar.edu/wrf/users/download/get_source.html) and is used to dynamically downscale ERA-Interim (ECMWF; now discontinued) global data into a regional climate model for the Great Lakes Region. The model results used in Hutson et al. are too big to be able to store or transfer easily, so this repository contains the information/methods needed to recreate said results. 

**Namelists**: this folder contains different namelist.input files (used to run wrf.exe in the WRF directory), one for each RCM version described in Table 1 and Table 2 in Hutson et al. (submitted). The extension after "namelist.input." for each file describes the RCM version for which it was used. 

**namelist.wps**: the namelist used for all WRF preprocessing (i.e., what is used to run geogrid.exe, ungrib.exe, and metgrid.exe in the WPS directory). All geography (geog) data was downloaded from https://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html, including the optional "lake depth" data required to run the one-dimensional WRF-Lake model (sf_lake=1). 

**module_sf_lake.F.4degC**: The default module_sf_lake.F physics module (the module used to run the one-dimesional WRF-Lake model while running WRF) is edited in Hutson et al. (submitted) so that the lake surface temperature for all five Great Lakes is equal to 4 degC (277 K). This was done to avoid anomalously cold LSTs (and anomalous ice all year long) that resulted from WRF's default LST initialization. This is only recommended for simulations initializing in the early winter (i.e., January 1st). The edits can be found on lines 5268-5274. For these edits to be used in WRF, this file must be copied to the "module_sf_lake.F" file in the WRF/phys/ directory before compiling. 

**METGRID.TBL**: The METGRID.TBL specifications used for all simulations. This file contains information on how ERA-I's land-sea mask, sea-surface temperatures, and skin temperatures were interpolated to the WRF grid. 

**interp.f** and **interpolation.py**: This program contains the remapping method used in Accadia et al. (2002). In Hutson et al. (submitted), method is used to compare precipitation output from WRF to other datasets with different resolutions. The Accadia method should always be used to remap precipitation from the finer grid to the coarser grid. Before first use, interp.f must be compiled with "f2py". Specifically, run "f2py -c -m interp interp.f" at the command line. More information is detailed in the doc string in interpolation.py. This program is provided by Russel Manser (Texas Tech University). 

**rcm_functions.py**: Contains different functions used to post-process the RCM output. Specifically, rcm_functions contains details used for both interpolation methods (Accadia remapping and Bilinear interpolation, the latter of which is used to interpolate a coarse grid to a more fine grid). 

**coodinateSystems.py**: Contains functions needed to conduct the linear interpolation contained within rcm_functions.py (specifically, the function used to convert to and from Earth-Center-Earth-Fixed (ECEF) coordinates). 
