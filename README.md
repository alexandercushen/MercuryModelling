# MercuryModelling
My work using SWMF to study mercury. All SWMF runfiles and plotting utilities are stored here.

-------------

Versions: contains subdirectories of the different versions I have worked on developing.

nightside_v1: Started in June 2024. Starting from the ground-up to build a stable, MHD-AEPIC implementation focused on the nightside. Previous development tests showed that the code is unstable, even without PC active; likely to do with boundary conditions or grid discontinuities at the pole.  

multifluid-sw_v1: Multifluid solar wind with seperate H+ and He2+ fluids.

-------------

Plotting: jupyter notebooks used for analysis from 2024 - late 2025. After this time, most development is done in Sypder, saved in "analysis". 

These files contain most of the plotting routines I have developed, as well as the post-processing scripts needed to translate .plt files into portable pickle files, to speed up plotting rates.

-------------

Analysis: cleaned up analysis tools used after post-processing, to generate plots.
