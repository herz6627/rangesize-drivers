# rangesize-drivers

# ADD publication info

This repository contains all code associated with Herzog XXXX and associated file structure. All code is within the ./R/ directory. The output, figs, and data directories are empty place holder directories to allow users to easily see how the data was structured when writing the code. Files for the ./data/ directory can be found at the following links:

## APT011
    Nippert, J. 2024. APT01 Daily precipitation amounts measured at multiple sites across konza prairie ver 21. Environmental Data Initiative. https://doi.org/10.6073/pasta/beceae506b1b975e418629c6c4c94813 (Accessed 2024-02-24).
## AWE012
    Nippert, J. 2024. AWE01 Meteorological data from the konza prairie headquarters weather station ver 24. Environmental Data Initiative. https://doi.org/10.6073/pasta/a267f9b0995f6fa91340ba5886ee2273 (Accessed 2024-02-24).
## CGR023
    Joern, A. 2024. CGR02 Sweep sampling of Grasshoppers on Konza Prairie LTER watersheds ver 23. Environmental Data Initiative. https://doi.org/10.6073/pasta/e349433ba41c72fe98b6347ec9e7dd91 (Accessed 2024-02-24).
## KFH011
    Blair, J. and P. O'Neal. 2024. KFH01 Konza prairie fire history ver 22. Environmental Data Initiative. https://doi.org/10.6073/pasta/0db54f512bba21337d0de4b4ec0951e9 (Accessed 2024-02-24).
## PBG011
    Blair, J. 2024. PBG01 Plant species composition in the Patch Burning-grazing Experiment at Konza Prairie ver 17. Environmental Data Initiative. https://doi.org/10.6073/pasta/cffd16d53830b097a5a229d39990ff59 (Accessed 2024-02-24).
## PPS011
    Nippert, J., J. Blair, and J. Taylor. 2024. PPS01 Konza prairie plant species list ver 4. Environmental Data Initiative. https://doi.org/10.6073/pasta/50f594d82bc0b385406662fb5dcba59f (Accessed 2024-07-31)
## PVC021
    Hartnett, D., S. Collins, and Z. Ratajczak. 2024. PVC02 Plant species composition on selected watersheds at Konza Prairie ver 24. Environmental Data Initiative. https://doi.org/10.6073/pasta/c2dde97352fb0e25ab749765967997b9 (Accessed 2024-02-24).
## grasshopper diets
Welti, E.A.R., Qiu, F., Tetreault, H.M., Ungerer, M., Blair, J., Joern, A., 2019. Fire, grazing and climate shape plant–grasshopper interactions in a tallgrass prairie. Functional Ecology 33, 735–745. https://doi.org/10.1111/1365-2435.13272

Campbell, J.B., Arnett, W.H., Lambley, J.D., Jantz, O.K., Knutson, H., 1974. Grasshoppers (Acrididae) of the Flint Hills native tallgrass prairie in Kansas. Agricultural Experiment Station, Kansas State University 19, 73–145.

Mulkern, G.B., 1969. Behavioral Influences on Food Selection in Grasshoppers (orthoptera: Acrididae). Entomologia Experimentalis et Applicata 12, 509–523. https://doi.org/10.1111/j.1570-7458.1969.tb02549.x

__________________________________________________________________________

The R scripts are numbered with assending numeric prefixes to indicate order that scripts should be run in. For example, 1_clean-PVC-data.R should be run before all other R scripts, followed by 2_range-size-estimation.R, etc..

Included R files are described below:

## 1_clean-PVC-data.R
Cleans and formats datasets.
## 2_range-size-estimation.R
Download and clean occurrence records, then estimate range size. Depends on functions in 2_func_mcp-sh.R
Specific iDigBio and GBIF records used to estimate range size can be found at XXXX
## 2_func_mcp-sh.R
Custom code largely based on the XXX R package. 
## 3_gh-range-cor.R
Estimate number of grasshopper consumers per plant species and correllate number of grasshopper consumers by plant range size.
## 4_analyses-sep-by-dat.R
Create linear models for colonization, change in cover, and extirpation using environmental drivers for each Konza plant species. Reliant on function in 4_functions.R
## 4_functions.R
Helper functions to estimate environmental effects on plant population dynamics.
## 5_additional-analyses.R
Using coefficients estimated in 4_ and range sizes estimated in 2_, find any linear trends in driver effect size and range size. Dependent on functions in 5_functions.R
## 5_functions.R
Helper functions to clean data, run linear models and make figures.
