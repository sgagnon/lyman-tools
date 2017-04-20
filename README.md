# lyman-tools
Tools for fmri analysis in the Lyman pipeline

This toolbox is intended to provide some helpful scripts for basic fMRI analyses, 
building off the [Lyman](https://github.com/mwaskom/lyman) pipeline. Some code is adapted
from Lyman, and [Ian Ballard](https://github.com/iancballard?tab=repositories).

## ROI Tools

### Create a spherical ROI from peak activation

`create_sphere_frompeak.py`: Take the peak MNI coordinate from group `zstat1_localmax.csv`, transform into avg152 coordinate space, create a sphere of `sphere_rad`mm around the peak, and then mask with thresholded map. Output both the sphere and masked sphere niftis. This is useful to run before extracting parameter estimates from an ROI of interest to visualize the pattern of effect. 

## Timeseries analyses

These analyses assume that your functional data has been processed with some variant of 
`run_fmri.py -w preproc reg -t -u -reg epi`, such that the unsmoothed raw data is preprocessed, 
and then coregistered into epi space (1st run) for each subject. An onset file, 
containing info about condition labels and onsets, is also necessary.

### FIR model

`run_fir.py`: Extract timeseries for each subject, and use Nitime's [EventRelatedAnalyzer](http://nipy.org/nitime/api/generated/nitime.analysis.event_related.html) 
to calculate the FIR event-related estimated of the HRFs for events of interest.


### Extract timeseries, locked to event onset

`run__extractraw.py`: Extract "raw" (but preprocessed/realigned) timeseries, relative to 
event onsets. Integration (using the [composite trapezoidal rule](https://docs.scipy.org/doc/scipy-0.10.1/reference/generated/scipy.integrate.trapz.html)) 
is also performed over a specified window. This analysis is done on a trial-by-trial basis, 
so can be helpful for trial-wise correlations between regions or multivariate measures.

## PPI analyses

### Generate design files for PPI analysis

`run_create_ppidesign.py`: Extract the coregistered preprocessed timeseries from a given 
ROI, take the principal eigenvariate, and z-score within run. Then convolve the 
psychological regressor with the HRF, center, and multiply with the timeseries to produce 
the interaction. Save out the timeseries and interaction as regressors, to be input as 
a `regressor_file` in the Lyman model parameters (in the experiment file).
