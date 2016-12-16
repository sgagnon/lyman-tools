# lyman-tools
Tools for fmri analysis in the Lyman pipeline

This toolbox is intended to provide some helpful scripts for basic fMRI analyses, 
building off the [Lyman](https://github.com/mwaskom/lyman) pipeline. Some code is adapted
from Lyman, and [Ian Ballard](https://github.com/iancballard?tab=repositories).

## Timeseries analyses

These analyses assume that your functional data has been processed with some variant of 
`run_fmri.py -w preproc reg -t -u -reg epi`, such that the unsmoothed raw data is preprocessed, 
and then aligned into epi space (1st run) for each subject. An onset file, 
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