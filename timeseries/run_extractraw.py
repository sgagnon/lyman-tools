#! /usr/bin/env python
"""
Script to extract raw timeseries, time locked to event onsets
Make sure timeseries are registered before running this 
(e.g., run_fmri.py -e ap_memory_raw -w reg -timeseries -regspace epi -n 20)

EXAMPLE: 
python run_extractraw.py -extract_info AP_mvpa_raw -mask_type mask -mask_name lh-hippocampus-tail
"""

import sys
import argparse
import shutil
import imp
import re
from textwrap import dedent
import time

import os.path as op

import nibabel as nib
import pandas as pd
import numpy as np
import scipy as sp
from nilearn.masking import compute_epi_mask
from nilearn.input_data import NiftiMasker
from nilearn import image

import lyman
from lyman import tools

# Define some functions

def extract_mean(d_m, subid, exp, args, save_trials=False):

    onsets = pd.read_csv(exp['onsetfile'].format(subid=subid)) # load in onset/event data
    onsets = onsets[onsets.condition.isin(exp['cond_list'])].reset_index() # subset for conditions

    if args.mask_type == 'func':
        # Specify path for func mask (run 1 if registered)
        mask_path = exp['func_maskfile'].format(subid=subid, run_id=1)
    elif args.mask_type == 'mask':
        mask_path = exp['maskfile'].format(subid=subid, mask_name=args.mask_name)

    # Define mask for timeseries data (to convert 4D to 2D)
    # apply smoothing and standardization of features here if necessary
    func_masker = NiftiMasker(mask_img=mask_path,
                              smoothing_fwhm=exp['smoothing_fwhm'],
                              standardize=exp['standardize_feat'])
    if exp['standardize_feat']:
        print 'Scaling each feature'

    num_voxels = np.sum(nib.load(mask_path).get_data(), axis=None).astype(int)

    run_list = list(set(onsets.run))
    for run in run_list:
        print 'run: ' + str(run)
        run = int(run)

        ts_path = exp['tsfile'].format(subid=subid, run_id=run)
        func_masked = func_masker.fit_transform(ts_path)
        n_trs = func_masked.shape[0]

        # Compute mean for each TR
        mean_ts = np.mean(func_masked, axis=1)

        if exp['standardize_roi']:
            print 'Scaling ROI mean across time'
            mean_ts = sp.stats.mstats.zscore(mean_ts)

        if exp['percentsig_roi']:
            print 'Converting to %sig change using ts mean'
            mean_oftimecourse = np.mean(mean_ts)
            mean_ts = (mean_ts/mean_oftimecourse * 100) - 100

        

        for tr_shift in exp['tr_shift']:
            
            # Get y TRs and labels
            #reset index so matches for adding in mean activity
            run_events = onsets[onsets.run == run].reset_index()

            print 'tr_shift: ' + str(tr_shift)

            # If calculating %sig change rel baseline, get baseline for each trial
            baseline_trs_run = np.ceil((run_events.onset)/exp['tr']).astype(int)
            baseline_trs_run = baseline_trs_run.replace(to_replace=0, value=1) # to deal w/indexing later on, for onsets of 0
            baseline_trs_ind = baseline_trs_run-1 # back for indexing w/0 as 1st TR
            # print baseline_trs_ind.shape
            # print baseline_trs_ind

            # Figure out the TRs, timeshifting
            ev_trs_run = np.ceil((run_events.onset + tr_shift)/exp['tr']).astype(int)
            ev_trs_run = ev_trs_run.replace(to_replace=0, value=1) # to deal w/indexing later on, for onsets of 0
            ev_trs_ind = ev_trs_run-1 # back for indexing w/0 as 1st TR
            # print ev_trs_ind.shape
            # print ev_trs_ind

            # only take TRs that are actually in the timeseries
            if max(ev_trs_ind) >= n_trs:
                print 'Timeseries is ' + str(func_masked.shape[0]) + ' trs...' + \
                      'but you want TR: ' + str(max(ev_trs_ind + 1))
                in_ts_ind = ev_trs_ind < n_trs
                print in_ts_ind

                ev_trs_ind = ev_trs_ind[in_ts_ind]
                baseline_trs_ind = baseline_trs_ind[in_ts_ind]
                run_events = run_events.loc[in_ts_ind]

            # print ev_trs_run

            # this is old code, where taking the mean over voxels, but replaced with above code
            # the DF functionality filled in NaNs for missing timepoints
            # ts_vals = pd.DataFrame(func_masked[ev_trs_ind, :]).mean(axis=1)

            # Convert to %sig change rel to basleine
            if exp['percentsig_roi_relbaseline']:
                print 'Converting to %sig change relative to baseline, using 0 timepoint'
                ev_values = (mean_ts[ev_trs_ind]/mean_ts[baseline_trs_ind] * 100) - 100
                print ev_values
            else:
                ev_values = mean_ts[ev_trs_ind]
                print ev_values

            # if TRs go beyond what's in timeseries, no values (run_Events is pruned)
            run_events.loc[:, 'mean_activity'] = pd.Series(ev_values,
                                                           index=run_events.index)
            run_events.loc[:, 'subid'] = subid
            run_events.loc[:, 'mask'] = args.mask_name
            run_events.loc[:, 'time'] = tr_shift
            d_m = d_m.append(run_events, ignore_index=True)


    return d_m


def main(arglist):

    args = parse_args(arglist)

    # Import design/experiment information
    exp_file = op.join('/share/awagner/sgagnon/scripts/lyman-tools/timeseries/extract_info',
                       args.extract_info + ".py")
    exp = imp.load_source(args.extract_info, exp_file)

    def keep(k):
        return not re.match("__.*__", k)

    exp = {k: v for k, v in exp.__dict__.items() if keep(k)}

    subjects = pd.Series(lyman.determine_subjects(), name="subj")

    # Create dataframe, add each subjects data
    d_m = pd.DataFrame()
    for subid in subjects:
        print subid + '-----------------------------------------------'

        # Figure out which conditions to examine
        cond_list = list(pd.read_csv(exp['onsetfile'].format(subid=subid)).condition.unique())

        if 'nuisance' in cond_list:
            cond_list.remove('nuisance') # remove nuisance trials
        exp['cond_list'] = cond_list

        d_m = extract_mean(d_m, subid, exp, args)

    # save csv of data
    filepath = op.join(exp['expdir'], 'group/roi', 'extractraw_' + args.extract_info +
                       '_' + args.mask_name + '.csv')
    d_m.to_csv(filepath, index=False)

    # Integrate over specified time points and save
    dm_wide = d_m.pivot_table(index=['subid', 'mask', 'run', 'onset', 'condition'],
                              values='mean_activity', columns=['time'],
                              fill_value=np.nan).reset_index()
    dm_wide['integrated_activity'] = sp.trapz(dm_wide[exp['tr_integrate']], axis=1)
    filepath = op.join(exp['expdir'], 'group/roi', 'extractraw_integrated_' + args.extract_info +
                       '_' + args.mask_name + '.csv')
    dm_wide.to_csv(filepath, index=False)

def parse_args(arglist):
    """Take an arglist and return an argparse Namespace."""
    help = dedent("""
    Run FIR model on subject data
    """)
    parser = tools.parser
    parser.description = help
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument("-extract_info", help="info for experiment to extract")
    parser.add_argument("-mask_type", help="mask or func?")
    parser.add_argument("-mask_name", help="name of mask in sub's mask directory")
    return parser.parse_args(arglist)

if __name__ == "__main__":
    main(sys.argv[1:])
