#! /usr/bin/env python
"""
Extract mean parameter estimates from ROI
"""

import os
import sys
import re
import os.path as op
from textwrap import dedent
import argparse
from subprocess import call

from scipy import stats
import scipy as sp
import numpy as np
import pandas as pd
import nibabel as nib

import lyman
from lyman import tools


def main(arglist):
    """Main function for workflow setup and execution."""
    args = parse_args(arglist)

    project = lyman.gather_project_info()
    exp = lyman.gather_experiment_info(args.experiment, args.altmodel)
    subject_list = lyman.determine_subjects(args.subjects)

    # Get group info
    if args.group_info:
        print 'Group counts:'
        lymandir = os.environ["LYMAN_DIR"]
        group_info = pd.read_csv(op.join(lymandir, args.group_info))
        print group_info.group.value_counts()

    # Determine some parameters
    if args.altmodel:
        exp['exp_name'] = "-".join([args.experiment, args.altmodel])
    else:
        exp['exp_name'] = args.experiment

    if args.mni_space:
        exp['regspace'] = 'mni'
        exp['smoothing'] = 'smoothed'
        exp['contrast_exp'] = args.contrast_exp
        exp['group'] = args.group
        exp['threshold'] = args.threshold

    else:
        exp['regspace'] = 'epi'
        exp['smoothing'] = 'unsmoothed'
        exp['threshold'] = args.threshold
        exp['contrast_exp'] = args.contrast_exp

    # Use cope for the input file, to extract parameter estimates (instead of, e.g., z-stats)
    exp['input_file'] = 'cope1.nii.gz'

    # Read in the masks of interest
    print '******* Reading in mask info... ********'
    masks = pd.read_csv(op.join(project['analysis_dir'], exp['exp_name'], 'group', args.masks))
    print masks

    # Conditions to extract
    conditions = exp['condition_names']

    # Now extract the means
    print '******* Extracting mean copes... ********'
    df = extract_means(masks, subject_list, exp, project, conditions)
    # print df.head()

    # Output the data
    if args.group_info: # add in group data if relevant
        df = df.merge(group_info, how='left')

    print '******* Writing out data... ********'
    rois = df.roi.unique()
    for roi in rois:
        print roi
        data = df[df.roi.isin([roi])]
        if isinstance(args.threshold, str):
            filepath = op.join(project['analysis_dir'], exp['exp_name'],
                               'group/roi', 'pe_' + roi + '_' +args.threshold[:-3] + '.csv')
        else:
            filepath = op.join(project['analysis_dir'], exp['exp_name'],
                               'group/roi', 'pe_' + roi + '.csv')

        print 'Filepath = ' + filepath
        data.to_csv(filepath)


def extract_means(masks, subject_list, exp, project, conditions):

    df = pd.DataFrame(columns=('subid', 'roi', 'regspace', 'smoothing', 'mask_vox', 'hemi', 'cond', 'value'))

    # remove nuisance from contrasts
    if 'nuisance' in conditions:
        conditions.remove('nuisance')

    print conditions
    contrast_list = conditions
    contrast_names = conditions

    # Setup counter for printing out info:
    output = True

    for subid in subject_list:
        print subid
        for contrast,label in zip(contrast_list, contrast_names):

            # subject data in ffx folder
            fmri_file = op.join(project['analysis_dir'], exp['exp_name'],
                                subid, 'ffx', exp['regspace'], exp['smoothing'],
                                contrast, exp['input_file'])

            # Read data using nibabel:
            fmri_data = nib.load(fmri_file)
            func_arr = fmri_data.get_data()

            for roi, hemi, roi_type in masks.values:

                # if something under hemi, add to roi name, otherwise just use roi
                if not isinstance(hemi, float):
                    if roi_type == 'anat':
                        mask_name = hemi + '-' + roi
                    else:
                        mask_name = hemi + '-' + roi
                else:
                    mask_name = roi

                # Read in the mask as bool using nibabel:
                if exp['regspace'] == 'mni':
                    # threshold a group level contrast
                    if type(exp['threshold']) in [float, int]:
                        if output:
                            print 'Mask: MNI space, threshold = ' + str(exp['threshold'])

                        mask_file = op.join(project['analysis_dir'],
                                            exp['contrast_exp'], exp['group'],
                                            'mni', roi, 'zstat1.nii.gz')
                        mask_data = nib.load(mask_file)
                        mask_arr = mask_data.get_data()
                        mask_thresh = sp.stats.threshold(mask_arr, threshmin=exp['threshold'])
                        mask_arr = mask_thresh.astype(bool)

                    # use pre-thresholded contrast
                    elif type(exp['threshold']) == bool:
                        if output:
                            print 'Mask: MNI space, corrected threshold'

                        mask_file = op.join(project['analysis_dir'],
                                            exp['contrast_exp'], exp['group'],
                                            'mni', roi, 'zstat1_threshold.nii.gz')
                        mask_data = nib.load(mask_file)
                        mask_arr = mask_data.get_data().astype(bool)

                    # threshold is name of file that's been prethresholded
                    else:
                        if output:
                            print 'Mask: MNI space, specified mask: ' + exp['threshold']

                        mask_file = op.join(project['analysis_dir'],
                                            exp['contrast_exp'], exp['group'],
                                            'mni', roi, exp['threshold'])
                        mask_data = nib.load(mask_file)
                        mask_arr = mask_data.get_data().astype(bool)

                # native space
                else:    

                    # threshold contrast in subjects ffx, epi
                    if exp['threshold']:
                        if output:
                            print 'Mask: Native space, thresholded: ' + exp['threshold']

                        mask_file = op.join(project['analysis_dir'],
                                            exp['contrast_exp'], subid, 'ffx',
                                            'epi', exp['smoothing'], roi, 'zstat1.nii.gz')
                        mask_data = nib.load(mask_file)
                        mask_arr = mask_data.get_data()
                        mask_thresh = sp.stats.threshold(mask_arr, threshmin=exp['threshold'])
                        mask_arr = mask_thresh.astype(bool)

                    # individual-defined anatomical mask
                    else:
                        if output:
                            print 'Mask: Native space, specified mask: ' + mask_name

                        mask_file = op.join(project['data_dir'], subid, 
                                            'masks', mask_name + '.nii.gz')
                        mask_data = nib.load(mask_file)
                        mask_arr = mask_data.get_data().astype(bool)

                num_voxels = mask_arr.sum()

                # Mask the data
                func_masked = func_arr[mask_arr]

                # Save to dataframe
                row = pd.DataFrame([dict(subid = subid, 
                                         roi = roi, 
                                         regspace = exp['regspace'],
                                         smoothing = exp['smoothing'],
                                         mask_vox = num_voxels, 
                                         hemi = hemi,
                                         cond = label, 
                                         value = func_masked.mean()), ])
                df = df.append(row, ignore_index = True)
                output = False
    return df


def parse_args(arglist):
    """Take an arglist and return an argparse Namespace."""
    help = dedent("""
    Example: 
        extract_copes.py -exp ap_memory_raw -masks hipp_masks.csv -group_info subjects_groups.csv

      Usage Details
    -------------
    """)

    parser = tools.parser
    parser.description = help
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument("-experiment", help="experimental paradigm")
    parser.add_argument("-altmodel", help="alternate model to fit")
    parser.add_argument("-masks", help="csv file to load into pandas df")
    parser.add_argument("-mni_space", action="store_true", help="Flag if true, leave out for any anat ROIs")
    parser.add_argument("-contrast_exp", default="None", help="String (regular experiment if no localizer)")
    parser.add_argument("-threshold", help="True (standard thresh), number to specify")
    parser.add_argument("-group", default="group", help="group directory")
    parser.add_argument("-group_info", help="csv file in scripts dir w/subid -> group mapping")

    return parser.parse_args(arglist)

if __name__ == "__main__":
    main(sys.argv[1:])