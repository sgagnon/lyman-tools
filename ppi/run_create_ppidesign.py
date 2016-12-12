
#! /usr/bin/env python
"""
Script to beging PPI analysis
For each subject, extract a timeseries from specified ROIs. Multiply by 
a specified contrast

Takes one argument, the design name (corresponding to .py file in ppi/design dir)
EXAMPLE: 
python run_create_ppidesign.py AP_mvpa_raw_hipp
"""

import os
import sys
import imp
import os.path as op
import numpy as np
import glob
from os.path import abspath
import csv
import re

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats, optimize
from pandas import DataFrame, Series
from moss import glm
import seaborn as sns
import random as rd
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import scipy.stats

from IPython.parallel import Client
from IPython.display import Image
import multiprocessing

import nibabel as nib
from nipype.pipeline.engine import Node, MapNode, Workflow
from nipype.interfaces.io import DataGrabber, DataFinder, DataSink
from nipype.interfaces import fsl
from nipype.interfaces.fsl import ImageMeants
from nipype.interfaces.fsl import ImageStats

def vector_rejection(a, b):
    return a - (np.dot(a, b)/np.dot(b, b) * b)

def extract_roi(in_tuple):
    sub, exp_name, run, mask = in_tuple
    
    sub_path = op.join(paths['analy_dir'].format(exp=exp_name), sub, 'preproc', 'run_'+run)
    ts_path = op.join(paths['analy_dir'].format(exp=exp_name), sub, 'reg/epi/unsmoothed', 'run_'+run)

    #make sure to get coregistered preproc data
    preproc_data = op.join(ts_path, 'timeseries_xfm.nii.gz')

    mask_dir = op.join(paths['data_dir'], sub, 'masks')
    out_dir = mask_dir + '/extractions/'
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    mask_file = op.join(mask_dir, mask + '.nii.gz')
    out_f = out_dir + ('').join(map(str, in_tuple))+ '.txt'

    if os.path.exists(sub_path):# and not os.path.exists(out_f):
        meants = ImageMeants(in_file=preproc_data, eig=True, order=1,
                                mask=mask_file, out_file=out_f)
        meants.run()
        
def extract_roi_prob(in_tuple):
    
    sub, exp_name, run, mask = in_tuple
    
    sub_path = op.join(paths['analy_dir'].format(exp=exp_name), sub, 'preproc', 'run_' + run)
    ts_path = op.join(paths['analy_dir'].format(exp=exp_name), sub, 'reg/epi/unsmoothed', 'run_'+run)

    #make sure to get coregistered preproc data
    preproc_data = op.join(ts_path, 'timeseries_xfm.nii.gz')

    mask_dir = op.join(paths['data_dir'], sub, 'masks')
    out_dir = mask_dir + '/extractions/'

    prob_file = mask_dir + exp_name + '_' + mask + '_func_space.nii.gz'
    mask_file = op.join(mask_dir, mask + '.nii.gz')
    out_f = out_dir + ('').join(map(str, in_tuple)) + '.txt'
    tmp_out = mask_dir + sub + exp_name + run + '.nii.gz'

    if os.path.exists(sub_path):# and not os.path.exists(out_f):
        cmd = ['fslmaths',preproc_data, '-mul', prob_file,tmp_out]
        cmd = ' '.join(cmd)
        os.system(cmd)
        
        meants = ImageMeants(in_file = tmp_out, eig = True, order = 1, 
                                mask = mask_file, out_file = out_f)
        meants.run()
        os.remove(tmp_out)

def write_design(in_tuple):
    
    sub,exp_name,mask,design_file,contrast = in_tuple
    
    #hrf params
    hrf = getattr(glm,'GammaDifferenceHRF')
    tr = 2.
    hrf = hrf(tr = tr)
    
    #set up filenames
    design_dir = op.join(paths['data_dir'], sub, 'design')
    out_f = op.join(design_dir, 
                    'ppi_regressors_{exp}_{mask}_{contrast}.csv').format(contrast=contrast, 
                                                                            exp=exp_name, mask=mask) #out file
    mask_dir = op.join(paths['data_dir'], sub, 'masks', 'extractions')

    #load design data for this subject
    design_data = pd.read_csv(op.join(design_dir, design_file))

#     #load in pre-existing noise regressors
#     reg_file = design_dir + 'noise_regressors_' + exp_name + '.csv'
#     regressors = pd.read_csv(reg_file)
    regressors = pd.DataFrame()

    #initialize vars to fill
    convolved_ev = []
    ts = []
    run_list = []
    for run in runs:
        if (sub == 'ap155') & (int(run) == 6):
            print 'skipping run 6 for sub ap155'
        else:
            sub_file = op.join(paths['analy_dir'].format(exp=exp_name), sub, 
                            'preproc/run_' + str(run), 'unsmoothed_timeseries.nii.gz')

            if os.path.exists(sub_file):
                ntp = nib.load(sub_file).shape[-1] #get number of time points
                design = design_data[design_data['run']==int(run)]

                model = glm.DesignMatrix(design = design, tr = tr, ntp = ntp, hrf_model = hrf, hpf_cutoff = 128)
                psy_reg = model.design_matrix[contrast].values
            
                ##centre convolved ev (see fsl docs), and append to data -- do this within run?!
                diff = max(psy_reg) - (max(psy_reg) - min(psy_reg))/2.0
                psy_reg = psy_reg - diff
                convolved_ev.extend(psy_reg) #get timeseries for regressor of interest

                #load ts data
                fid = (sub,exp_name,run,mask)
                mask_f = op.join(mask_dir, ('').join(map(str,fid))+ '.txt')
                roi_ts = np.loadtxt(mask_f)
                roi_ts = roi_ts - np.mean(roi_ts) #mean center
                ts.extend(roi_ts)
            
                # add this on if dont have preexisting noise regressors
                run_list.extend([int(run)] * len(roi_ts))

    #update regressors dataframe
    ts = scipy.stats.zscore(ts) #add ts to the regressors DF
    #orthogonalize noise regressors to speed up computation
#         ts = vector_rejection(ts,regressors['ventricles'])
#         ts = vector_rejection(ts,regressors['wm'])
    regressors[mask] = ts
    regressors['interaction'] = convolved_ev * ts #interaction regressor
    regressors['run'] = run_list
    
    #write output
    regressors.to_csv(out_f, header=True,index = False, 
                    columns = [mask,'interaction','run'])
    print sub, mask

def main(design_name):
    
    # Import design/experiment information
    exp_file = op.join('/share/awagner/sgagnon/scripts/lyman-tools/ppi/design',
                       design_name + ".py")
    exp = imp.load_source(design_name, exp_file)

    def keep(k):
        return not re.match("__.*__", k)

    exp = {k: v for k, v in exp.__dict__.items() if keep(k)}

    # Extract timeseries
    sub_list = list(np.loadtxt(exp['subj_file'], 'string'))

    os.chdir(exp['home_dir'])
    runs = map(str, range(1,exp['num_runs']+1))

    global paths
    paths = dict(home_dir = exp['home_dir'], 
                data_dir = exp['data_dir'],
                analy_dir = exp['analy_dir'])
    
    in_tuples = []
    for sub in sub_list:
        for exp_name in exp['exps']:
            for run in runs:
                if (sub == 'ap155') & (int(run) > 5):
                    print 'Skipping run 6 for ap165!'
                else:
                    for mask in exp['masks']:
                        in_tuples.append((sub,exp_name,run,mask))

    pool = multiprocessing.Pool(processes = exp['n_cores'])
    pool.map(extract_roi,in_tuples)
    pool.terminate()
    pool.join()

    # Write PPI design

    in_tuples = []
    for sub in sub_list:
        for exp_name in exp['exps']:
            for mask in exp['masks']:
                in_tuples.append((sub,exp_name,mask,exp['design_file'],exp['contrast']))

    ## Run in parallel
    pool = multiprocessing.Pool(processes = exp['n_cores'])
    pool.map(write_design,in_tuples)
    pool.terminate()
    pool.join()

    ## Run in serial (for trouble shooting)
    # for tuple in in_tuples:
    #     print tuple
    #     write_design(tuple)


if __name__ == "__main__":
    main(sys.argv[1])
