#! /usr/bin/env python
"""
This script finds clusters in group FFX (MNI space, smoothed), and then outputs 
the peaks within a specified ROI as a csv file. 
"""

import numpy as np
import glob
import os
import os.path as op
from scipy import stats
import nibabel as nib
import pandas as pd

from moss import locator

from nipype import IdentityInterface, Node, Workflow
from nipype.interfaces import fsl, freesurfer
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as niu           # Data i/o

z_thresh = 2.3
grf_pthresh = 0.05
peak_dist = 30

exp = 'localizer-1respmod'
contrast = 'eye-hand'

hemi = 'lh'
rois = ['Angular G', 
        'Sup Par Lobule',
        'Supramarg G, ant',
        'Supramarg G, post',
        'Lat Occ Ctx, sup']

subjects = np.loadtxt('/Volumes/group/awagner/sgagnon/RM/scripts/subjects.txt', str)

filedir = '/Volumes/group/awagner/sgagnon/RM/analysis/{exp}/{{subid}}/ffx/mni/smoothed/{contrast}'
basedir = filedir.format(exp=exp,contrast=contrast)
peaks_outdir = '/Volumes/group/awagner/sgagnon/RM/analysis/'+exp+'/peaks_'+ contrast
outdir = op.join('/Volumes/group/awagner/sgagnon/RM/analysis',
                 exp, 'maxclusters_' + contrast + '.csv')

#####################
# Set up nodes
#####################
infosource = Node(IdentityInterface(fields=['subid']),
                  name="infosource")
infosource.iterables = [('subid', subjects)]

templates = {'func': op.join(basedir, 'zstat1.nii.gz'),
             'mask': op.join(basedir, 'mask.nii.gz')}
selectfiles = Node(nio.SelectFiles(templates),
                   name="selectfiles")

smoothest = Node(fsl.SmoothEstimate(), name="smoothest")

cluster = Node(fsl.Cluster(threshold=z_thresh,
                           pthreshold=grf_pthresh,
                           peak_distance=peak_dist,
                           use_mm=True), name="cluster")

cluster_nogrf = Node(fsl.Cluster(threshold=z_thresh,
                                 peak_distance=peak_dist,
                                 use_mm=True,
                                 out_localmax_txt_file=True), name="cluster")

datasink = Node(nio.DataSink(), name='sinker')
datasink.inputs.base_directory = peaks_outdir

#####################
# Construct workflow, and run
#####################
wf = Workflow(name=exp+'_' + contrast)

wf.connect([(infosource, selectfiles, [('subid', 'subid')]),
            (selectfiles, cluster_nogrf, [('func', 'in_file')]),
            (infosource, datasink, [('subid', 'container')]),
            (cluster_nogrf, datasink, [('localmax_txt_file', 'peaks')])
            ])

wf.run()


#####################
# Just extract ROIs
#####################

def reformat_cluster_table(cluster_dir, rois, hemi):
    """Add some info to an FSL cluster file and format it properly."""
    df = pd.read_table(cluster_dir, delimiter="\t")
    df = df[["Cluster Index", "Value", "x", "y", "z"]]
    df.columns = ["Cluster", "Value", "x", "y", "z"]
    df.index.name = "Peak"

    # Find out where the peaks most likely are
    if len(df):
        coords = df[["x", "y", "z"]].values
        loc_df = locator.locate_peaks(coords)
        df = pd.concat([df, loc_df], axis=1)
        mni_coords = locator.vox_to_mni(coords).T
        for i, ax in enumerate(["x", "y", "z"]):
            df[ax] = mni_coords[i]
    
    df_trim = df.loc[df['MaxProb Region'].isin(rois)]
    
    if hemi == 'lh':
        df_trim = df_trim.loc[df_trim.x < 0]
    else:
        df_trim = df_trim.loc[df_trim.x > 0]

    return df_trim.reset_index()
    


df = pd.DataFrame()
for subid in subjects:
    print subid

    cluster_dir = '/Volumes/group/awagner/sgagnon/RM/analysis/'+ exp +'/peaks_' + contrast+'/{subid}/peaks/_subid_{subid}/zstat1_localmax.txt'.format(subid=subid)

    # reformat clusters (note that xyz from FSL is in diff voxel coordinates, not MNI)
    df_sub = reformat_cluster_table(cluster_dir, rois, hemi)
    df_sub['subid'] = subid
    df_sub['exp'] = exp
    df_sub['contrast'] = contrast
    
    df = df.append(df_sub)


df.to_csv(outdir)