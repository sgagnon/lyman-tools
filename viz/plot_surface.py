#! /usr/bin/env python
"""
This script includes a bunch of functions for surface plotting with pysurfer
Make sure you run %gui qt and %matplotlib inline in a nb before running
"""

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from surfer import Brain, io
import surfer as sf
from surfer import utils, io

import numpy as np
import glob
import os.path as op
from scipy import stats
import nibabel as nib
import pandas as pd
import time

import lyman
from lyman.tools.plotting import multi_panel_brain_figure, crop, add_colorbars

def calculate_sat_point(template, contrast, sign, sig_to_z=True, subj=None):
    """Calculate the point at which the colormap should saturate."""
    data = []
    for hemi in ["lh", "rh"]:
        hemi_file = template.format(contrast=contrast, subj=subj, hemi=hemi)
        hemi_data = nib.load(hemi_file).get_data()

        if sig_to_z:
            stat_sign = np.sign(hemi_data)
            p_data = 10 ** -np.abs(hemi_data)
            z_data = stats.norm.ppf(p_data)
            z_data[np.sign(z_data) != stat_sign] *= -1
            hemi_data = z_data

        data.append(hemi_data)

    data = np.concatenate(data)


    if sign == "pos":
        z_max = max(3.71, np.percentile(data, 98))
    elif sign == "neg":
        z_max = max(3.71, np.percentile(-data, 98))
    elif sign == "abs":
        z_max = max(3.71, np.percentile(np.abs(data), 98))
    return z_max


def add_mask_overlay(b, mask_file):
    """Gray-out vertices outside of the common-space mask."""
    mask_data = nib.load(mask_file).get_data()

    # Plot the mask
    mask_data = np.logical_not(mask_data.astype(bool)).squeeze()
    if mask_data.any():
        b.add_data(mask_data, min=0, max=10, thresh=.5,
                   colormap="bone", alpha=.6, colorbar=False)


def add_stat_overlay(b, stat_file, thresh, max, sign, sig_to_z=False, color="Reds_r", alpha=1, output=False, colorbar=False):
    """Plot a surface-encoded statistical overlay."""
    stat_data = nib.load(stat_file).get_data()

    # Possibly convert -log10(p) images to z stats
    if sig_to_z:
        stat_sign = np.sign(stat_data)
        p_data = 10 ** -np.abs(stat_data)
        z_data = stats.norm.ppf(p_data)
        z_data[np.sign(z_data) != stat_sign] *= -1
        stat_data = z_data

    # Plot the statistical data
    stat_data = stat_data.squeeze()
    if sign in ["pos", "abs"] and (stat_data > thresh).any():
        stat_data[stat_data<0] = 0
        b.add_data(stat_data, thresh, max, thresh,
                   colormap=color, colorbar=colorbar, alpha=alpha, )

    if sign in ["neg", "abs"] and (stat_data < -thresh).any():
        stat_data[stat_data>0] = 0
        b.add_data(-stat_data, thresh, max, thresh,
                   colormap="Blues_r", colorbar=colorbar)
    
    if output:
        return stat_data
    
def plot_contrasts(subj, hemi, exps, contrasts, colors, 
                   z_thresh, save_name, alpha=1, save_views=['lateral', 'medial'], save_file=True,
                   snap_views = ['lat', 'med'],
                   plot_conjunction=False, conjunct_color='Purples_r', colorbar=False, sign='pos', corrected=True, 
                   contour=False, n_contours=7):
    
    if contour:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'], config_opts={"cortex": "low_contrast"})
    else:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'])

    
    # determine exp mapping
    if len(exps) == 1:
        exps = exps * len(contrasts)
    
    sig_thresh = -np.log10(stats.norm.sf(z_thresh)); 
    sig_thresh = np.round(sig_thresh) * 10; 
    sign = sign
    
    if corrected:
        sig_name = "cache.th%d.%s.sig.masked.mgh" % (sig_thresh, sign)
    else:
        sig_name = "sig.mgh"

    # Plot the mask
    temp_base = '/Volumes/group/awagner/sgagnon/RM/analysis/%s/group/fsaverage/%s' % (exps[0], contrasts[0])
    mask_temp = op.join(temp_base, "{hemi}/mask.mgh")
    mask_file = mask_temp.format(contrast=contrasts[0],
                                 hemi=hemi, subj=subj)
    add_mask_overlay(b, mask_file)

    z_max = {}; stat_file = {}; sig_data={}
    for exp, contrast in zip(exps, contrasts):
        print exp, contrast
        temp_base = '/Volumes/group/awagner/sgagnon/RM/analysis/%s/group/fsaverage/%s' % (exp, contrast)
        stat_temp = op.join(temp_base, "{hemi}/osgm", sig_name)
        z_max[contrast] = calculate_sat_point(stat_temp, contrast, sign, subj=subj)
        print z_max[contrast]

        stat_file[contrast] = stat_temp.format(contrast=contrast,
                                               hemi=hemi, subj=subj)

    z_max_max = z_max[max(z_max, key=z_max.get)]
    for contrast, color in zip(contrasts, colors):
        print contrast, color
        
        if contour:
            b.add_contour_overlay(stat_file[contrast], colormap=color, colorbar=colorbar,
                                  line_width=4, n_contours=7, min=z_thresh, 
                                  remove_existing=False)
        else:
            sig_data[contrast] = add_stat_overlay(b, stat_file[contrast], z_thresh, z_max_max, sign,
                                                  sig_to_z=True, color=color, alpha=alpha, output=True, colorbar=colorbar)

    if plot_conjunction & (sign == 'pos'):
        print 'plotting conjunction'
        conjunct = np.min(np.vstack(sig_data.values()), axis=0)

        b.add_data(conjunct, min=z_thresh, thresh=z_thresh, max=z_max_max, 
                   colormap=conjunct_color, colorbar=False, alpha=alpha)

    if save_file:
        b.save_imageset(save_name, save_views)
    else:
        view_dict = dict(lat=dict(lh=[160, 50],
                                  rh=[20, 50]),
                         fro=dict(lh=[135, 80],
                                  rh=[45, 80]),
                         par=dict(lh=[230, 55],
                                  rh=[310, 55]),
                         med=dict(lh=[325, 90],
                                  rh=[215, 90]))


        snapshots = dict()
        for view in snap_views:
            a, e = view_dict[view][hemi]
            b.show_view(dict(azimuth=a, elevation=e))
            time.sleep(0.5)
            
            # crop white space
            arr = b.screenshot()
            x,y = np.argwhere((arr != 255).any(axis=-1)).T
            cropped = arr[x.min() - 5:x.max() + 5, y.min() - 5:y.max() + 5, :]
            
            snapshots[view] = cropped
        return b, snapshots

    
def plot_contrasts_sublevel(subj, hemi, exps, contrasts, colors, 
                            z_thresh, save_name,regspace='epi', alpha=1, 
                            save_views=['lateral', 'medial'], save_file=True,
                            plot_conjunction=False, conjunct_color='Purples_r', colorbar=False, sign='pos', contour=False, n_contours=7):
    
    if contour:
        b = Brain(subj, hemi, 'inflated', background="white", views=['parietal'], config_opts={"cortex": "low_contrast"})
        # brighten image
        par_light = b.brains[0]._f.scene.light_manager.lights[-1]
        par_light.intensity = .50
        par_light.elevation = 60
        par_light.azimuth = 30 if hemi == "rh" else -30
        par_light.activate = True
    else:
        b = Brain(subj, hemi, 'inflated', background="white", views=['parietal'])
        
    # determine exp mapping
    if len(exps) == 1:
        exps = exps * len(contrasts)
    
    # Plot the mask
    temp_base = '/Volumes/group/awagner/sgagnon/RM/analysis/%s/%s/ffx/%s/smoothed/%s' % (exps[0], subj, regspace, contrasts[0])
    mask_temp = op.join(temp_base, "{hemi}.mask.mgz")
    mask_file = mask_temp.format(contrast=contrasts[0],hemi=hemi)
    add_mask_overlay(b, mask_file)

    z_max = {}; stat_file = {}; sig_data={}
    for exp, contrast in zip(exps, contrasts):
        print exp, contrast
        temp_base = '/Volumes/group/awagner/sgagnon/rm/analysis/{exp}/{subj}/ffx/{regspace}/smoothed/{contrast}'
        stat_temp = op.join(temp_base, "{hemi}.zstat1.mgz").format(hemi=hemi, exp=exp, subj=subj, 
                                                                   regspace=regspace, contrast=contrast)
        z_max[contrast] = calculate_sat_point(stat_temp, contrast, sign, subj=subj, sig_to_z=False)

        stat_file[contrast] = stat_temp.format(contrast=contrast,
                                               hemi=hemi, subj=subj)

    z_max_max = z_max[max(z_max, key=z_max.get)]
    for contrast, color in zip(contrasts, colors):
        print contrast
        
        if contour:
            b.add_contour_overlay(stat_file[contrast], colormap=color, colorbar=colorbar,
                                  line_width=4, n_contours=7, min=z_thresh,
                                  remove_existing=False)
        else:    
            sig_data[contrast] = add_stat_overlay(b, stat_file[contrast], z_thresh, z_max_max, sign,
                                                  sig_to_z=False, color=color, alpha=alpha, output=True, colorbar=colorbar)

            if plot_conjunction & (sign == 'pos'):
                print 'plotting conjunction'
                conjunct = np.min(np.vstack(sig_data.values()), axis=0)
                print 'Overlapping voxels:' + str(sum(conjunct > z_thresh))

                if sum(conjunct > z_thresh) > 0:
                    b.add_data(conjunct, min=z_thresh, thresh=z_thresh, max=z_max_max, 
                               colormap=conjunct_color, colorbar=False, alpha=alpha)

    if save_file:
        b.save_imageset(save_name, save_views)
    else:
        return b
    