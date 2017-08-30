#! /usr/bin/env python
"""
This script includes a bunch of functions for surface plotting with pysurfer
Make sure you run %gui qt and %matplotlib inline in a nb before running
"""

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from surfer import Brain, io, project_volume_data
import surfer as sf
from surfer import utils, io

import numpy as np
import glob
import os.path as op
import os
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
        hemi_file = template.format(hemi=hemi)
        
        if sig_to_z:
			hemi_data = nib.load(hemi_file).get_data()

			stat_sign = np.sign(hemi_data)
			p_data = 10 ** -np.abs(hemi_data)
			z_data = stats.norm.ppf(p_data)
			z_data[np.sign(z_data) != stat_sign] *= -1
			hemi_data = z_data

        else:
			reg_file = op.join(os.environ['FREESURFER_HOME'], 'average/mni152.register.dat')
			hemi_data = project_volume_data(hemi_file, hemi, reg_file=reg_file, smooth_fwhm=5)

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


def add_stat_overlay(b, stat_file, thresh, max, sign, hemi='lh', sig_to_z=False, color="Reds_r", alpha=1, output=False, colorbar=False):
    """Plot a surface-encoded statistical overlay."""
    
    # Possibly convert -log10(p) images to z stats
    if sig_to_z:
		stat_data = nib.load(stat_file).get_data()
		
		stat_sign = np.sign(stat_data)
		p_data = 10 ** -np.abs(stat_data)
		z_data = stats.norm.ppf(p_data)
		z_data[np.sign(z_data) != stat_sign] *= -1
		stat_data = z_data
    else:
		reg_file = op.join(os.environ['FREESURFER_HOME'], 'average/mni152.register.dat')
		stat_data = project_volume_data(stat_file, hemi, reg_file=reg_file, smooth_fwhm=5, verbose=False)

    # Plot the statistical data
    stat_data = stat_data.squeeze()
    if sign in ["pos", "abs"] and (stat_data > thresh).any():
        stat_data[stat_data<0] = 0
        b.add_data(stat_data, thresh, max, thresh,
                   colormap=color, colorbar=colorbar, alpha=alpha)

    if sign in ["neg", "abs"] and (stat_data < -thresh).any():
        stat_data[stat_data>0] = 0
        b.add_data(-stat_data, thresh, max, thresh,
                   colormap="Blues_r", colorbar=colorbar)
    
    if output:
        return stat_data
    
def plot_contrasts(subj, hemi, exps, contrasts, colors, 
                   z_thresh, save_name, alpha=1, save_views=['lateral', 'medial'], save_file=True,
                   base_exp='/Volumes/group/awagner/sgagnon/RM',
                   group='group', regspace='fsaverage',
                   skip_reg_dir=False,
                   contrast_num=1,
                   snap_views = ['lat', 'med'], sig_to_z=True,
                   plot_conjunction=False, conjunct_color='Purples_r', colorbar=False, sign='pos', corrected=True, 
                   contour=False, n_contours=7, view_dict=None):
                   
    if view_dict is None:
        view_dict = dict(lat_rot=dict(lh=[160, 50],
                                      rh=[20, 50]),
                         lat=dict(lh=[180, 90],
                                  rh=[180, -90]),
                         fro=dict(lh=[135, 80],
                                  rh=[45, 80]),
                         par=dict(lh=[230, 55],
                                  rh=[310, 55]),
                         med_rot=dict(lh=[325, 90],
                                      rh=[215, 90]),
                         med=dict(lh=[0,90],
                                  rh=[0,-90]))
    
    if contour:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'], config_opts={"cortex": "low_contrast"})
    else:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'])

    
    # determine exp mapping
    if len(exps) == 1:
        exps = exps * len(contrasts)
    
    # figure out threshold and stats files depending on surf analysis or MNI
    if sig_to_z:
		sig_thresh = -np.log10(stats.norm.sf(z_thresh)); 
		sig_thresh = np.round(sig_thresh) * 10; 
		
		if corrected:
			sig_name = "cache.th%d.%s.sig.masked.mgh" % (sig_thresh, sign)
		else:
			sig_name = "sig.mgh"
    else:
		sig_thresh = z_thresh
		
		if corrected:
			sig_name = "zstat"+str(contrast_num)+"_threshold.nii.gz"
		else:
			sig_name = "zstat"+str(contrast_num)+".nii.gz"


    # Plot the mask
    if skip_reg_dir:
        temp_base = op.join(base_exp, 'analysis/%s', group, '%s') % (exps[0], contrasts[0])
    else:
        temp_base = op.join(base_exp, 'analysis/%s', group, regspace,'%s') % (exps[0], contrasts[0])
    
    if regspace == 'fsaverage':
		mask_temp = op.join(temp_base, "{hemi}/mask.mgh")
		mask_file = mask_temp.format(contrast=contrasts[0],
									 hemi=hemi, subj=subj)
    else:
		mask_temp = op.join(temp_base, "{hemi}.group_mask.mgz")
		mask_file = mask_temp.format(hemi=hemi)
	
    add_mask_overlay(b, mask_file)

    z_max = {}; stat_file = {}; sig_data={}
    for exp, contrast in zip(exps, contrasts):
        print exp, contrast
        
        if skip_reg_dir:
            if corrected:
                temp_base = op.join(base_exp, 'analysis/%s', group, '%s') % (exp, contrast)
            else:
                temp_base = op.join(base_exp, 'analysis/%s', group, '%s', 'stats') % (exp, contrast)
        else:
            temp_base = op.join(base_exp, 'analysis/%s', group, regspace,'%s') % (exp, contrast)
        
    	if regspace == 'fsaverage':
			stat_temp = op.join(temp_base, "{hemi}/osgm", sig_name)
        else:
			stat_temp = op.join(temp_base, sig_name)
			
        z_max[contrast] = calculate_sat_point(stat_temp, contrast, sign, subj=subj, sig_to_z=sig_to_z)
        print z_max[contrast]

        stat_file[contrast] = stat_temp.format(hemi=hemi)

    z_max_max = z_max[max(z_max, key=z_max.get)]
    for contrast, color in zip(contrasts, colors):
        print contrast, color
        
        if contour:
            b.add_contour_overlay(stat_file[contrast], colormap=color, colorbar=colorbar,
                                  line_width=4, n_contours=7, min=z_thresh, 
                                  remove_existing=False)
        else:
            sig_data[contrast] = add_stat_overlay(b, stat_file[contrast], z_thresh, z_max_max, sign,
                                                  hemi=hemi, sig_to_z=sig_to_z, color=color, alpha=alpha, 
                                                  output=True, colorbar=colorbar)

    if plot_conjunction & (sign == 'pos'):
        print 'plotting conjunction'
        conjunct = np.min(np.vstack(sig_data.values()), axis=0)

        b.add_data(conjunct, min=z_thresh, thresh=z_thresh, max=z_max_max, 
                   colormap=conjunct_color, colorbar=False, alpha=alpha)

    if save_file:
        b.save_imageset(save_name, save_views)
    else:

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

def plot_groups(subj, hemi, exps, contrast, colors, 
                z_thresh, save_name, alpha=1, save_views=['lateral', 'medial'], save_file=True,
                base_exp='/Volumes/group/awagner/sgagnon/RM',
                groups=['group_control', 'group_stress'], 
                group_mask_dir='group_control-stress', regspace='fsaverage',
                skip_reg_dir=False,
                contrast_num=1,
                snap_views = ['lat', 'med'], sig_to_z=True,
                plot_conjunction=False, conjunct_color='Purples_r', colorbar=False, sign='pos', corrected=True, 
                contour=False, n_contours=7, view_dict=None, add_border=False, border_max=False):

# plot one contrast for 2 diff groups
    
    if view_dict is None:
        view_dict = dict(lat_rot=dict(lh=[160, 50],
                                      rh=[20, 50]),
                         lat=dict(lh=[180, 90],
                                  rh=[180, -90]),
                         fro=dict(lh=[135, 80],
                                  rh=[45, 80]),
                         par=dict(lh=[230, 55],
                                  rh=[310, 55]),
                         med_rot=dict(lh=[325, 90],
                                      rh=[215, 90]),
                         med=dict(lh=[0,90],
                                  rh=[0,-90]))
    
    if contour:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'], config_opts={"cortex": "low_contrast"})
    else:
        b = Brain(subj, hemi, 'semi7', background="white", views=['parietal'])

    
    # determine exp mapping
    if len(exps) == 1:
        exps = exps * len(groups)
    
    # figure out threshold and stats files depending on surf analysis or MNI
    if sig_to_z:
        sig_thresh = -np.log10(stats.norm.sf(z_thresh)); 
        sig_thresh = np.round(sig_thresh) * 10; 
        
        if corrected:
            sig_name = "cache.th%d.%s.sig.masked.mgh" % (sig_thresh, sign)
        else:
            sig_name = "sig.mgh"
    else:
        sig_thresh = z_thresh

        if corrected:
            sig_name = "zstat"+str(contrast_num)+"_threshold.nii.gz"
        else:
            sig_name = "zstat"+str(contrast_num)+".nii.gz"


    # Plot the mask
    if skip_reg_dir:
        temp_base = op.join(base_exp, 'analysis/%s', group_mask_dir, '%s') % (exps[0], contrast)
    else:
        temp_base = op.join(base_exp, 'analysis/%s', group_mask_dir, regspace,'%s') % (exps[0], contrast)
    
    if regspace == 'fsaverage':
        print 'figure out mask'
# 		mask_temp = op.join(temp_base, "{hemi}/mask.mgh")
# 		mask_file = mask_temp.format(contrast=contrast,
# 									 hemi=hemi, subj=subj)
    else:
		mask_temp = op.join(temp_base, "{hemi}.group_mask.mgz")
		mask_file = mask_temp.format(hemi=hemi)
	
    add_mask_overlay(b, mask_file)

    z_max = {}; stat_file = {}; sig_data={}
    for exp, group in zip(exps, groups):
        print exp, contrast, group
        
        if skip_reg_dir:
            if corrected:
                temp_base = op.join(base_exp, 'analysis/%s', group, '%s') % (exp, contrast)
            else:
                temp_base = op.join(base_exp, 'analysis/%s', group, '%s', 'stats') % (exp, contrast)
        else:
            temp_base = op.join(base_exp, 'analysis/%s', group, regspace,'%s') % (exp, contrast)
        
    	if regspace == 'fsaverage':
			stat_temp = op.join(temp_base, "{hemi}/osgm", sig_name)
        else:
			stat_temp = op.join(temp_base, sig_name)
			
        z_max[group] = calculate_sat_point(stat_temp, contrast, sign, subj=subj, sig_to_z=sig_to_z)
        print z_max[group]

        stat_file[group] = stat_temp.format(hemi=hemi)

    z_max_max = z_max[max(z_max, key=z_max.get)]
    print 'Colorbar max: z = ' + str(z_max_max)
    for group, color in zip(groups, colors):
        print contrast, color
        
        if contour:
            b.add_contour_overlay(stat_file[group], colormap=color, colorbar=colorbar,
                                  line_width=4, n_contours=7, min=z_thresh, 
                                  remove_existing=False)
        else:
            print stat_file[group]
            sig_data[group] = add_stat_overlay(b, stat_file[group], z_thresh, z_max_max, sign,
                                                  hemi=hemi, sig_to_z=sig_to_z, color=color, alpha=alpha, 
                                                  output=True, colorbar=colorbar)

    if plot_conjunction & (sign == 'pos'):
        print 'plotting conjunction'
        conjunct = np.min(np.vstack(sig_data.values()), axis=0)

        b.add_data(conjunct, min=z_thresh, thresh=z_thresh, max=z_max_max, 
                   colormap=conjunct_color, colorbar=False, alpha=alpha)

    # add some overlay (filename should be add_border)
    if add_border:
        print 'adding border'

        if border_max:
            border_max = calculate_sat_point(add_border, '', sign, subj=subj, sig_to_z=sig_to_z)    
        else:
            border_max=z_max_max
        
        print border_max
        print z_thresh

        reg_file = op.join(os.environ['FREESURFER_HOME'], 'average/mni152.register.dat')
        border_data = project_volume_data(add_border, hemi, reg_file=reg_file, smooth_fwhm=5, verbose=False)
        
        if len(border_data[border_data > 2.3]) > 0:
            
            b.add_contour_overlay(border_data, hemi=hemi, n_contours=4, line_width=4,
                                  min=z_thresh, max=border_max, colormap='binary_r', colorbar=False)

    if save_file:
        b.save_imageset(save_name, save_views)
    else:

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
                            save_views=['lateral', 'medial'], save_file=True, sig_to_z=False,
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
        z_max[contrast] = calculate_sat_point(stat_temp, contrast, sign, subj=subj, sig_to_z=sig_to_z)

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
                                                  sig_to_z=sig_to_z, color=color, alpha=alpha, output=True, colorbar=colorbar)

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
    