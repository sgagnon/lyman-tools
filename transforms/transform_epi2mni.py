import os

import sys
import shutil
import os.path as op
from textwrap import dedent
import argparse

import subprocess as sp

import nipype
from nipype import Node, SelectFiles, DataSink, IdentityInterface
from nipype.interfaces import fsl

import lyman
import lyman.workflows as wf
from lyman import tools
from lyman.tools import add_suffix, submit_cmdline


##########################################
# Set up some info for this experiment
##########################################
exp_name = 'mvpa_raw'
altmodel = None
space = 'mni'
smoothing = 'unsmoothed'
subjects = None
regtype = 'model'
interpolation = "trilinear"

# localizer
# time_list = [4.5]
# searchlight_path = "searchlight/localizer_acc_{subject_id}.nii.gz"

# reinstatement
time_list = [0,2,4,6,8,10,12]
searchlight_dir = 'searchlight_test'
searchlight_path = searchlight_dir + "/sourcehit_time{time}_acc_{{subject_id}}.nii.gz"

##########################################
# Pull info from project
##########################################
project = lyman.gather_project_info()
exp = lyman.gather_experiment_info(exp_name, altmodel)
warp_method = project['normalization']

os.environ["SUBJECTS_DIR"] = project["data_dir"]
subject_list = lyman.determine_subjects(subjects)
subj_source = tools.make_subject_source(subject_list)

exp_base = exp_name
if altmodel is not None:
    exp_name = "-".join([exp_base, altmodel])
    
data_dir = project["data_dir"]
analysis_dir = op.join(project["analysis_dir"], exp_name)
working_dir = op.join(project["working_dir"], exp_name)

##########################################
# Set up paths to files
##########################################
reg_templates = dict(
    masks="{subject_id}/preproc/run_{run}/functional_mask.nii.gz",
    means="{subject_id}/preproc/run_{run}/mean_func.nii.gz",
                     )
reg_templates.update(dict(
        searchlight=op.join(searchlight_path)))
reg_lists = reg_templates.keys()
print reg_lists

aff_ext = "mat" if warp_method == "fsl" else "txt"
reg_templates["warpfield"] = op.join(data_dir, "{subject_id}",
                                     "normalization/warpfield.nii.gz")
reg_templates["affine"] = op.join(data_dir, "{subject_id}",
                                  "normalization/affine." + aff_ext)

# Rigid (6dof) functional-to-anatomical matrices
rigid_stem = op.join(analysis_dir, 
                     "{subject_id}/preproc/run_{run}/func2anat_")
if warp_method == "ants" and space == "mni":
    reg_templates["rigids"] = rigid_stem + "tkreg.dat"
else:
    reg_templates["rigids"] = rigid_stem + "flirt.mat"
    
ref_file = fsl.Info.standard_image("avg152T1_brain.nii.gz")

in_file = op.join(analysis_dir, reg_templates['searchlight'])
out_fname = op.basename(add_suffix(in_file, "warp"))
out_file =  op.join(analysis_dir, searchlight_dir, out_fname)
out_rigid = op.join(analysis_dir, searchlight_dir, op.basename(add_suffix(out_file, "anat")))


##########################################
# Warp images
##########################################

print subject_list

for time in time_list:

    # update accordingly for time
    if len(time_list) > 1:
        in_file_spec = in_file.format(time=str(time))
        out_file_spec = out_file.format(time=str(time))
        out_rigid_spec = out_rigid.format(time=str(time))
    else:
        in_file_spec = in_file
        out_file_spec = out_file
        out_rigid_spec = out_rigid

    for subid in subject_list:
        print subid

        continuous_interp = dict(trilinear="trilin",
                                 spline="cubic")[interpolation]
        interp = "nearest" if "mask" in in_file_spec else continuous_interp
        cmdline_rigid = ["mri_vol2vol",
                         "--mov", in_file_spec.format(subject_id=subid),
                         "--reg", reg_templates['rigids'].format(run=1, subject_id=subid),
                         "--fstarg",
                         "--" + interp,
                         "--o", out_rigid_spec.format(subject_id=subid),
                         "--no-save-reg"]
        cmdline = " ".join(cmdline_rigid)
        print cmdline
        os.system(cmdline)

        continuous_interp = dict(trilinear="trilin",
                                 spline="BSpline")[interpolation]
        interp = "NN" if "mask" in in_file else continuous_interp
        cmdline_warp = ["WarpImageMultiTransform",
                        "3",
                        out_rigid_spec.format(subject_id=subid),
                        out_file_spec.format(subject_id=subid),
                        reg_templates['warpfield'].format(subject_id=subid),
                        reg_templates['affine'].format(subject_id=subid),
                        "-R", ref_file]
        if interp != "trilin":
            cmdline_warp.append("--use-" + interp)
        cmdline = " ".join(cmdline_warp)
        print cmdline
        os.system(cmdline)

##########################################
# Combine across subjects into 4d group image
##########################################

##### LOCALIZER
#################################
# cmdline_merge = ["fslmerge",
#                  "-t",
#                  "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_4D",
#                  "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_*_warp.nii.gz"]
# cmdline = " ".join(cmdline_merge)
# print cmdline
# os.system(cmdline)

# fslmaths /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_4D -sub 0.3333 -Tmean /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_mean
# fslmaths /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_4D -sub 0.3333 /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_4D

##### REINSTATEMENT
#################################

for time in time_list:

    # Create a 4D image from all the subjects data
    cmdline_merge = ["fslmerge",
                     "-t",
                     "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_4D".format(time=str(time)),
                     "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_*_warp.nii.gz".format(time=str(time))]
    cmdline = " ".join(cmdline_merge)
    print cmdline
    os.system(cmdline)

    # subtract chance, and take the mean across subjects
    cmdline_mean = ["fslmaths",
                    "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_4D".format(time=str(time)),
                    "-sub 0.3333 -Tmean",
                    "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_mean".format(time=str(time))]
    cmdline = " ".join(cmdline_mean)
    print cmdline
    os.system(cmdline)

    # subtract chance (0.33) from 4D (for subsequent 1-sample t-testing against chance)
    cmdline_mean = ["fslmaths",
                    "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_4D".format(time=str(time)),
                    "-sub 0.3333",
                    "/share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_4D".format(time=str(time))]
    cmdline = " ".join(cmdline_mean)
    print cmdline
    os.system(cmdline)


##########################################
# Plot group mean accuracy on surface
##########################################

# from surfer import Brain, project_volume_data

##### LOCALIZER
#################################
# brain = Brain("fsaverage", "split", "inflated",  views=['lat', 'med', 'ven'], background="white")
# volume_file = "/Volumes/group/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_t_tstat1.nii.gz"

# reg_file = os.path.join(os.environ['FREESURFER_HOME'], "average/mni152.register.dat")

# for hemi in ['lh', 'rh']:
#     zstat = project_volume_data(volume_file, hemi, reg_file, smooth_fwhm=0.5)
#     brain.add_overlay(zstat, hemi=hemi, min=2.3)
# brain.save_image('/Volumes/group/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_t_tstat1_t23_surf.png')


##### REINSTATEMENT
#################################

# for time in time_list:
#     brain = Brain("fsaverage", "split", "inflated",  views=['lat', 'med', 'ven'], background="white")
#     volume_file = "/Volumes/group/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_mean.nii.gz".format(time=str(time))

#     reg_file = os.path.join(os.environ['FREESURFER_HOME'], "average/mni152.register.dat")

#     for hemi in ['lh', 'rh']:
#         zstat = project_volume_data(volume_file, hemi, reg_file, smooth_fwhm=0.5)
#         brain.add_overlay(zstat, hemi=hemi, min=0.05)
#     brain.save_image('/Volumes/group/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight_test/sourcehit_time{time}_acc_mean_surf.png'.format(time=str(time)))

##########################################
# 1-samp t-test
##########################################

# cmdline = 'randomise_parallel -i /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_4D -o /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_t -1 -T'
# print cmdline
# os.system(cmdline)

#cluster -i /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_t_tfce_corrp_tstat1 -t 0.95 -c /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_t_tfce_tstat1 --scalarname="1-p" > /share/awagner/sgagnon/AP/analysis/mvpa_raw/searchlight/localizer_acc_corrp1.txt

