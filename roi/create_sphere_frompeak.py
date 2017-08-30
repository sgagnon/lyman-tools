#! /usr/bin/env python
"""
From peak voxel in MNI space, create a sphere for an ROI
"""

import os
import sys
import re
import os.path as op
from textwrap import dedent
import argparse
from subprocess import call

import numpy as np
import pandas as pd
import nibabel as nib

import lyman
from lyman import tools


def main(arglist):
    """Main function for workflow setup and execution."""
    args = parse_args(arglist)

    if isinstance(args.peak_num, list):
        peak_num = args.peak_num
    else:
        peak_num = int(args.peak_num)

    def mni_to_vox(mni_coords, output=True):
        """Given xyz mni coordinates, return ijk voxel coordinates using image affine.
        Example: 
            mni_coords = np.array([[18, -36, 0], [34, 0, -24]]); 
            mni_to_vox(mni_coords)
        """
        try:
            fsldir = os.environ["FSLDIR"]
        except KeyError:
            raise RuntimeError("mni_to_vox requires FSLDIR to be defined.")
        mni_file = op.join(fsldir, "data/standard/avg152T1.nii.gz")
        aff = nib.load(mni_file).affine
        vox_coords = np.zeros_like(mni_coords)
        for i, coord in enumerate(mni_coords):
            coord = coord.astype(float)
            vox_coords[i] = np.dot(np.linalg.inv(aff), np.r_[coord, 1])[:3].astype(int)
            if output:
                print 'MNI ' + str(coord.astype(int)) + ' to: ' + str(vox_coords[i])
        return vox_coords

    # Get some project information
    project = lyman.gather_project_info()

    # If peak number, read in csv and get voxel coordinates in avg152 coordinate space
    if isinstance(peak_num, int):
        in_fname = op.join(project['analysis_dir'], args.experiment, args.group, 
                        'mni', args.contrast, 'zstat1_localmax.csv')
        out_fname = in_fname[:-12] + 'peak' + str(peak_num) + '_'+ args.sphere_rad +'mm_sphere.nii.gz'
        masked_fname = in_fname[:-12] + 'peak' + str(peak_num) + '_'+ args.sphere_rad +'mm_sphere_masked.nii.gz'
        threshold_fname = in_fname[:-12] + 'threshold.nii.gz'
        df = pd.read_csv(in_fname)

        if len(df):
            coords = df[["x", "y", "z"]].values
            vox_coords = mni_to_vox(coords)

    # Select appropriate coordinate, and create the sphere w/FSL
    # Code borrowed from: http://www.jonaskaplan.com/lab/files/make_roi_sphere.py
    if vox_coords.shape[0] > peak_num:
        vox_coord = vox_coords[peak_num]

        command = "fslmaths $FSLDIR/data/standard/avg152T1.nii.gz -mul 0 -add 1 -roi %s 1 %s 1 %s 1 0 1 tmp" % (vox_coord[0],vox_coord[1],vox_coord[2])
        print command
        call(command,shell=True)

        command = "fslmaths tmp -kernel sphere %s -fmean tmp" % args.sphere_rad
        print command
        call(command,shell=True)

        command = "fslmaths tmp -thr .00001 -bin %s" % out_fname
        print command
        call(command,shell=True)

        # Mask the sphere with activation
        command = "fslmaths tmp -mas %s -bin %s" % (threshold_fname, masked_fname)
        print command
        call(command,shell=True)

        command = "rm tmp.nii.gz"
        call(command,shell=True)

        

    else: print('Peak unavailable, try another one')


def parse_args(arglist):
    """Take an arglist and return an argparse Namespace."""
    help = dedent("""
      Usage Details
    -------------
    """)

    parser = tools.parser
    parser.description = help
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument("-experiment", help="experimental paradigm")
    parser.add_argument("-contrast", help="contrast")
    parser.add_argument("-group", default="group", help="group directory")
    parser.add_argument("-peak_num", help="Number of peak in localmax csv")
    parser.add_argument("-sphere_rad", default="5", help="radius of peak sphere")

    return parser.parse_args(arglist)

if __name__ == "__main__":
    main(sys.argv[1:])