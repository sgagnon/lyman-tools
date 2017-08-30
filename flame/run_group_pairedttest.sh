
basedir=/share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/safe-threat_sourcehit-cr

####################################################################################################
# 1) re-run run_group.py -s subjects_shock_byshockCondSH.txt for group_stress, using update 
# in mixedfx workflow that also outputs merged dofs
# Note that this is changed on steph-backup desktop under group_analysis branch
####################################################################################################

####################################################################################################
# 2) merge FFX copes, dofs, varcopes across subjects into 4d image

# each dim4 21
fslinfo /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/safe_sourcehit-cr/cope_merged.nii.gz
fslinfo /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/threat_sourcehit-cr/cope_merged.nii.gz

fslmerge -t $basedir/cope_merged.nii.gz /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/safe_sourcehit-cr/cope_merged.nii.gz \
                                        /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/threat_sourcehit-cr/cope_merged.nii.gz
fslmerge -t $basedir/dof_merged.nii.gz /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/safe_sourcehit-cr/dof_merged.nii.gz \
                                       /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/threat_sourcehit-cr/dof_merged.nii.gz
fslmerge -t $basedir/varcope_merged.nii.gz /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/safe_sourcehit-cr/varcope_merged.nii.gz \
                                           /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/threat_sourcehit-cr/varcope_merged.nii.gz
fslinfo $basedir/cope_merged.nii.gz # dim4 42

# move group mask over
cp /share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/mni/safe_sourcehit-cr/group_mask.nii.gz \
/share/awagner/sgagnon/AP/analysis/ap_memory_raw-byshockCond/group_stress/safe-threat_sourcehit-cr/.
####################################################################################################


####################################################################################################
# Make a paired t-test design for n subjs
# Glm_gui -> paired t-test default ()
####################################################################################################


####################################################################################################
# Fit the mixed effects model
model=safe-threat_sourcehit-cr
flameo --copefile=$basedir/cope_merged.nii.gz \
--covsplitfile=$basedir/$model.grp \
--designfile=$basedir/$model.mat \
--dofvarcopefile=$basedir/dof_merged.nii.gz \
--ld=stats \
--maskfile=$basedir/group_mask.nii.gz \
--runmode=flame1 \
--tcontrastsfile=$basedir/$model.con \
--varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################


####################################################################################################
# Estimate the smoothness of the data + Correct for multiple comparisons
# Code in Finish paired t-test analysis.ipynb
####################################################################################################


