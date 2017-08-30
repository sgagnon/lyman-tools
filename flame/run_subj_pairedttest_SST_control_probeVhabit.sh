# group=control
# group=stress
# 
# model_name=probe-habit
# indir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group_$group/mni
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group_$group/$model_name
# mkdir $basedir
# contrast_name1=ASSIGNED_shortcut_good_1
# contrast_name2=ASSIGNED_habit_good_1

# model_name=nav_probe-habit
# indir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group_$group/mni
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group_$group/$model_name
# mkdir $basedir
# contrast_name1=NAVIGATE_shortcut_good_1
# contrast_name2=NAVIGATE_habit_good_1

####################################################################################################
# 2) merge FFX copes, dofs, varcopes across subjects into 4d image
# 
# # each dim4 18
# cd $basedir
# fslinfo $indir/$contrast_name1/cope_merged.nii.gz
# fslinfo $indir/$contrast_name2/cope_merged.nii.gz
# 
# fslmerge -t $basedir/cope_merged.nii.gz $indir/$contrast_name1/cope_merged.nii.gz \
#                                         $indir/$contrast_name2/cope_merged.nii.gz
# fslmerge -t $basedir/dof_merged.nii.gz $indir/$contrast_name1/dof_merged.nii.gz \
#                                        $indir/$contrast_name2/dof_merged.nii.gz
# fslmerge -t $basedir/varcope_merged.nii.gz $indir/$contrast_name1/varcope_merged.nii.gz \
#                                            $indir/$contrast_name2/varcope_merged.nii.gz
# fslinfo $basedir/cope_merged.nii.gz # dim4 36
# 
# # move group mask over
# cp $indir/$contrast_name1/group_mask.nii.gz $basedir/.
# cp $indir/$contrast_name1/rh.group_mask.mgz $basedir/.
# cp $indir/$contrast_name1/lh.group_mask.mgz $basedir/.
####################################################################################################


####################################################################################################
# Make a paired t-test design for n subjs
# Glm_gui -> paired t-test default ()
# 1) Higher level design: inputs = 18* 2 = 36, press enter
# 2) Wizard
# 3) 2 groups, paired
# 4) Save to directory (basedir), takes a while
####################################################################################################


####################################################################################################
# Fit the mixed effects model
# cd $basedir
# model=$model_name
# flameo --copefile=$basedir/cope_merged.nii.gz \
# --covsplitfile=$basedir/$model.grp \
# --designfile=$basedir/$model.mat \
# --dofvarcopefile=$basedir/dof_merged.nii.gz \
# --ld=stats \
# --maskfile=$basedir/group_mask.nii.gz \
# --runmode=flame1 \
# --tcontrastsfile=$basedir/$model.con \
# --varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################

####################################################################################################
# Estimate the smoothness of the data + Correct for multiple comparisons
# Code in Finish paired t-test analysis.ipynb
####################################################################################################

####################################################################################################
# Planning (PROBE > HABIT) x (CONTROL > STRESS)
####################################################################################################

# model_name=probe-habit_control-stress
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/$model_name
# mkdir $basedir
# cd $basedir
# contrast_name1=ASSIGNED_shortcut_good_1
# contrast_name2=ASSIGNED_habit_good_1
# 
# indirbase=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt

####################################################################################################
# 2) merge FFX copes, dofs, varcopes across subjects into 4d image

# cd $basedir
# fslinfo $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# 
# fslmerge -t $basedir/cope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# fslmerge -t $basedir/dof_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/dof_merged.nii.gz
# fslmerge -t $basedir/varcope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/varcope_merged.nii.gz
# 
# fslinfo $basedir/cope_merged.nii.gz # 76 long
# fslinfo $basedir/dof_merged.nii.gz # 76 long
# fslinfo $basedir/varcope_merged.nii.gz # 76 long
# 
# # move group mask over
# cp $indirbase/group_control-stress/mni/$contrast_name1/group_mask.nii.gz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/rh.group_mask.mgz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/lh.group_mask.mgz $basedir/.
####################################################################################################


####################################################################################################
# Make a ANOVA: 2-groups, 2-levels per subject (2-way Mixed Effect ANOVA)
# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM#ANOVA:_2-groups.2C_2-levels_per_subject_.282-way_Mixed_Effect_ANOVA.29
# Glm_gui -> paired t-test default ()
# 1) Higher level design: inputs = (18 + 20) * 2 = 76, press enter
# 2) Wizard
# 3) 2 groups, paired
# 4) main EVs: 38 + 2 (run effect, group x run type interaction)
# 5) move s1 to last ev, and replace s1 with g x run interaction
# 6) change group to 1s and 2s for cov split file
# 6) set up contrasts: run effect, interaction
# 7) View design mat
# 8) Save to directory (basedir), takes a while
# 9) confirm that saved stuff looks ok
####################################################################################################

####################################################################################################
# Fit the mixed effects model
# cd $basedir
# model=$model_name
# flameo --copefile=$basedir/cope_merged.nii.gz \
# --covsplitfile=$basedir/$model.grp \
# --designfile=$basedir/$model.mat \
# --dofvarcopefile=$basedir/dof_merged.nii.gz \
# --ld=stats \
# --maskfile=$basedir/group_mask.nii.gz \
# --runmode=flame1 \
# --tcontrastsfile=$basedir/$model.con \
# --varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################

####################################################################################################
# Estimate the smoothness of the data + Correct for multiple comparisons
# Code in Finish paired t-test analysis.ipynb
####################################################################################################

# jeanette version
####################################################################################################
# model_name=probe-habit_control-stress_mumford
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/$model_name
# cd $basedir

# echo $contrast_name1
# echo $contrast_name2

# fslinfo $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# 
# fslmerge -t $basedir/cope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# fslmerge -t $basedir/dof_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/dof_merged.nii.gz
# fslmerge -t $basedir/varcope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/varcope_merged.nii.gz
# 
# fslinfo $basedir/cope_merged.nii.gz # 76 long
# fslinfo $basedir/dof_merged.nii.gz # 76 long
# fslinfo $basedir/varcope_merged.nii.gz # 76 long
# 
# # move group mask over
# cp $indirbase/group_control-stress/mni/$contrast_name1/group_mask.nii.gz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/rh.group_mask.mgz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/lh.group_mask.mgz $basedir/.


# Fit the mixed effects model
# model_name=probe-habit_control-stress_mumford
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/$model_name
# cd $basedir
# model=$model_name
# flameo --copefile=$basedir/cope_merged.nii.gz \
# --covsplitfile=$basedir/$model.grp \
# --designfile=$basedir/$model.mat \
# --dofvarcopefile=$basedir/dof_merged.nii.gz \
# --ld=stats \
# --maskfile=$basedir/group_mask.nii.gz \
# --runmode=flame1 \
# --tcontrastsfile=$basedir/$model.con \
# --varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################


####################################################################################################
##### Jeanette version on restricted mask (combined ROIs)
# Fit the mixed effects model
# model_name=probe-habit_control-stress_mumford
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/probe-habit_control-stress_mumford_restrictedROI
# cd $basedir
# 
# # make group restricted mask
# maskfile=/share/awagner/sgagnon/SST/data/fsaverage/masks/SST_ROIs_combined.nii.gz
# fslmaths $maskfile -mul $basedir/group_mask.nii.gz -bin $basedir/group_mask_restrictedROI.nii.gz
# 
# model=$model_name
# flameo --copefile=$basedir/cope_merged.nii.gz \
# --covsplitfile=$basedir/$model.grp \
# --designfile=$basedir/$model.mat \
# --dofvarcopefile=$basedir/dof_merged.nii.gz \
# --ld=stats \
# --maskfile=$basedir/group_mask_restrictedROI.nii.gz \
# --runmode=flame1 \
# --tcontrastsfile=$basedir/$model.con \
# --varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################

# FDR corrected
fslmaths $basedir/stats/zstat3.nii.gz -ztop $basedir/stats/pstat3.nii.gz
fdr -i $basedir/stats/pstat3.nii.gz -m $basedir/stats/mask.nii.gz -q 0.05 -v


# compare probe 1 vs 2 by group
####################################################################################################
# model_name=probe1v2_control-stress_mumford
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/$model_name
# cd $basedir
# 
# contrast_name1=ASSIGNED_shortcut_good_1
# contrast_name2=ASSIGNED_shortcut_good_2
# echo $contrast_name1
# echo $contrast_name2

# fslinfo $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz
# fslinfo $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz
# fslinfo $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# 
# fslmerge -t $basedir/cope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/cope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/cope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/cope_merged.nii.gz
# fslmerge -t $basedir/dof_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/dof_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/dof_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/dof_merged.nii.gz
# fslmerge -t $basedir/varcope_merged.nii.gz $indirbase/group_control/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name1/varcope_merged.nii.gz \
# $indirbase/group_control/mni/$contrast_name2/varcope_merged.nii.gz \
# $indirbase/group_stress/mni/$contrast_name2/varcope_merged.nii.gz
# 
# fslinfo $basedir/cope_merged.nii.gz # 76 long
# fslinfo $basedir/dof_merged.nii.gz # 76 long
# fslinfo $basedir/varcope_merged.nii.gz # 76 long

# move group mask over
# cp $indirbase/group_control-stress/mni/$contrast_name1/group_mask.nii.gz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/rh.group_mask.mgz $basedir/.
# cp $indirbase/group_control-stress/mni/$contrast_name1/lh.group_mask.mgz $basedir/.


# Fit the mixed effects model
# model_name=probe1v2_control-stress_mumford
# basedir=/share/awagner/sgagnon/SST/analysis/mvpa-goodruntype_filt/group/$model_name
# cd $basedir
# model=$model_name
# flameo --copefile=$basedir/cope_merged.nii.gz \
# --covsplitfile=$basedir/$model.grp \
# --designfile=$basedir/$model.mat \
# --dofvarcopefile=$basedir/dof_merged.nii.gz \
# --ld=stats \
# --maskfile=$basedir/group_mask.nii.gz \
# --runmode=flame1 \
# --tcontrastsfile=$basedir/$model.con \
# --varcopefile=$basedir/varcope_merged.nii.gz
####################################################################################################
