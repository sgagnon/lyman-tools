import os.path as op

smoothing = 'unsmoothed'
regspace = 'epi'
project = 'SST'
design = 'navcond_byruntype.csv'
func_exp = 'mvpa'
onset_exp = 'mvpa-runtype'
smoothing_fwhm = 0
standardize_feat = False
standardize_roi = False
percentsig_roi = True
percentsig_roi_relbaseline = False
tr = float(2)
tr_shift = [0, 2, 4, 6, 8, 10, 12] # in seconds
tr_integrate = [0, 2, 4, 6, 8, 10] # in seconds
n_runs = 12

basedir = op.join('/share/awagner/sgagnon', project)
analydir = op.join(basedir, 'analysis', func_exp)
expdir = op.join(basedir, 'analysis', onset_exp)
subjfile = op.join(basedir, 'scripts/subjects_14runs.txt')

# Filepath templates
tsfile = op.join(analydir, "{subid}", 'reg', regspace,
                 smoothing, "run_{run_id}", 'timeseries_xfm.nii.gz')
func_maskfile = op.join(analydir, "{subid}", 'reg', regspace,
                        smoothing, "run_{run_id}", 'functional_mask_xfm.nii.gz')
maskfile = op.join(basedir, 'data', "{subid}", 'masks',
                   "{mask_name}.nii.gz")
meanfile = op.join(analydir, "{subid}", 'preproc',
                   "run_{run_id}", 'mean_func.nii.gz')
onsetfile = op.join(basedir, 'data', "{subid}", 'design', design)
