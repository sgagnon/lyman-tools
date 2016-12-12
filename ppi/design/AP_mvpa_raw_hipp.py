import os.path as op

home_dir = '/share/awagner/sgagnon/AP'
exps = ['mvpa_raw']
design_file = 'AP_sourcehits-vs-CR.csv'
contrast= 'sourcehit-vs-CR' # contrast to use in PPI analysis
data_dir = op.join(home_dir, 'data') 
analy_dir = op.join(home_dir, 'analysis', '{exp}')
num_runs = 6
masks = ['lh-hippocampus-tail', 'rh-hippocampus-tail']
subj_file = op.join(home_dir, 'scripts', 'subjects.txt')
n_cores = 16
