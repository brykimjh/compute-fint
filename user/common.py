# md executeable
COMM_charmm = '/N/project/RP-FM-CV/br3/1_research-data/git/brykim/strg_fm/project/charmm/c42a2_gnu_xxlarge_ifort_sqm_resd_ml_mollif/exec/gnu/charmm'

# qm description
COMM_cv_ind = [1,2,3] # cv indices
COMM_near_cv_ind = [4,5,6,7,8,9] # indices of nearest neighbor atoms to cvs (covalent + 1)
COMM_add_qm_ind = [] # additional qm atoms

COMM_fm_ind = []
COMM_fm_ind.extend(COMM_cv_ind)
COMM_fm_ind.extend(COMM_near_cv_ind)

COMM_qm_ind = []
COMM_qm_ind.extend(COMM_fm_ind)
COMM_qm_ind.extend(COMM_add_qm_ind)

# g16_list
COMM_meth0 = 'AM1' # lo level method
COMM_meth1 = 'B3LYP/6-31+G(d,p)' # hi level method
COMM_nproc = 2 # number of processors for g16
COMM_nchar = 0 # charge of qm system
COMM_nmult = 1 # multiplicity
COMM_cutof = 12.0 # cutoff distance for mm atoms in g16

# g16_slurm
COMM_walltime = '0-3:59:00' # walltime of g16 calculation
COMM_ncores = 24 # max number of cores per node
COMM_scratch = '/N/scratch/brykim' # directory where the gaussian calculation is executed

# cv info (using atoms in COMM_cv_ind)
COMM_natoms = 3 # number of unique cv atoms
COMM_cvs = [
    [[1,2], 1000.0], # cv_ind, force constant
    [[1,3], 1000.0],
    [[2,3], 1000.0]
]

# string parameters
COMM_nskip = 10 # number of MD steps to exclude for force averaging
COMM_delt = 0.001 # time step
COMM_smth = 0.01 # smooth scale

# ann hidden layer
COMM_nnode = 9
COMM_hidden = [
        [f'{COMM_nnode}', 'tanh'],
        ['linear']
    ]
# ann activation function for each layer
# 1 = tanh, 0 = linear
COMM_act = '1 0'
COMM_x1tare = 5 # mollifer %
