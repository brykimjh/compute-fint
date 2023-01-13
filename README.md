# compute_fint
Python code which computes internal force corrections between semiempirical and ab initio methods.

$ python3 run.py to run script.

| File/Directory| Description   |
| ------------- | ------------- |
| 1_string_mfep/| contains coordinate files from QM/MM molecular dynamics   |
| 2_method/     | contains semiempirical and ab initio forces along with internal force corrections (target.csv)  |
| lib/          | scripts for computing internal forces (geom_opt.py contains B_Matrix transformatiion)  |
| scratch/      | scratch directory  |
| user/         | contains custom input files for differenct chemical systems, and job submission scripts for different clusters  |
| job.slurm     | job submission script (python3 run.py to run code)  |
| run.py        | driver file for RP_FM_CV (enhanced sampling and machine learning components are omitted)  |
| step.txt      | indicates current iteration of string QM/MM_MD  |
