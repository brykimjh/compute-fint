# compute_fint
Python code which computes internal force corrections between semiempirical and ab initio methods.

"Reaction Path-Force Matching in Collective Variables: Determining Ab Initio QM/MM Free Energy Profiles by Fitting Mean Force," Kim, B.; Snyder, R.; Nagaraju, M.; Zhou, Y.; Ojeda-May, P.; Keeton, S.; Hege, M.; Shao, Y.; Pu, J. J. Chem. Theory Comput. 2021, 17, 4961-4980 (doi: 10.1021/acs.jctc.1c00245; PMID:34283604).

![alt text](https://pubs.acs.org/cms/10.1021/acs.jctc.1c00245/asset/images/medium/ct1c00245_0014.gif)

Author: Bryant Kim  
Email: brykimjh@gmail.com  
Website: https://wikipugr.sitehost.iu.edu/wiki/index.php?title=The_Pu_Research_Group  

To run python script:
```
$ python3 run.py
```

| File/Directory| Description   |
| ------------- | ------------- |
| 1_string_mfep/| contains coordinate files from QM/MM molecular dynamics   |
| 2_method/     | contains semiempirical and ab initio forces along with internal force corrections (target.csv)  |
| lib/          | scripts for computing internal forces (geom_opt.py contains B-Matrix transformatiion)  |
| scratch/      | scratch directory  |
| user/         | contains custom input files for different chemical systems, and job submission scripts for different clusters  |
| job.slurm     | job submission script (python3 run.py to run code)  |
| run.py        | driver file for RP-FM-CV (enhanced sampling and machine learning components are omitted)  |
| step.txt      | indicates current iteration of string QM/MM-MD  |
