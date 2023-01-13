import os,sys
sys.path.append('user')
import numpy as np

from common import *
#from lib import xtrap
from lib import data
#from lib import ml_train
#from lib import predict
#from lib import ml_string

#################
## SUBROUTINES ##
#################

def read_file(file_name):
    try:
        file = open(file_name, 'r')
    except IOError:
        print('Error: file (%s) not found!\n' % (file_name))
        sys.exit()
    lines = file.readlines()
    file.close()
    array = []
    for line in lines:
        array.append(line.split())
    return array

################
## MAIN BLOCK ##
################

def main(argv1,argv2):
    dn0 = f'0_slurm/{argv1}strg'
    dn1 = f'1_string_mfep/{argv1}strg'
    dn2 = f'2_method/{argv1}strg'

# INITIAL FORCE CALC ONLY STEP
# --------------------------------------------------------------------
#    if (argv1 == 0):
#        fn = f'{dn1}/next.txt'
#        f = read_file(fn)
#        last = int(f[0][0])
#        nmd = int(f[0][1])

# 1. make the g16 input files        
#        xtrap.main(dn1,dn2,last,nmd)

# 2. execute the qm calculation        
#        xtrap.slurm(dn0,dn1,dn2,last,nmd)

# 3. update step
#        os.system(f'echo "{argv1+1} {argv2}" > step.txt')

# ML/STRING AND FORCE CALC STEP
# --------------------------------------------------------------------
    if (1 <= argv1 < argv2):
        ii = argv1-1 # previous step
        pd1 = f'1_string_mfep/{ii}strg' # prev dir1
        pd2 = f'2_method/{ii}strg' # prev dir2

        fn = f'{pd1}/next.txt'
        f = read_file(fn)
        last = int(f[0][0])
        nmd = int(f[0][1])

# 1. get training data
        data.main(pd1,pd2,last,nmd)

# 2. train cvs to the corrections with ml
#        ml_train.main(pd2)

# 3. predictions with the ml model
#        predict.main(pd2)

# 4. run string-md with the ml correction
#        ml_string.main(dn1,last,nmd,pd1,pd2)

# FINAL ML/STRING ONLY STEP
# --------------------------------------------------------------------
    if (argv1 == argv2):
        pass

f = read_file('step.txt')
argv1 = int(f[0][0]) # current step
argv2 = int(f[0][1]) # max step

# run main
main(argv1,argv2)

##############
## END MAIN ##
##############

