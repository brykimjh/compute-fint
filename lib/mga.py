import os,sys
import random
import numpy as np
import statistics
np.set_printoptions(suppress=True)

#################
## SUBROUTINES ##
#################

def fitness(nfiles,cv_avg,cv_fint0,cv_fint1,s):
    ilist0 = [0 for i in range(nfiles)]
    ilist1 = [0 for i in range(nfiles)]
    ilist2 = [0.0 for i in range(nfiles)]
    for i in range(len(s)):
        ii = s[i]
        ifint0 = cv_fint0[i][ii]
        ifint1 = cv_fint1[i][ii]

        ilist0[i] = ifint0
        ilist1[i] = ifint1
        ilist2[i] = ifint1-ifint0

#    imean0 = statistics.mean(ilist0)
#    ivar0 = statistics.stdev(ilist0)
#    ivar1 = statistics.stdev(ilist1)
#    ivar2 = statistics.stdev(ilist2)
#    diff = abs(cv_avg-imean0)

#    ierr = diff+ivar0+ivar1+ivar2
#    return ierr
    return 1
    
################
## MAIN BLOCK ##
################

def main(avg,fint0,fint1):
    nfiles = len(fint0)
    nint = len(fint0[0])

# 1. generate solutions
    nsolutions = 1000
    solutions = []
    for s in range(nsolutions):
        samp = [random.randint(0,nint-1) for _ in range(nfiles)]
        solutions.append(samp)
#    print('\nsol',solutions[:5])

    niter = 10000
#    niter = 1
    ibest = 0
    iconverge = 0
    for i in range(niter):
        iconverge += 1
        rankedsolutions = []
        for s in solutions:
            rankedsolutions.append((fitness(nfiles,avg,fint0,fint1,s),s))
        rankedsolutions.sort()
        
#        print(f'=== Gen {i} best solutions ===')
        itop = rankedsolutions[0]
        itest = 1/itop[0]
        if (itest > ibest):
            iconverge = 0
            ibest = itest
#            print(i,avg,itop)
            idx = itop[1]
            imean = 0.0
            ilist0 = [0 for n in range(nfiles)]
            ilist1 = [0 for n in range(nfiles)]
            ilist2 = [0.0 for n in range(nfiles)]
            for n in range(nfiles):
                ind = idx[n]
                ifint0 = fint0[n][ind]
                ifint1 = fint1[n][ind]

                ilist0[n] = ifint0
                ilist1[n] = ifint1
                ilist2[n] = ifint1-ifint0
            
#            imean0 = statistics.mean(ilist0)
#            imean1 = statistics.mean(ilist1)
#            ivar0 = statistics.stdev(ilist0)
#            ivar1 = statistics.stdev(ilist1)
#            ivar2 = statistics.stdev(ilist2)
#            iout = [avg,imean0,imean1,ivar0,ivar1,ivar2]
        '''
            print('\n',iout)
            print('\t',ilist0)
            print('\t',ilist1)

        if (i % 10 == 0):
            print(i,iconverge)
        '''

        if (iconverge == 100):
            return ilist0,ilist1
            break
#        if (iconverge == 100):
#            return ilist0,ilist1,iout
#            break
    
        bestsolutions = rankedsolutions[:100]
#        bestsolutions = rankedsolutions[:10]
#        print('\nbest',bestsolutions)

        elements = []
        for s in bestsolutions:
#            print(s[1][0])
            for j in range(nfiles):
                elements.append(s[1][j])
#        print('\nele',elements)    

        newGen = []
        for _ in range(nsolutions):
            tmp = []
            for j in range(nfiles):
                jj = random.choice(elements) + random.randint(-2,2)
                if (jj < 0):
                    jj = 0
                if (jj > nint-1):
                    jj = nint-1
                tmp.append(jj)
#            print('\ntmp',tmp)
            newGen.append(tmp)
#        print('\nnGen',newGen)
        solutions = newGen

##############
## END MAIN ##
##############

