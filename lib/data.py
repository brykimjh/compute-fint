import os,sys
import time
import numpy as np
import concurrent.futures

from common import *
#from lib import xtrap
from lib import geom_opt
from lib import mga
from scipy.spatial import KDTree

np.set_printoptions(suppress=True)

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

def get_zcv(pd1,last,nmd):
    fn = f'{pd1}/data/z{last}.txt'
    f = read_file(fn)
#    nimg = len(f)
    nimg = 1
    ncv = len(COMM_cvs)

    zdat = np.zeros((nimg,nmd,ncv))
    fdz = np.zeros((nimg,nmd,ncv))
    for i in range(nmd):
        md = last-i
        fn = f'{pd1}/data/z{md}.txt'
        f1 = read_file(fn)
        fn = f'{pd1}/data/fdz_{md}.txt'
        f2 = read_file(fn)

        for j in range(nimg):
            for k in range(ncv):
                zdat[j][i][k] = float(f1[j][k])
                fdz[j][i][k] = float(f2[j][k])
    return zdat,fdz

def print_header(fout):
    fn = fout[0]
    col = fout[1]

    if (col == 1):
        natm = len(COMM_cv_ind)
    if (col == 3):
        natm = len(COMM_fm_ind)
    if (col == 5):
        natm = len(COMM_cv_ind)

    fp = open(fn, "w")
    for i in range(natm):
        ii = i+1
        if (col > 1 and i == 0):
            fp.write(f'ind')
        if (col == 1):
            if (i == 0):
                text = f'cv{ii}'
            else:
                text = f',cv{ii}'
        if (col == 3):
            text = f',f{ii}x,f{ii}y,f{ii}z'
        if (col == 5):
            text = f',nlam{ii},rdist{ii},iflo{ii},ifhi{ii},fcorr{ii}'
        fp.write(f'{text}')
    fp.write(f'\n')
    fp.close()

def get_fm_coord(ind,cfile):
    natm = len(COMM_fm_ind)
    fn = f'scratch/crd_tmp{ind}'

    CC = f'''\
    grep -A {natm+1} 'EXT' {cfile} | tail -n +2 | awk '{{print $5,$6,$7}}' > {fn}
    '''
#    print(CC)
    os.system(CC)

# read the temporary crd file
    f = read_file(fn)
    os.system(f'rm -rf {fn}')

# save the coordinates to an array
    crd = np.zeros((natm,3))
    for i in range(natm):
        ii = int(COMM_fm_ind[i])-1
        for j in range(3):
            crd[i][j] = float(f[ii][j])

    return crd

def get_lines(ind,ffile):
    fn = f'{ffile}'
    lines_tmp = f'scratch/lines_tmp{ind}'
    SYS = f'''
    grep -n 'Forces (' {fn} | cut -f1 -d: > {lines_tmp}
    '''
#    print(SYS)
    os.system(SYS)
    lines = read_file(f'{lines_tmp}')
    iflo = int(lines[0][0])
    ifhi = int(lines[1][0])
    os.system(f'rm -rf {lines_tmp}')
    return iflo,ifhi

def grep_forces(ind,ffile,iflo,ifhi,natm):
    fn = f'{ffile}'
    flo_tmp = f'scratch/flo_tmp{ind}'
    fhi_tmp = f'scratch/fhi_tmp{ind}'
    a0 = iflo+3
    a1 = a0+natm-1
    b0 = ifhi+3
    b1 = b0+natm-1
    SYS = f'''
sed '{a0},{a1}!d' {fn} | awk '{{ print $3,$4,$5 }}' > {flo_tmp}
sed '{b0},{b1}!d' {fn} | awk '{{ print $3,$4,$5 }}' > {fhi_tmp}
    '''
    #print(SYS)
    os.system(SYS)
    flo = read_file(f'{flo_tmp}')
    fhi = read_file(f'{fhi_tmp}')
    os.system(f'rm -rf {flo_tmp} {fhi_tmp}')
    return flo,fhi

def fconvert(flo,fhi):
    natm = len(flo)
# constants used in gaussian16 https://gaussian.com/constants/
    hartree = 627.5095
    bohr = 0.52917721092

    lof = np.zeros((natm,3))
    hif = np.zeros((natm,3))
    for i in range(natm):
        for j in range(3):
            lo = float(flo[i][j])
            hi = float(fhi[i][j])
            lof[i][j] = -lo * hartree/bohr
            hif[i][j] = -hi * hartree/bohr
    return lof,hif

def get_force(ind,ffile):
    natm = len(COMM_fm_ind)
# 1. get the line numbers for flo and fhi in *.log
    iflo,ifhi = get_lines(ind,ffile)

# 2. get flo and fhi
    flo,fhi = grep_forces(ind,ffile,iflo,ifhi,natm)

# 3. convert the forces to kcal/mol/a
    lof,hif = fconvert(flo,fhi)

    return lof,hif

def get_ic(lof):
    ncv = len(COMM_cvs)
    natm = len(lof)
    nic = int(natm*(natm-1)/2)

    ic_list = [[0 for j in range(2)] for i in range(nic)]
    cantor_list = [0 for i in range(nic)]
    ind = 0
    for i in range(natm-1):
        k1 = i+1
        for j in range(i+1,natm):
            k2 = j+1
            ic_list[ind] = [k1,k2]
            cantor_list[ind] = int((k1 + k2)*(k1 + k2 + 1)/2 + k2)
            ind += 1

#    print(ic_list)
#    print(cantor_list)

    ic = []
    for i in range(ncv):
        k1 = COMM_cvs[i][0][0]
        k2 = COMM_cvs[i][0][1]
        cantor = int((k1 + k2)*(k1 + k2 + 1)/2 + k2)
        ind = cantor_list.index(cantor)
        rm1 = cantor_list.pop(ind)
        rm2 = ic_list.pop(ind)
#        print(cantor,ind,rm1,rm2)
        ic.append(rm2)
    ic.extend(ic_list)
#    print(ic_list)
#    print(ic)
    return ic

def print_forces(ind,lof,hif,f4,f5):
    natm = len(lof)
    o1 = open(f4, 'a')
    o2 = open(f5, 'a')
    for i in range(natm):
        if (i == 0):
            o1.write(f'{ind}')
            o2.write(f'{ind}')
        for j in range(3):
            f0 = lof[i][j]
            f1 = hif[i][j]
            o1.write(f',{f0}')
            o2.write(f',{f1}')
    o1.write('\n')
    o2.write('\n')
    o1.close()
    o2.close()

def get_g16_lines(ind,ffile):
    fn = f'{ffile}'
    lines_tmp = f'scratch/crd_tmp{ind}'
    SYS = f'''
    grep -n 'No Z-Matrix' {fn} | cut -f1 -d: > {lines_tmp}
    grep -n 'Recover connectivity' {fn} | cut -f1 -d: >> {lines_tmp}
    '''
#    print(SYS)
    os.system(SYS)
    lines = read_file(f'{lines_tmp}')
    lc0 = int(lines[0][0])
    lc1 = int(lines[1][0])
    os.system(f'rm -rf {lines_tmp}')
    return lc0,lc1

def grep_coord(ind,ffile,lc0,lc1):
    ncv = len(COMM_cv_ind)
    fn = f'{ffile}'
    crd_tmp1 = f'scratch/gcrd_cv_tmp{ind}'
    crd_tmp2 = f'scratch/gcrd_noncv_tmp{ind}'
    a0 = lc0+1
    a1 = a0+ncv-1
    b0 = lc0+1+ncv
    b1 = lc1-1
    SYS = f'''
sed '{a0},{a1}!d' {fn} | awk '{{ print $2,$3,$4 }}' > {crd_tmp1}
sed '{b0},{b1}!d' {fn} | awk '{{ print $2,$3,$4 }}' > {crd_tmp2}
    '''
#    print(SYS)
    os.system(SYS)

# cv coords    
    f1 = read_file(f'{crd_tmp1}')
    gcrd_cv = np.zeros((ncv,3))
    for i in range(ncv):
        for j in range(3):
            gcrd_cv[i][j] = float(f1[i][j])

# non-cv coords    
    f2 = read_file(f'{crd_tmp2}')
    natm = len(f2)
    ng16 = ncv+natm
    gcrd_noncv = np.zeros((natm,3))
    for i in range(natm):
        for j in range(3):
            gcrd_noncv[i][j] = float(f2[i][j])
    os.system(f'rm -rf {crd_tmp1} {crd_tmp2}')

# get coordinates of nearby non-cv atoms
    np.set_printoptions(suppress=True)
    pts = gcrd_noncv
    gidx = []
    for i in range(ncv):
        cen = gcrd_cv[i]
        T = KDTree(gcrd_noncv)
        idx = T.query_ball_point(cen,r=3.0)
        gidx.extend(idx)
#        print(idx)
    gidx = sorted(list(set(gidx)))
#    print(gidx)
    gcrd = gcrd_cv
    gcrd = np.append(gcrd,pts[gidx],axis=0)
#    print(gcrd)
    cidx = [i for i in range(ncv)]
    cidx1 = [gidx[i]+ncv for i in range(len(gidx))]
    cidx.extend(cidx1)
#    print(cidx)
    return ng16,gcrd,cidx

def get_stnd_crd_forces(ind,cfile,ffile):
# 1. coordinates of the atoms in pulays procedure
    crd = get_fm_coord(ind,cfile)

# 2. cartesian forces of the atoms in pulays procecure
    lof,hif = get_force(ind,ffile)

# 3. get the list of internal coordinates
    ic = get_ic(lof)
    return crd,lof,hif,ic

def update_forces(cidx,flo0,fhi0):
    natm = len(cidx)
#    print(cidx)
    flo = np.zeros((natm,3))
    fhi = np.zeros((natm,3))
    for i in range(natm):
        ii = cidx[i]
        flo[i] = flo0[ii]
        fhi[i] = fhi0[ii]
#    print(flo)
    return flo,fhi

def get_g16_crd_forces(ind,ffile):
# 1. get the line numbers for the coords in *.log
    lc0,lc1 = get_g16_lines(ind,ffile)

# 2. get atoms within the cutoff to the cv atoms
    ng16,gcrd,cidx = grep_coord(ind,ffile,lc0,lc1)

# 3. get the line numbers for flo and fhi in *.log
    iflo,ifhi = get_lines(ind,ffile)

# 4. get flo and fhi
    flo0,fhi0 = grep_forces(ind,ffile,iflo,ifhi,ng16)

# 5. get the forces of atoms within the cutoff
    flo,fhi = update_forces(cidx,flo0,fhi0)

# 6. convert the forces to kcal/mol/a
    lof,hif = fconvert(flo,fhi)

# 7. get the list of internal coordinates
    ic = get_ic(lof)
#    print(ic)
#    print(len(ic))
    return gcrd,lof,hif,ic
    
def read_ic():
    fn = 'user/geom_opt/ic.txt'
    f = read_file(fn)
    print(f)

def reorder_fint(fd,fint0,fint1):
    ncv = len(COMM_cvs)
    nfint = len(fint0[0])

    ifint0 = np.zeros((ncv,nfint))
    ifint1 = np.zeros((ncv,nfint))
    for i in range(ncv):
        kref = fd[i]

# internal forces to np array        
        klo0 = np.array(fint0[i])
        khi0 = np.array(fint1[i])

# get indices of sorted forces
        kdiff = abs(kref-klo0)
        kind = np.argsort(kdiff)

# reordered internal forces
        ifint0[i] = klo0[kind]
        ifint1[i] = khi0[kind]
    return ifint0,ifint1

def ml_prep(ind,iimg,iconfig,ffile,cfile,fd):

# 1. get the coordinates and forces and ic list
    crd,lof,hif,ic = get_stnd_crd_forces(ind,cfile,ffile)
#    crd,lof,hif,ic = get_g16_crd_forces(ind,ffile)
#    read_ic()

# 2. get cvs and internal forces
    rdist,fint0,fint1 = geom_opt.main(ic,crd,lof,hif)

# 3. reorder internal forces
    ifint0,ifint1 = reorder_fint(fd,fint0,fint1)

    nic = len(fint0[0])
    return nic,ind,iimg,iconfig,rdist,ifint0,ifint1

def get_ml_inputs(pd1,pd2,zdat,fdz,last,nfrm):
    fdir = f'{pd2}/log'
    nimg = len(zdat)
    nmd = len(zdat[0])
    ncv = len(COMM_cvs)

# 1 . prepare the output files
    fout = [
#        [f'{pd2}/data/input0.csv',1],
#        [f'{pd2}/data/target0.csv',1],
        [f'{pd2}/data/geom_opt.csv',5],
#        [f'{pd2}/data/cflo.csv',3],
#        [f'{pd2}/data/cfhi.csv',3]
    ]

# 2. print out the header
    for i in range(len(fout)):
        ii = fout[i]
        print_header(ii)

#    nimg,nmd,nfrm = 1,1,1
# 3. get cvs and calculate the internal forces
    start = time.perf_counter()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        ind = 0
        runs = []
        for j in range(nimg):
            iimg = j+1
            iconfig = 0
            for i in range(nmd):
                md = last-i
                zd = zdat[j][i]
                fd = fdz[j][i]

                for k in range(nfrm):
                    kk = k+1
                    ind += 1
                    iconfig += 1
                    cfile = f'{pd1}/img{iimg}/crd/md{md}/{kk}.crd'
                    ffile = f'{pd2}/log/{ind}.log'
#                    ml_prep(ind,iimg,iconfig,ffile,cfile)
                    cmd = executor.submit(ml_prep,ind,iimg,iconfig,ffile,cfile,fd)
                    runs.append(cmd)

# organize the concurrent results
        df = []
        for f in concurrent.futures.as_completed(runs):
            results = f.result()
            df.append(results)

# get the minimum number of ic's
        ic_min = 0
        for i in range(len(df)):
            ii = int(df[i][0])
            if (1/ii > ic_min):
                ic_min = ii
        
# store cv's and fint's for the ordered lambdas within the min
        rdist = np.zeros((ind,ncv))
        nsamp = nmd*nfrm
        jfint0 = np.zeros((nimg,ncv,nsamp,ic_min))
        jfint1 = np.zeros((nimg,ncv,nsamp,ic_min))
#        print(df[0])
        for i in range(len(df)):
            i0 = int(df[i][1])-1 # ind
            i1 = int(df[i][2])-1 # iimg
            i2 = int(df[i][3])-1 # iconfig

            # cv's
            for j in range(ncv):
                rdist[i0][j] = df[i][4][j]

                # fint's
                for k in range(ic_min):
                    jfint0[i1][j][i2][k] = df[i][5][j][k]
                    jfint1[i1][j][i2][k] = df[i][6][j][k]

    finish = time.perf_counter()
    print(f'{finish-start:.2f} second(s)')

# return cv's and fint's    
    return rdist,jfint0,jfint1

def genetic_algo(i,j,avg,fint0,fint1):
    #ilist0,ilist1,iout = mga.main(avg,fint0,fint1)
    ilist0,ilist1 = mga.main(avg,fint0,fint1)
    #print(i+1,j+1,iout)
    print('config:',i+1,'cv:',j+1)
    print('SE force:',ilist0,'AI force:',ilist1)
    return i,j,ilist0,ilist1

def fint_opt(jfint0,jfint1,fdz):
    nimg = len(jfint0)
    ncv = len(jfint0[0])
    nfrm = len(jfint0[0][0])
    print(jfint0.shape)
    nmd = len(fdz[0])

    avg_fd = np.zeros((nimg,ncv))
    for i in range(nimg):
        for j in range(nmd):
            for k in range(ncv):
                fd = fdz[i][j][k]
                avg_fd[i][k] += fd/nmd
    print(avg_fd)

    start = time.perf_counter()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        ind = 0
        runs = []
        for i in range(nimg):
            for j in range(ncv):
                avg = avg_fd[i][j]
                fint0 = jfint0[i][j]
                fint1 = jfint1[i][j]
#                print(avg,fint0.shape)
#                genetic_algo(i,j,avg,fint0,fint1)
                cmd = executor.submit(genetic_algo,i,j,avg,fint0,fint1)
                runs.append(cmd)

# organize the concurrent results
    df = []
    for f in concurrent.futures.as_completed(runs):
        results = f.result()
        df.append(results)

# organize the results to a 2d array
    kfint0 = np.zeros((nimg*nfrm,ncv))
    kfint1 = np.zeros((nimg*nfrm,ncv))

    # iterate over the concurrent results
    for i in range(len(df)):
        i0 = int(df[i][0]) # nimg
        i1 = int(df[i][1]) # ncv
        ii = i0*nfrm

        # assign fint's to each sample/cv
        for j in range(nfrm):
            jj = ii+j
            kfint0[jj][i1] = df[i][2][j]
            kfint1[jj][i1] = df[i][3][j]
            
    finish = time.perf_counter()
    print(f'{finish-start:.2f} second(s)')

# return the optimized fint's for all configurations
    return kfint0,kfint1

def ml_print(pd2,rdist,kfint0,kfint1):
    nconfig = len(kfint0)
    ncv = len(kfint0[0])
    fcorr = kfint1-kfint0

    o0 = f'{pd2}/data/geom_opt.csv'
    o1 = f'{pd2}/data/input0.csv'
    o2 = f'{pd2}/data/target.csv'

    f0 = open(o0, 'w')
    f1 = open(o1, 'w')
    f2 = open(o2, 'w')

    f0.write(f'ind')
    for j in range(ncv):
        jj = j+1
        j0 = f',r{jj},flo{jj},fhi{jj},fcorr{jj}'
        j1 = f'i{jj}'
        j2 = f't{jj}'

        f0.write(f'{j0}')
        f1.write(f'{j1}') if (j == 0) else f1.write(f',{j1}')
        f2.write(f'{j2}') if (j == 0) else f2.write(f',{j2}')
    f0.write('\n')
    f1.write('\n')
    f2.write('\n')

    for i in range(nconfig):
        ii = i+1
        f0.write(f'{ii}')
        for j in range(ncv):
            j1 = rdist[i][j]
            flo = kfint0[i][j]
            fhi = kfint1[i][j]
            j2 = fhi-flo

            f0.write(f',{j1},{flo},{fhi},{j2}')
            f1.write(f'{j1}') if (j == 0) else f1.write(f',{j1}')
            f2.write(f'{j2}') if (j == 0) else f2.write(f',{j2}')
        f0.write('\n')
        f1.write('\n')
        f2.write('\n')
    f0.close()
    f1.close()
    f2.close()

################
## MAIN BLOCK ##
################

def main(pd1,pd2,last,nmd):
    dir2 = f'{pd2}/data'
    CC = f'''\
    rm -rf {dir2}
    mkdir {dir2}
    '''
    print(CC)
    os.system(CC)

# 1. get nfrm
#    nsavc,nfrm = xtrap.get_isim(pd1)
    nsavc,nfrm = 1,1

# 2. get resd distances and mean forces for cvs
    zdat,fdz = get_zcv(pd1,last,nmd)

# 3. get the cvs and calculate internal forces
    rdist,jfint0,jfint1 = get_ml_inputs(pd1,pd2,zdat,fdz,last,nfrm)

# 4. optimize the internal forces to mean forces in cvs
    kfint0,kfint1 = fint_opt(jfint0,jfint1,fdz)
#    print(kfint0,kfint1)

# 5. print out results for ml
    ml_print(pd2,rdist,kfint0,kfint1)
    print('Target force corrections in:\n\t2_method/0strg/data/target.csv')
    
##############
## END MAIN ##
##############

