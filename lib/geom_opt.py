import numpy as np
from numpy import linalg as LA
from common import *

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

def read_ic():
    fn = 'user/geom_opt/ic.txt'
    f = read_file(fn)
    nint = len(f)

    ic_type = [0 for i in range(nint)]
    ic = [[0 for j in range(4)] for i in range(nint)]
    for i in range(nint):
        ic_type[i] = int(f[i][0])
        for j in range(4):
            jj = j+1
            ic[i][j] = int(f[i][jj])
    return ic_type,ic

def calc_dr(c1,c2):
    rsum = 0.0
    dr12 = [0.0 for i in range(3)]
    for i in range(3):
        dr12[i] = c1[i]-c2[i]
        rsum += (dr12[i])**2
    r12 = rsum**(1/2)
    return dr12,r12

def bmat_bond(dr12,r12):
    dr21 = [0.0 for i in range(3)]
    evec = [0.0 for i in range(6)]

    for i in range(3):
        dr21[i] = -dr12[i]
    for i in range(3):
        ii = i
        evec[ii] = dr12[i]/r12
    for i in range(3):
        ii = i+3
        evec[ii] = dr21[i]/r12
    return evec

def calc_ang(dr12,r12,dr32,r23):
    dp2 = 0
    for i in range(3):
        i1 = dr12[i]
        i2 = dr32[i]
        dp2 += i1*i2
    cos2 = dp2 / (r12*r23)
    sin2 = (1-(cos2**2))**(1/2)
    return cos2,sin2

def bmat_ang(dr12,r12,dr23,r23,cos2,sin2):
    dr32 = [0.0 for i in range(3)]
    for i in range(3):
        dr32[i] = -dr23[i]

    evec = [0.0 for i in range(9)]
    for i in range(3):
        i1 = i
        evec[i1] = (cos2*dr12[i]/r12 - dr32[i]/r23) / r12*sin2
    for i in range(3):
        i3 = i+6
        evec[i3] = (cos2*dr32[i]/r23 - dr12[i]/r12) / r23*sin2
    for i in range(3):
        i1 = i
        i2 = i+3
        i3 = i+6
        evec[i2] = -(evec[i1] + evec[i3])
    return evec

def calc_cp(dr21,r12,dr32,r23):
    e12 = [0.0 for i in range(3)]
    e23 = [0.0 for i in range(3)]
    cp2 = [0.0 for i in range(3)]
    for i in range(3):
        e12[i] = dr21[i]/r12
        e23[i] = dr32[i]/r23

    cp2[0] = (e12[1]*e23[2]-e12[2]*e23[1])
    cp2[1] = (e12[0]*e23[2]-e12[2]*e23[0])*-1
    cp2[2] = (e12[0]*e23[1]-e12[1]*e23[0])
    return cp2

def bmat_dihe(cp2,cp3,cos2,sin2,cos3,sin3,r12,r23,r34):
    r32 = r23
    r43 = r34
    evec = [0.0 for i in range(12)]
    for i in range(3):
        ii = i
        evec[ii] = -cp2[i]/r12*sin2**2
    for i in range(3):
        ii = i+3
        evec[ii] = ((r23 - r12*cos2) / r23*r12*sin2) * (cp2[i]/sin2) - (cos3 / r23*sin3) * (cp3[i]/sin3)
    for i in range(3):
        ii = i+6
        evec[ii] = ((r32 - r43*cos3) / r32*r43*sin3) * (cp3[i]/sin3) - (cos2 / r32*sin2) * (cp2[i]/sin2)
    for i in range(3):
        ii = i+9
        evec[ii] = -cp3[i]/r34*sin3**2
    return evec

def get_bmat0(ic,crd,ic_type):
    nint = len(ic)
    natm = len(crd)
    n3 = natm*3
    ncv = len(COMM_cvs)

    bmat = np.zeros((nint,n3))
#    rdist = np.zeros(ncv)
    rdist = [0.0 for i in range(ncv)]
    for n in range(nint):
        iic = ic_type[n]

# atom indices
        a1 = ic[n][0]-1
        a2 = ic[n][1]-1
        a3 = ic[n][2]-1
        a4 = ic[n][3]-1
# column indices
        ind1,ind2,ind3,ind4 = a1*3,a2*3,a3*3,a4*3

        # bond stretch
        if (iic >= 2):
            c1 = [0.0 for i in range(3)]
            c2 = [0.0 for i in range(3)]
            for j in range(3):
                c1[j] = crd[a1][j]
                c2[j] = crd[a2][j]

            dr12,r12 = calc_dr(c1,c2)
            dr21 = [-dr12[i] for i in range(3)]

        # cv distance
        if (n < ncv):
            rdist[n] = r12

        # b-matrix for bond stretch
        if (iic == 2):
            evec = bmat_bond(dr12,r12)
            for j in range(3):
                j1,j2 = ind1+j,ind2+j
                bmat[n][j1] = evec[j+0]
                bmat[n][j2] = evec[j+3]

        # angle bend
        if (iic >= 3):
            c3 = [0.0 for i in range(3)]
            for j in range(3):
                c3[j] = crd[a3][j]

            dr23,r23 = calc_dr(c2,c3)
            dr32 = [-dr23[i] for i in range(3)]
            cos2,sin2 = calc_ang(dr12,r12,dr32,r23)

        # b-matrix for angle bend
        if (iic == 3):
            evec = bmat_ang(dr12,r12,dr23,r23,cos2,sin2)
            for j in range(3):
                j1,j2,j3 = ind1+j,ind2+j,ind3+j
                bmat[n][j1] = evec[j+0]
                bmat[n][j2] = evec[j+3]
                bmat[n][j3] = evec[j+6]

        # torsion
        if (iic == 4):
            c4 = [0.0 for i in range(3)]
            for j in range(3):
                c4[j] = crd[a4][j]

            dr34,r34 = calc_dr(c3,c4)
            dr43 = [-dr34[i] for i in range(3)]
            cos3,sin3 = calc_ang(dr43,r34,dr23,r23)
            cp2 = calc_cp(dr21,r12,dr32,r23)
            cp3 = calc_cp(dr34,r34,dr23,r23)

        # b-matrix for torsion
        if (iic == 4):
            evec = bmat_dihe(cp2,cp3,cos2,sin2,cos3,sin3,r12,r23,r34)
            for j in range(3):
                j1,j2,j3,j4 = ind1+j,ind2+j,ind3+j,ind4+j
                bmat[n][j1] = evec[j+0]
                bmat[n][j2] = evec[j+3]
                bmat[n][j3] = evec[j+6]
                bmat[n][j4] = evec[j+9]
    return bmat,rdist

def get_bmat(ic,crd):
    nint = len(ic)
    natm = len(crd)
    n3 = natm*3
    ncv = len(COMM_cvs)

    bmat = np.zeros((nint,n3))
    rdist = np.zeros(ncv)
    for i in range(nint):
# atom indices
        a1 = ic[i][0]-1
        a2 = ic[i][1]-1

# bmat column indices
        i1 = a1*3
        i2 = a2*3

# cartesian displacement and bond distance for a1 and a2
        disp = crd[a1]-crd[a2]
        sum_sq = np.dot(disp.T, disp)
        dist = np.sqrt(sum_sq)

        if (i < ncv):
            rdist[i] = dist

# construct the b-matrix elements for bond stretches
        for j in range(3):
            j1 = i1+j
            j2 = i2+j
            bmat[i][j1] = disp[j]/dist
            bmat[i][j2] = -disp[j]/dist
    return bmat,rdist

def get_gmat(bmat):
    bmat_t = bmat.transpose()
    gmat = bmat.dot(bmat_t)
    return gmat

def get_linalg(gmat):
# solve for the eigenvalues (w) and eigenvectors (v)
    w,v = LA.eig(gmat)

# save only the real eigenvalues and eigenvectors for a symmetric matrix
    lam = w.real
    klmat = v.real
    klmat_t = klmat.transpose()

# inverse the eigenvalues
    rlam = np.reciprocal(lam)
    ilam = rlam.tolist()
    return ilam,klmat

def get_ginv(gmat,ilam,klmat,nlam):
    nint = len(gmat)

# select the nonzero eigenvalues        
    nonz = np.zeros(nint)
    for i in range(nlam):
        nonz[i] = ilam[i]

# diagonalize the nonzero eigenvalues        
    nonz_diag = np.diag(nonz)
        
# construct the g-inverse matrix
    klmat_nonz_diag = klmat.dot(nonz_diag)
    klmat_t = klmat.transpose()
    ginv = klmat_nonz_diag.dot(klmat_t)
    return ginv

def get_fint(bmat,ginv,lof,hif):
    nint = len(bmat)
    natm = len(lof)
    n3 = natm*3

# reshape cartesian forces
    flo = lof.reshape((n3,1))
    fhi = hif.reshape((n3,1))

# form the l.h.s. of the gradient vector
    ginv_bmat = ginv.dot(bmat)

# get the internal forces
    iflo = ginv_bmat.dot(flo).flatten()
    ifhi = ginv_bmat.dot(fhi).flatten()
    return iflo,ifhi

def geom_opt(bmat,gmat,ilam,klmat,lof,hif):
    nint = len(bmat)
    natm = len(lof)
    ncv = len(COMM_cvs)
#    imax = 3*natm-6-6+1 # 3 rot + 3 trans
    imax = nint

    fint0 = np.zeros((imax,ncv))
    fint1 = np.zeros((imax,ncv))

# calculate the internal force with different g-inverses
    nzero = [0 for i in range(ncv)]
    for i in range(imax):
#        nlam = i+6 # 3 rot + 3 trans
        nlam = i+1
        ginv = get_ginv(gmat,ilam,klmat,nlam)
        iflo,ifhi = get_fint(bmat,ginv,lof,hif)
        for j in range(ncv):
            jj = nzero[j]
            j0 = iflo[j]
            j1 = ifhi[j]
            if (abs(j0) <= 3 and jj == 0):
                j0 = 1e6
                fint0[i][j] = j0
                fint1[i][j] = j1
            else:
                nzero[j] = 1
                fint0[i][j] = j0
                fint1[i][j] = j1
    fint0_t = fint0.T
    fint1_t = fint1.T
    fint0 = fint0_t.tolist()
    fint1 = fint1_t.tolist()
    return fint0,fint1

################
## MAIN BLOCK ##
################

def main(ic,crd,lof,hif):
    ic_type,tmp = read_ic()
    bmat,rdist = get_bmat0(tmp,crd,ic_type)

# 1. get the b-matrix
#    bmat,rdist = get_bmat(ic,crd)

# 2. get the g-matrix
    gmat = get_gmat(bmat)

# 3. get the inverse eigenvalues (ilam) and eigenvectors (klmat)
    ilam,klmat = get_linalg(gmat)

# 4. get the internal forces from g-inverse
    fint0,fint1 = geom_opt(bmat,gmat,ilam,klmat,lof,hif)
    return rdist,fint0,fint1

##############
## END MAIN ##
##############

