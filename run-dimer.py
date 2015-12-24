import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit


import logging, sys

def point_charge_energy(dimer):
    m1, m2 = dimer.molecules
    e = 0.
    for s1 in m1.sites:
        for s2 in m2.sites:
            r = s1.distance_to(s2)
            e += s1.charge*s2.charge/r
    return e

def m_to_xyz(m):
    if m == -1:
        return 1 # y
    if m == 0:
        return 2 # z
    if m == 1:
        return 0 # x

def T(l1, l2, m1, m2, r1, r2, c):
    if l1 == l2 == 0:
        return 1.
    if l1 == 0 and l2 == 1:
        alpha = m_to_xyz(m2)
        return r2[alpha]
    if l1 == 1 and l2 == 0:
        alpha = m_to_xyz(m1)
        return r1[alpha]
    if l1 ==1 and l2 == 1:
        alpha = m_to_xyz(m1)
        beta = m_to_xyz(m2)
        return 3*r1[alpha]*r2[beta] + c[alpha, beta]

def multipole_energy(dimer, rank):
    mol1, mol2 = dimer.molecules
    # pre-compute stuff for interaction function T
    R_vec = mol1.center_of_mass() - mol2.center_of_mass()
    R = np.linalg.norm(R_vec)
    e12 = R_vec/R
    r1 = []
    r2 = []
    for e in mol1.frame.local_axes.T:
        r1.append(np.dot(e, e12))
    for e in mol2.frame.local_axes.T:
        r2.append(-np.dot(e, e12))
    c = np.zeros((3,3))
    for i, e1 in enumerate(mol1.frame.local_axes.T):
        for j, e2 in enumerate(mol2.frame.local_axes.T):
            c[i,j] = np.dot(e1, e2)
    # compute energy
    U = 0.
    for l1 in xrange(rank+1):
        for m1 in xrange(-l1, l1+1):
            for l2 in xrange(rank+1):
                for m2 in xrange(-l2, l2+1):
                    Qlm_1 = mol1.get_reference_Qlm(l1, m1)
                    Qlm_2 = mol2.get_reference_Qlm(l2, m2)
                    U += Qlm_1*Qlm_2*pow(R, -l1 -l2 -1)*T(l1,l2,m1,m2,r1,r2,c)
    return U


logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

molecule_name = 'water'
rank = 1
radius = 0.5

data = {
        'name': molecule_name,
        'theory': 'b3lyp_augccpvdz',
        'sphere params': (rank, radius),
        }

pc_energies = []
mult_energies = []

for i in xrange(1, 11):
    parser = fftb.XYZ()
    parser.read_file(filename='data/dimer/dimer-%i.xyz' % i)

    parser.write_to_data('multipoles', fftb.GDMA(data=data).multipoles)

    data1, data2 = parser.split_data([1,2,3], [4,5,6])
    data1.update(data)
    data2.update(data)

    print data1
    print data2

    molecule1 = fftb.LebedevMolecule(data1)
    frame1 = [molecule1.get_atom_by_index(i).coordinates for i in [2, 1, 3]]
    molecule1.align_with_frame(frame1)

    molecule2 = fftb.LebedevMolecule(data2)
    frame2 = [molecule2.get_atom_by_index(i).coordinates for i in [5, 4, 6]]
    molecule2.align_with_frame(frame2)

    print molecule1.charges
    for atom in molecule1.atoms:
        print atom
    print molecule2.charges
    for atom in molecule2.atoms:
        print atom

    dimer = fftb.TwoLebedevMolecules('water_dimer', molecule1, molecule2)
    dimer.write_xyz(filename='dimer_spheres.xyz', here=True)


    epc = dimer.point_charge_energy()
    emult = dimer.multipole_energy(rank)

    pc_energies.append(epc)
    mult_energies.append(emult)

x = pc_energies
y =  mult_energies

import matplotlib.pyplot as plt
#import seaborn

plt.ylabel('Multipolar Energies, au')
plt.xlabel('Point Charge Energies, au')
plt.title('Point charge vs multipoles energies of water dimer')
xlim = [min(x), max(x)]

correlation = np.corrcoef(x,y)[0,1]
R2 = correlation**2
popt, pcov = curve_fit(lambda xdata,m,n: m*xdata+n, np.array(x), np.array(y))
a, b = popt
da, db = np.sqrt(np.diag(pcov))
xx = np.arange(xlim[0], xlim[1], 0.01)
yy = a*xx + b

plt.plot(xx, yy, '-', label='y = %.3f x + %.3f\nR2 = %.2f' % (a, b, R2))
plt.plot(x, y, 'o')
plt.legend(loc='upper right')

plt.savefig('dimer-energies.pdf')
plt.show()

