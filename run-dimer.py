import fftoolbox as fftb
import numpy as np
from scipy.optimize import curve_fit


import logging, sys


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

    data1['frame'] = [d['coordinates'] for d in data1['atoms']]
    molecule1 = fftb.LebedevMolecule(data1)

    data2['frame'] = [d['coordinates'] for d in data2['atoms']]
    molecule2 = fftb.LebedevMolecule(data2)

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
#plt.show()

