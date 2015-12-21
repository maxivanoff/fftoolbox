import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

molecule_name = 'water'

data = {
        'name': molecule_name,
        'theory': 'b3lyp_augccpvdz',
        'sphere params': (2, 0.5),
        }

parser = fftb.Gaussian()
parser.read_file(filename='data/dimer/dimer.log')

parser.write_to_data('multipoles', fftb.GDMA(data=data).multipoles)

data1, data2 = parser.split_data([1,2,3], [4,5,6])
data1.update(data)
data2.update(data)

molecule1 = fftb.LebedevMolecule(data1)
molecule2 = fftb.LebedevMolecule(data2)

print molecule1.charges
for atom in molecule1.atoms:
    print atom
print molecule2.charges
for atom in molecule2.atoms:
    print atom

dimer = fftb.Complex('water_dimer', molecule1, molecule2)
dimer.write_xyz(filename='dimer_spheres.xyz', here=True)

