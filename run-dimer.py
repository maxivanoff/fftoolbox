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
frame1 = [molecule1.get_atom_by_index(i).coordinates for i in [2, 1, 3]]
molecule1.align_with_frame(frame1)

molecule2 = fftb.LebedevMolecule(data2)
frame2 = [molecule2.get_atom_by_index(i).coordinates for i in [5, 4, 6]]
molecule2.align_with_frame(frame1)

print molecule1.charges
for atom in molecule1.atoms:
    print atom
print molecule2.charges
for atom in molecule2.atoms:
    print atom

dimer = fftb.Complex('water_dimer', molecule1, molecule2)
dimer.write_xyz(filename='dimer_spheres.xyz', here=True)

print dimer.point_charge_energy()

