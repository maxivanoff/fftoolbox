import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

amber = {
        'O-sp3':('sp3', 0.35, 109.5/2.),
        'S-sp3':('sp3', 0.70, 90.),
        'O-sp2': ('sp2', 0.35, 60.0),
        'N-sp2': ('sp2', 0.35, 0.0)
        }

data = {
        'name': 'cf3sh',
        'theory': 'mp2_augccpvtz',
        'density': 1.5,
        'symmetry': True,
        'hybridizations': {'S': amber['S-sp3']},
        'representation': ('cartesian', 2)
        }

results = fftb.Results(fftb.Gaussian(data=data).data['multipoles'])

parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())

grid = fftb.vdwGrid(data)
molecule = fftb.HybridMolecule(data)

ls = fftb.LeastSquaresCharges(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)

print results

# save
molecule.write_xyz(filename='amber-%s_%s.xyz' % (data['name'], data['theory']))
molecule.write_mol2(filename='amber-%s_%s.mol2' % (data['name'], data['theory']))

ff = fftb.ForceFieldXML()
ff.write_file(molecule=molecule, xmlfilename='data/forcefields/amber-%s_%s.xml' % (data['name'], data['theory']))
ff.load_forcefields(filename='amber-%s_%s.xml' % (data['name'], data['theory']), molecule=molecule)


