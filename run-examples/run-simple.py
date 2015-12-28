import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'cf3sh',
        'theory': 'mp2_augccpvtz',
        'density': 1.5,
        'symmetry': [4,5,6],
        'exclude': [],
        'representation': ('cartesian', 1)
        }

results = fftb.Results(fftb.Gaussian(data=data).data['multipoles'])

parser = fftb.GaussianCube(data=data)
data.update(parser.data.copy())

grid = fftb.vdwGrid(data)
molecule = fftb.HM(data)

ls = fftb.LSC(grid=grid, molecule=molecule) 
ls.solve()
results.add(ls)

# save
molecule.write_xyz(filename='%s_%s.xyz' % (data['name'], data['theory']))
molecule.write_mol2(filename='%s_%s.mol2' % (data['name'], data['theory']))

#ff = fftb.ForceFieldXML()
#ff.write_file(molecule=molecule)
#ff.load_forcefields(filename='%s_%s.xml' % (data['name'], data['theory']), molecule=molecule)

print results
print fftb.Report(data['name'], ls)

for at_i, at_grid in grid.atomic_grids.items():
    ls.update_grid(at_grid)
    print fftb.Report('%s, atom %i %s' % (data['name'], at_i, at_grid.element), ls)

