import fftoolbox as fftb
from glob import iglob

reports = fftb.Reports()

for filename in iglob('./data/cub/*mp2_augccpvtz_d1.5.cub'):
    molecule_name = filename.split('/')[-1].split('_mp2_augccpvtz_d1.5.')[0]

    data = {
            'name': molecule_name,
            'theory': 'mp2_augccpvtz',
            'density': 1.5,
            'symmetry': False,
            'representation': ('cartesian', 1)
            }

    parser = fftb.GaussianCube(data=data)
    data.update(parser.data.copy())

    grid = fftb.vdwGrid(data)
    molecule = fftb.HM(data)

    ls = fftb.LSC(molecule=molecule, grid=grid)
    ls.solve()

    reports.add(ls)

print reports
