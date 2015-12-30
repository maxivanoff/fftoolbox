import fftoolbox as fftb

import logging, sys

logger = logging.getLogger(__name__)
lf = '%(levelname)s: %(funcName)s at %(filename)s +%(lineno)s\n%(message)s\n'
logging.basicConfig(level=logging.DEBUG, format=lf)

data = {
        'name': 'water',
        'theory': 'mp2_augccpvtz',
        'hybridizations': {'O': ('sp3', 0.70, 109.47/2.)},
        'representation': ('cartesian', 2)
        }


parser = fftb.XYZ(data=data)
data.update(parser.data.copy())

molecule = fftb.HybridMolecule(data)

for site in molecule.get_sites_by_element('O'):
    site.set_charge(0.0)
    site.set_epsilon(0.16)
    r0 = pow(2, 1./6.)*3.12/2.
    site.set_r0(r0)
for site in molecule.get_sites_by_element('H'):
    site.set_charge(0.241)
    site.set_epsilon(0.0)
    site.set_r0(0.0)
for site in molecule.get_sites_by_element('EP'):
    site.set_charge(-0.241)
    site.set_epsilon(0.0)
    site.set_r0(0.0)
print molecule

molecule.write_mol2(filename='tip5p.mol2', here=True)

ff = fftb.ForceFieldXML()
ff.write_file(molecule=molecule, xmlfilename='tip5p.xml')
ff.load_forcefields(filename='tip5p.xml', molecule=molecule)


