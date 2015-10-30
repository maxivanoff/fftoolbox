from sphere import LebedevSphere
from fftoolbox.molecule import Molecule
from fftoolbox.atom import Atom

import logging

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class LebedevAtom(Atom, LebedevSphere):

    def __init__(self, element=None, coordinates=None, index=None, ref_multipoles=None, rank=None, radius=None):
        Atom.__init__(self, index=index, element=element, coordinates=coordinates, representation=None)
        LebedevSphere.__init__(self, name=self.name, rank=rank, radius=radius, \
                origin=coordinates, ref_multipoles=ref_multipoles)

class LebedevMolecule(LebedevSphere):

    def __init__(self, data):
        rank, radis = data['sphere params']
        representation = ('spherical', rank)
        LebedevSphere.__init__(self, name=data['name'], rank=rank, radius=radius, \
                origin=None, ref_multipoles=data['multipoles'])

class DistributedLebedevMolecule(Molecule):

    def __init__(self, data):
        self.sphere_params = data['sphere params']
        Molecule.__init__(self, data)

    def add_atoms(self, atoms):
        for a in atoms:
            index, element, crds, multipoles = a
            rank, radius = self.sphere_params[element]
            atom = LebedevAtom(element=element, coordinates=crds, ref_multipoles=multipoles, \
                    index=index, rank=rank, radius=radius)
            self.add_atom(atom)
        for atom in self.atoms:
            pass
            


