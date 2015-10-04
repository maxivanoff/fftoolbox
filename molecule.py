from atom import Atom, Site
from groups import BuriedGroup, CarbonylGroup
from multipole import Multipole
from multipole import GroupOfAtoms

import numpy as np
import logging
import os

molLogger = logging.getLogger('molecule')
WORKDIR = os.path.dirname(__file__)

class Molecule(Multipole):
    
    def __init__(self, data):
        self.data = data
        try:
            self.mults = data['multipoles']
        except KeyError:
            self.mults = ('cartesian', 0)
        try:
            self.sym = data['symmetry']
        except KeyError:
            self.sym = False
        try:
            self.hybrids = data['hybridizations']
        except KeyError:
            self.hybrids = []
        molLogger.info('Start assembling %s Molecule\nsymmetry: %r\nmultipoles: %r\nhybridizations: %r' % (data['name'], self.sym, self.mults, self.hybrids))
        Multipole.__init__(self, name=data['name'], multipoles=self.mults)
        self.set_atoms(data['atoms'])
        self.set_groups()
        self.set_sym_sites()
        self.set_multipole_matrix()
        molLogger.info('%s Molecule is created' % (self.name))

    def copy(self):
        for i, a in enumerate(self.atoms):
            element, crds, charge = self.data['atoms'][i]
            a_updated = (element, a.coordinates, charge)
            self.data['atoms'][i] = a_updated
        molecule = Molecule(data=self.data)
        return molecule

    def get_ep_data(self):
        data = {}
        for a in self.atoms:
            if len(a.sites) > 1:
                data[a.element] = a.get_ep_data()
        return data
    
    def set_atoms(self, atoms):
        self.atoms = []
        for i, a in enumerate(atoms):
            element, crds, charge = a
            if element[0] == 'X' or element[0] == 'E': # extra point connected to previous atom
                s = Site(coordinates=crds, name='EP_%s' % atom.element, charge=charge)
                atom.sites.append(s)
            else: # regular atom
                atom = Atom(name=element, element=element, coordinates=crds, index=i+1, charge=charge)
                self.atoms.append(atom)

    def set_groups(self):
        """
        define bonded neighboring atoms for each atom
        assign new names to atoms according to the functional group: 
            carbonyl, methyl groups are implemented
        """
        # identify functional groups
        buried_atoms = []
        carbonyl_carbons = []
        for atom in self:
            atom.neighbors = filter(atom.bonded_to, self.atoms)
            if len(filter(lambda a: a.element=='H', atom.neighbors))>1:
                buried_atoms.append(atom)
            if len(filter(lambda a: a.element=='O', atom.neighbors))==1 and len(atom.neighbors)==3:
                carbonyl_carbons.append(atom)
            molLogger.debug("Neighbors are assigned to %s:\n%s" % (atom.name, atom.neighbors))

        # set extra points
        for atom in self.atoms:
            if atom.element in self.hybrids:
                atom.set_hybridization(self.hybrids[atom.element])
        
        # create groups
        exclude = ['methane', 'benzene', 'methane', 'tip5p', 'tip3p', 'ammonia', 'water']
        groups = []
        if not self.name in exclude:
            for count, atom in enumerate(carbonyl_carbons):
                groups.append(CarbonylGroup(center=atom, count=count, sym=self.sym))
            for count, atom in enumerate(buried_atoms):
                groups.append(BuriedGroup(center=atom, count=count, sym=self.sym))

        
        # symmetrical atoms have identical names
         
        try:
            template, symmetrical = self.sym[0], self.sym[1:] 
            t = self.get_atom_by_index(template)
            for i in symmetrical:
               a = self.get_atom_by_index(i)
               a.name = t.name
        except TypeError:
            pass

        molLogger.info("Names of equivalent sites: %s" % self.sites_names_eq)
        molLogger.info("Names of non-equivalent sites: %s" % self.sites_names_noneq)

    def __add__(self, molecule):
        complex = Complex()
        complex.add_molecule(self.copy())
        complex.add_molecule(molecule.copy())
        return complex

class Complex(GroupOfAtoms):

    def __init__(self):
        GroupOfAtoms.__init__(self)
        self.molecules = []

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    @property
    def atoms(self):
        atoms = []
        for mol in self.molecules:
            atoms += mol.atoms
        return atoms

    @property
    def sites(self):
        sites = []
        for mol in self.molecules:
            sites += mol.sites
        return sites

