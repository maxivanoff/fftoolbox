import numpy as np
import logging
import os

import units
from bonds import Bond
from atom import HybridAtom, FFSite, Atom
from multipole import Multipole
from parser import ForceFieldXML
from multipole import GroupOfSites
from groups import BuriedGroup, CarbonOxygenGroup

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)
WORKDIR = os.path.dirname(__file__)

class AtomsInMolecule(object):

    masses = {'H': 1.0079, 'C': 12.0107, 'O': 15.9994, 'S': 32.065}
    
    def __init__(self, name=None, atoms=None, sym=False):
        self._name = name
        self._atoms = []
        self.sym = sym
        if atoms:
            self.add_atoms(atoms)
            self.set_groups()
            logger.info("Names of equivalent atoms: %s" % self.atoms_names_eq)
            logger.info("Names of non-equivalent atoms: %s" % self.atoms_names_noneq)

    @property
    def name(self):
        return self._name

    def add_atom(self, atom):
        self._atoms.append(atom)

    def add_atoms(self, atoms):
        extra_points = []
        for i, atom in enumerate(atoms):
            index = atom['index']
            element = atom['element']
            crds = atom['coordinates']
            multipoles = atom['multipoles']
            atom = HybridAtom(index=index, element=element, coordinates=crds, 
                    multipoles=multipoles, sym=self.sym)
            self.add_atom(atom)

    @property
    def num_atoms(self):
        return len(self._atoms)

    @property
    def atoms(self):
        return iter(sorted(self._atoms, key=lambda a: a.index))

    def get_atoms_by_name(self, name):
        return iter(filter(lambda s: s.name == name, self.atoms))

    def get_atoms_by_element(self, element):
        return iter(filter(lambda s: s.element == element, self.atoms))

    @property
    def atoms_noneq(self):
        atoms = []
        for s in self.atoms:
            if not s.name in [ss.name for ss in atoms]:
                atoms.append(s)
        return iter(atoms)

    def get_atom_by_index(self, index):
        return next(a for a in self.atoms if a.index==index)

    @property
    def atoms_names_noneq(self):
        return [a.name for a in self.atoms_noneq]

    @property
    def atoms_names_eq(self):
        return [a.name for a in self.atoms]

    def center_of_mass(self):
        M = sum([self.masses[atom.element] for atom in self.atoms])
        MR = np.zeros(3)
        for atom in self.atoms:
            MR += self.masses[atom.element]*atom.coordinates
        return MR/M

    def set_groups(self):
        """
        define bonded neighboring atoms for each atom
        assign new names to atoms according to their functional group
        """
        # identify functional groups
        buried_atoms = []
        carbons_with_oxygen = []
        for atom in self.atoms:
            atom.set_neighbors(self.atoms)
            if len(filter(lambda a: a.element=='H', atom.neighbors))>1:
                buried_atoms.append(atom)
            if not len(filter(lambda a: a.element=='O', atom.neighbors))==0 and atom.element=='C' and len(atom.neighbors)==3:
                carbons_with_oxygen.append(atom)
            logger.debug("Neighbors are assigned to %s:\n%s" % (atom.name, atom.neighbors))
        # create groups
        # group instances just change the atoms names
        exclude = ['methane', 'benzene', 'methane', 'tip5p', 'tip3p', 'ammonia', 'water']
        self.buried_groups = []
        if not self.name in exclude:
            for count, atom in enumerate(carbons_with_oxygen):
                CarbonOxygenGroup(center=atom, count=count, sym=self.sym)
            for count, atom in enumerate(buried_atoms):
                g = BuriedGroup(center=atom, count=count, sym=self.sym)
                self.buried_groups.append(g)
        # symmetrical atoms have identical names
        try:
            template, symmetrical = self.sym[0], self.sym[1:] 
            t = self.get_atom_by_index(template)
            for i in symmetrical:
               a = self.get_atom_by_index(i)
               a.set_name(t.name)
        except TypeError:
            pass

    def copy(self):
        # take coordinates of the current molecule
        # and create a new one
        data = self.data.copy()
        data['atoms'] = []
        for i, a in enumerate(self.atoms):
            index, element, crds, mult = self.data['atoms'][i]
            a_updated = (index, element, a.coordinates, mult)
            self.data['atoms'].append(a_updated)
        molecule = Molecule(data=data)
        return molecule

    @property
    def bonds(self):
        bonds = []
        for atom in self.atoms:
            for nghbr in atom.neighbors:
                bond = Bond(atom, nghbr)
                if not bond in bonds:
                    bonds.append(bond)
        return bonds

class MoleculeWithFrames(AtomsInMolecule):

    def __init__(self, name=None, atoms=None, sym=None, framed_atoms=None):
        if framed_atoms is None: framed_atoms = []
        AtomsInMolecule.__init__(self, name=name, atoms=atoms, sym=sym)
        self.frames = []
        self.set_frames(framed_atoms)

    def set_frames(self, framed_atoms):
        self.frames = []
        for atom in self.atoms:
            if atom.element in framed_atoms:
                atom.set_frame()
                self.frames.append(atom.frame)


class HybridMolecule(Multipole, MoleculeWithFrames):
    
    def __init__(self, data):
        self.data = data
        try:
            name = data['name']
        except KeyError:
            name = None
        try:
            self.theory = data['theory']
        except KeyError:
            self.theory = None
        try:
            representation = data['representation']
        except KeyError:
            representation = None
        try:
            sym = data['symmetry']
        except KeyError:
            sym = False
        try:
            self.energy = data['energy']
        except KeyError:
            self.energy = None
        try:
            atoms = data['atoms']
        except KeyError:
            atoms = None
        try:
            self.hybrid_atoms = data['hybridizations']
        except KeyError:
            self.hybrid_atoms = {}
        logger.info('Start assembling %s Molecule\ntheory: %r\nsymmetry: %r\nrepresentation: %r' % (name, self.theory, sym, representation))
        Multipole.__init__(self, name=name, representation=representation)
        MoleculeWithFrames.__init__(self, name=name, atoms=atoms, sym=sym, framed_atoms=self.hybrid_atoms.keys())
        self.set_hybridizations()
        self.set_frames_from_sites() # in case extra points are loaded from the file
        self.set_sym_sites()
        if representation is not None:
            self.set_multipole_matrix()
        logger.info("Names of equivalent sites: %s" % self.sites_names_eq)
        logger.info("Names of non-equivalent sites: %s" % self.sites_names_noneq)
        logger.info('%s Molecule is created' % (self.name))

    @property
    def sites(self):
        sites = []
        for atom in self.atoms:
            for site in atom.sites:
                sites.append(site)
        return iter(sorted(sites, key=lambda s: s.index))

    @property
    def num_sites(self):
        sites = []
        for atom in self.atoms:
            for site in atom.sites:
                sites.append(site)
        return len(sites)


    def get_ep_data(self):
        data = {}
        for a in self.atoms:
            if a.num_sites > 1:
                data[a.element] = a.get_ep_data()
        return data
    
    def add_atoms(self, atoms):
        extra_points = []
        for i, atom in enumerate(atoms):
            index = atom['index']
            element = atom['element']
            crds = atom['coordinates']
            multipoles = atom['multipoles']
            if element.startswith('EP'):
                charge = multipoles['charge']
                s = FFSite(element='EP', coordinates=crds, index=index, charge=charge)
                extra_points.append(s)
            else: # regular atom
                atom = HybridAtom(index=index, element=element, coordinates=crds, 
                        multipoles=multipoles, representation=self.representation, sym=self.sym)
                self.add_atom(atom)
        # find atom the extra point is connected to
        for s in extra_points:
            closest = sorted(self.atoms, key=lambda a: a.distance_to(s))
            atom = closest[0]
            atom.add_site(s)
            s.set_attachment(atom)
            logger.debug('Site %s is appended to sites of atom %s' % (s.name, atom.name))

    def set_hybridizations(self):
        # set extra points
        for frame in self.frames:
            h, d, a = self.hybrid_atoms[frame.center.element]
            index = self.get_max_index()
            ep = (index+1, h, d, a)
            frame.center.set_hybridization(ep)

    def set_frames_from_sites(self):
        for a in self.atoms:
            a.set_frame_from_sites()

    def __repr__(self):
        o = '%s\n' % self.name
        for s in self.sites:
            o += '%s %s xyz = %.3f %.3f %.3f charge = %s r0 = %s epsilon = %s\n' % (s.element, s.name, s.x, s.y, s.z, s.charge, s.r0, s.epsilon)
        return o

