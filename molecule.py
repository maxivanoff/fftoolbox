import numpy as np
import logging
import os

from atom import Atom, Site
from groups import BuriedGroup, CarbonOxygenGroup
from multipole import Multipole
from multipole import GroupOfAtoms
from parser import ForceFieldXML
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger('molecule')
WORKDIR = os.path.dirname(__file__)

class Molecule(Multipole):
    
    def __init__(self, data):
        self.data = data
        try:
            self.theory = data['theory']
        except KeyError:
            self.theory = None
        try:
            self.mults = data['representation']
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
        try:
            self.energy = data['energy']
        except KeyError:
            self.energy = None
        logger.info('Start assembling %s Molecule\nsymmetry: %r\nmultipoles: %r\nhybridizations: %r' % (data['name'], self.sym, self.mults, self.hybrids))
        Multipole.__init__(self, name=data['name'], multipoles=self.mults)
        self.set_atoms(data['atoms'])
        self.set_groups()
        self.set_frames() # in case extra points are loaded from the file
        self.set_sym_sites()
        self.set_multipole_matrix()
        logger.info('%s Molecule is created' % (self.name))

    def copy(self):
        data = self.data.copy()
        data['atoms'] = []
        for i, a in enumerate(self.atoms):
            element, crds, charge = self.data['atoms'][i]
            a_updated = (element, a.coordinates, charge)
            self.data['atoms'].append(a_updated)
        molecule = Molecule(data=data)
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
                logger.debug('Site %s is appended to sites of atom %s' % (s.name, aton.name))
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
        carbons_with_oxygen = []
        for atom in self:
            atom.neighbors = filter(atom.bonded_to, self.atoms)
            if len(filter(lambda a: a.element=='H', atom.neighbors))>1:
                buried_atoms.append(atom)
            if not len(filter(lambda a: a.element=='O', atom.neighbors))==0 and atom.element=='C':
                carbons_with_oxygen.append(atom)
            logger.debug("Neighbors are assigned to %s:\n%s" % (atom.name, atom.neighbors))

        # set extra points
        for atom in self.atoms:
            if atom.element in self.hybrids:
                atom.set_hybridization(self.hybrids[atom.element])
        
        # create groups
        exclude = ['methane', 'benzene', 'methane', 'tip5p', 'tip3p', 'ammonia', 'water']
        self.groups = []
        if not self.name in exclude:
            for count, atom in enumerate(carbons_with_oxygen):
                g = CarbonOxygenGroup(center=atom, count=count, sym=self.sym)
                self.groups.append(g)
            for count, atom in enumerate(buried_atoms):
                g = BuriedGroup(center=atom, count=count, sym=self.sym)
                self.groups.append(g)

        
        # symmetrical atoms have identical names
         
        try:
            template, symmetrical = self.sym[0], self.sym[1:] 
            t = self.get_atom_by_index(template)
            for i in symmetrical:
               a = self.get_atom_by_index(i)
               a.name = t.name
        except TypeError:
            pass

        logger.info("Names of equivalent sites: %s" % self.sites_names_eq)
        logger.info("Names of non-equivalent sites: %s" % self.sites_names_noneq)

    def set_frames(self):
        for a in self.atoms:
            a.set_frame()

    def __add__(self, molecule):
        c = Complex()
        c.add_molecule(self.copy())
        c.add_molecule(molecule.copy())
        return c

    def __repr__(self):
        o = '%s\n' % self.name
        for s in self.sites:
            o += '%s %s xyz = %.3f %.3f %.3f charge=%s r0=%s epsilon=%s\n' % (s.element, s.name, s.x, s.y, s.z, s.charge, s.r0, s.epsilon)
        return o

class Complex(GroupOfAtoms):

    def __init__(self, data=None, molecules_data=None):
        GroupOfAtoms.__init__(self)
        self.qm_energy = None
        self.molecules = []
        # load molecules
        x = 0
        for mol_data in molecules_data:
            natoms = mol_data['num atoms']
            ffname = mol_data['xml']
            mol_name = mol_data['name']
            sym = mol_data['symmetry']

            atoms = data['atoms'][x:x+natoms]
            x += natoms
            d = {
                    'name': mol_name,
                    'atoms':atoms,
                    'symmetry': sym,
                }

            m = Molecule(data=d)
            ff = ForceFieldXML()
            ff.load_forcefields(filename=ffname, molecule=m, here=True)
            self.add_molecule(m)

    def ff_energy(self):
        m1, m2 = self.molecules
        e = 0.
        for s1 in m1.sites:
            for s2 in m2.sites:
                # electrostatic
                r = s1.distance_to(s2)
                e += s1.charge*s2.charge/r*units.au_to_kcal
                # Lennard-Jones
                eps = np.sqrt(s1.epsilon*s2.epsilon)
                r0 = (s1.r0 + s2.r0)*units.angst_to_au
                r6 = pow(r0/r, 6)
                r12 = pow(r6, 2)
                e += eps*(r12 - 2*r6)
        return e

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

