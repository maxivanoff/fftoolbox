import numpy as np
import logging
import os

from atom import Atom, Site
from bonds import Bond
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
        self.set_frames_from_sites() # in case extra points are loaded from the file
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
        extra_points = []
        for i, a in enumerate(atoms):
            index, element, crds, charge = a
            if element[0] == 'X' or element[0] == 'E': 
                s = Site(coordinates=crds, name='EP', index=index, charge=charge)
                extra_points.append(s)
            else: # regular atom
                atom = Atom(name=element, element=element, coordinates=crds, index=index, charge=charge)
                self.atoms.append(atom)
        # find atom the extra point is connected to
        for s in extra_points:
            closest = sorted(self.atoms, key=lambda a: a.distance_to(s))
            atom = closest[0]
            atom.sites.append(s)
            s.name = 'EP_%s' % atom.element
            logger.debug('Site %s is appended to sites of atom %s' % (s.name, atom.name))

    @property
    def bonds(self):
        bonds = []
        for atom in self.atoms:
            for nghbr in atom.neighbors:
                bond = Bond(atom, nghbr)
                if not bond in bonds:
                    bonds.append(bond)
            for s in atom.sites[1:]:
                bond = Bond(atom, s)
                if not bond in bonds:
                    bonds.append(bond)
        return bonds

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
            if not len(filter(lambda a: a.element=='O', atom.neighbors))==0 and atom.element=='C' and len(atom.neighbors)==3:
                carbons_with_oxygen.append(atom)
            logger.debug("Neighbors are assigned to %s:\n%s" % (atom.name, atom.neighbors))

        # set extra points
        for atom in self.atoms:
            if atom.element in self.hybrids:
                h, d, a = self.hybrids[atom.element]
                index = self.get_max_index()
                ep = (index+1, h, d, a)
                atom.set_hybridization(ep)
        
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

    def set_frames_from_sites(self):
        for a in self.atoms:
            a.set_frame_from_sites()

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

    def __init__(self, complex_data=None, molecules_data=None):
        if complex_data is None: complex_data = {'name': None}
        if molecules_data is None: molecules_data = []
        GroupOfAtoms.__init__(self, name=complex_data['name'])
        self.qm_energy = None
        self.molecules = []
        # load molecules
        num_extra = 0
        num_atoms = 0
        for mol_data in molecules_data:
            natoms = mol_data['num atoms']
            mol_name = mol_data['name']
            sym = mol_data['symmetry']

            atoms = []
            try:
                atoms_order = mol_data['atoms order']
            except KeyError:
                atoms_order = [a+1+num_atoms for a in range(natoms)]
            for i in atoms_order:
                index, elem, crds, q = complex_data['atoms'][i-1] # indices start with 1
                index += num_extra
                atom = (index, elem, crds, q)
                atoms.append(atom)

            d = {
                    'name': mol_name,
                    'atoms': atoms,
                    'symmetry': sym,
                }

            m = Molecule(data=d)
            self.add_molecule(m)
            num_extra += len(m.extra_sites)
            num_atoms += len(m.atoms)

        self.decomposition = {
                'electrostatic': None,
                'dispersion': None,
                'exchange': None,
                }

    def load_forcefields(self, ffnames):
        num_extra = 0
        for m, ffname in zip(self.molecules, ffnames):
            m.reindexate(shift=num_extra)
            ff = ForceFieldXML()
            ff.load_forcefields(filename=ffname, molecule=m, here=True)
            num_extra += len(m.extra_sites)

    def reindexate(self):
        n = 0
        for m in self.molecules:
            m.reindexate(shift=n)
            n += len(m.sites)

    def ff_energy(self):
        m1, m2 = self.molecules
        elec, disp, rep = [0.]*3
        for s1 in m1.sites:
            for s2 in m2.sites:
                # electrostatic
                r = s1.distance_to(s2)
                elec += s1.charge*s2.charge/r*units.au_to_kcal
                # Lennard-Jones
                eps = np.sqrt(s1.epsilon*s2.epsilon)
                r0 = (s1.r0 + s2.r0)*units.angst_to_au
                r6 = pow(r0/r, 6)
                r12 = pow(r6, 2)
                disp += -eps*2*r6
                rep += eps*r12
        e = elec + disp + rep
        self.decomposition = {
                'electrostatic': elec,
                'dispersion': disp,
                'repulsion': rep,
                'total': e,
                'vdw': disp+rep,
                }
        return e

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    @property
    def bonds(self):
        bonds = []
        for m in self.molecules:
            bonds += m.bonds
        return bonds

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

