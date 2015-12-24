import numpy as np
import logging
import os

import units
from multipole import GroupOfSites
from molecule import AtomsInMolecule

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)
WORKDIR = os.path.dirname(__file__)

def m_to_xyz(m):
    if m == -1:
        return 1 # y
    if m == 0:
        return 2 # z
    if m == 1:
        return 0 # x

def T(l1, l2, m1, m2, r1, r2, c):
    if l1 == l2 == 0:
        return 1.
    if l1 == 0 and l2 == 1:
        alpha = m_to_xyz(m2)
        return r2[alpha]
    if l1 == 1 and l2 == 0:
        alpha = m_to_xyz(m1)
        return r1[alpha]
    if l1 ==1 and l2 == 1:
        alpha = m_to_xyz(m1)
        beta = m_to_xyz(m2)
        return 3*r1[alpha]*r2[beta] + c[alpha, beta]

class Complex(GroupOfSites, AtomsInMolecule):

    def __init__(self, name, mol1, mol2):
        GroupOfSites.__init__(self, name)
        AtomsInMolecule.__init__(self, name)
        self.mol1 = mol1
        self.mol2 = mol2
        self.molecules = (mol1, mol2)
        i = 1
        for m in self.molecules:
            for site in m.sites:
                self.add_site(site)
                site.set_index(i)
                i += 1
            for atom in m.atoms:
                self.add_atom(atom)

    def point_charge_energy(self):
        e = 0.
        for s1 in self.mol1.sites:
            for s2 in self.mol2.sites:
                r = s1.distance_to(s2)
                e += s1.charge*s2.charge/r
        return e

    @property
    def bonds(self):
        bonds = []
        for m in self.molecules:
            bonds += m.bonds
        return bonds

class TwoLebedevMolecules(Complex):

    def __init__(self, name, mol1, mol2):
        Complex.__init__(self, name, mol1, mol2)

    def multipole_energy(self, rank):
        # pre-compute stuff for interaction function T
        R_vec = self.mol1.center_of_mass() - self.mol2.center_of_mass()
        R = np.linalg.norm(R_vec)
        e12 = R_vec/R
        r1 = []
        r2 = []
        for e in self.mol1.frame.local_axes.T:
            r1.append(np.dot(e, e12))
        for e in self.mol2.frame.local_axes.T:
            r2.append(-np.dot(e, e12))
        c = np.zeros((3,3))
        for i, e1 in enumerate(self.mol1.frame.local_axes.T):
            for j, e2 in enumerate(self.mol2.frame.local_axes.T):
                c[i,j] = np.dot(e1, e2)
        # compute energy
        U = 0.
        for l1 in xrange(rank+1):
            for m1 in xrange(-l1, l1+1):
                for l2 in xrange(rank+1):
                    for m2 in xrange(-l2, l2+1):
                        Qlm_1 = self.mol1.get_reference_Qlm(l1, m1)
                        Qlm_2 = self.mol2.get_reference_Qlm(l2, m2)
                        U += Qlm_1*Qlm_2*pow(R, -l1 -l2 -1)*T(l1,l2,m1,m2,r1,r2,c)
        return U

        
class OldComplex(GroupOfSites):

    def __init__(self, complex_data=None, molecules_data=None):
        if complex_data is None: complex_data = {'name': None}
        if molecules_data is None: molecules_data = []
        GroupOfSites.__init__(self, name=complex_data['name'])
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


