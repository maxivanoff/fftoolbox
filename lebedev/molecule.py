from atom import LebedevAtom
from scipy.special import sph_harm as Y
from fftoolbox.multipole import Ylmc, Ylms, Multipole
from fftoolbox.molecule import Molecule
from fftoolbox.atom import Site
import lebedev_write
import numpy as np
import logging
logger = logging.getLogger(__name__)

class LebedevMolecule(Multipole):

    def __init__(self, data):
        Multipole.__init__(self, data['name'], representation=data['representation'])
        self.num_sites = data['points']
        self.radius = data['radius']
        self.set_sites()
        self.set_sym_sites()
        self.set_multipole_matrix()
        logger.info("Lebedev molecule is created.\nNumber of charged sites: %i\nRadius: %.1f" % (len(self.sites), self.radius))

    def set_sites(self):
        sites = lebedev_write.Lebedev(self.num_sites)
        self.atoms = []
        self.weights = np.zeros(len(sites))
        for i, site in enumerate(sites):
            xyz = np.array(site[0:3])*self.radius
            w = site[3]*np.sqrt(4*np.pi)
            self.weights[i] = w
            s = Site(name='%i' % i, index=i+1, coordinates=xyz)
            self.atoms.append(s)

    def get_basis(self, nodes=None, l=None):
        N = len(nodes)
        B = np.zeros((l*l, N))
        i = 0
        for ll in range(l):
            for mm in range(-ll,ll+1):
                if mm < 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Ylms(abs(mm),ll,  p.theta, p.phi).real
                if mm == 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Y(0, ll, p.theta, p.phi).real
                if mm > 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Ylmc(mm, ll, p.theta, p.phi).real
                i += 1
        return B

class DistributedLebedevMolecule(Molecule, Multipole):

    def __init__(self, data):
        try:
            self.rep = data['representation']
        except KeyError:
            self.rep = ('spherical', 1)
        try:
            self.sym = data['symmetry']
        except KeyError:
            self.sym = False
        self.distrib = data['distributions']
        Multipole.__init__(self, data['name'], representation=self.rep)
        # set atoms
        self.set_atoms(data['atoms'])
        self.set_groups()
        # distribute point charges
        self.set_spheres()
        self.set_sym_sites()
        self.set_multipole_matrix()
        logger.info('Distributed Lebedev %s Molecule is created' % (self.name))

    def set_spheres(self):
        for a in self.atoms:
            rank, radius = self.distrib[a.element]
            a.set_sphere(rank, radius)

    def set_atoms(self, atoms):
        self.atoms = []
        for a in atoms:
            index, elem, crds, multipoles = a
            atom = LebedevAtom(name=elem, element=elem, coordinates=crds, ref_multipoles=multipoles, index=index, representation=self.rep)
            self.atoms.append(atom)

