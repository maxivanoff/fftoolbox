from fftoolbox.multipole import Ylmc, Ylms, Multipole
from fftoolbox.atom import FFSite
from fftoolbox.molecule import Molecule
from fftoolbox.atom import Atom
import lebedev_write

from scipy.special import sph_harm as Y
import numpy as np
import logging

logger = logging.getLogger(__name__)

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"


class LebedevSphere(Multipole):

    rank_to_num = {
                    1:6,
                    2:14,
                    3:26,
                    4:38,
                    5:50,
                    6:74,
                    }

    def __init__(self, index=None, name=None, rank=None, radius=None, origin=None, ref_multipoles=None):
        self.reference_multipoles = ref_multipoles 
        if rank > self.reference_multipoles['rank']:
            raise ValueError('Multipoles only up to rank %i are available' % self.reference_multipoles['rank'])
        representation = ('spherical', rank)
        self.rank = rank
        Multipole.__init__(self, name, origin=origin, representation=representation)
        if rank == 0:
            charge = self.reference_multipoles['00']
            s = FFSite(index=index, name=self.name, coordinates=self.origin, charge=charge, attachment=self)
            self.add_site(s)
        else:
            points = lebedev_write.Lebedev(self.rank_to_num[rank])
            for i, point in enumerate(points):
                # coordinates and weights
                xyz = np.array(point[0:3])*radius
                w = point[3]*4*np.pi
                # create site and compute charge
                name = '%s-%i' % (self.name, i)
                s = FFSite(index=index, element='EP', coordinates=xyz, attachment=self)
                q = w*self.compute_charge(rank, s.r, s.theta, s.phi)
                s.set_charge(q)
                # shift site relative to the atom center
                shifted = self.origin + s.coordinates
                s.set_coordinates(shifted)
                self.add_site(s)
        self.set_sym_sites()
        self.set_multipole_matrix()
        logger.info("LebedevSphere %s is created.\nNumber of charged sites: %i\nRadius: %.1f" 
                % (self.name, self.num_sites, radius))

    def compute_charge(self, n, r, theta, phi):
        q = 0.
        for l in xrange(n+1):
            k = np.sqrt( (2. * l + 1.) / 4. / np.pi) / r**l
            for m in xrange(-l,l+1):
                if m < 0:
                    mname = '%i%is' % (l, abs(m))
                    YY = Ylms
                if m == 0:
                    mname = '%i0' % l
                    YY = Y
                if m > 0:
                    mname = '%i%ic' % (l, m)
                    YY = Ylmc
                try:
                    Q = self.reference_multipoles[mname]
                except KeyError:
                    Q = 0.
                q += k*Q*YY(abs(m), l, theta, phi).real
        return q

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

class LebedevAtom(Atom, LebedevSphere):

    def __init__(self, element=None, coordinates=None, index=None, ref_multipoles=None, rank=None, radius=None):
        Atom.__init__(self, index=index, element=element, coordinates=coordinates, representation=None)
        LebedevSphere.__init__(self, name=self.name, rank=rank, radius=radius, \
                origin=coordinates, ref_multipoles=ref_multipoles)

class LebedevMolecule(LebedevSphere):

    def __init__(self, data):
        self.data = data
        rank, radius = data['sphere params']
        self.rank = rank
        representation = ('spherical', rank)
        self.molecule = Molecule(data)
        LebedevSphere.__init__(self, name=data['name'], rank=rank, radius=radius, \
                origin=None, ref_multipoles=data['multipoles'])
    
    def multipoles_data(self):
        data = {
                self.name: (self.origin, self.rank, self.reference_multipoles), 
                }
        return data


class DistributedLebedevMolecule(Molecule):

    def __init__(self, data):
        self.sphere_params = data['sphere params']
        Molecule.__init__(self, data)

    def multipoles_data(self):
        data = {}
        for a in self.atoms:
            d = (a.origin, a.rank, a.reference_multipoles)
            data[a.name] = d
        return data

    def add_atoms(self, atoms):
        for a in atoms:
            index, element, crds, multipoles = a
            try:
                rank, radius = self.sphere_params[element]
            except:
                rank, radius = self.sphere_params
            atom = LebedevAtom(element=element, coordinates=crds, ref_multipoles=multipoles, \
                    index=index, rank=rank, radius=radius)
            self.add_atom(atom)
            

