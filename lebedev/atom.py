import logging
import numpy as np

import lebedev_write
from fftoolbox.atom import Site, Atom
from scipy.special import sph_harm as Y
from fftoolbox.multipole import Ylmc, Ylms, Multipole

logger = logging.getLogger(__name__)


class LebedevAtom(Atom, Multipole):

    def __init__(self, name=None, element=None, coordinates=None, index=None, ref_multipoles=None, representation=None):
        Atom.__init__(self, name=name, element=element, coordinates=coordinates, index=index, charge=0.0)
        Multipole.__init__(self, name=name, representation=representation)
        self.origin = coordinates
        self.ref_multipoles = ref_multipoles
        self.rank = ref_multipoles['rank']
        self.atoms = [self]

    def set_sphere(self, rank=1, radius=1.0):
        if rank is None:
            rank = self.rank
        sites = lebedev_write.Lebedev(rank_to_num[rank])
        self._extra_sites = []
        for i, site in enumerate(sites):
            # coordinates and weights
            xyz = np.array(site[0:3])*radius
            w = site[3]*4*np.pi
            # create site and compute charge
            s = Site(name='%s-X%i' % (self.name,i), coordinates=xyz)
            q = self.compute_charge(rank, s.r, s.theta, s.phi)
            s.charge = q*w
            # shift site relative to the atom center
            shifted = self.coordinates + s.coordinates
            s.set_coordinates(shifted)
            self._extra_sites.append(s)
        self.set_sym_sites()
        self.set_multipole_matrix()

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
                    Q = self.ref_multipoles[mname]
                except KeyError:
                    Q = 0.
                q += k*Q*YY(abs(m), l, theta, phi).real
        return q

rank_to_num ={
        1:6,
        2:14,
        3:26,
        4:38
        }


