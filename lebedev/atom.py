import logging
import numpy as np

import lebedev_write
from fftoolbox.atom import Atom, FFSite
from scipy.special import sph_harm as Y
from fftoolbox.multipole import Ylmc, Ylms

logger = logging.getLogger(__name__)

class LebedevAtom(Atom):

    def __init__(self, element=None, coordinates=None, index=None, ref_multipoles=None, representation=None):
        Atom.__init__(self, index=index, element=element, coordinates=coordinates, representation=representation)
        self.ref_multipoles = ref_multipoles
        self.rank = ref_multipoles['rank']
        logger.info('LebedevAtom instance is created for %s' % self.name)

    def set_sphere(self, rank=1, radius=1.0):
        if rank is None:
            rank = self.rank
        sites = lebedev_write.Lebedev(rank_to_num[rank])
        for i, site in enumerate(sites):
            # coordinates and weights
            xyz = np.array(site[0:3])*radius
            w = site[3]*4*np.pi
            # create site and compute charge
            s = FFSite(name='EP_%s-%i' % (self.name, i), element='EP', coordinates=xyz)
            q = w*self.compute_charge(rank, s.r, s.theta, s.phi)
            s.set_charge(q)
            # shift site relative to the atom center
            shifted = self.coordinates + s.coordinates
            s.set_coordinates(shifted)
            self.add_site(s)
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


