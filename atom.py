from __future__ import division
from numpy.linalg import norm
import numpy as np
import logging

import units
from frame import Frame
from multipole import Multipole

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)

class Coordinates(object):
    def __init__(self, coordinates):
        self._coordinates = np.float64(coordinates)
        self._spherical = self.cartesian_to_spherical()

    def cartesian_to_spherical(self):
        spherical = np.zeros(3)
        xy2 = self.x**2 + self.y**2                     # x2 + y2
        spherical[0] = np.sqrt(xy2 + self.z**2)         # r2 = x2 + y2 + z2
        spherical[1] = np.arctan2(self.y, self.x)       # theta = arctan(y/x)
        spherical[2] = np.arctan2(np.sqrt(xy2), self.z) # phi = arctan(xy/z)
        return spherical
     
    def set_coordinates(self, new_crds):
        self._coordinates = np.copy(new_crds)
        self._spherical = self.cartesian_to_spherical()
        logger.debug("Coordinates were changed to %r" % (new_crds))

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def x(self):
        return self._coordinates[0]
    
    @property
    def y(self):
        return self._coordinates[1]
    
    @property
    def z(self):
        return self._coordinates[2]

    @property
    def r(self):
        return self._spherical[0]

    @property
    def theta(self):
        return self._spherical[1]

    @property
    def phi(self):
        return self._spherical[2]
    
    def distance_to(self, site):
        d = self.coordinates - site.coordinates
        return np.linalg.norm(d)

    def __eq__(self, a):
        for i in xrange(3):
            if not self.coordinates[i] == a.coordinates[i]:
                return False
        return True

class FFSite(Coordinates):

    def __init__(self, index=None, element=None, coordinates=None, charge=None, r0=None, epsilon=None, atom=None):
        Coordinates.__init__(self, coordinates)
        self.index = index
        self.element = element
        self.charge = charge
        self.r0 = r0
        self.epsilon = epsilon
        self.set_atom(atom)
        logger.debug("Creating FFSite instance:\n%s" % self)

    @property
    def name(self):
        if self.element == 'EP':
            return 'EP_%s' % self._atom.element
        else:
            return self._atom.name

    def set_atom(self, atom):
        self._atom = atom
    
    def __repr__(self):
        return '%s %s %f %f %f %s' % (self.name, self.index, self.coordinates[0], self.coordinates[1], self.coordinates[2], self.charge)

class Atom(Coordinates, Multipole):

    """
    An atom in a molecule. 
    Its index is based on the atom's order in the molecule.
    Has neighbors to which its covalently bonded based on the vdW radii.
    vdW radii are taken from wikipedia.
    """
    vdw_radius_angst = {'O':1.52, 'N':1.55, 'S':1.8, 'C':1.7, 'H':1.2, 'Na': 2.27, 'F':1.47, 'Cl':1.88, 'CL':1.88, 'K': 2.75}
    vdw_radius_bohr = dict((name, radius/0.52917721092) for name, radius in vdw_radius_angst.iteritems())

    def __init__(self, index=None, element=None, coordinates=None, representation=None):
        Multipole.__init__(self, name=element, origin=coordinates, representation=representation)
        Coordinates.__init__(self, coordinates)
        self.index = index
        self.element = element
        self.neighbors = None

    def set_neighbors(self, atoms):
        self.neighbors = []
        for atom in atoms:
            if not atom is self and self.bonded_to(atom):
                self.neighbors.append(atom)

    @property
    def radius(self):
        return self.vdw_radius_bohr[self.element]

    def bonded_to(self, atom):
        d = self.distance_to(atom)
        radius = atom.radius + self.radius
        if d <= radius*0.65 and d > 1e-4:
            return True
        else: return False

class HybridAtom(Atom):

    def __init__(self, element=None, coordinates=None, index=None, multipoles=None, representation=None):
        Atom.__init__(self, index=index, element=element, coordinates=coordinates, representation=representation)
        charge = multipoles['charge']
        center = FFSite(index=index, element=element, coordinates=coordinates, charge=charge, atom=self)
        self.add_site(center)
        self.frame = None

    @property
    def center(self):
        return self._sites[0]

    @property
    def extra_sites(self):
        return iter(self._sites[1:])

    @property
    def num_extra_sites(self):
        return len(self._sites[1:])

    def free_extra(self):
        self._sites = [self._sites[0]]

    def set_hybridization(self, ep_property):
        index, hybridization, distance, angle = ep_property
        if not self.frame: self.frame = Frame(self, hybridization)
        self.free_extra()
        ep_coordinates = self.frame.ep_crds(distance, angle)
        for i, crds in enumerate(ep_coordinates):
            s = FFSite(element='EP', coordinates=crds, index=index+i, atom=self)
            self.add_site(s)

    def set_frame_from_sites(self):
        if self.frame is None and self.num_extra_sites > 0:
            num_domains = len(self.sites) + len(self.neighbors) - 1
            h = 'sp%i' % (num_domains - 1)
            self.frame = Frame(self, h)

    def get_ep_data(self):
        data = {}
        if self.num_extra_sites == 2:
            a, b = self.extra_sites
            v1 = a.coordinates - self.center.coordinates
            v2 = b.coordinates - self.center.coordinates
            data['distance'] = norm(v1)*units.au_to_angst
            v1 = v1/norm(v1)
            v2 = v2/norm(v2)
            cosa = np.dot(v2, v1)
            data['angle'] = np.arccos(cosa)/np.pi*180./2.
        if self.num_extra_sites == 1:
            # from EP to atom
            a = self.extra_sites[0]
            v = self.a.coordinates - self.center.coordinates
            # sum of two neighbors
            n1 = self.center.coordinates - self.neighbors[0].coordinates
            n2 = self.center.coordinates - self.neighbors[1].coordinates
            n = n1/norm(n1) + n2/norm(n2)
            cosa = np.dot(v, n)/norm(v)/norm(n)
            if abs(abs(cosa) - 1.0) < 0.00001: cosa = 1.0*np.sign(cosa)
            data['angle'] = np.arccos(cosa)/np.pi*180
            data['distance'] = norm(v)*units.au_to_angst
        return data
    
