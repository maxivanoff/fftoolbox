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

class Site(object):

    def __init__(self, coordinates=None, name=None, element=None, index=None, charge=None):
        self.coordinates = np.float64(coordinates)
        self.spherical = self.cart_to_sphere()
        self.name = name
        self.element = element
        self.index = index
        self.charge = charge
        self.r0 = None
        self.epsilon = None
        if not self.name == 'H+':
            logger.debug("Creating Site instance: %s" % self)
    
    def cart_to_sphere(self):
        spherical = np.zeros(3)
        xy2 = self.x**2 + self.y**2                     # x2 + y2
        spherical[0] = np.sqrt(xy2 + self.z**2)         # r2 = x2 + y2 + z2
        spherical[1] = np.arctan2(self.y, self.x)       # theta = arctan(y/x)
        spherical[2] = np.arctan2(np.sqrt(xy2), self.z) # phi = arctan(xy/z)
        return spherical
     
    def set_coordinates(self, new_crds):
        self.coordinates = np.copy(new_crds)
        self.spherical = self.cart_to_sphere()
        #logger.debug("%s site's coordinates are changed to %r" % (self.name, new_crds))

    @property
    def x(self):
        return self.coordinates[0]
    
    @property
    def y(self):
        return self.coordinates[1]
    
    @property
    def z(self):
        return self.coordinates[2]

    @property
    def r(self):
        return self.spherical[0]

    @property
    def theta(self):
        return self.spherical[1]

    @property
    def phi(self):
        return self.spherical[2]
    
    
    def distance_to(self, site):
        d = self.coordinates - site.coordinates
        return np.linalg.norm(d)

    def __eq__(self, a):
        for i in xrange(3):
            if not self.coordinates[i] == a.coordinates[i]:
                return False
        return True
    
    def __repr__(self):
        return '%s %s %f %f %f %s' % (self.name, self.index, self.coordinates[0], self.coordinates[1], self.coordinates[2], self.charge)


class Atom(Multipole):

    # vdW radii are taken from wikipedia:
    vdw_radius_angst = {'O':1.52, 'N':1.55, 'S':1.8, 'C':1.7, 'H':1.2, 'Na': 2.27, 'F':1.47, 'Cl':1.88, 'CL':1.88, 'K': 2.75}
    vdw_radius_bohr = dict((name, radius/0.52917721092) for name, radius in vdw_radius_angst.iteritems())
    
    def __init__(self, name=None, element=None, coordinates=None, index=None, multipoles=None, representation=None):
        charge = multipoles['charge']
        if name is None: name = element
        Multipole.__init__(self, name=name, origin=coordinates, representation=representation)
        self.center = Site(coordinates=coordinates, name=name, element=element, index=index, charge=charge)
        self._extra = []
        self.index = index
        self.element = element
        self.neighbors = None
        self.frame = None

    def set_name(self, name):
        self.name = name
        self.center.name = name

    @property
    def sites(self):
        sites = [self.center] + self._extra
        return iter(sites)

    @property
    def extra_sites(self):
        return iter(self._extra)

    @property
    def num_sites(self):
        return len(self._extra) + 1 

    def add_extra_site(self, site):
        self._extra.append(site)

    def free_extra(self):
        self._extra = []

    def set_hybridization(self, ep_property):
        index, hybridization, distance, angle = ep_property
        if not self.frame: self.frame = Frame(self, hybridization)
        self.free_extra()
        ep_coordinates = self.frame.ep_crds(distance, angle)
        for i, crds in enumerate(ep_coordinates):
            s = Site(coordinates=crds, name = 'EP_%s' % self.element, index=index+i)
            self.add_extra_site(s)

    def set_frame_from_sites(self):
        if self.frame is None and len(self._extra) > 0:
            num_domains = len(self.sites) + len(self.neighbors) - 1
            h = 'sp%i' % (num_domains - 1)
            self.frame = Frame(self, h)

    def get_ep_data(self):
        data = {}
        if self.num_sites == 3:
            a, b = self._extra
            v1 = a.coordinates - self.center.coordinates
            v2 = b.coordinates - self.center.coordinates
            data['distance'] = norm(v1)*units.au_to_angst
            v1 = v1/norm(v1)
            v2 = v2/norm(v2)
            cosa = np.dot(v2, v1)
            data['angle'] = np.arccos(cosa)/np.pi*180./2.
        if self.num_sites == 2:
            # from EP to atom
            a = self._extra[0]
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
    
    def bonded_to(self, atom):
        d = self.center.distance_to(atom)
        radius = atom.radius + self.radius
        if d <= radius*0.65 and d > 1e-4:
            return True
        else: return False

    @property
    def coordinates(self):
        return self.center.coordinates

    @property
    def x(self):
        return self.center.x
    
    @property
    def y(self):
        return self.center.y
    
    @property
    def z(self):
        return self.center.z
    
    @property
    def radius(self):
        return self.vdw_radius_bohr[self.element]
    
    def __repr__(self):
        return self.center.__repr__()

