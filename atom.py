from __future__ import division
import numpy as np
from numpy.linalg import norm
import logging

from frame import Frame
import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

atom_logger = logging.getLogger('atom')

class Site(object):

    def __init__(self, coordinates=None, name=None, charge=None):
        self.coordinates = np.float64(coordinates)
        self.spherical = self.cart_to_sphere()
        self.name = name
        self.element = name
        self.charge = charge
        self.r0 = None
        self.epsilon = None
        self.sites = [self]
        self.i = -1
        if not self.name == 'H+':
            atom_logger.debug("Creating Site instance: %s" % self)
    
    def cart_to_sphere(self):
        spherical = np.zeros(3)
        xy2 = self.x**2 + self.y**2                     # x2 + y2
        spherical[0] = np.sqrt(xy2 + self.z**2)         # r2 = x2 + y2 + z2
        spherical[1] = np.arctan2(self.y, self.x)       # theta = arctan(y/x)
        spherical[2] = np.arctan2(np.sqrt(xy2), self.z) # phi = arctan(xy/z)
        return spherical
     
    def set_coordinates(self, new_crds):
        self.coordinates = np.copy(new_crds)

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
    
    def __iter__(self):
        return self
    
    def next(self):
        if self.i < len(self.sites)-1:
            self.i += 1
            return self.sites[self.i]
        else:
            self.i = -1
            raise StopIteration
    
    def distance_to(self, site):
        d = self.coordinates - site.coordinates
        return np.linalg.norm(d)
    
    def __repr__(self):
        return '%s %f %f %f %s' % (self.name, self.coordinates[0], self.coordinates[1], self.coordinates[2], self.charge)


class Atom(Site):
    vdw_radius_angst = {'O':1.52, 'N':1.55, 'S':1.8, 'C':1.7, 'H':1.2, 'Na': 2.27, 'F':1.47, 'Cl':1.88, 'CL':1.88, 'K': 2.75 , 'XX':0.0 }#taken from wikipedia
    vdw_radius_bohr = dict((name, radius/0.52917721092) for name, radius in vdw_radius_angst.iteritems())
 
    
    def __init__(self, name=None, element=None, coordinates=None, index=None, charge=None):
        if not name: name = element
        Site.__init__(self, coordinates=coordinates, name=name, charge=charge)
        self.element = element
        self.index = index
        self.neighbors = None
        self.frame = None

    def set_hybridization(self, ep_property):
        hybridization, distance, angle = ep_property
        if not self.frame: self.frame = Frame(self, hybridization)
        self.sites = (self,)
        for i, crds in enumerate(self.frame.ep_crds(distance, angle)):
            s = Site(crds, 'EP_' + self.name, None)
            self.sites += (s,)
    
    def get_ep_data(self):
        data = {}
        if len(self.sites) == 3:
            v1 = self.sites[1].coordinates - self.sites[0].coordinates
            v2 = self.sites[2].coordinates - self.sites[0].coordinates
            data['distance'] = norm(v1)*units.au_to_angst
            v1 = v1/norm(v1)
            v2 = v2/norm(v2)
            cosa = np.dot(v2, v1)
            data['angle'] = np.arccos(cosa)/np.pi*180./2.
        if len(self.sites) == 2:
            n1 = self.coordinates - self.neighbors[0].coordinates
            n2 = self.coordinates - self.neighbors[1].coordinates
            n = n1/norm(n1) + n2/norm(n2)
            v = self.sites[1].coordinates - self.sites[0].coordinates
            cosa = np.dot(v, n)/norm(v)/norm(n)
            if abs(abs(cosa) - 1.0) < 0.00001: cosa = 1.0*np.sign(cosa)
            data['angle'] = np.arccos(cosa)/np.pi*180
            data['distance'] = norm(v)*units.au_to_angst
        return data
    
    def bonded_to(self, atom):
        d = self.distance_to(atom)
        radius = atom.radius + self.radius
        if d <= radius*0.65 and d > 1e-4:
            return True
        else: return False
    
    @property
    def radius(self):
        return self.vdw_radius_bohr[self.element]

