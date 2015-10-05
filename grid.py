import numpy as np
import os
import logging
import parser

import units

grid_logger = logging.getLogger('grid')

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

WORKDIR = os.path.dirname(__file__)

def ortho_axis(plane):
    for i in [0, 1, 2]:
        if not i in plane: return i

def ortho_plane(axis):
    if axis==0: return 'yz'
    if axis==1: return 'xz'
    if axis==2: return 'xy'

class Property(object):

    def __init__(self, value=None,property_name=None, theory=None):
        self.name = property_name
        self.value = value
        self.theory = theory

    def set_property(self, value=None, property_name=None, theory=None):
        self.value = value
        self.name = property_name
        self.theory = theory

    def set_value(self, value=None):
        self.value = value

    def __str__(self):
        return '%.10f %s %s' % (self.value, self.name, self.theory)

    
class GridPoint(object):
    """
    A grid point that contains coordinates, properties at the given coordinate

    self.labels are needed for the grid manipulation
    self.signature is a string with coordinates, used in the creating of potential energy surface
    """
    def __init__(self, coordinates=None, signature=None, value=None):
        self.value = value
        self.labels = set()
        if type(coordinates)==np.ndarray and signature is None:
            self.coordinates = np.float64(coordinates)
            self.signature = '%.3f_%.3f_%.3f' % (coordinates[0], coordinates[1], coordinates[2]) 
        if type(signature)==str and coordinates is None:
            self.signature = signature
            self.coordinates = np.array([np.float64(s) for s in self.signature.split('_')])
        for i in range(3):
            if abs(self.coordinates[i])<1e-4: self.coordinates[i]=0.
            if self.coordinates[i]<0.: self.labels.add('<'+ortho_plane(i))
            if self.coordinates[i]>0.: self.labels.add('>'+ortho_plane(i))
            if self.coordinates[i]==0.: 
                self.labels.add(ortho_plane(i))
        self.spherical = self.cart_to_sphere()
    
    def cart_to_sphere(self):
        spherical = np.zeros(3)
        xy2 = self.x**2 + self.y**2                     # x2 + y2
        spherical[0] = np.sqrt(xy2 + self.z**2)         # r2 = x2 + y2 + z2
        spherical[1] = np.arctan2(self.y, self.x)       # theta = arctan(y/x)
        spherical[2] = np.arctan2(np.sqrt(xy2), self.z) # phi = arctan(xy/z)
        return spherical
 

    def add_label(self, label): 
        self.labels.add(label)

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

class Grid(object):
    
    def __init__(self, data):
        
        self.molecule_name = data['name']
        try:
            self.theory = data['theory']
        except KeyError: 
            self.theory = None
        try:
            self.exclude = data['exclude']
        except KeyError:
            self.exclude = []
        try:
            self.include = data['include']
        except KeyError:
            self.include = []
        try:
            self.origin = data['origin']
            self.num_points = data['num_points']
            self.vectors = data['vectors']
            grid_logger.info("Grid instance is to be created:\nNumber of points: %s\nvectors:\n%s\nOrigin: %s"\
                    % (self.num_points, self.vectors, self.origin))
        except KeyError:
            self.origin = None
            self.num_points = None
            self.vectors = None
        try:
            self.atoms = data['atoms']
        except KeyError:
            self.atoms = []

        self.points = []

    def create_from_list(self, coordinates=None, values=None):
        self.points = []
        nc = len(coordinates)
        nv = len(values)
        if not nc == nv:
            raise ValueError('Number of grid points %i does not match number of values %i' % (nc, nv))
        for xyz, value in zip(coordinates, values):
            point = GridPoint(coordinates=xyz, value=value)
            self.points.append(point)
        grid_logger.info("grid is created: %i points" % len(self.points))

    def create_grid(self, values):
        self.points = []
        count = 0
        for i in xrange(self.num_points[0]):
            for j in xrange(self.num_points[1]):
                for k in xrange(self.num_points[2]):
                    xyz = self.origin + i*self.vectors[0] + j*self.vectors[1] + k*self.vectors[2]
                    point = GridPoint(coordinates=xyz, value=values[count])
                    self.points.append(point)
                    count += 1
        grid_logger.info("grid is created: %i points" % len(self.points))

    def proper(self, point):
        for label in self.exclude:
            if label in point.labels: 
                return False
        for label in self.include:
            if not label in point.labels:
                return False
        return True

    def get_coordinates(self, gridformat=False, shift = None):
        coords = np.array([p.coordinates for p in self.points])
        if not gridformat:
            return coords
        else:
            x = np.copy(coords[:,0])
            y = np.copy(coords[:,1])
            z = np.copy(coords[:,2])
            if shift is not None:
                x += -shift[0]
                y += -shift[1]
                z += -shift[2]
            num_points = np.sqrt(len(self.points))
            x.resize(num_points, num_points)
            y.resize(num_points, num_points)
            z.resize(num_points, num_points)
            return x, y, z
    
    def get_values(self):
        return np.array([p.value for p in self.points])
    
    def build_LEM(self, filename, path2xyz=False,vmax=0.001, r_sphere=0.2):
        if path2xyz is False:
            path2xyz = '%s/data/xyz/%s_%s.xyz' % (WORKDIR, self.molecule_name, self.theory)
        s = 'from pymol.cgo import *\nfrom pymol import cmd\ncmd.load("%s")\nobj = [ BEGIN, LINES, ]\n' % (path2xyz)
        for p in self.points:
            crds = p.coordinates*units.au_to_angst
            if p.value is None: 
                s_color = 'x = 0.0\ncolor = [COLOR, 1-x, 1-x, 1]\n'
            elif p.value >= 0:
                s_color = 'x = %f\ncolor = [COLOR, 1, 1-x, 1-x]\n' % (p.value/vmax)
            elif p.value < 0:
                s_color = 'x = %f\ncolor = [COLOR, 1-x, 1-x, 1]\n' % (-p.value/vmax)
            s_sphere = 'sphere = [ SPHERE, %f, %f, %f,%f]\n' % (crds[0], crds[1], crds[2], r_sphere)
            s = s + s_color + s_sphere + 'obj += color+sphere\n'
        s = s + 'obj.append(END)\ncmd.load_cgo(obj,"cgo01")\n'

        file = open(filename,'w')
        file.write(s)
        file.close()

    def __str__(self):
        s=''
        for p in self.points:
            s += '%.3f %.3f %.3f' % (p.coordinates[0], p.coordinates[1], p.coordinates[2])
            for l in p.labels:
                s += ' %s' % (l)
            s += '\n'
        return s
 
class vdwGrid(Grid):
    vdw_radius_angst = {'O':1.52, 'N':1.55, 'S':1.8, 'C':1.7, 'H':1.2, 'Na': 2.27, 'F':1.47, 'Cl':1.88, 'K': 2.75 }#taken from wikipedia
    vdw_radius_bohr = dict((name, radius/0.52917721092) for name, radius in vdw_radius_angst.iteritems())
    
    def __init__(self, data):
        Grid.__init__(self, data)

        self.scale = [1.4, 2.]
        self.create_vdw_grid(values=data['values'])        

    def set_vdw_labels(self, atoms=None):
        if atoms is None or len(atoms)==0:
            return
        if self.vdw_range is not None:
            grid_logger.info('Setting vdW labels: vdW ranges are user-defined')
        else:
            grid_logger.info('Setting vdW labels: vdW ranges are default, scale=%s' % self.scale)
        grid_logger.debug("vdW surface is built based on %i atoms:\n%r" % (len(atoms), atoms))
        for p in self.points:
            if inside_range(atoms, p.coordinates): p.add_label('inside')
            else: p.add_label('outside')

    def inside_range(self, xyz):
        inside_small = False
        inside_large = False
        for atom_name, atom_crds, mass in self.atoms:
            d = np.linalg.norm(atom_crds - xyz)
            low_lim = self.vdw_radius_bohr[atom_name]*self.scale[0]
            high_lim = self.vdw_radius_bohr[atom_name]*self.scale[1]
            if d < low_lim: inside_small = True
            if d < high_lim: inside_large = True
        if not inside_small and inside_large: return True
        else: return False


    def create_vdw_grid(self, values):
        self.points = []
        count = 0
        for i in xrange(self.num_points[0]):
            for j in xrange(self.num_points[1]):
                for k in xrange(self.num_points[2]):
                    xyz = self.origin + i*self.vectors[0] + j*self.vectors[1] + k*self.vectors[2]
                    point = GridPoint(coordinates=xyz, value=values[count])
                    if self.inside_range(xyz) and self.proper(point):
                        self.points.append(point)
                    count += 1
        grid_logger.info("vdW grid is created: %i points" % len(self.points))


    

