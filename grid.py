import numpy as np
import os
import logging
import parser
from time import time

import units

logger = logging.getLogger(__name__)

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

class GridPoint(object):
    """
    A grid point that contains coordinates, properties at the given coordinate

    self.labels are needed for the grid manipulation
    self.signature is a string with coordinates, used in the creating of potential energy surface
    """
    def __init__(self, coordinates=None, value=None):
        self.value = value
        self.coordinates = np.float64(coordinates)
        self.labels = set()
        for i in range(3):
            if abs(self.coordinates[i]) < 1e-4: self.coordinates[i] = 0.
            plane = ortho_plane(i)
            if self.coordinates[i] < 0.: self.labels.add('<' + plane)
            if self.coordinates[i] > 0.: self.labels.add('>'+ plane)
            if self.coordinates[i] == 0.: self.labels.add(plane)
        self.spherical = self.cartesian_to_spherical()
    
    def cartesian_to_spherical(self, origin=None):
        if origin is None:
            origin = np.zeros(3)
        x, y, z = self.coordinates - origin
        spherical = np.zeros(3)
        xy2 = x**2 + y**2                     # x2 + y2
        spherical[0] = np.sqrt(xy2 + z**2)         # r2 = x2 + y2 + z2
        spherical[1] = np.arctan2(y, x)       # theta = arctan(y/x)
        spherical[2] = np.arctan2(np.sqrt(xy2), z) # phi = arctan(xy/z)
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
    
    def __init__(self, data=None):
        self.clear_points()
        if data is None: data = {}
        try:
            self.molecule_name = data['name']
        except:
            self.molecule_name = None
        try:
            self.theory = data['theory']
        except KeyError: 
            self.theory = None
        try:
            self.atoms = data['atoms']
        except KeyError:
            self.atoms = []
        # points to include/exclude
        try:
            self.exclude = data['exclude']
        except KeyError:
            self.exclude = []
        try:
            self.include = data['include']
        except KeyError:
            self.include = []
        # cubic grid parameters
        try:
            self.origin = data['origin']
            self.num_cubic_points = data['num_points']
            self.vectors = data['vectors']
            logger.info("Grid instance is to be created:\nNumber of points: %s\nvectors:\n%s\nOrigin: %s"\
                    % (self.num_cubic_points, self.vectors, self.origin))
        except KeyError:
            self.origin = None
            self.num_cubic_points = None
            self.vectors = None

    def clear_points(self):
        self._points = []

    def add_point(self, point):
        self._points.append(point)

    @property
    def points(self):
        return self._points

    @property
    def num_points(self):
        return len(self.points)

    def create_from_list(self, coordinates=None, values=None):
        self.clear_points()
        nc = len(coordinates)
        nv = len(values)
        if not nc == nv:
            raise ValueError('Number of grid points %i does not match number of values %i' % (nc, nv))
        for xyz, value in zip(coordinates, values):
            point = GridPoint(coordinates=xyz, value=value)
            self.add_point(point)
        logger.info("grid is created: %i points" % len(self.points))

    def create_grid(self, values=None):
        self.clear_points()
        if values is None:
            values = [0.]*np.prod(self.num_cubic_points)
        count = 0
        for i in xrange(self.num_cubic_points[0]):
            for j in xrange(self.num_cubic_points[1]):
                for k in xrange(self.num_cubic_points[2]):
                    xyz = self.origin + i*self.vectors[0] + j*self.vectors[1] + k*self.vectors[2]
                    point = GridPoint(coordinates=xyz, value=values[count])
                    self.add_point(point)
                    count += 1
        logger.info("grid is created: %i points" % len(self.points))

    def contains_proper_labels(self, point):
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
    vdw_radius_angst = {'O':1.52, 'N':1.55, 'S':1.8, 'C':1.7, 'H':1.2, 'Na': 2.27, 'F':1.47, 'Cl':1.88, 'Br': 1.9, 'K': 2.75 }#taken from wikipedia
    vdw_radius_bohr = dict((name, radius/0.52917721092) for name, radius in vdw_radius_angst.iteritems())
    
    def __init__(self, data):
        Grid.__init__(self, data)
        self.scale = [1.66, 2.2]
        self.create_vdw_grid(values=data['values'])
        self.create_atomic_grids()

    def inside_vdw_range(self, xyz):
        inside_small = False
        inside_large = False
        for atom in self.atoms:
            atom_crds = atom['coordinates']
            element = atom['element']
            d = np.linalg.norm(atom_crds - xyz)
            low_lim = self.vdw_radius_bohr[element]*self.scale[0]
            high_lim = self.vdw_radius_bohr[element]*self.scale[1]
            if d < low_lim: inside_small = True
            if d < high_lim: inside_large = True
        if not inside_small and inside_large:
            return True
        else: 
            return False

    def create_atomic_grids(self):
        self.atomic_grids = {}
        for point in self.points:
            d_min = 100
            for atom in self.atoms:
                index = atom['index']
                at_xyz = atom['coordinates']
                element = atom['element']
                d = np.linalg.norm(at_xyz - point.coordinates)
                if d < d_min:
                    d_min = d
                    element_min = element
                    i_min = index
            try:
                atomic_grid = self.atomic_grids[i_min]
            except KeyError:
                atomic_grid = Grid()
                self.atomic_grids[i_min] = atomic_grid
            atomic_grid.element = element_min
            atomic_grid.add_point(point)
        check_sum = 0.
        for index, grid in self.atomic_grids.items():
            check_sum += grid.num_points
            logger.info("vdW around atom with index %i and element %s is created with %i points" % (index, grid.element, grid.num_points))
        if not check_sum == self.num_points:
            raise ValueError('Number of atomic grid points (%i) does not match total number of points in the grid (%i)' % (check_sum, self.num_points))
        
    def create_vdw_grid(self, values):
        self.clear_points()
        count = 0
        for i in xrange(self.num_cubic_points[0]):
            for j in xrange(self.num_cubic_points[1]):
                for k in xrange(self.num_cubic_points[2]):
                    xyz = self.origin + i*self.vectors[0] + j*self.vectors[1] + k*self.vectors[2]
                    point = GridPoint(coordinates=xyz, value=values[count])
                    if self.inside_vdw_range(xyz) and self.contains_proper_labels(point):
                        self.add_point(point)
                    count += 1
        logger.info("vdW grid is created: %i points" % len(self.points))


    
class LebedevGrid(Grid):

    def __init__(self, data):
        self.molecule_name = data['name']
        self.theory = data['theory']
        self.points = [] 
        self.radius = data['grid radius']
        self.create_grid(points=data['grid points'], values=data['values'])
        self.weights = data['weights']

    def create_grid(self, points, values):
        self.points = []
        logger.info("Creating grid from list of values.\nNumber of grid points: %i\nNumber of grid values: %i\nRadius: %.1f" % (len(points), len(values), self.radius))
        for p, v in zip(points, values):
            point = GridPoint(coordinates=p, value=v)
            self.points.append(point)

    def get_values(self):
        logger.info("Grid point values were multiplied by the quadrature weights\nNumber of weights:%i\nNumber of sites:%i" % (len(self.weights), len(self.points)))
        return np.array([p.value*w for p, w in zip(self.points, self.weights)])

