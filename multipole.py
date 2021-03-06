import numpy as np
from numpy.linalg import norm
from scipy.special import sph_harm as Y
import logging
import os

from units import au_to_angst
import parser

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

WORKDIR = os.path.dirname(__file__)
logger = logging.getLogger(__name__)

def rodrigues_rotation_matrix(rot_axis, angle):
    rot_axis = rot_axis/norm(rot_axis)
    sinA = np.sin(angle)
    cosA = np.cos(angle)
    I = np.identity(3)
    W = np.array([[0., -rot_axis[2], rot_axis[1]],\
                  [rot_axis[2], 0., -rot_axis[0]],\
                  [-rot_axis[1], rot_axis[0], 0.]])
    R = I + W*sinA + np.dot(np.dot(W, W),I*(1 - cosA))
    return R

def Rlm(l, m, r, theta, phi):
    return r**l * np.sqrt(4 * np.pi / (2 * l + 1)) * Y(m, l, theta, phi)

def Rlmc(l, m, r, theta, phi):
    return r**l * np.sqrt(4 * np.pi / (2 * l + 1)) * Ylmc(m, l, theta, phi)

def Rlms(l, m, r, theta, phi):
    return r**l * np.sqrt(4 * np.pi / (2 * l + 1)) * Ylms(m, l, theta, phi)

def Ylmc(m, l, theta, phi):
    #v = np.sqrt(0.5) * (np.conj(Y(m, l, theta, phi)) + Y(m, l, theta, phi))
    v = np.sqrt(0.5) * (Y(-m, l, theta, phi) + (-1)**m*Y(m, l, theta, phi))
    #v = np.sqrt(0.5) * ((-1)**m*Y(-m, l, theta, phi) + Y(m, l, theta, phi))
    if abs(v.imag) > 0.0001: raise ValueError("Non-zero imaginary part in Ylmc")
    return v.real

def Ylms(m, l, theta, phi):
    #v = 1j * np.sqrt(0.5) * (np.conj(Y(m, l, theta, phi)) - Y(m, l, theta, phi))
    v = 1j * np.sqrt(0.5) * (Y(-m, l, theta, phi) - (-1)**m*Y(m, l, theta, phi))
    #v = 1j * np.sqrt(0.5) * (-(-1)**m*Y(-m, l, theta, phi) + Y(m, l, theta, phi))
    if abs(v.imag) > 0.0001: raise ValueError("Non-zero imaginary part in Ylms")
    return v.real


class GroupOfSites(object):

    def __init__(self, name=None):
        self._name = name
        self._sites = []

    def free_sites(self):
        self._sites = []

    @property
    def name(self):
        return self._name

    def set_name(self, name):
        self._name = name

    @property
    def sites(self):
        return iter(self._sites)
    
    def add_site(self, site):
        self._sites.append(site)

    @property
    def num_sites(self):
        return len(self._sites)
    
    @property
    def num_sites_noneq(self):
        sites = []
        for s in self.sites:
            if not s.name in [ss.name for ss in sites]:
                sites.append(s)
        return len(sites)

    def reindexate(self, shift):
        for s in self.sites:
            s.index += shift

    def set_sym_sites(self):
        sites = {}
        sym_sites = []
        for i, name in enumerate(self.sites_names_eq):
            if not name in sites:
                sites[name] = len(sites.keys())
            sym_sites.append(sites[name])
        self.sym_sites = np.array(sym_sites)
        logger.debug("Symmetry indeces in %s:\n%r" % (self.name, self.sym_sites))
                
    def get_coordinates(self):
        crds = np.zeros((self.num_sites, 3))
        for i, s in enumerate(self.sites):
            crds[i][:] = s.coordinates[:]
        return crds

    def get_max_index(self):
        try:
            return max([s.index for s in self.sites])
        except ValueError:
            return 0
    
    def get_sites_by_name(self, name):
        return iter(filter(lambda s: s.name == name, self.sites))

    def get_sites_by_element(self, element):
        return iter(filter(lambda s: s.element == element, self.sites))
    
    @property
    def sites_names_noneq(self):
        return [s.name for s in self.sites_noneq]

    @property
    def sites_names_eq(self):
        return [s.name for s in self.sites]

    @property
    def sites_names(self):
        return [s.name for s in self.sites]
    
    @property
    def sites_noneq(self):
        sites = []
        for s in self.sites:
            if not s.name in [ss.name for ss in sites]:
                sites.append(s)
        return sites

    def write_xyz(self, filename, here=False):
        out = '%i\n' % self.num_sites
        try:
            out += 'e = %s\n' % self.energy
        except AttributeError:
            out += '\n'
        for s in self.sites:
            out += '%s %.10f %.10f %.10f\n' % (s.element, s.x*au_to_angst, s.y*au_to_angst, s.z*au_to_angst)
        if here:
            file = open('./%s' % (filename), 'w')
        else:
            file = open('%s/data/xyz/%s' % (WORKDIR, filename), 'w')
        file.write(out)
        file.close()

    def write_mol2(self, filename, here=False):
        mol2 = parser.Mol2()
        mol2.write_file(filename, here, self)

    def color_charges(self, filename, xyzname=False, here=False, vmax=None, r_sphere=0.08):
        if here is False:
            path2xyz = '%s/data/xyz/%s' % (WORKDIR, xyzname)
        else:
            path2xyz = './%s' % xyzname
        s = 'from pymol.cgo import *\nfrom pymol import cmd\ncmd.load("%s")\nobj = [ BEGIN, LINES, ]\n' % (path2xyz)
        for atom in self.atoms:
            vmax = max([abs(ss.charge) for ss in atom.sites]) + 0.3
            for site in atom.sites:
                crds = site.coordinates*au_to_angst
                if site.charge is None: 
                    s_color = 'x = 0.0\ncolor = [COLOR, 1-x, 1-x, 1]\n'
                elif site.charge <= 0:
                    s_color = 'x = %f\ncolor = [COLOR, 1, 1-x, 1-x]\n' % (-site.charge/vmax)
                elif site.charge > 0:
                    s_color = 'x = %f\ncolor = [COLOR, 1-x, 1-x, 1]\n' % (site.charge/vmax)
                s_sphere = 'sphere = [ SPHERE, %f, %f, %f,%f]\n' % (crds[0], crds[1], crds[2], r_sphere)
                s = s + s_color + s_sphere + 'obj += color+sphere\n'
        s = s + 'obj.append(END)\ncmd.load_cgo(obj,"cgo01")\n'
        file = open(filename, 'w')
        file.write(s)
        file.close()

    def charges_to_sites(self, charges):
        for site in self.sites:
            q = charges[site.name]
            site.set_charge(q)
        logger.info('Charges were transfered to sites')

class Multipole(GroupOfSites):
    """
    This is Multipole
    """
    def __init__(self, name=None, origin=None, representation=None, ref_multipoles=None):
        if ref_multipoles is None:
            ref_multipoles = None
        self.reference_multipoles = ref_multipoles
        GroupOfSites.__init__(self, name)
        if origin is None:
            self.origin = np.zeros(3)
        else:
            self.origin = origin
        self.representation = representation
        logger.info('%s Multipole is created.\nMultipoles representation: %s' % (self.name, self.representation))

    def get_reference_Qlm(self, l, m):
        if m < 0:
            label = '%i%is' % (l, abs(m))
        if m > 0:
            label = '%i%ic' % (l, abs(m))
        if m == 0:
            label = '%i0' % l
        try:
            return self.reference_multipoles[label]
        except KeyError:
            return 0.

    def rotate(self, rot_axis, angle):
        R = rodrigues_rotation_matrix(rot_axis, angle)
        for atom in self.atoms:
            atom.rotate(R)
        if self.representation:
            self.set_multipole_matrix()

    def translate(self, vector):
        for atom in self.atoms:
            atom.translate(vector)
        if self.representation:
            self.set_multipole_matrix()

    def move_to_center(self, center=None):
        if center is None:
            center = self.center_of_mass()
        self.translate(-center)
        if self.representation:
            self.set_multipole_matrix()

    def set_multipole_matrix(self):
        rep, l = self.representation
        if rep == 'cartesian':
            multipole = Cartesian(l, self.get_coordinates(), self.sym_sites, self.origin, self.name)
        elif rep == 'spherical':
            multipole = Spherical(l, self.get_coordinates(), self.sym_sites, self.origin, self.name)
        self.multipoles_names = multipole.names
        self.QtoM = multipole.rotation_matrix_direct
        self.QtoM_normed = np.zeros(self.QtoM.shape)
        for i, u in enumerate(self.QtoM):
            n = np.linalg.norm(u)
            if n < 0.001:
                self.QtoM_normed[i,:] = np.zeros(len(u))
            else:
                self.QtoM_normed[i,:] = u/n
        self.MtoQ = multipole.rotation_matrix_inverse

    def charges_to_multipoles(self, charges):
        try:
            self.QtoM
        except AttributeError:
            raise ValueError('Transformation matrix from charges to multipoles was not set up. Include "representation" key in the data dictionary.')
        Q = np.array([])
        for name in self.sites_names_noneq:
            Q = np.append(Q, charges[name])
        M = np.dot(self.QtoM, Q)
        multipoles = {}
        for multipole, m_value in zip(self.multipoles_names, M):
            multipoles[multipole] = m_value
        return multipoles
    
    def multipoles_to_charges(self, multipoles):
        if self.MtoQ is None:
            raise ValueError('Cannot convert multipoles to charges')
        M = np.array([])
        for multipole in self.multipoles_names:
            M = np.append(M, multipoles[multipole])
        Q = np.dot(self.MtoQ, M)
        charges = {}
        for name, q_value in zip(self.sites_names_noneq, Q):
            charges[name] = q_value
        return charges

    @property
    def charges(self):
        charges = {}
        for s in self.sites:
            charges[s.name] = s.charge
        return charges

    @property
    def multipoles(self):
        return self.charges_to_multipoles(self.charges)


class MultipoleMatrix(object):

    def __init__(self, sym_sites=None, formula=None):
        # build matrix
        rotation_matrix = np.zeros((len(self.names), len(sym_sites)))
        for i, m_name in enumerate(self.names):
            rotation_matrix[i][:] = formula.u(m_name).real
        # reduce matrix
        self.rotation_matrix_direct = np.zeros((len(self.names), max(sym_sites)+1))
        for i, _ in enumerate(self.names):
            self.rotation_matrix_direct[i] = np.bincount(sym_sites, weights=rotation_matrix[i])
        try:
            self.rotation_matrix_inverse = np.linalg.inv(self.rotation_matrix_direct) 
        except np.linalg.LinAlgError:
            self.rotation_matrix_inverse = None
        logger.debug("Multipole conversion matrix is set up.\nmultipoles = %s\ntotal number of components: %i \nQ to M matrix: %s" % (self.names, len(self.names), self.rotation_matrix_direct.shape)) 
        #logger.critical("Charges to multipoles matrix:\n%r" % self.rotation_matrix_direct)

class Spherical(MultipoleMatrix):

    def __init__(self, l=None, coordinates=None, sym_sites=None, origin=None, name=None):
        try:
            self.names = []
            for ll in xrange(l+1):
                for mm in xrange(ll+1):
                    if mm==0:
                        self.names.append('%i%i' % (ll, mm))
                    else:
                        self.names.append('%i%ic' % (ll, mm))
                        self.names.append('%i%is' % (ll, mm))
        except TypeError:
            self.names = l
        #cartesian to spherical (r, theta, phi) = (r, azimuth, polar)
        def arctan(a,b):
            if a==b==0:
                return 0.
            if b==0:
                return (-1)*np.pi*np.sign(a)/2
            else:
                return np.arctan(a/b)
        spherical = np.zeros(coordinates.shape)
        x, y, z = coordinates[:,0], coordinates[:,1], coordinates[:,2]
        x += -origin[0]
        y += -origin[1]
        z += -origin[2]
        xy2 = x**2 + y**2               # x2 + y2
        spherical[:,0] = np.sqrt(xy2 + z**2)           # r2 = x2 + y2 + z2
        spherical[:,1] = np.arctan2(y, x) # theta = arctan(y/x)
        spherical[:,2] = np.arctan2(np.sqrt(xy2), z)     # phi = arctan(xy/z)
        formula = SphericalFormulas(spherical, origin)
        MultipoleMatrix.__init__(self, sym_sites, formula)
        logger.debug('Spherical representation for %s is set up.\nOrigin = %r' % (name, origin))

class Cartesian(MultipoleMatrix):

    def __init__(self, l=None, coordinates=None, sym_sites=None, origin=None, name=None):
        self.names = []
        for i in xrange(l+1):
            self.names += self.l_to_names(i)
        formula = CartesianFormulas(coordinates, origin)
        MultipoleMatrix.__init__(self, sym_sites, formula)
        logger.debug('Cartesian representation for %s is set up.\nOrigin = %r' % (name, origin))

    def l_to_names(self, l):
        if l == 0: return ['charge']
        if l == 1: return 'X Y Z'.split()
        if l == 2: return 'XX YY ZZ XY XZ YZ'.split()

class Formulas(dict):

    def __init__(self, coordinates=None, origin=None):
        self.origin = origin
        self.coordinates = coordinates
        dict.__init__(self)


class SphericalFormulas(Formulas):

    def __init__(self, coordinates=None, origin=None):
        Formulas.__init__(self, coordinates, origin)
        self[0] = Rlm
        self['c'] = Rlmc
        self['s'] = Rlms

    def u(self, m_name):
        l, m = [int(t) for t in m_name[:2]]
        try:
            x = m_name[2] 
        except IndexError:
            x = 0
        u = np.array([])
        for crds in self.coordinates:
            r, theta, phi = crds
            u = np.append(u, self[x](l, m, r, theta, phi))
        return u


class CartesianFormulas(Formulas):

    def __init__(self, coordinates=None, origin=None):
        Formulas.__init__(self, coordinates, origin)
        self[0] = self.total_charge
        self[1] = self.dipole
        self[2] = self.quadrupole
        self[3] = self.hexadecapole

    def name_to_num(self, m_name):
        def convert(a):
            if a == 'X': return 0
            if a == 'Y': return 1
            if a == 'Z': return 2
        if m_name == 'charge':
            return
        else:
            return [convert(a) for a in m_name]
    
    def u(self, m_name):
        components = self.name_to_num(m_name)
        if m_name == 'charge': c = 0
        else: c = len(m_name) 
        u = np.array([])
        for crds in self.coordinates:
            u = np.append(u, self[c](crds, components))
        return u
    
    def total_charge(self, crds, components):
        return 1.

    def dipole(self, crds, components):
        c = components[0]
        return crds[c] - self.origin[c]
                
    def quadrupole(self, crds, components):
        a2 = np.sum(crds**2)
        m, n = components
        am = crds[m] - self.origin[m]
        an = crds[n] - self.origin[n]
        return 3.0 / 2.0 * am * an - 0.5 * a2 * self.delta(m,n)
                    
    def octapole(self, crds, components):
        m, n, k = components
        a2 = np.sum(crds**2)
        am = crds[m] - self.origin[m]
        an = crds[n] - self.origin[n]
        ak = crds[k] - self.origin[k]
        return 5. / 2. * am * an * ak - 0.5 * a2 * (am * self.delta(n,k) + an * self.delta(m,n) + ak * self.delta(m,n))
        
    def hexadecapole(self, crds, components):
        m, n, k, l = components
        am = crds[m] - self.origin[m]
        an = crds[n] - self.origin[n]
        ak = crds[k] - self.origin[k]
        al = crds[l] - self.origin[l]
        return 1. / (1. * 2. * 3. * 4.) * (105. * am * an * ak * al - 15. * a2 * (am * an * self.delta(k,l) + am * ak * self.delta(n,l) + am * al * self.delta(n,k) + an * ak * self.delta(m,l) + an * al * self.delta(m,k) + ak * al * self.delta(m,n)) + 3. * a2**2 * (self.delta(m,n) * self.delta(k,l) + self.delta(m,k) * self.delta(n,l) + self.delta(m,l) * self.delta(n,k)))

    def delta(self, i, j):
        if i==j: return 1
        else: return 0


