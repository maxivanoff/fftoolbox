import numpy as np
from numpy.linalg import norm
from scipy.special import sph_harm as Y
import openbabel
import logging
import os

from units import au_to_angst
import parser

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

WORKDIR = os.path.dirname(__file__)
logger = logging.getLogger('multipole')

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
    return r**l * np.sqrt(4 * np.pi / (2 * l + 1)) * Ylmc(l, m, theta, phi)

def Rlms(l, m, r, theta, phi):
    return r**l * np.sqrt(4 * np.pi / (2 * l + 1)) * Ylms(l, m, theta, phi)

def Ylmc(l, m, theta, phi):
    #v = np.sqrt(0.5) * (np.conj(Y(m, l, theta, phi)) + Y(m, l, theta, phi))
    v = np.sqrt(0.5) * (Y(-m, l, theta, phi) + (-1)**m*Y(m, l, theta, phi))
    #v = np.sqrt(0.5) * ((-1)**m*Y(-m, l, theta, phi) + Y(m, l, theta, phi))
    if abs(v.imag) > 0.0001: raise ValueError("Non-zero imaginary part in Ylmc")
    return v.real

def Ylms(l, m, theta, phi):
    #v = 1j * np.sqrt(0.5) * (np.conj(Y(m, l, theta, phi)) - Y(m, l, theta, phi))
    #v = 1j * np.sqrt(0.5) * (Y(-m, l, theta, phi) - (-1)**m*Y(m, l, theta, phi))
    v = 1j * np.sqrt(0.5) * (-(-1)**m*Y(-m, l, theta, phi) + Y(m, l, theta, phi))
    if abs(v.imag) > 0.0001: raise ValueError("Non-zero imaginary part in Ylms")
    return v.real


class GroupOfAtoms(object):

    def __init__(self, name=None):
        self.name = name
        self.i = -1

    def set_sym_sites(self):
        sites = {}
        sym_sites = []
        for i, name in enumerate(self.sites_names_eq):
            if not name in sites:
                sites[name] = len(sites.keys())
            sym_sites.append(sites[name])
        self.sym_sites = np.array(sym_sites)
        logger.debug("Symmetry indeces:\n%r" % self.sym_sites)
                
    def get_coordinates(self):
        crds = np.zeros((len(self.sites), 3))
        for i, s in enumerate(self.sites):
            crds[i][:] = s.coordinates[:]
        return crds

    def get_max_index(self):
        i = [s.index for s in self.sites]
        return max(i)
    
    def get_sites(self, name):
        return filter(lambda s: s.name==name, self.sites)
    
    def get_atoms_by_element(self, element):
        return filter(lambda a: a.element==element, self.atoms)

    def get_atom_by_index(self, index):
        return filter(lambda a: a.index==index, self.atoms)[0]

    def get_sites_by_name(self, name):
        return filter(lambda s: s.name == name, self.sites)

    def get_atom(self, index):
        return next(a for a in self.atoms if a.index==index)

    @property
    def atoms_names_noneq(self):
        return [a.name for a in self.atoms_noneq]

    @property
    def atoms_names_eq(self):
        return [a.name for a in self.atoms]
    
    @property
    def sites_names_noneq(self):
        return [s.name for s in self.sites_noneq]

    @property
    def sites_names(self):
        return self.sites_names_noneq

    @property
    def sites_names_eq(self):
        return [s.name for s in self.sites]
    
    @property
    def sites(self):
        sites = []
        for atom in self:
            sites += atom.sites
        return sorted(sites, key=lambda s: s.index)

    @property
    def sites_noneq(self):
        sites = []
        for s in self.sites:
            if not s.name in [ss.name for ss in sites]:
                sites.append(s)
        return sites
    
    @property
    def atoms_noneq(self):
        atoms = []
        for a in self.atoms:
            if not a.name in [aa.name for aa in atoms]:
                atoms.append(a)
        return atoms

    def write_xyz(self, filename, here=False):
        out = '%i\n' % len(self.sites)
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
        """
        xyzname = '%s_tmp.xyz' % filename.split('.mol2')[0] 
        fullxyzname = '%s/%s' % (path2xyz, xyzname)
        self.write_xyz(xyzname, here=here)
        obmol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("xyz", "mol2")
        conv.ReadFile(obmol, fullxyzname)
        a = obmol.GetAtom(1)
        a.GetPartialCharge()
        for i, site in enumerate(self.sites):
            a = obmol.GetAtom(i+1)
            print a.GetType(), site.charge
            a.SetPartialCharge(site.charge)
        conv.WriteFile(obmol, '%s/%s' % (path2mol2, filename))
        os.system('rm %s' % (fullxyzname))
        """
        

    def __iter__(self):
        return self

    def next(self):
        if self.i < len(self.atoms)-1:
            self.i += 1
            return self.atoms[self.i]
        else:
            self.i = -1
            raise StopIteration

class Multipole(GroupOfAtoms):
    """
    This is Multipole
    """
    def __init__(self, name=None, origin=None, multipoles=('cartesian', 0)):
        GroupOfAtoms.__init__(self, name)
        self.origin = origin
        self.m_rep = multipoles[0]
        self.l = multipoles[1]
        logger.info('%s Multipole is created.\nMultipoles representation: %s\nl = %i' % (self.name, self.m_rep, self.l))

    def rotate(self, rot_axis, angle):
        R = rodrigues_rotation_matrix(rot_axis, angle)
        for site in self.sites:
            site.set_coordinates(np.dot(R, site.coordinates))
        self.set_multipole_matrix()

    def translate(self, vector):
        for site in self.sites:
            site.set_coordinates(site.coordinates + vector)
        self.set_multipole_matrix()

    def set_multipole_matrix(self):
        if self.m_rep == 'cartesian':
            multipole = Cartesian(self.l, self.get_coordinates(), self.sym_sites, self.origin)
        elif self.m_rep == 'spherical':
            multipole = Spherical(self.l, self.get_coordinates(), self.sym_sites, self.origin)
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
        logger.debug("Multipole conversion matrix is set up.\nmultipoles = %s; total number of components: %i \nQ to M matrix: %s" % (self.names, len(self.names), self.rotation_matrix_direct.shape)) 

class Spherical(MultipoleMatrix):

    def __init__(self, l=None, coordinates=None, sym_sites=None, origin=None):
        try:
            self.names = []
            for ll in xrange(l):
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
        #r = np.sqrt(x**2 + y**2 + z**2)
        #phi = np.arccos(z/r)
        #theta = np.array([])
        #for xx, yy in zip(x,y):
        #    if yy>=0 and xx>0:
        #        s = 0
        #    if xx<=0:
        #        s = np.pi
        #    if xx>0 and yy<0:
        #        s = 2*np.pi
        #    if xx==0 and yy==0:
        #        s = 0
        #    theta = np.append(theta, arctan(yy,xx) + s)
        #spherical[:,0] = r 
        #spherical[:,1] = theta 
        #spherical[:,2] = phi
        xy2 = x**2 + y**2               # x2 + y2
        spherical[:,0] = np.sqrt(xy2 + z**2)           # r2 = x2 + y2 + z2
        spherical[:,1] = np.arctan2(y, x) # theta = arctan(y/x)
        spherical[:,2] = np.arctan2(np.sqrt(xy2), z)     # phi = arctan(xy/z)
        formula = SphericalFormulas(spherical, origin)
        MultipoleMatrix.__init__(self, sym_sites, formula)

class Cartesian(MultipoleMatrix):

    def __init__(self, l=None, coordinates=None, sym_sites=None, origin=None):
        self.names = []
        for i in xrange(l+1):
            self.names += self.l_to_names(i)
        formula = CartesianFormulas(coordinates, origin)
        MultipoleMatrix.__init__(self, sym_sites, formula)

    def l_to_names(self, l):
        if l == 0: return ['charge']
        if l == 1: return 'X Y Z'.split()
        if l == 2: return 'XX YY ZZ XY XZ YZ'.split()

class Formulas(dict):

    def __init__(self, coordinates=None, origin=None):
        self.coordinates = coordinates
        if origin == None:
            self.origin = np.zeros(3)
        else:
            self.origin = origin
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


