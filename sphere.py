from fftoolbox.multipole import Ylmc, Ylms, Multipole, GroupOfSites
from fftoolbox.atom import FFSite
from fftoolbox.molecule import AtomsInMolecule, MoleculeWithFrames
from fftoolbox.atom import Atom, MultipolarAtom
from fftoolbox.frame import BasicFrame
from units import au_to_angst
import lebedev_write

from scipy.special import sph_harm as Y
import numpy as np
import logging
import os

WORKDIR = os.path.dirname(__file__)
logger = logging.getLogger(__name__)

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"


class LebedevSphere(GroupOfSites):

    rank_to_num = {
                    1:6,
                    2:14,
                    3:26,
                    4:38,
                    5:50,
                    6:74,
                    }

    def __init__(self, index=None, name=None, rank=None, radius=None, origin=None, ref_multipoles=None, frame=None):
        GroupOfSites.__init__(self, name)
        self.set_local_frame(frame)
        if index is None:
            index = 0
        if origin is None:
            origin = np.zeros(3)
        self.origin_of_sphere = origin
        self.reference_multipoles = ref_multipoles 
        #if rank > self.reference_multipoles['rank']:
        #    raise ValueError('%s multipoles only up to rank %i are available' % (name, self.reference_multipoles['rank']))
        self.rank = rank
        self.sphere_radius = radius
        if rank == 0:
            charge = self.reference_multipoles['00']
            s = FFSite(index=index, name=self.name, coordinates=self.origin_of_sphere, charge=charge, attachment=self)
            self.add_site(s)
        else:
            s = FFSite(index=index, name=self.name, coordinates=self.origin_of_sphere, attachment=self)
            s.set_charge(0.0)
            self.add_site(s)
            points = lebedev_write.Lebedev(self.rank_to_num[rank])
            self.weights = []
            for i, point in enumerate(points):
                # coordinates and weights
                xyz = np.array(point[0:3])*radius
                w = point[3]*4*np.pi
                self.weights.append(w)
                # create site and compute charge
                name = '%s-%i' % (self.name, i)
                s = FFSite(index=index+i+1, element='EP', coordinates=xyz, attachment=self)
                q = w*self.compute_charge(rank, s.r, s.theta, s.phi)
                s.set_charge(q)
                s.set_r0(0.0)
                s.set_epsilon(0.0)
                # shift site relative to the atom center
                shifted = self.origin_of_sphere + s.coordinates
                s.set_coordinates(shifted)
                self.add_site(s)
        logger.info("LebedevSphere %s is created.\nNumber of charged sites: %i\nRadius: %.1f" % (self.name, self.num_sites, radius))

    def set_local_frame(self, frame=None):
        if not frame is None:
            self.frame = BasicFrame(frame)
        else:
            self.frame = None

    def align_with_frame(self):
        for site in self.sites:
            local_xyz = site.coordinates - self.origin_of_sphere
            global_xyz = np.dot(self.frame.local_axes, local_xyz) + self.origin_of_sphere
            site.set_coordinates(global_xyz)
        logger.info('Sphere %s was aligned with local frame' % self.name)

    def align_with_frame_and_recompute_2(self):
        points = lebedev_write.Lebedev(self.rank_to_num[self.rank])
        for j, site in enumerate(self.sites):
            if j == 0:
                continue
            i = j - 1
            point = points[i]
            # coordinates and weights
            xyz = np.array(point[0:3])*self.sphere_radius
            w = point[3]*4*np.pi
            rotated_local_xyz = np.dot(self.frame.local_axes, xyz)
            site.set_coordinates(rotated_local_xyz)
            # compute charge
            q = w*self.compute_charge(self.rank, site.r, site.theta, site.phi)
            site.set_charge(q)
            site.set_r0(0.0)
            site.set_epsilon(0.0)
            # shift site relative to the atom center
            shifted = self.origin_of_sphere + site.coordinates
            site.set_coordinates(shifted)
        logger.info('Sphere %s was aligned with local frame and charges recomputed' % self.name)

    def align_with_frame_and_recompute(self):
        for j, site in enumerate(self.sites):
            if j == 0:
                continue
            i = j - 1
            original_local_xyz = site.coordinates - self.origin_of_sphere
            rotated_local_xyz = np.dot(self.frame.local_axes, original_local_xyz)
            site.set_coordinates(rotated_local_xyz)
            w = self.weights[i]
            q = w*self.compute_charge(self.rank, site.r, site.theta, site.phi)
            site.set_charge(q)
            new_global_xyz = self.origin_of_sphere + site.coordinates
            site.set_coordinates(new_global_xyz)
        logger.info('Sphere %s was aligned with local frame and charges recomputed' % self.name)

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

class LebedevAtom(MultipolarAtom, LebedevSphere):

    def __init__(self, element=None, coordinates=None, index=None, ref_multipoles=None, rank=None, radius=None, sym=False):
        representation = ('spherical', rank)
        MultipolarAtom.__init__(self, index=index, element=element, coordinates=coordinates, representation=representation, sym=sym)
        LebedevSphere.__init__(self, index=index, name=self.name, rank=rank, radius=radius, \
                origin=coordinates, ref_multipoles=ref_multipoles)
        self.set_sym_sites()
        self.set_multipole_matrix()

class LebedevMolecule(LebedevSphere, AtomsInMolecule, Multipole):

    def __init__(self, data):
        self.data = data
        rank, radius = data['sphere params']
        self.rank = rank
        representation = ('spherical', rank)
        atoms = data['atoms']
        name = data['name']
        try:
            sym = data['symmetry']
        except KeyError:
            sym = False
        try:
            frame = data['frame']
        except KeyError:
            frame = None
        AtomsInMolecule.__init__(self, name=name, atoms=atoms, sym=sym)
        origin = self.center_of_mass()
        Multipole.__init__(self, name=name, origin=origin, representation=representation)
        LebedevSphere.__init__(self, name=name, rank=rank, radius=radius, \
                origin=origin, ref_multipoles=data['multipoles'], frame=frame)
        if not self.frame is None:
            self.align_with_frame()
        self.set_sym_sites()
        self.set_multipole_matrix()

    def multipoles_data(self):
        data = {
                self.name: (self.origin, self.rank, self.reference_multipoles), 
                }
        return data

    def color_charges(self, filename, xyzname=False,vmax=None, r_sphere=0.08):
        path2xyz = '%s/data/xyz/%s' % (WORKDIR, xyzname)
        s = 'from pymol.cgo import *\nfrom pymol import cmd\ncmd.load("%s")\nobj = [ BEGIN, LINES, ]\n' % (path2xyz)
        vmax = max([abs(ss.charge) for ss in self.sites]) + 0.3
        for site in self.sites:
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


class DistributedLebedevMolecule(MoleculeWithFrames, Multipole):

    def __init__(self, data):
        self.sphere_params = data['sphere params']
        atoms = data['atoms']
        name = data['name']
        try:
            self.theory = data['theory']
        except KeyError:
            self.theory = None
        try:
            sym = data['symmetry']
        except KeyError:
            sym = False
        try:
            representation = data['representation']
        except KeyError:
            representation = None
        Multipole.__init__(self, name=name, representation=representation)
        MoleculeWithFrames.__init__(self, name, atoms, sym)
        self.set_frames(self.framed_atoms)
        for frame in self.frames:
            frame.center.align_with_frame_and_recompute()
        self.set_sym_sites()
        if representation is not None:
            self.set_multipole_matrix()

    def multipoles_data(self):
        data = {}
        for a in self.atoms:
            d = (a.origin, a.rank, a.reference_multipoles)
            data[a.name] = d
        return data

    def add_atoms(self, atoms):
        self.framed_atoms = []
        for atom in atoms:
            element = atom['element']
            crds = atom['coordinates']
            multipoles = atom['multipoles']
            try:
                rank, radius = self.sphere_params[element]
            except:
                rank, radius = self.sphere_params
            index = self.get_max_index() + 1
            atom = LebedevAtom(element=element, coordinates=crds, ref_multipoles=multipoles, \
                    index=index, rank=rank, radius=radius, sym=self.sym)
            self.add_atom(atom)
            # only atoms with rank > 0 require frame
            if rank > 0:
                self.framed_atoms.append(atom.element)
            
    @property
    def sites(self):
        sites = []
        for atom in self.atoms:
            for site in atom.sites:
                sites.append(site)
        return iter(sorted(sites, key=lambda s: s.index))

    @property
    def num_sites(self):
        sites = []
        for atom in self.atoms:
            for site in atom.sites:
                sites.append(site)
        return len(sites)
    def __repr__(self):
        o = '%s\n' % self.name
        for s in self.sites:
            o += '%s %s xyz = %.3f %.3f %.3f charge = %s r0 = %s epsilon = %s\n' % (s.element, s.name, s.x, s.y, s.z, s.charge, s.r0, s.epsilon)
        return o

