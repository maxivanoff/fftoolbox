from __future__ import division
import os
import re
import numpy as np
import xml.etree.ElementTree as ET
import units
import logging

parserLogger = logging.getLogger('parser')

WORKDIR = os.path.dirname(__file__)

class GaussianCube(object):
    
    number_to_name = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 17:'Cl', 16:'S', 0:'X'}
        
    def __init__(self, cubefile=None, filename=None, data=None):
        "Cubefile is in atomic units"
        self.property = 'esp'
        self.origin = None
        self.vectors = np.zeros((3,3), dtype=np.float64)
        self.num_points = np.zeros(3, dtype=np.int)
        self.values = ()
        self.atoms = ()
        self.density = None
        if data:
            filename = '%s/data/cub/%s_%s_d%s.cub' % (WORKDIR, data['name'], data['theory'], data['density'])
        self.read_file(filename=filename, cubefile=cubefile)

    def set_grids_density(self):
        cube_sides = np.array([np.linalg.norm(v) for v in self.vectors])
        self.density = np.prod(self.num_points)/np.prod((self.num_points - 1)*cube_sides)


    def read_file(self, filename=None, cubefile=None):
        if filename:
            cubefile = open(filename, 'r')
        cubefile.readline()
        cubefile.readline()

        # Read number of atoms and origin coordinate
        tmp = np.array([float(t) for t in cubefile.readline().split()])
        self.num_atoms = int(tmp[0])
        self.origin = np.array(tmp[1:4])

        # Read cell geometry
        for i in xrange(3):
            tmp = np.array([float(t) for t in cubefile.readline().split()])
            self.num_points[i] = tmp[0]
            self.vectors[i][:] = np.array(tmp[1:4])

        #Read atoms coordinates
        atoms = list()
        for ati in xrange(self.num_atoms):
            line = cubefile.readline().split()
            crds = np.array([float(t) for t in line[2:]])
            atomic_n = int(line[0])
            element = self.number_to_name[atomic_n]
            atoms.append((element, crds, None))
        self.atoms = tuple(atoms)

        # Read values
        self.values = []
        while True:
            line = cubefile.readline()
            if line == '':
                break
            else:
                tmp = [float(t) for t in line.split()]
                self.values += tmp
        if filename:
            cubefile.close()

        self.data = {
                'origin': self.origin,
                'vectors': self.vectors,
                'num_points': self.num_points,
                'atoms': self.atoms,
                'values': self.values,
                'property': self.property,
                }

    def __str__(self):
        return 'Gaussian cube parser'

class QChem(object):

    def __init__(self, filename=None, data=None):
        self.atoms = ()
        self.multipoles = {}
        if data:
            filename = '%s/data/qcout/%s_%s.qcout' % (WORKDIR, data['name'], data['theory'])
        self.read_file(filename)

    def read_file(self, filename):
        txt = open(filename, 'r').read()

        # coordinates
        s_coords = txt.split('Orientation (Angstroms)')[-1].split('Nuclear Repulsion Energy')[0]
        atoms = []
        for s_xyz in s_coords.split('\n')[2:-1]:
            try:
                tmp = s_xyz.split()
                element = tmp[1]
                crds = np.array([float(t) for t in tmp[2:5]])*units.angst_to_au
                atoms.append((element, crds, None))
            except IndexError:
                pass
        self.atoms = (atoms)

        #multipoles
        s_multipoles = txt.split('The traceless molecular multipole moments')[-1].split('Multipole moments of the Stewart Charges')[0]
        s_charge = s_multipoles.split('Charge (ESU x 10^10)')[1].split('Dipole Moment (Debye)')[0].split()[0]
        multipoles = [float(s_charge)]
        s_dipole = s_multipoles.split('Dipole Moment (Debye)')[1].split('Quadrupole Moments (Debye-Ang)')[0].split()
        multipoles += [float(s_dipole[1])*units.debye_to_au, \
                       float(s_dipole[3])*units.debye_to_au, \
                       float(s_dipole[5])*units.debye_to_au]
        s_q = s_multipoles.split('Quadrupole Moments (Debye-Ang)')[1].split('Octopole Moments (Debye-Ang^2)')[0].split()
        for i in range(5):
            multipoles.append(float(s_q[2*i+1])*units.debye_to_au*units.angst_to_au)
        """
        s_o = s_multipoles.split('Octopole Moments (Debye-Ang^2)')[1].split('Hexadecapole Moments (Debye-Ang^3)')[0].split()
        for i in range(7):
            multipoles.append(float(s_o[2*i+1])*units.debye_to_au*units.angst_to_au**2)
        s_h = s_multipoles.split('Hexadecapole Moments (Debye-Ang^3)')[1].split('5th-Order Moments (Debye-Ang^4)')[0].split()
        for i in range(9):
            multipoles.append(float(s_h[2*i+1])*units.debye_to_au*units.angst_to_au**3)
        s_fifth = s_multipoles.split('5th-Order Moments (Debye-Ang^4)')[1].split('6th-Order Moments (Debye-Ang^5)')[0].split()
        for i in range(11):
            multipoles.append(float(s_fifth[2*i+1])*units.debye_to_au*units.angst_to_au**4)
        """
        names = []
        for l in xrange(6):
            for m in xrange(l+1):
                if m > 0:
                    names.append('%i%ic' % (l,m))
                    names.append('%i%is' % (l,m))
                else:
                    names.append('%i%i' % (l,m))
        self.multipoles = dict((a,b) for a,b in zip(names, multipoles))

        self.data = {
                'atoms': self.atoms,
                'multipoles': self.multipoles,
                }

    def __str__(self):
        return 'QChem parser'

class Gaussian(object):
    atom_name = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 17:'Cl', 0:'X'}

    def __init__(self, filename=None, data=None):
        self.energy = None
        self.atoms = ()
        self.multipoles = {}
        if data:
            filename = '%s/data/log/%s_%s.log' % (WORKDIR, data['name'], data['theory'])
        self.read_file(filename)

    def read_file(self, filename):
        self.geometries_input = []
        self.geometries_standard = []
        self.energies = []
        self.pcm_energy = None
        self.E2 = None
        logfile = open('%s' % (filename), 'r')
        while True:
            line = logfile.readline()
            # end of file
            if re.search(r'Normal termination',line) or re.search(r'Error termination',line): break
            # level of theory
            m = re.search(r' \#P (\w+\/.+)', line)
            if m:
                self.theory_level = m.group(1).split()[0]
            # geometries in input orientation
            if re.search(r'Input orientation:', line):
                atoms = list()
                for i in range(4): logfile.readline()
                while True:
                    line = logfile.readline()
                    m = re.search(r'(\d) *(\d) *(\d) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*)', line)
                    if m is None: break
                    atomic_num = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    crds = np.array([x, y, z])*units.angst_to_au
                    atoms.append((atomic_num, crds))
                self.geometries_input.append(atoms)
            # electronic energy
            scf = re.search(r'SCF Done:  E\(.*\) = *(-?\d*\.\d*)', line)
            pcm = re.search(r'After PCM corrections, the energy is *(-?\d*\.\d*)', line)
            E2 = re.search(r'E2 = *(-?\d*\.\d*)', line)
            if pcm:
                self.pcm_energy = float(pcm.group(1))
            if scf:
                self.energies.append(float(scf.group(1)))
            if E2:
                self.E2 = float(E2.group(1)) 
            # geometries in standard orientation
            if re.search(r'Standard orientation:', line):
                atoms = list()
                for i in range(4): logfile.readline()
                while True:
                    line = logfile.readline()
                    m = re.search(r'(\d) *(\d) *(\d) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*)', line)
                    if m is None: break
                    atomic_num = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    crds = np.array([x, y, z])*units.angst_to_au
                    element = self.atom_name[atomic_num]
                    atoms.append((element, crds, None))
                self.geometries_standard.append(atoms)

            monopole = re.search(r'Charge= *(-?\d\.\d*) electrons', line)
            if monopole:
                self.multipoles['charge'] = float(monopole.group(1))
            dipole = re.search(r'X= *(-?\d\.\d*) *Y= *(-?\d\.\d*) *Z= *(-?\d\.\d*)', line)
            if dipole:
                d = np.array([float(dipole.group(i+1)) for i in xrange(3)])*units.debye_to_au
                self.multipoles['X'] = d[0]
                self.multipoles['Y'] = d[1]
                self.multipoles['Z'] = d[2]
            quadrupole1 = re.search(r'XX= *(-?\d\.\d*) *YY= *(-?\d\.\d*) *ZZ= *(-?\d\.\d*)', line)
            if quadrupole1:
                Q1 = np.array([float(quadrupole1.group(i+1)) for i in xrange(3)])*units.debye_to_au*units.angst_to_au
                self.multipoles['XX'] = Q1[0]
                self.multipoles['YY'] = Q1[1]
                self.multipoles['ZZ'] = Q1[2]
            quadrupole2 = re.search(r'XY= *(-?\d\.\d*) *XZ= *(-?\d\.\d*) *YZ= *(-?\d\.\d*)', line)
            if quadrupole2:
                Q2 = np.array([float(quadrupole2.group(i+1)) for i in xrange(3)])*units.debye_to_au*units.angst_to_au
                self.multipoles['XY'] = Q2[0]
                self.multipoles['XZ'] = Q2[1]
                self.multipoles['YZ'] = Q2[2]
            """
            octapole1 = re.search(r'XXX= *(-?\d\.\d*) *YYY= *(-?\d\.\d*) *ZZZ= *(-?\d\.\d*) *XYY= *(-?\d\.\d*)', line)
            if octapole1:
                O1 = np.array([float(octapole1.group(i+1)) for i in xrange(4)])*units.debye_to_au*units.angst_to_au**2
                self.multipoles['XXX'] = O1[0]
                self.multipoles['YYY'] = O1[1]
                self.multipoles['ZZZ'] = O1[2]
                self.multipoles['XYY'] = O1[3]
            octapole2 = re.search(r'XXY= *(-?\d\.\d*) *XXZ= *(-?\d\.\d*) *XZZ= *(-?\d\.\d*) *YZZ= *(-?\d\.\d*)', line)
            if octapole2:
                O2 = np.array([float(octapole2.group(i+1)) for i in xrange(4)])*units.debye_to_au*units.angst_to_au**2
                self.multipoles['XXY'] = O2[0]
                self.multipoles['XXZ'] = O2[1]
                self.multipoles['XZZ'] = O2[2]
                self.multipoles['YZZ'] = O2[3]
            octapole3 = re.search(r'YYZ= *(-?\d\.\d*) *XYZ= *(-?\d\.\d*)', line)
            if octapole3:
                O3 = np.array([float(octapole3.group(i+1)) for i in xrange(2)])*units.debye_to_au*units.angst_to_au**2
                self.multipoles['YYZ'] = O3[0]
                self.multipoles['XYZ'] = O3[1]
            """

        self.atoms = self.geometries_standard[-1]
        if self.pcm_energy:
            self.energy = self.pcm_energy
        elif self.E2:
            self.energy = self.energies[-1] + self.E2
        else:
            self.energy = self.energies[-1]

        self.data = {
                'atoms': self.atoms,
                'energy': self.energy,
                'multipoles': self.multipoles,
                }
        

    def __str__(self):
        return 'Gaussian log parser'

class Mol2(object):

    def __init__(self, filename=None, data=None):
        if data:
            filename = '%s/data/mol2/%s_%s.mol2' % (WORKDIR, data['name'], data['theory'])
        if filename:
            filename = '%s/data/mol2/%s' % (WORKDIR, filename)
        self.atoms = ()
        self.read_file(filename)

    def read_file(self, filename):
        mol2file = open(filename, 'r')
        s = mol2file.read()
        s = s.split('@<TRIPOS>')
        self.num_atoms = int(s[1].split()[2])
        geom = s[2].split('\n')
        atoms = []
        for line in geom[1:self.num_atoms+1]:
            tmp = line.split()
            atom_name = tmp[1]
            charge = float(tmp[8])
            crds = np.array([float(t)*units.angst_to_au for t in tmp[2:5]])
            atoms.append((atom_name, crds, charge))
        self.atoms = (atoms)
        mol2file.close()
        self.data = {'atoms': self.atoms}

    def __str__(self):
        return 'mol2 parser'

class XYZ(object):

    def __init__(self, filename=None, data=None, here=True):
        if here == True:
            path2xyz = '.'
        else:
            path2xyz = '%s/data/xyz' % WORKDIR
        if data:
            filename = '%s/%s_%s.xyz' % (path2xyz, data['name'], data['theory'])
        else:
            filename = '%s/%s' % (path2xyz, filename)
        self.read_file(filename)

    def read_file(self, filename):
        xyzfile = open(filename, 'r')
        s = xyzfile.readline().split()
        num_atoms = int(s[0])
        xyzfile.readline()
        atoms = []
        for i in xrange(num_atoms):
            tmp = xyzfile.readline().split()
            atom_name = tmp[0]
            crds = np.array([float(t)*units.angst_to_au for t in tmp[1:4]])
            atoms.append((atom_name, crds, None))
        xyzfile.close()
        self.data = {'atoms': atoms}

    def __str__(self):
        return 'xyz parser'

