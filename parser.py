from __future__ import division
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.etree import ElementTree
from xml.dom import minidom
import numpy as np
import logging
import os
import re

import units

__author__ = "Maxim Ivanov"
__email__ = "maxim.ivanov@marquette.edu"

logger = logging.getLogger(__name__)
WORKDIR = os.path.dirname(__file__)

class Parser(object):

    def __init__(self, filename=None, data=None, here=False, dname=None, suffix=None):
        if here == True:
            path2xyz = '.'
        else:
            path2xyz = '%s/data/%s' % (WORKDIR, dname)
        if data:
            self.filename = '%s/%s_%s%s' % (path2xyz, data['name'], data['theory'], suffix)
        else:
            self.filename = '%s/%s' % (path2xyz, filename)
        self.atoms = []

    def add_atom(self, index=None, element=None, crds=None, charge=None, multipoles=None):
        if multipoles is None:
            multipoles = {}
            multipoles['charge'] = charge
        atom = (index, element, crds, multipoles)
        self.atoms.append(atom)

    def read_file(self, filename):
        logger.info('%s will read data from %s' % (self, filename))

    def __repr__(self):
        return 'Just Parser'

class GaussianCube(Parser):
    
    number_to_name = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 17:'Cl', 16:'S', 35:'Br'}
    name_to_number = dict(zip(number_to_name.values(), number_to_name.keys()))
        
    def __init__(self, filename=None, data=None, here=False):
        "Cubefile is in atomic units"
        self.property = 'esp'
        self.origin = None
        self.vectors = np.zeros((3,3), dtype=np.float64)
        self.num_points = np.zeros(3, dtype=np.int)
        self.values = ()
        self.atoms = ()
        self.density = None
        if data:
            try:
                suffix = '_d%s.cub' % data['density']
            except KeyError:
                suffix = '.cub'
            Parser.__init__(self, filename=filename, data=data, dname='cub', suffix=suffix, here=here)
            self.read_file(filename=self.filename)
        if filename:
            Parser.__init__(self, filename=filename, here=here)
            self.read_file(filename=self.filename)


    def set_grids_density(self):
        cube_sides = np.array([np.linalg.norm(v) for v in self.vectors])
        self.density = np.prod(self.num_points)/np.prod((self.num_points - 1)*cube_sides)

    def write_file(self, grid=None, molecule=None, values=None, filename=None, here=False):
        if here is False:
            filename = '%s/data/cub/local/%s' % (WORKDIR, filename)
        s = ' %s Potential\n Electrostatic potential\n' % molecule.name
        s += '%5s' % molecule.num_atoms
        a = '%.6f' % grid.origin[0]
        b = '%.6f' % grid.origin[1]
        c = '%.6f' % grid.origin[2]
        s += '%12s%12s%12s    1\n' % (a, b, c)
        for i in xrange(3):
            a = '%.6f' % grid.vectors[i][0]
            b = '%.6f' % grid.vectors[i][1]
            c = '%.6f' % grid.vectors[i][2]
            s += '%5s%12s%12s%12s\n' % (grid.num_cubic_points[i], a, b, c)
        for a in molecule.atoms:
            atomn = self.name_to_number[a.element]
            ns = '%.6f' % atomn
            xs = '%.6f' % a.x
            ys = '%.6f' % a.y
            zs = '%.6f' % a.z
            s += '%5s%12s%12s%12s%12s\n' % (atomn, ns, xs, ys, zs)
        # write values
        nz = grid.num_cubic_points[2]  
        countz = nz   
        countc = 1
        for i, v in enumerate(values):
            x = '%.5E' % v
            s += '%13s' % x 
            if countc == 6:
                s += '\n' 
                countc = 0
            if i == countz - 1:
                if not s[-1] == '\n': s += '\n'
                countz += nz
                countc = 0
            countc += 1
        file = open(filename, 'w')
        file.write(s)
        file.close()
        logger.info('Cubefile was written to %s' % filename)


    def read_file(self, filename=None):
        Parser.read_file(self, filename)
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
        for ati in xrange(self.num_atoms):
            line = cubefile.readline().split()
            crds = np.array([float(t) for t in line[2:]])
            atomic_n = int(line[0])
            element = self.number_to_name[atomic_n]
            index = ati + 1
            self.add_atom(index, element, crds)
        # Read values
        values = []
        while True:
            line = cubefile.readline()
            if line == '':
                break
            else:
                tmp = [float(t) for t in line.split()]
                values += tmp
        self.values = np.array(values)
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

class QChem(Parser):

    def __init__(self, filename=None, data=None, here=False):
        self.atoms = ()
        self.multipoles = {}
        Parser.__init__(self, filename=filename, data=data, dname='qcout', suffix='.qcout', here=here)
        self.read_file(filename=self.filename)

    def read_file(self, filename):
        Parser.read_file(self, filename)
        txt = open(filename, 'r').read()

        # coordinates
        s_coords = txt.split('Orientation (Angstroms)')[-1].split('Nuclear Repulsion Energy')[0]
        atoms = []
        index=1
        for s_xyz in s_coords.split('\n')[2:-1]:
            try:
                tmp = s_xyz.split()
                element = tmp[1]
                crds = np.array([float(t) for t in tmp[2:5]])*units.angst_to_au
                self.add_atom(index, element, crds)
                index+=1
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

class Gaussian(Parser):
    atom_name = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 17:'Cl', 0:'X', 16:'S', 35:'Br'}

    def __init__(self, filename=None, data=None, here=False, orientation='standard'):
        self.energy = None
        self.multipoles = {}
        self.orientation=orientation
        Parser.__init__(self, filename=filename, data=data, dname='log', suffix='.log', here=here)
        self.read_file(filename=self.filename)

    def read_file(self, filename, orientation=None):
        Parser.read_file(self, filename)
        if orientation is None:
            orientation = self.orientation
        self.geometries_input = []
        self.geometries_standard = []
        self.energies = []
        self.pcm_energy = None
        self.s1_energies = []
        self.E2 = None
        logfile = open('%s' % (filename), 'r')
        while True:
            line = logfile.readline()
            if line == '': break
            # end of file
            # if re.search(r'Normal termination',line) or re.search(r'Error termination',line): break
            # level of theory
            m = re.search(r' \#P (\w+\/.+)', line)
            if m:
                self.theory_level = m.group(1).split()[0]
            # geometries in input orientation
            if re.search(r'Input orientation:', line):
                self.atoms = []
                for i in range(4): logfile.readline()
                index=1
                while True:
                    line = logfile.readline()
                    m = re.search(r'(\d+) *(\d+) *(\d+) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*)', line)
                    if m is None: break
                    atomic_num = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    crds = np.array([x, y, z])*units.angst_to_au
                    element = self.atom_name[atomic_num]
                    self.add_atom(index, element, crds)
                    index+=1
                self.geometries_input.append(self.atoms)
            # electronic energy
            scf = re.search(r'SCF Done:  E\(.*\) = *(-?\d*\.\d*)', line)
            S1_energy = re.search(r' Total Energy, E\(TD-HF/TD-KS\) = *(-?\d*\.\d*)', line)
            pcm = re.search(r'After PCM corrections, the energy is *(-?\d*\.\d*)', line)
            E2 = re.search(r'E2 = *(-?\d*\.\d*)', line)
            if pcm:
                self.pcm_energy = float(pcm.group(1))
            if scf:
                self.energies.append(float(scf.group(1)))
            if E2:
                self.E2 = float(E2.group(1)) 
            if S1_energy:
                self.s1_energies.append(float(S1_energy.group(1)))
            # geometries in standard orientation
            if re.search(r'Standard orientation:', line):
                self.atoms = []
                for i in range(4): logfile.readline()
                index=1
                while True:
                    line = logfile.readline()
                    m = re.search(r'(\d+) *(\d+) *(\d+) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*)', line)
                    if m is None: break
                    atomic_num = int(m.group(2))
                    x = float(m.group(4))
                    y = float(m.group(5))
                    z = float(m.group(6))
                    crds = np.array([x, y, z])*units.angst_to_au
                    element = self.atom_name[atomic_num]
                    self.add_atom(index, element, crds)
                    index+=1
                self.geometries_standard.append(self.atoms)

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

        if orientation == 'standard':
            self.atoms = self.geometries_standard[-1]
        if orientation == 'input':
            self.atoms = self.geometries_input[-1]
        if self.pcm_energy:
            self.energy = self.pcm_energy
        elif self.E2:
            self.energy = self.energies[-1] + self.E2
        else:
            self.energy = self.energies[-1]
        if self.s1_energies:
            self.S1 = self.s1_energies[-1]
        else:
            self.S1 = None
        self.data = {
                'atoms': self.atoms,
                'ground energy': self.energy,
                'S1 energy': self.S1,
                'multipoles': self.multipoles,
                }
        

    def __str__(self):
        return 'Gaussian log parser'

class Mol2(Parser):

    def __init__(self, filename=None, data=None, here=False):
        Parser.__init__(self, filename=filename, data=data, dname='mol2', suffix='.mol2', here=here)

    def read_file(self, filename=None):
        if filename is None:
            filename = self.filename
        Parser.read_file(self, filename)
        mol2file = open(filename, 'r')
        s = mol2file.read()
        s = s.split('@<TRIPOS>')
        self.num_atoms = int(s[1].split()[2])
        geom = s[2].split('\n')
        atoms = []
        index = 1
        for line in geom[1:self.num_atoms+1]:
            tmp = line.split()
            element = tmp[1]
            charge = float(tmp[8])
            crds = np.array([float(t)*units.angst_to_au for t in tmp[2:5]])
            self.add_atom(index, element, crds, charge=charge)
            index += 1
        mol2file.close()
        self.data = {'atoms': self.atoms}

    def write_file(self, filename, here=False, molecule=None):
        if here:
            path2mol2 = '.'
        else:
            path2mol2 = '%s/data/mol2' % WORKDIR
        s = '@<TRIPOS>MOLECULE\n%s\n%i %i 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n' % (molecule.name, molecule.num_sites, len(molecule.bonds))
        try:
            molecules = molecule.molecules
        except:
            molecules = [molecule]
        for mi, mol in enumerate(molecules):
            for site in mol.sites:
                sid = '%i' % site.index
                aname = '%s' % amberType[site.name.split('-')[0]]
                x = '% .3f' % (site.x*units.au_to_angst)
                y = '% .3f' % (site.y*units.au_to_angst)
                z = '% .3f' % (site.z*units.au_to_angst)
                smi = '%i' % (mi+1)
                myname = '%s' % site.name
                charge = '% .4f' % site.charge
                mname = '%-3.3s' % mol.name
                s += ' %s %s %s %s %s %s %s %s %s\n' %(sid.ljust(3), aname.ljust(3), x.ljust(8), y.ljust(8), z.ljust(8), aname.ljust(3), smi.ljust(2), mname, charge.ljust(7))
        s += '@<TRIPOS>BOND\n'
        for i, bond in enumerate(molecule.bonds):
            bid = '%i' % (i+1)
            a1 = '%i' % bond.a1.index
            a2 = '%i' % bond.a2.index
            s += ' %s %s %s 1\n' % (bid.ljust(3), a1.ljust(3), a2.ljust(3))
        with open('%s/%s' % (path2mol2, filename), 'w') as f:
            f.write(s)

    def __str__(self):
        return 'mol2 parser'

class XYZ(Parser):

    def __init__(self, filename=None, data=None, here=False):
        Parser.__init__(self, filename=filename, data=data, dname='xyz', suffix='.xyz')

    def read_file(self, filename):
        Parser.read_file(self, filename)
        xyzfile = open(filename, 'r')
        s = xyzfile.readline().split()
        num_atoms = int(s[0])
        comment = xyzfile.readline()
        atoms = []
        for i in xrange(num_atoms):
            tmp = xyzfile.readline().split()
            element = tmp[0]
            crds = np.array([float(t)*units.angst_to_au for t in tmp[1:4]])
            index = i + 1
            self.add_atom(index, element, crds)
        xyzfile.close()
        self.data = {'atoms': atoms,
                     'comment': comment,
                     }

    def __str__(self):
        return 'xyz parser'

class GDMA(Parser):

    def __init__(self, filename=None, data=None, here=False):
        Parser.__init__(self, filename=filename, data=data, here=here, dname='gdma', suffix='.out')
        self.read_file(filename=self.filename)

    def read_file(self, filename):
        Parser.read_file(self, filename)
        with open(filename, 'r') as f:
            i=1
            multipoles = {}
            max_rank = 0
            while True:
                line = f.readline()
                if line == '': break
                m = re.search(r'([A-Z][a-z]?) *x = +(-?\d*\.\d*) *y = +(-?\d*\.\d*) *z = +(-?\d*\.\d*)', line)
                if m:
                    # add atom
                    elem =  m.group(1)
                    x, y, z = map(lambda v: float(v), [m.group(2), m.group(3), m.group(4)])
                    crds = np.array([x, y, z])
                    multipoles = {}
                    self.add_atom(i, elem, crds, multipoles=multipoles)
                    i+=1
                m = re.search(r'Maximum rank = +(\d+)', line)
                if m:
                    rank = int(m.group(1))
                    multipoles['rank'] = rank
                    if rank > max_rank:
                        max_rank = rank
                    while True:
                        tmp = f.readline().split()
                        if len(tmp) == 0: 
                            break
                        if tmp[0].startswith('|'):
                            shift = 3
                        else:
                            shift = 0
                        keys = tmp[shift:][::3]
                        values = tmp[shift+2:][::3]
                        for key, value in zip(keys, values):
                            multipoles[key[1:]] = float(value)
                m = re.search(r' *x = +0.0+, *y = +0.0+, *z = +0.0+', line)
                if m:
                    total_multipoles = {}
                    while True:
                        tmp = f.readline().split()
                        if len(tmp) == 0: break
                        if tmp[0].startswith('|'):
                            shift = 3
                        else:
                            shift = 0
                        keys = []
                        values =[]
                        for word in tmp[shift:]:
                            if word == '=':
                                continue
                            elif word.startswith('Q'):
                                keys.append(word)
                            elif word.startswith('='):
                                values.append(float(word[1:]))
                            else:
                                values.append(float(word))
                        for key, value in zip(keys, values):
                            total_multipoles[key[1:]] = value
            total_multipoles['rank'] = max_rank
            self.data = {'atoms': self.atoms,
                         'multipoles': total_multipoles,
                         }
    def __repr__(self):
        return "GDMA parser"
           
            

class ForceFieldXML(object):

    LJ = {
            'CT': (1.9080, 0.1094), # methyl carbon
            'C' : (1.9080, 0.0860), # sp2: CO2, C=O
            'HC': (1.4870, 0.0157), # connected to CT / methyl hydrogen,
            'H' : (0.6000, 0.0157), # connected to N
            'HO': (0.0000, 0.0000), # connected to OH / hydroxyl
            'N' : (1.8240, 0.1700), # sp2
            'N3': (1.8750, 0.1700), # sp3
            'OH': (1.7210, 0.2104),
            'O' : (1.6612, 0.2100), # C=O
            'O2': (1.6612, 0.2100), # CO2
            'S' : (2.0000, 0.2500),
            }


    def __init__(self):
        pass

    def load_forcefields(self, filename=None, here=False, molecule=None):
        if here is False:
            filename = '%s/data/forcefields/%s' % (WORKDIR, filename)
        tree = ElementTree.parse(filename)
        root = tree.getroot()
        for atom in root.findall('atom'):
            element = atom.get('element')
            name = atom.get('name')
            a_charge = float(atom.find('charge').text)
            r0 = float(atom.find('r0').text)
            epsilon = float(atom.find('epsilon').text)
            # get extra points data
            try:
                h = atom.find('hybridization').text
                site = atom.find('site')
                distance = float(site.find('distance').text)
                angle = float(site.find('angle').text)
                e_charge = float(site.find('charge').text)
                extra_exists = True
            except AttributeError:
                extra_exists = False
            # load force fields to atoms
            for a in molecule.get_sites_by_name(name):
                a.set_charge(a_charge)
                a.r0 = r0
                a.epsilon = epsilon
                logger.debug('Load forcefields to %s: charge = %.4f; epsilon = %.4f; r0 = %.4f' % (a.name, a.charge, a.epsilon, a.r0))
                # load force fields to extra points
                if extra_exists:
                    index = molecule.get_max_index()
                    ep = (index+1, h, distance, angle)
                    atom = a.get_attachment()
                    atom.set_hybridization(ep)
                    for s in atom.extra_sites:
                        s.set_charge(e_charge)
                        s.r0 = 0.0
                        s.epsilon = 0.0
                    logger.debug('Force fields for extra points are loaded at %s %s\nNumber of extra points: %i\nCharge: %.3f\nDistance: %.3f\nAngle: %.3f' % (h, atom.name, atom.num_extra_sites, s.charge, distance, angle))
        
    def write_file(self, molecule=None, xmlfilename=None):
        top = Element('forcefield', name=molecule.name)
        comment = Comment('distances in Angstroms, energies in kcal/mol, angles in degrees')
        top.append(comment)
        for atom in molecule.atoms_noneq:
            trunc_name = atom.name.split('-')[0]
            amber_name = amberType[trunc_name]
            r0, e = self.LJ[amber_name]
            atomElem = SubElement(top, 'atom', name=atom.name, element=atom.element, amber=amber_name)
            SubElement(atomElem, 'charge').text = '%.4f' % atom.center.charge
            SubElement(atomElem, 'r0').text = str(r0)
            SubElement(atomElem, 'epsilon').text = str(e)
            try:
                SubElement(atomElem, 'hybridization').text = atom.frame.hybridization
            except AttributeError:
                pass
            # extra points
            ep_data = atom.get_ep_data()
            for site in atom.extra_sites: 
                siteElem = SubElement(atomElem, 'site', name=site.name)
                SubElement(siteElem, 'charge').text = '%.4f' % site.charge
                SubElement(siteElem, 'r0').text = '0.0'
                SubElement(siteElem, 'epsilon').text = '0.0'
                SubElement(siteElem, 'distance').text = '%.4f' % ep_data['distance']
                SubElement(siteElem, 'angle').text = '%.4f' % ep_data['angle']
                break # extra points are symmetrical
        s = prettify(top)
        if xmlfilename is None:
            xmlfilename = '%s/data/forcefields/%s_%s.xml' % (WORKDIR, molecule.name, molecule.theory)
        with open(xmlfilename, 'w') as f:
            f.write(s)

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

amberType = {
        'C_CH3': 'CT',
        'H_CH3': 'HC',
        'O_CO2': 'O2',
        'C_CO2': 'C',
        'N_NH4': 'N3',
        'H_NH4': 'H',
        'N': 'N',
        'S': 'S',
        'O': 'O',
        'H': 'HO',
        'EP_S':'EP',
        'EP_N':'EP',
        'EP_O':'EP',
        }

class Sphere(object):
    
    def __init__(self, data=None):
        filename = '%s/data/sphere/%s_%s_p%ir%.1f.qesp' % (WORKDIR, data['name'], data['theory'], data['grid points'], data['grid radius'])
        self.read_file(filename)

    def read_file(self, filename):
        file = open(filename, 'r')
        self.xyz = []
        self.values = []
        self.atoms = []
        self.property = 'esp'
        for line in file:
            tmp = [float(t) for t in line.split()]
            xyz = np.array(tmp[0:3])
            self.xyz.append(xyz)
            self.values.append(tmp[3])
        file.close()

        num_points = len(self.values)
        self.weights = np.array([np.sqrt(4*np.pi*w) for x, y, z, w in Lebedev(num_points)])

        self.data = {
                    'values': self.values,
                    'grid points': self.xyz,
                    'weights': self.weights,
                }

    def __str__(self):
        return "Sphere parser"



