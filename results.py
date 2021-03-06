from collections import defaultdict
import units

class Reports(object):

    def __init__(self):
        self.data = defaultdict(list)
        self.keys = self.data.keys()
        self.n = 0

    def add(self, charges):
        #self.data['Name'].append(charges.name)
        self.data['RMSD'.ljust(10, ' ')].append(charges.rmsd)
        self.data['RMAE'.ljust(10, ' ')].append(charges.rmad)
        self.data['R2'.ljust(10, ' ')].append(charges.R2)
        self.data['Max. Error'.ljust(10, ' ')].append(charges.max_error)
        a, da, b, db = charges.ab
        self.data['alpha'.ljust(10, ' ')].append(a)
        self.n += 1

    def add_whatever(self, key, value):
        self.data[key].append(value)

    def __str__(self):
        num = len(self.data.keys())
        s = ''
        j = 0
        for key in self.data.keys():
            j += 1
            if not j == num:
                s += '%s & ' % key
            else:
                s += '%s\\\\ \n\\hline \n' % key
        for i in xrange(self.n):
            j = 0
            for key, values in self.data.items():
                j += 1
                if key == 'name':
                    s += '%s & ' % values[i].ljust(20, ' ')
                    continue
                try:
                    value_str = '%.3f' % values[i]
                except TypeError:
                    value_str = '%s' % values[i]
                value_str = value_str.ljust(10, ' ')
                if not j == num:
                    s += '%s & ' % value_str
                else:
                    s += '%s' % value_str
            s += '\\\\ \n'
        return s


class Report(object):

    def __init__(self, name=None, charges=None):
        self.name = name
        self.rmsd = charges.rmsd
        self.max_error = charges.max_error
        self.rmad = charges.rmad
        self.R2 = charges.R2
        self.a, self.da, self.b, self.db = charges.ab
        self.data = {
                'rmsd': self.rmsd,
                'max error': self.max_error,
                'rmad': self.rmad,
                'R2': self.R2,
                'alpha': self.a,
                }

    def __repr__(self):
        s = '%s\nrmsd = %.3f kcal/mol\nmax error = %.3f kcal/mol\nrmad = %.3f\ny = (%.3f +/- %.3f) * b + %.3f +/- %.3f\nR2 = %.3f\n' \
                % (self.name, self.rmsd, self.max_error, self.rmad, self.a, self.da, self.b, self.db, self.R2)
        return s

class Results(object):

    def __init__(self, ref_multipoles=None):
        self.charges = []
        self.multipoles = []
        self.rmsd = []
        if ref_multipoles:
            self.add_reference(ref_multipoles)
        self.extra = []

    def add_reference(self, ref_multipoles):
        self.multipoles.append(ref_multipoles)
    
    def add(self, o, name=None, qnames=None):
        try:
            rmsd = o.rmsd
        except:
            rmsd = 0.
        charges = o.charges.copy()
        charges['RMSD'] = rmsd
        if name is None:
            charges['name'] = o.molecule.name
        else:
            charges['name'] = name
        multipoles = o.multipoles.copy()
        multipoles['RMSD'] = rmsd
        if name is None:
            multipoles['name'] = o.molecule.name
        else:
            multipoles['name'] = name
        if len(self.charges) == 0:
            self.qnames = charges.keys()
            self.mnames = sorted(multipoles.keys())
        if qnames:
            self.qnames = qnames
        self.charges.append(charges)
        self.multipoles.append(multipoles)
        self.rmsd.append(rmsd)

    def add_extra(self, ep, name):
        extra = {}
        extra['a name'] = name
        for atom, geom in ep.items():
            extra['%s-d' % atom] = geom['distance']
            extra['%s-a' % atom] = geom['angle']
        self.extra.append(extra)

    def tabular(self, names, dicts):
        a = ''
        for name in names:
            a += '%s' % name
            for d in dicts:
                try:
                    a += ' %.3f' % d[name]
                except KeyError:
                    a += ' 0.0'
            a += '\n'
        return a

    def latex(self, names, dicts):
        s = ''
        rec = ''
        for d in dicts:
            title = ''
            #for n in sorted(names):
            for n in names:
                nn = n
                for _ in range(10-len(n)): nn += ' '
                title += '& %s ' % nn
            title += '\\\\ \n'
            if not title == rec:
                s += title
                rec = title
            #for n in sorted(names):
            for n in names:
                try:
                    ss = '%.4f' % d[n]
                except KeyError:
                    ss = '--'
                except TypeError:
                    ss = '%.10s' % d[n]
                for _ in range(10-len(ss)): ss += ' '
                s += '& %s ' % ss
            s += '\\\\ \n'
        return s

    def __str__(self):
        start = '\\begin{table}\n\\begin{tabular}{ccccccccccccccc}\n'
        end = '\\end{tabular}\n\\caption{PBE/def2-SV(P)+d}\n\end{table}\n'
        s = ''
        s += 'charges:\n'
        #s += start
        s += self.latex(self.qnames, self.charges)
        #s += end
        #s += start
        s += 'multipoles:\n'
        s += self.latex(self.mnames, self.multipoles)
        #s += end
        if self.extra:
            s += 'extra points:\n'
            #s += start
            s += self.latex(sorted(self.extra[0].keys()), self.extra)
            #s += end
        return s

