class Energy(object):

    def __init__(self, terms=['electrostatic', 'LJ']):
        self.terms = list()
        for name in terms:
            if name == 'electrostatic':
                self.terms.append(Electrostatic(forcefields))
            if name == 'LJ':
                self.terms.append(LennardJones(forcefields))
    
    def between_sites(self, site1, site2):
        energy = 0.
        for term in self.terms:
            energy += term.get_energy(site1, site2)
        return energy

    def between_molecules(self, molecules):
        mol1, mol2 = molecules
        energy = 0.
        for site1 in mol1.sites:
            for site2 in mol2.sites:
                energy += self.between_sites(site1, site2)
        return energy

