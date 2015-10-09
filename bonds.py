

class Bond(object):

    def __init__(self, a1, a2):
        self.name = '%s--%s' % (a1.name, a2.name) 
        self.a1 = a1
        self.a2 = a2
        self.atoms = [a1, a2]

    def __eq__(self, b):
        if b.a1 in self.atoms and b.a2 in self.atoms:
            return True
        else:
            return False

    def __repr__(self):
        return self.name

