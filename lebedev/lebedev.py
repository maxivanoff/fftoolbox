    def get_basis(self, nodes=None, l=None):
        N = len(nodes)
        B = np.zeros((l*l, N))
        i = 0
        for ll in range(l):
            for mm in range(-ll,ll+1):
                if mm < 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Ylms(ll, abs(mm), p.theta, p.phi).real * np.sqrt(p.w)
                if mm == 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Y(0, ll, p.theta, p.phi).real * np.sqrt(p.w)
                if mm > 0:
                    for j, p in enumerate(nodes):
                        B[i, j] = Ylmc(ll, mm, p.theta, p.phi).real * np.sqrt(p.w)
                i += 1
        return B

