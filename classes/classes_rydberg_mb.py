
class AtomicNLJMIterator:

    def __init__(self, nmin:int, nmax:int, lmax:int, iter_mj=True):
        """
        An iterator through valid atomic states,
        returned as a tuple:
        (n,l,j,m_j)

        """
        self.nmin = nmin
        self.nmax = nmax
        self.lmax = lmax
        self.iter_mj = iter_mj
        # Note: nljm starts at unphysical value due to
        #  iterator incrementing it prior to first output
        self.nljm = [nmin, 0, .5, -1.5]
        self.end = [nmax, lmax, lmax+.5, lmax+.5]

    def __iter__(self):
        return self  # The iterator object itself

    def __next__(self):
        if self.nljm == self.end:
            raise StopIteration  # End of iteration
        
        if self.nljm[3] < self.nljm[2] : self.nljm[3] += 1.
        elif self.nljm[2] < self.nljm[1] + .5:
            self.nljm[2] += 1.
            self.nljm[3] = -self.nljm[2]
        elif self.nljm[1] < self.lmax and self.nljm[1] < self.nljm[0] - 1:
            self.nljm[1] += 1
            self.nljm[2] = self.nljm[1] - .5
            self.nljm[3] = -self.nljm[2]
        elif self.nljm[0] < self.nmax:
            self.nljm[0] += 1
            self.nljm[1] = 0
            self.nljm[2] = .5
            self.nljm[3] = -.5

        if not self.iter_mj : self.nljm[3] = self.nljm[2]
        
        if self.iter_mj : return tuple(self.nljm)
        else : return tuple(self.nljm[:3])
