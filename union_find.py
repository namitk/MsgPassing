class union_find:
    """ 
    A simple Union find implementation : Adopted from the book 
    'Python Algorithms - Mastering Basic Algorithms in the Python Language' by Magnus Lie Hetland
    """
    def find(self, C, u):
        if C[u] != u:
            C[u] = self.find(C, C[u])
        return C[u]

    def union(self, C, R, u, v):
        u, v = self.find(C, u), self.find(C, v)
        if R[u] > R[v]:
            C[v] = u
        else:
            C[u] = v
        if R[u] == R[v]:
            R[v] += 1
