import decimal as dc
from itertools import product
info = \
'''
    This class was initially developed to accomodate fast, high-precision sparse
    matrix multiplications and powers. WARNING! This class does not covers all
    matrix operations, it only cover the basic operations used by PyNUCTRAN, i.e.
    Multiplication and Powers.
    This class uses the basic Python dictionaries to store data. Sparse matrix
    elements are accessed at an incredible speed via the use of hash table.

    SPARSE STORAGE. Only the non-zero elements are stored in smatrix.data dictio-
    nary. The keys of smatrix.data are the tuple specifying the position of the
    element in the dense matrix version. smatrix.common_column and smatrix.common_rows
    are dictionaries that stores the collection (also a dict.) of position tuples
    with common column or row indices, respectively. The keys are the common column/row
    indices.

    SPARSE MULTIPLICATION. Consider sparse matrices A and B. We want to evaluate A*B.
    Firstly, the cartesian products (more or less like a possible combinations) of 
    A.common_column and B.common_row are evaluated for all common index, x. The product
    of these elements are evaluated, and the value of A[i][x] x B[x][j] will contributes
    to AB[i][j]. For example,

    AB[i][j] = A[i][1]*B[1][j] + A[i][2]*B[2][j] + A[i][3]*B[3][j] + ... 

    The above summations can be accumulated using the dictionary's get().
    For a more comprehensive understanding, consider reading the code below. Good luck!

    SPARSE POWER. Suppose we want to evaluate the power of a sparse matrix, i.e. A^n.
    Let n be a large integer number. A naive method is given by,

    A^n = A x A x A x .... (n times)

    Fortunately, this process can be accelerated using the binary decomposition method,
    for instance,

    let C = A x A (power raised to 2)
    C = C x C     (power raised to 4)
    C = C x C     (power raised to 8)
    :
    :
    until... 
    C = C x C     (power raised to n)

    This algorithm has a complexity of O(log n).

    Prepared by M.R.Omar, 22/10/2021.
'''
class smatrix:
    __one__ = dc.Decimal('1.0')
    __zero__ = dc.Decimal('0.0')

    def __init__(self, shape: tuple):
        self.shape = shape
        self.data = {}
        self.common_column = {}
        self.common_row = {}
        return

    # Inserts a non-zero element of the sparse matrix. (row, col) is the position
    # of the element in the dense matrix version.
    def addelement(self, row, col, value):
        self.common_column.setdefault(col,{}).setdefault((row, col), None)
        self.common_row.setdefault(row,{}).setdefault((row, col), None)
        self.data[(row,col)] = self.data.get((row, col), smatrix.__zero__) + value
        return

    # Initializes smatrix from a python list.
    @classmethod
    def fromlist(cls, A: list) -> 'smatrix':
        result = cls(shape=(len(A), len(A[0])))

        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                if not A[i][j] == smatrix.__zero__:
                    result.addelement(i, j, A[i][j])
        return result

    # Converts smatrix into its dense version.
    def todense(self):
        result = [[smatrix.__zero__ for _ in range(self.shape[0])] for _ in range(self.shape[1])]
        for key in self.data.keys():
            result[key[0]][key[1]] = self.data[key]
        return result

    # Overrides the multiplication operator for class smatrix.
    # This method defines the sprase matrix multiplication.
    def __mul__(self, other: 'smatrix'):
        
        result = smatrix(shape=(self.shape[0], self.shape[1]))
        for j in self.common_column.keys():
            L1 = self.common_column.get(j,{})
            L2 = other.common_row.get(j,{})
            # Prepare the cartesian product of L1 an L2
            comb = [(xA, xB) for xA in L1.keys() for xB in L2.keys()]
            for x in comb:
                    result.addelement(x[0][0], x[1][1],
                      self.data[x[0]] * other.data[x[1]])
        return result

    # Creates a sparse IDENTITY matrix.
    @classmethod
    def identity(cls, n) -> 'smatrix':
        result = cls((n,n))
        for i in range(n):
            result.addelement(i, i, smatrix.__one__)
        return result

    # Overrides the matrix power operator **. Implements 
    # the binary decomposition method for matrix power.
    def __pow__(self, power: int) -> 'smatrix':
        c = self
        res = smatrix.identity(self.shape[0])
        while power > 0:
            if power & 1:
                res = res * c
            c = c * c
            power >>= 1
        return res
