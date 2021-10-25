import decimal as dc
from collections import defaultdict
from itertools import product


'''
    This class was initially developed to accomodate fast, high-precision sparse
    matrix multiplications and powers. WARNING! This class does not covers all
    matrix operations, it only cover the basic operations used by PyNUCTRAN, i.e.
    Multiplication and Powers.
   
    This class uses the basic Python dictionaries to store data. Sparse matrix
    elements are accessed at an incredible speed via the use of hash table.

    Of course, the existing SCIPY library can be used, however, it does not allow
    dtype=Decimal, and this is frustrating. Therefore, I must endure writing a new
    specialized class handling sparse matrix power to preserve the accuracy.

    SPARSE STORAGE. Only the non-zero elements are stored in smatrix.data dictio-
    nary. The keys of smatrix.data are the tuple specifying the position of the
    non-zero elements in the dense matrix version. smatrix.common_column (cc) and 
    smatrix.common_rows (cr) are dictionaries that stores the collection (also a dict.)
    of position tuples with common column or row indices, respectively. The keys are
    the common column/row indices.

    SPARSE MULTIPLICATION. Consider sparse matrices A and B. We want to evaluate A*B.
    Firstly, the cartesian products (more or less like a possible combinations) of 
    A.common_column and B.common_row are evaluated for all common index, x. The product
    of these elements are evaluated, and the value of A[i][x] x B[x][j] will contributes
    to AB[i][j]. For example,

    AB[i][j] = A[i][1]*B[1][j] + A[i][2]*B[2][j] + A[i][3]*B[3][j] + ... 

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
        self.data = defaultdict(lambda: smatrix.__zero__)
        self.common_column = defaultdict(defaultdict)
        self.common_row = defaultdict(defaultdict)

        return

    # Initializes smatrix from a python list.
    @classmethod
    def fromlist(cls, A: list) -> 'smatrix':
        result = cls(shape=(len(A), len(A[0])))

        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                if not A[i][j] == smatrix.__zero__:
                    #result.addelement(i, j, A[i][j])
                    result.common_column[j][(i,j)] = None
                    result.common_row[i][(i,j)] = None
                    result.data[(i,j)] = A[i][j]
        return result

    # Overrides the multiplication operator for class smatrix.
    # This method defines the sprase matrix multiplication.
    def __mul__(self, other: 'smatrix'):
        s = self
        result = smatrix(shape=(s.shape[0], s.shape[1]))
        scc, ocr, rcc, rcr, sd, od, rd = self.common_column, other.common_row, \
                                    result.common_column, result.common_row, \
                                    self.data, other.data, result.data

        for j in scc:
            # Prepare the cartesian product of L1 an L2
            #comb = [(xA, xB) for xA in scc[j] for xB in ocr[j]]
            comb = product(scc[j], ocr[j])

            for x in comb:
                rd[(x[0][0],x[1][1])] += sd[x[0]] * od[x[1]]

        for x in rd:
            rcc[x[1]][x], rcr[x[0]][x] = None, None

        return result

    # Overrides the matrix power operator **. Implements 
    # the binary decomposition method for matrix power.
    def __pow__(self, power: int) -> 'smatrix':
        c = self
        n = c.shape[0]
        result = smatrix((n,n))
        # Prepare identity matrix.
        for x in range(n):
            i, j = x, x
            result.common_column[j][(i,j)] = None
            result.common_row[i][(i,j)] = None
            result.data[(i,j)] += smatrix.__one__

        while power > 0:
            if power & 1:
                result = result * c
            c = c * c
            power >>= 1
        return result
