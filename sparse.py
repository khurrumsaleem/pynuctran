import decimal as dc
from collections import defaultdict
import time as tm
import psutil
from functools import lru_cache, partial
from multiprocessing.pool import ThreadPool as Pool

'''
    This class was initially developed to accomodate fast, high-precision sparse
    matrix multiplications and powers. WARNING! This class does not covers all
    matrix operations, it only covers the basic operations used by PyNUCTRAN, i.e.
    Multiplication and Powers.
   
    This class uses the basic Python dictionaries to store data. Sparse matrix
    elements are accessed at an incredible speed via the use of hash table.

    Of course, the existing SCIPY library can be used, however, it does not allow
    dtype=Decimal, and this is frustrating. Therefore, we must endure writing a new
    specialized class handling sparse matrix power to preserve the accuracy.

    SPARSE STORAGE. Only the non-zero elements are stored in smatrix.data dictio-
    nary. A modified compressed sparse row format is used. Here, rowptr is not used,
    instead, the non-zero values are kept in a 2D dictionary.

    SPARSE MULTIPLICATION. Consider sparse matrices A and B. We want to evaluate A*B.
    Here we implement the row-wise multiplication algorithm. Each row can be vectorized
    into multiple concurrent threads. This sparse multiplication technique is originated
    from Gustafson (1978): 
    
    Fred G. Gustavson. 1978. Two Fast Algorithms for Sparse Matrices: Multiplication and
    Permuted Transposition. ACM Trans. Math. Softw. 4, 3 (Sept. 1978), 250–269. 
    DOI:https://doi.org/10.1145/355791.355796

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
# Evaluates the cartesian product of self.common_columns and other.common_rows.

def row_operation(irow, sd, od, rd):
    rdr, sdd = rd[irow], sd[irow]
    for icol in sdd:
        odd, x = od[icol], sdd[icol]
        for ocol in odd:
            rdr[ocol] += x * odd[ocol]

class smatrix:
    __one__ = dc.Decimal('1.0')
    __zero__ = dc.Decimal('0.0')
    __pool_obj__ = Pool(len(psutil.Process().cpu_affinity()))
    def __init__(self, shape: tuple):


        self.shape = shape
        self.data = defaultdict(lambda: defaultdict(lambda: dc.Decimal('0.0')))
        
        return

    # Initializes smatrix from a python list.
    @classmethod
    def fromlist(cls, A: list) -> 'smatrix':
        result = cls(shape=(len(A), len(A[0])))

        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                a = A[i][j]
                if a != smatrix.__zero__:
                    result.data[i][j] = a

        return result


    # Overrides the multiplication operator for class smatrix.
    # This method defines the sprase matrix multiplication.

    def __mul__(self, other: 'smatrix'):
        sx = self.shape[0]
        sy = other.shape[1]
        result = smatrix(shape=(sx, sy))
        sd, od, rd = self.data, other.data, result.data


        #ncpu = len(psutil.Process().cpu_affinity())
        
        smatrix.__pool_obj__.map(partial(row_operation, rd=rd, sd=sd, od=od), range(sx))
        return result

    # Overrides the matrix power operator **. Implements 
    # the binary decomposition method for matrix power.
    @lru_cache(maxsize = 200)
    def bpow(self, bit: int):
        result = self
        if (bit == 0):
            return result
        result = self.bpow(bit - 1)
        return result * result
