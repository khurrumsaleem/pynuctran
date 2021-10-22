import decimal as dc
from itertools import product

class smatrix:
    __one__ = dc.Decimal('1.0')
    __zero__ = dc.Decimal('0.0')

    def __init__(self, shape: tuple):
        self.shape = shape
        self.data = {}
        self.common_column = {}
        self.common_row = {}
        return

    def addelement(self, row, col, value):
        self.common_column.setdefault(col,{}).setdefault((row, col), None)
        self.common_row.setdefault(row,{}).setdefault((row, col), None)
        self.data[(row,col)] = self.data.get((row, col), smatrix.__zero__) + value
        return

    @classmethod
    def fromlist(cls, A: list) -> 'smatrix':
        result = cls(shape=(len(A), len(A[0])))

        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                if not A[i][j] == smatrix.__zero__:
                    result.addelement(i, j, A[i][j])
        return result

    def todense(self):
        result = [[smatrix.__zero__ for _ in range(self.shape[0])] for _ in range(self.shape[1])]
        for key in self.data.keys():
            result[key[0]][key[1]] = self.data[key]
        return result


    def __mul__(self, other: 'smatrix'):
        
        result = smatrix(shape=(self.shape[0], self.shape[1]))
        for j in self.common_column.keys():
            L1 = self.common_column.get(j,{})
            L2 = other.common_row.get(j,{})
            comb = list(product(L1,L2))
 
            for comb_pair in comb:
                    result.addelement(comb_pair[0][0], comb_pair[1][1],
                      self.data[comb_pair[0]] * other.data[comb_pair[1]])
        return result

    @classmethod
    def identity(cls, n) -> 'smatrix':
        result = cls((n,n))
        for i in range(n):
            result.addelement(i, i, smatrix.__one__)
        return result


    def __pow__(self, power: int) -> 'smatrix':
        if self.shape[0] != self.shape[1]:
            print('Fatal error: Matrix power of non-square matrix is encountered.')
            return None
        if power == 0:
            return smatrix.identity(self.shape[0])
        tmp = self.__pow__(power//2)
        if (power % 2):
            return self * tmp * tmp
        else:
            return tmp * tmp
