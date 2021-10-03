'''
...............................................................................
A PYTHON LIBRARY FOR NUCLEAR TRANSMUTATION SOLVER (PyNUCTRAN)
License: MIT

Initially developed, designed  and  proposed  by M. R. Omar for the purpose of 
simulating various nuclear transmutations such as decays,  fissions as well as 
neutron  absorptions.  PYNUCTRAN was  initially  developed to avoid cumbersome 
numerical issues of solving the nuclide depletion equations.

This code does not directly solve  Bateman's  equations.  Instead, it uses the 
pi-distribution to  estimate  the  evolution  of  species  concentrations in a 
nuclide depletion problem. The pi-distribution is given by

    pi(i,l) = c * product(j=1 to J_i) d(j,l) + (-1)**d(j,l) * exp(-rate(j)*dt)

c is the normalization factor of the distribution.
pi(i,0) is the probability of no removal happens.
rate(j) is the rate of transmutation event-j.
d(j,l) is the kronecker delta.
dt is the time substep interval.


define w as an array consists of the current weight of all isotopes.
define I as the total number of isotopes
define J(i) as the total number of transmutation events for isotope-i.

The calculation is based on the following iteration,

    w(t) = A^(t/dt) w(0)

where w(0) is the initial concentration of all species, A is the transfer matrix 
which is defined as follows:
        _                      _
        |   p(1->1) ... p(I->1)  |
    A = |     :     '.     :     |
        |_  p(1->I) ... p(I->I) _|

and p(i->k) is the transfer probability which can be derived using
pi-distribution using the following formula,

    p(k->i) = sum of pi(k,j) for all events j that mutates 
              species k into i.

note also that matrix A  is a square matrix (IxI) with its columns as the parent
species and rows as the daughter species. Also, w and w(0) are Lx1 column matrix.

.................................................................................
Created on 3-10-21.
(c) M. R. Omar, School of Physics, Universiti Sains Malaysia, 11800 Penang, MY.

'''

import numpy             as np
import matplotlib.pyplot as plt
import copy
import matplotlib.pyplot as plt
import time              as tm
import decimal           as dc
class solver:

    # shared private constants
    __const_no_product__        = -1
    __const_default_steps__     = 100

    def __init__(self, isotope_names: list):

        self.isotope_names = isotope_names
        self.__I__         = len(self.isotope_names)
        self.lambdas       = [ []    for i in range(self.__I__)]
        self.G             = [[[-1]] for i in range(self.__I__)]
        self.P             = [ []    for i in range(self.__I__)]
        
    def add_removal(self, isotope_index: int, 
                          rate         : float, 
                          products     : list = [-1]):
        d_rate = dc.Decimal(rate)
        i = isotope_index
        self.lambdas[i].append(d_rate)
        self.G[i]      .append(products)

    def prepare_transfer_matrix(self, dt: float, consolidate: bool = False) -> np.ndarray:

        # A list of tuples storing the position of short lived species in the matrix.
        sl_positions = []

        # Initialize the sparse matrix.
        A = [ [dc.Decimal(0.0) for i in range(self.__I__)] for i in range(self.__I__)]
        for i in range(self.__I__):

            # Total number of events. This includes the no-removal event.
            n_events = len(self.G[i])

            # Initialize the normalization factor, c. Norm is the sum of all probs.
            norm = dc.Decimal(0.0)

            # Compute the probability of removals.
            for j in range(n_events):
                self.P[i].append(dc.Decimal(1.0))
                for l in range(1, n_events):
                    kronecker = dc.Decimal(l == j)
                    self.P[i][j] = self.P[i][j] * (kronecker + (dc.Decimal(-1))**kronecker * \
                                      np.exp(-self.lambdas[i][l-1] * dc.Decimal(dt)))
                norm = norm + self.P[i][j]

            # Construct the sparse transfer matrix.
            for j in range(n_events):
                self.P[i][j] = self.P[i][j] / norm
                for k in self.G[i][j]:
                    if not k == -1:
                        A[k][i] += self.P[i][j]
                        if A[k][i] == dc.Decimal(1.0):
                            sl_positions.append([k,i])
                if j == 0:
                        A[i][i] += self.P[i][j]

        if not consolidate:
            return np.matrix(A)

        # Consolidates short lived species...
        for pos in sl_positions:
            A[pos[0]] = [x + y for x, y in zip(A[pos[0]], copy.deepcopy(A[pos[1]]))]
            A[pos[1]] = [dc.Decimal(0.0) for i in range(self.__I__)]
            A[pos[0]][pos[1]] = dc.Decimal(0.0)
        return np.matrix(A)

    def solve(self, n0, t: float, steps: int, consolidate: bool = False) -> np.ndarray:
        long_n0 = np.transpose(np.matrix([dc.Decimal(x) for x in n0]))
        dt = t / steps
        A = self.prepare_transfer_matrix(dt, consolidate)
        return np.matmul(np.linalg.matrix_power(A, steps), long_n0)

############ OBSOLETE METHOD. THE DIRECT SIMULATION METHOD. #############
#    def process(self, k, r_w, dw):
#        for l in range(1, len(self.G[k])):
#            if self.P[k][l] == dc.Decimal(1.0):
#                for g in self.G[k][l]:
#                    self.process(g, r_w, dw)
#                return     
#        r_w[k] = r_w[k] + dw
#        return
#
#    def run(self, w_in: list):
#
#        w_vs_steps = []
#
#        # Prepare the initial weights forall isotopes.
#        w0 = [dc.Decimal(w_in[i]) for i in range(len(w_in))]
#        w_vs_steps.append(copy.deepcopy(w0))
#
#
#        # Set the first timestep w to w0
#        w = copy.deepcopy(w0)
#        for t in range(self.steps):
#            
#            r_w = [dc.Decimal(0.0) for i in range(self.__I__)]
#
#            for i in range(self.__I__):
#                for j in range(len(self.G[i])):
#                    dw = self.P[i][j] * w[i]
#                    for g in self.G[i][j]:
#                        if not g == self.__const_no_product__:
#                            self.process(g, r_w, dw)
#                r_w[i] = r_w[i] + self.P[i][0] * w[i]
#
#            # Record isotope weights data...
#            w_vs_steps.append(copy.deepcopy(r_w))
#
#            # update w...
#            w = copy.deepcopy(r_w)
#
#
#        end_time = tm.process_time()
#        return w_vs_steps
#
######### OBSOLETE METHOD. THE DIRECT SIMULATION METHOD. ##############
