import numpy as np
import decimal as dc
from scipy.sparse import csr_matrix
import time as tm
import copy
from pynuctran.sparse import *

'''
SECTION II: PyNUCTRAN SOLVER MODULE............................................ SEC. II
A PYTHON LIBRARY FOR NUCLEAR TRANSMUTATION SOLVER (PyNUCTRAN)
License: MIT

Initially developed, designed  and  proposed  by M. R. Omar for the purpose of 
simulating various nuclear transmutations such as decays,  fissions as well as 
neutron  absorptions.  PYNUCTRAN was developed to avoid cumbersome 
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


class solver:

    # shared private constants
    __no_product__        = -1
    __zero__ = dc.Decimal(0.0)
    __one__ = dc.Decimal(1.0)
    __negone__ = dc.Decimal(-1.0)

    def __init__(self, species_names: list):

        # species_names stores the list of species defined by the user.
        # __I__ stores the total number of species defined by the user.
        # lambdas is a 2D array storing the rates of removal events, indexed by
        # species_id (i) and next by removal event index (j).
        # G is a 3D array storing the isotope_id of daughter species of each 
        # removal events. It is indexed by parent's species_id, removal event_id,
        # and lastly the daughter list.
        # P is a 2D array that stores the per-calculated event probabilities, π(i,j).
        # A is a square matrix which is used in CRAM. It is used for PyNUCTRAN's 
        # verification.
        self.species_names = species_names        
        self.__I__         = len(self.species_names)
        self.lambdas       = [ []    for i in range(self.__I__)]
        self.G             = [[[solver.__no_product__]] for i in range(self.__I__)]
        self.P             = [ []    for i in range(self.__I__)]
        self.A             = [ [0.0 for i in range(self.__I__)] for i in range(self.__I__) ]
        self.fission_yields = [ []    for i in range(self.__I__)]

    # *************************************************************************
    # add_removal(...) Adds a removal event to the solver.
    #
    # Parameters:
    # isotope_index - The ID of the isotope species based on the species list 
    #                 given during the initialization of solver class (refer 
    #                 to __init__(...) class constructor.)
    # rate          - The rate of the event, for decay, this is equivalent to 
    #                 branching_ratio*decay_rate. For fission, rate is equals 
    #                 to the total fission rate.
    # products      - The removal event product(s). For reactions other than 
    #                 fission, only one product is allowed. Here the product is 
    #                 a python list with one element. For fission reactions, 
    #                 the number of products must be >1. Products that are not 
    #                 tracked must be set to -1, i.e., [2,-1], [-1].. etc.
    # fission_yield - A list of fission yield. The length of the list must 
    #                 equal to the number of products.
    # *************************************************************************
    def add_removal(self, species_index: int, 
                          rate         : float, 
                          products     : list = [-1],
                          fission_yields: list = None):
        d_rate = dc.Decimal('%g' % rate)
        i = species_index
        self.lambdas[i].append(d_rate)
        self.G[i]      .append(products)

        # If fission_yield is supplied and the product is >= 1. Here, we know 
        # that the removal is a fission reaction.
        if not fission_yields is None and len(products) > 1:
            # First we check if the fission yield size is the same with the 
            # number of products.
            if len(fission_yields) >= len(products):

                # Update the fission yield table. It must be in Decimal, since it will be used
                # by PyNUCTRAN. 
                self.fission_yields[i] = \
                    [dc.Decimal('%g' % y) for y in fission_yields]

                # Update the transmutation matrix A elements (for CRAM use),
                # accounting the new removal event.
                self.A[i][i] -= rate
                for k in range(len(products)):       
                    if not products[k] <= solver.__no_product__:
                        self.A[products[k]][i] += rate * np.longfloat(fission_yields[k])
            else:
                print('Fatal Error: Insufficient fission yields given for species ' \
                      + self.species_names[i] + ' products.')
                exit()

        # For non-fission case (decay, (n,2n),(n,3n),(n,a),(n,p))... the case if fission_yield is not supplied.
        # Of course, other 
        elif fission_yields is None and len(products) == 1:
            # Update the transmutation matrix A elements (for CRAM use),
            # accounting the new removal event.
            for product in products:
                self.A[i][i] -= np.longfloat(rate)
                if not product <= solver.__no_product__:
                    self.A[product][i] += np.longfloat(rate)
        else:
            print('Fatal Error: Invalid removal definition for isotope ' + self.species_names[i])
            print('Non-fission events MUST only have ONE daughter product.')
            print('Fission events must have >1 products to track.')
            exit()
                
    # Prepare the transmutation matrix A for CRAM.
    def prepare_transmutation_matrix(self) -> csr_matrix:
        t0 = tm.process_time()
        l = csr_matrix(np.array(self.A))
        t1 = tm.process_time()
        print('Done building transmutation matrix. Size = %s CPU time = %f.' % (l.shape, t1-t0))
        return l
      
    '''
        ***********************************************************************************
        THIS SUB-SECTION IS THE CORE OF THE PI-DISTRIBUTION METHOD DEPLETION CALCULATION
        IMPLEMENTED IN PyNUCTRAN. 
        ***********************************************************************************
        
        prepare_transfer_matrix(dt) is a function that constructs the transfer matrix,
        based on the provided removal events parameters specified via add_removal(...)
        method. dt=time_step/substeps is the substep interval.
        
        TODO: To further clean-up the code for fast and efficient computation of the 
        transfer matrix.

        Update-1: Of course, understanding the math of preparing the transfer matrix
        is relatively easy and straightforward. Unfortunately, the matrix preparation
        requires high presicion calculation. Even a small binary operation float error
        will affect the accuracy of pi-distribution. Therefore, I preserved high
        presicion calculation for the calculation of pi-distribution. Once the distri-
        bution is computed, it is converted into np.longfloat and the transfer matrix
        is saved using the Compressed Sparse Row (CSR) format.
        
    '''
   
    def prepare_transfer_matrix(self, dt: np.float64, consolidate: bool = False) -> smatrix:
        __zero__ = dc.Decimal('0.0')
        __one__ =   dc.Decimal('1.0')
        __negone__ = dc.Decimal('-1.0')
        sl_positions = []
        # Initialize the sparse matrix.
        A = [ [__zero__ for _ in range(self.__I__)] for _ in range(self.__I__)]
        long_dt = dc.Decimal('%.15g' % dt)

        for i in range(self.__I__):

            n_events = len(self.G[i])
            norm = __zero__

            # Compute the probability of removals... π(i,j).
            E = [(-self.lambdas[i][l-1]*long_dt).exp() for l in range(1,n_events)]
            for j in range(n_events):
                self.P[i].append(__one__)
                for l in range(1, n_events):
                    kron = l == j
                    self.P[i][j] = self.P[i][j] * \
                        (kron + (__negone__)**kron * E[l-1])
                norm = norm + self.P[i][j]

            if norm == __zero__:
                continue
            # Construct the sparse transfer matrix.
            for j in range(n_events):
                self.P[i][j] = self.P[i][j] / norm
               
                n_daughters = len(self.G[i][j])
                for l in range(n_daughters):
                    # For fission case, we need to multiply the probability with the fission yield.
                    # Sidenote: fission reaction will always have more than one daughters,
                    k = self.G[i][j][l]
                    if not k == solver.__no_product__:
                        if n_daughters > 1:
                            A[k][i] += (self.P[i][j] * self.fission_yields[i][l])
                        else:
                            A[k][i] += self.P[i][j]
                        if A[k][i] == __one__ and consolidate:
                            sl_positions.append([k,i])
                # Add a removal event.
                if j == 0:
                    A[i][i] += self.P[i][j]
        if not consolidate:
            return smatrix.fromlist(A)

        # Consolidates short lived species...
        for pos in sl_positions:
            A[pos[0]] = [x + y for x, y in zip(A[pos[0]], copy.deepcopy(A[pos[1]]))]
            A[pos[1]] = [__zero__ for i in range(self.__I__)]
            A[pos[0]][pos[1]] = __zero__
        return smatrix.fromlist(A)

    '''
        ***********************************************************************************
        THIS SUB-SECTION IS THE CORE OF THE PI-DISTRIBUTION METHOD DEPLETION CALCULATION
        IMPLEMENTED IN PyNUCTRAN. 
        ***********************************************************************************
        
        solve(n0, t, steps) returns the species concentrations after t seconds. n0 is the
        initial species concentrations. t is to total time step. substeps is the total
        number of substeps. consolidate=True removes all short-lived species from the trans-
        fer matrix.

    '''
    def solve(self, w0: dict, t: np.float64, substeps: int = 9999999999, consolidate: bool = False) -> dict:
        # Prepare the sparse version of w0 column matrix.
        w0_matrix = [[dc.Decimal('0.0')] for i in range(self.__I__)]
        for i in range(self.__I__):
            if self.species_names[i] in w0.keys():
                w0_matrix[i][0] = dc.Decimal.from_float(w0[self.species_names[i]])
        converted_w = smatrix.fromlist(w0_matrix)

        # Converts all necessary parameters into decimal.
        t_long = dc.Decimal.from_float(t)
        dt = t_long / dc.Decimal.from_float(substeps)

        # Prepare the transfer matrix.
        t0 = tm.process_time()
        A = self.prepare_transfer_matrix(dt, consolidate)
        t1 = tm.process_time()
        print('Done building transfer matrix. Size = %s CPU time = %f secs.' % (A.shape, t1-t0))

        # Compute the matrix power.
        t0 = tm.process_time()
        An = A**(substeps)
        t1 = tm.process_time()
        print('Done computing sparse matrix power. CPU time = %f secs.' % (t1-t0))

        # Sparse multiplication of w(t) = A^n * w(0).
        t0 = tm.process_time()
        w = An * converted_w
        t1 = tm.process_time()
        print('Done computing concentrations. CPU time = %f secs.' % (t1-t0))

        # Convert sparse matrix n into a python dictionary.
        output = {}
        for i in range(self.__I__):
            output[self.species_names[i]] = w.data.get((i,0), solver.__zero__)
        
        return output

