'''
A Python Library for Pseudo-MC Nuclear Transmutation Solver (PyNUCTRAN)
License: MIT

Initially developed, designed and proposed by Dr M. R. Omar for the purpose of 
simulating various nuclear transmutations such as decays,  fissions as well as 
neutron  absorptions.  PYNUCTRAN was  initially  developed to avoid cumbersome 
numerical issues of solving the nuclide depletion equations.

This code does not directly solve Bateman's  equations.  Instead, it simulates
the actual transmutation processes via Poisson statistics.  However, this code 
does  not  incorporate  pseudo-random number  generations,  thus,  making it a
"Quasi-Monte Carlo" technique. 

The algorithm used is embeded in solver.advance_time_step() method. The pseudo
code of the algorithm is given as follows:
                
define p(i,l) = product(j=1 to J_i) d(j,l) + (-1)**d(j,l) * exp(-lambda(j)*dt)
**l=0 indicates that p(i,0) is the probability of no removal happens.

define w as an array consists of the current weight of all isotopes.
define I as the total number of isotopes
define e as the precision level, say, 1E-20.
define J(i) as the total number of removals of isotope-i.

PSEUDO-CODE...................................... COMMENTS....................

routine advance_time_steps(w)
    input w                                       current weight 
    for i ... I                                   loop over I isotopes 
        w[i] < e ? continue for
        for j ... J[i]                            loop over J removals 
            for k in K(i,j)                       loop over K products 
            w[k] = w[k] + w[i] * P(i,j)
        w[i] = w[i] * P(i,0)
    return w                                      return advanced weight)

end
..............................................................................

To iterate over time steps, the following procedure is recommended:

PSEUDO-CODE...................................... COMMENTS ...................

w = [w0_1, w0_2, ..., w0I]                        assign initial weights
for t ... total_steps                             loop over time steps
   w = advance_time_steps(w)
   
end

..............................................................................


Created on 3-09-21.
(c) M. R. Omar, School of Physics, Universiti Sains Malaysia, 11800 Penang, MY.

'''

import numpy             as np
import matplotlib.pyplot as plt
import copy
import matplotlib.pyplot as plt
import time              as tm

class solver:

    # shared private constants
    __const_no_product__        = -1
    __const_default_steps__     = 100
    __const_default_precision__ = 1E-20
    __const_default_dt__        = 1.0

    def __init__(self, isotope_names: list, 
                       time_interval: float = __const_default_dt__,
                       steps        : int   = __const_default_steps__,
                       precision    : float = __const_default_precision__):

        self.dt            = time_interval
        self.steps         = steps
        self.isotope_names = isotope_names
        self.__I__         = len(self.isotope_names)
        self.lambdas       = [ []    for i in range(self.__I__)]
        self.G             = [[[-1]] for i in range(self.__I__)]
        self.__J__         = [  0    for i in range(self.__I__)]
        self.__K__         = [ []    for i in range(self.__I__)]
        self.P             = [ []    for i in range(self.__I__)]
        self.presicion = precision
        
    def add_removal(self, isotope_index: int, 
                          rate         : float, 
                          products     : list = [-1]):
        i = isotope_index
        self.lambdas[i].append(rate)
        self.G[i]      .append(products)
        self.__K__[i]  .append(len(products))

    def preprocess_p_table(self):

        for i in range(self.__I__):

            n_removals = len(self.lambdas[i])
            # Store J values... add by one to include the no-removal at the 
            # first element.
            self.__J__[i] = n_removals + 1

            # Compute the probability of no-removal.
            self.P[i].append(1.0)
            for l in range(len(self.lambdas[i])):
                self.P[i][0] *= np.exp(-self.lambdas[i][l]*self.dt)

            # Compute the probability of removals.
            for j in range(n_removals):
                self.P[i].append(1.0)
                for l in range(n_removals):
                    kronecker = float(l == j)
                    self.P[i][j+1] *= kronecker + (-1)**kronecker * \
                                      np.exp(-self.lambdas[i][l] * self.dt)


    def advance_time_step(self, w: list):
        # Loop over the isotopes...
        for i in range(self.__I__):

            if w[i] < self.presicion:
                continue
            
            # Loop over the removals...
            for j in range(1, self.__J__[i]):
                # Loop over the daughters...
                for k in range(self.__K__[i][j-1]):
                    # If the number of 
                    if self.G[i][j][0] != solver.__const_no_product__:
                        w[self.G[i][j][k]] += w[i] * self.P[i][j]
            w[i] *= self.P[i][0]
        return
    
    def run(self, w: list):
        # record the start time.
        start_time = tm.process_time()

        w_vs_steps = []
        w0 = copy.deepcopy(w)
        w_vs_steps.append(w0)

        for t in range(self.steps):
            self.advance_time_step(w)
            w_vs_steps.append(copy.deepcopy(w))

        end_time = tm.process_time()
        return w_vs_steps, end_time-start_time


    def plot_concentrations(self, w                : list, 
                                  isotopes_to_plot : list, 
                                  colors           : list = []):

        new_w = np.transpose(np.matrix(w)).tolist()

        for i in range(len(isotopes_to_plot)):
            if not colors == []:
                plt.plot(new_w[isotopes_to_plot[i]], colors[i])
            else:
                plt.plot(new_w[isotopes_to_plot[i]])
        names_to_plot =[]
        for i in isotopes_to_plot:
            names_to_plot.append(self.isotope_names[i])
        plt.xlabel('Step')
        plt.ylabel('Isotope Concentrations')
        plt.legend(names_to_plot, loc='upper right')
        plt.show()

        return
    
