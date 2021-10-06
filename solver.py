import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import decimal as dc
from scipy.sparse import csr_matrix
import xml.etree.ElementTree as ET

'''
    SECTION I: CHEBYSHEV RATIONAL APPRIXIMATION METHOD OF ORDER-16/48........... SEC. I

    This code is obtained from MIT's CRPG group.
    This code is only used for the purpose of PyNUCTRAN's verification. 

    Not to be publicly disclosed to PyNUCTRAN user. 

'''
class cram:

    @staticmethod
    def order16(A, n0, dt):
        """ Chebyshev Rational Approximation Method, order 16
        Algorithm is the 16th order Chebyshev Rational Approximation Method,
        implemented in the more stable incomplete partial fraction (IPF) form
        [cram16]_.
        .. [cram16]
            Pusa, Maria. "Higher-Order Chebyshev Rational Approximation Method and
            Application to Burnup Equations." Nuclear Science and Engineering 182.3
            (2016).
        Parameters
        ----------
        A : scipy.linalg.csr_matrix
            Matrix to take exponent of.
        n0 : numpy.array
            Vector to operate a matrix exponent on.
        dt : float
            Time to integrate to.
        Returns
        -------
        numpy.array
            Results of the matrix exponent.
        """

        alpha = np.array([+2.124853710495224e-16,
                        +5.464930576870210e+3 - 3.797983575308356e+4j,
                        +9.045112476907548e+1 - 1.115537522430261e+3j,
                        +2.344818070467641e+2 - 4.228020157070496e+2j,
                        +9.453304067358312e+1 - 2.951294291446048e+2j,
                        +7.283792954673409e+2 - 1.205646080220011e+5j,
                        +3.648229059594851e+1 - 1.155509621409682e+2j,
                        +2.547321630156819e+1 - 2.639500283021502e+1j,
                        +2.394538338734709e+1 - 5.650522971778156e+0j],
                        dtype=np.complex128)
        theta = np.array([+0.0,
                        +3.509103608414918 + 8.436198985884374j,
                        +5.948152268951177 + 3.587457362018322j,
                        -5.264971343442647 + 16.22022147316793j,
                        +1.419375897185666 + 10.92536348449672j,
                        +6.416177699099435 + 1.194122393370139j,
                        +4.993174737717997 + 5.996881713603942j,
                        -1.413928462488886 + 13.49772569889275j,
                        -10.84391707869699 + 19.27744616718165j],
                        dtype=np.complex128)

        n = A.shape[0]

        alpha0 = 2.124853710495224e-16

        k = 8

        y = np.array(n0, dtype=np.float64)
        for l in range(1, k+1):
            y = 2.0*np.real(alpha[l]*sla.spsolve(A*dt - theta[l]*sp.eye(n), y)) + y

        y *= alpha0
        return y

    @staticmethod
    def order48(A, n0, dt):
        """ Chebyshev Rational Approximation Method, order 48
        Algorithm is the 48th order Chebyshev Rational Approximation Method,
        implemented in the more stable incomplete partial fraction (IPF) form
        [cram48]_.
        .. [cram48]
            Pusa, Maria. "Higher-Order Chebyshev Rational Approximation Method and
            Application to Burnup Equations." Nuclear Science and Engineering 182.3
            (2016).
        Parameters
        ----------
        A : scipy.linalg.csr_matrix
            Matrix to take exponent of.
        n0 : numpy.array
            Vector to operate a matrix exponent on.
        dt : float
            Time to integrate to.
        Returns
        -------
        numpy.array
            Results of the matrix exponent.
        """

        theta_r = np.array([-4.465731934165702e+1, -5.284616241568964e+0,
                            -8.867715667624458e+0, +3.493013124279215e+0,
                            +1.564102508858634e+1, +1.742097597385893e+1,
                            -2.834466755180654e+1, +1.661569367939544e+1,
                            +8.011836167974721e+0, -2.056267541998229e+0,
                            +1.449208170441839e+1, +1.853807176907916e+1,
                            +9.932562704505182e+0, -2.244223871767187e+1,
                            +8.590014121680897e-1, -1.286192925744479e+1,
                            +1.164596909542055e+1, +1.806076684783089e+1,
                            +5.870672154659249e+0, -3.542938819659747e+1,
                            +1.901323489060250e+1, +1.885508331552577e+1,
                            -1.734689708174982e+1, +1.316284237125190e+1])
        theta_i = np.array([+6.233225190695437e+1, +4.057499381311059e+1,
                            +4.325515754166724e+1, +3.281615453173585e+1,
                            +1.558061616372237e+1, +1.076629305714420e+1,
                            +5.492841024648724e+1, +1.316994930024688e+1,
                            +2.780232111309410e+1, +3.794824788914354e+1,
                            +1.799988210051809e+1, +5.974332563100539e+0,
                            +2.532823409972962e+1, +5.179633600312162e+1,
                            +3.536456194294350e+1, +4.600304902833652e+1,
                            +2.287153304140217e+1, +8.368200580099821e+0,
                            +3.029700159040121e+1, +5.834381701800013e+1,
                            +1.194282058271408e+0, +3.583428564427879e+0,
                            +4.883941101108207e+1, +2.042951874827759e+1])
        theta = np.array(theta_r + theta_i * 1j, dtype=np.complex128)

        alpha_r = np.array([+6.387380733878774e+2, +1.909896179065730e+2,
                            +4.236195226571914e+2, +4.645770595258726e+2,
                            +7.765163276752433e+2, +1.907115136768522e+3,
                            +2.909892685603256e+3, +1.944772206620450e+2,
                            +1.382799786972332e+5, +5.628442079602433e+3,
                            +2.151681283794220e+2, +1.324720240514420e+3,
                            +1.617548476343347e+4, +1.112729040439685e+2,
                            +1.074624783191125e+2, +8.835727765158191e+1,
                            +9.354078136054179e+1, +9.418142823531573e+1,
                            +1.040012390717851e+2, +6.861882624343235e+1,
                            +8.766654491283722e+1, +1.056007619389650e+2,
                            +7.738987569039419e+1, +1.041366366475571e+2])
        alpha_i = np.array([-6.743912502859256e+2, -3.973203432721332e+2,
                            -2.041233768918671e+3, -1.652917287299683e+3,
                            -1.783617639907328e+4, -5.887068595142284e+4,
                            -9.953255345514560e+3, -1.427131226068449e+3,
                            -3.256885197214938e+6, -2.924284515884309e+4,
                            -1.121774011188224e+3, -6.370088443140973e+4,
                            -1.008798413156542e+6, -8.837109731680418e+1,
                            -1.457246116408180e+2, -6.388286188419360e+1,
                            -2.195424319460237e+2, -6.719055740098035e+2,
                            -1.693747595553868e+2, -1.177598523430493e+1,
                            -4.596464999363902e+3, -1.738294585524067e+3,
                            -4.311715386228984e+1, -2.777743732451969e+2])
        alpha = np.array(alpha_r + alpha_i * 1j, dtype=np.complex128)
        n = A.shape[0]

        alpha0 = 2.258038182743983e-47

        k = 24

        y = np.array(n0, dtype=np.float64)
        for l in range(k):
            y = 2.0*np.real(alpha[l]*sla.spsolve(A*dt - theta[l]*sp.eye(n), y)) + y

        y *= alpha0
        return y

'''
    SECTION II: PyNUCTRAN SOLVER MODULE........................................ SEC. II
'''
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

        # isotope_names stores the list of species defined by the user.
        self.isotope_names = isotope_names
        # __I__ stores the total number of species defined by the user.
        self.__I__         = len(self.isotope_names)
        # lambdas is a 2D array storing the rates of removal events, indexed by
        # species_id (i) and next by removal event index (j).
        self.lambdas       = [ []    for i in range(self.__I__)]
        # G is a 3D array storing the isotope_id of daughter species of each 
        # removal events. It is indexed by parent's species_id, removal event_id,
        # and lastly the daughter list.
        self.G             = [[[-1]] for i in range(self.__I__)]
        self.P             = [ []    for i in range(self.__I__)]
        self.A             = [ [0.0 for i in range(self.__I__)] for i in range(self.__I__) ]
        self.fission_yield = [ []    for i in range(self.__I__)]

    '''
        Adds a removal event to the solver.

        Parameters:
        isotope_index - The ID of the isotope species based on the species list given
                        given during the initialization of solver class (refer to __init__(...)
                        class constructor.)
        rate          - The rate of the event, for decay, this is equivalent to branching_ratio*decay_rate.
                        For fission, rate is equals to the total fission rate.
        products      - The removal event product(s). For reactions other than fission, only one
                        product is allowed. Here the product is a python list with one element.
                        For fission reactions, the number of products must be >1. Products that
                        are not tracked must be set to -1. For instance, [2,-1], [-1].. etc.
        fission_yield - A list of fission yield. The length of the list must equal to the number of
                        products. 
        
    '''
    def add_removal(self, isotope_index: int, 
                          rate         : float, 
                          products     : list = [-1],
                          fission_yield: list = None):
        d_rate = dc.Decimal(rate)
        i = isotope_index
        self.lambdas[i].append(d_rate)
        self.G[i]      .append(products)

        # If fission_yield is supplied and the product is >= 1. Here, we know that
        # the removal is a fission reaction.
        if not fission_yield is None and len(products) > 1:
            # First we check if the fission yield size is the same with the 
            # number of products.
            if len(fission_yield) == len(products):
                # Update the fission yield table.
                self.fission_yield[isotope_index] = fission_yield

                # Update the transmutation matrix A elements (for CRAM use),
                # accounting the new removal event.
                self.A[isotope_index][isotope_index] -= rate
                for k in range(len(products)):       
                    if not products[k] <= -1:
                        self.A[products[k]][isotope_index] += rate * fission_yield[k]
            else:
                print('Fatal Error: Insufficient fission yields provided for isotope ' + self.isotope_names[isotope_index])
                exit()

        # For non-fission case (decay, (n,2n),(n,3n),(n,a),(n,p))... the case if fission_yield is not supplied.
        # Of course, other 
        elif fission_yield is None and len(products) == 1:
            # Update the transmutation matrix A elements (for CRAM use),
            # accounting the new removal event.
            for product in products:
                self.A[isotope_index][isotope_index] -= rate
                if not product <= -1:
                    self.A[product][isotope_index] += rate
        else:
            print('Fatal Error: Invalid removal definition for isotope ' + self.isotope_names[isotope_index])
            print('Non-fission events must have only one daughter product.')
            print('Fission events must have >1 products to track.')
            exit()
                

    def prepare_transmutation_matrix(self) -> csr_matrix:
        return csr_matrix(np.array(self.A))

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
                no_of_daughters = len(self.G[i][j])
                for l in range(no_of_daughters):
                    # For fission case, we need to multiply the probability with the fission yield.
                    # Sidenote: fission reaction will always have more than one daughters,
                    # len(G[i][j]) > 1
                    if no_of_daughters > 1:
                        k = self.G[i][j][l]
                        if not k < 0:
                            A[k][i] += self.P[i][j] * dc.Decimal(self.fission_yield[i][l])
                            if A[k][i] == dc.Decimal(1.0):
                                sl_positions.append([k,i])
                    
                    # Here we process removal events with only one product, i.e., (n,a),(n,p),(n,2n),
                    # and so on.
                    else:
                        k = self.G[i][j][l]
                        if not k < 0:
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

    def solve(self, n0: np.array, t: float, steps: int, consolidate: bool = False) -> np.ndarray:
        long_n0 = np.transpose(np.matrix([dc.Decimal(x) for x in n0]))
        dt = t / steps
        A = self.prepare_transfer_matrix(dt, consolidate)
        n = np.matmul(np.linalg.matrix_power(A, steps), long_n0)
        return np.array(n, dtype=np.float64)

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

    def plot_concentrations(self, w                : list, 
                                  isotopes_to_plot : list, 
                                  colors           : list = []):

        new_w = np.transpose(np.matrix(w)).tolist()
        names_to_plot =[]
        for i in range(len(isotopes_to_plot)):
            names_to_plot.append(self.isotope_names[i])
            if not colors == []:
                plt.scatter([j for j in range(0,len(new_w[isotopes_to_plot[i]]))],new_w[isotopes_to_plot[i]], color=colors[i], s=4)
            else:
                plt.scatter([j for j in range(0,len(new_w[isotopes_to_plot[i]]))],new_w[isotopes_to_plot[i]], s=4)
#        names_to_plot =[]
#        for i in isotopes_to_plot:
#            names_to_plot.append(self.isotope_names[i])
        plt.xlabel('Sub step')
        plt.ylabel('Isotope Concentrations')
        plt.legend(names_to_plot, loc='upper right')

        return


'''
    SECTION III: DEPLETION DATA PRE-PROCESSING ................................................ SEC. III
'''
class depletion_scheme:

    '''
        Defines the depletion scheme in the code, based on ENDFB17 data. The nuclide data are
        stored in an xml file 'chains_endfb71.xml'. Here, the depletion chains are created
        based on the user specified reaction rates and species.

        Parameters:

        xml_data_location: A string specifying the location of chains_endfb71.xml on the disk.
        species_names    : A python list containing the names of the species involved in the
                           depletion scheme. The species name must be equal to the nuclide ids
                           defined in chains_endfb71.xml. For example, [U235, Pu239, Am241,...].
        rxn_rates        : A 2D python dictionary containing the reaction rates of various
                           removal events. For example,

                           rxn_rates = {
                                'U238' : {'(n,gamma)': 1E-4, 'fission': 1E-5},
                                'Pu239': {'(n,gamma)': 1E-5},
                           }

    '''
    @staticmethod
    def build_chain(xml_data_location: str, species_names: list, solver: solver, rxn_rates):

        t0 = tm.process_time()

        tree = ET.parse(xml_data_location)
        root = tree.getroot()

        for species in root:
            species_name = species.attrib['name']
            if not species_name in species_names:
                continue
            if 'half_life' in species.attrib:
                decay_rate = np.log(2) / float(species.attrib['half_life'])
            else:
                decay_rate = 0.0
            
            removals = list(species)

            for removal in removals:
                if removal.tag == 'decay_type':
                    decay_rate_adjusted = float(removal.attrib['branching_ratio']) * decay_rate
                    parent = species_name
                    daughter  = removal.attrib['target']
                    parent_id = species_names.index(parent)
                    if daughter in species_names:
                        daughter_id = species_names.index(daughter)
                        solver.add_removal(parent_id, decay_rate_adjusted, [daughter_id])
                    else:
                        solver.add_removal(parent_id, decay_rate_adjusted, [-1])
                        # print('Added removal of type decay, parent = %8s\t daughter = %8s\t rate = %.15g.' % (parent, daughter, decay_rate_adjusted))
                
                # If reaction rates are not provided then we skip this.
                if not rxn_rates is None:
                    if species_name in rxn_rates.keys():

                        # Process all absorption reactions, except fission.
                        if removal.tag == 'reaction_type' and 'target' in removal.attrib:
                            parent = species_name
                            parent_id = species_names.index(parent)
                            if removal.attrib['type'] in rxn_rates[parent].keys() and \
                               not removal.attrib['type'] == 'fission':
                                daughter = removal.attrib['target']
                                removal_rate = float(rxn_rates[parent][removal.attrib['type']])
                                if daughter in species_names:
                                    daughter_id = species_names.index(daughter)
                                    solver.add_removal(parent_id, removal_rate, [daughter_id])
                                else:
                                    solver.add_removal(parent_id, removal_rate, [-1])

                        # Process fission reaction.
                        if removal.tag == 'neutron_fission_yields':
                            parent = species_name
                            parent_id = species_names.index(parent)
                            yield_data = list(removal)
                            energy = 0.0
                            products = []
                            yields = []
                            if 'fission' in rxn_rates[parent].keys():
                                for data in yield_data:
                                    if data.tag == 'energies':
                                        energy = sorted([float(e) for e in data.text.split()])[0]
                                    if data.tag == 'fission_yields':
                                        if float(data.attrib['energy']) == energy:
                                            for param in list(data):
                                                if param.tag == 'products':
                                                    products = param.text.split()
                                                if param.tag == 'data':
                                                    yields = [float(y) for y in param.text.split()]
                             
                                total_fission_rate = rxn_rates[parent]['fission']
                                yields_to_add = []
                                daughters_id_to_add = []
                                for product in products:
                                    # For fission, we separate the products into (1) defined-products and
                                    # (2) undefined-products. A defined-products are a list of fission products
                                    # registered in the species_names list. Undefined-products are combined together
                                    # and not tracked, however, the production of these products still affects 
                                    # the parent species concentration. 
                                    if product in species_names:
                                        daughters_id_to_add.append(species_names.index(product))
                                        yields_to_add.append(yields[products.index(product)])
                                parent_id = species_names.index(species_name)
                                solver.add_removal(parent_id, total_fission_rate, daughters_id_to_add, yields_to_add)
                                # Creating a dummy removal events to account for the undefined-products.
                               
        # Report the data processing time.
        t1 = tm.process_time()
        print('Done building chains. CPU time = %.10g secs' % (t1-t0))
        return

    '''
        Gets the list of species available in the nuclides data.
    '''
    @staticmethod
    def get_all_species_names(xml_data_location: str) -> list:
        tree = ET.parse(xml_data_location)
        root = tree.getroot()

        species_names = []
        for species in root:
            species_names.append(species.attrib['name'])
        return species_names
