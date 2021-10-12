"""
  Full U238 Depletion Problem [Using ENDFB71 nuclides library]
  
  Time step: 2.84E+4 (2 days).
   
"""

from pynuctran.solver import *
import numpy as np
from scipy.sparse.csr import csr_matrix

# Create a list of species involved in the depletion problem.
# It is possible to include all nuclides, but it is memory expensive.
# To include all nuclides replace the definition of isotopes with the
# following:
#
# To include the entire 3000++ nuclides (this requires calculation on a legit machine.
# isotopes = depletion_scheme.get_all_species_names('path\\to\\chains_endfb71.xml')
# 
# To include range of nuclides, with mass number between AMin and AMax:
# isotopes = depletion_scheme.get_all_species_names_range('path\\to\\chains_endfb71.xml', AMin, AMax)
#
#

# Create a list of species involved in the depletion problem.
# It is possible to include all nuclides, but it is memory expensive.
# To include 3821 nuclides:
isotopes = depletion_scheme.get_all_species_names_range('E:\\chain_endfb71.xml',1,400)

# Initialize the PyNUCTRAN solver.
sim = solver(species_names=isotopes)

# Register the rates of neutron-induced reaction events.
# The pre-defined reaction ids are:
# [(n,gamma), (n,p), (n,a), (n,2n), (n,3n), (n,4n), fission]
rxn_rates = {
    'U238': {
        '(n,gamma)' : 1E-4,
        'fission'   : 1E-5
    }
}

# Build the depletion scheme based on the nuclides data stored in chains_endfb71.xml.
depletion_scheme.build_chains(sim, rxn_rates)

# Setup the initial concentration.
w0 = {
    'U238': 1.0
}

# Runs the calculation.
total_time = 2.84E+4
n_final = sim.solve(w0, total_time)



#------------------ OBTAINING REFERENCE SOLUTION [CRAM48]----------------------------
# The CRAM solver was derived from MIT's CRPG codes repository:
# https://github.com/mit-crpg/opendeplete/blob/master/opendeplete/integrator/cram.py
#
# CRAM method is included in this library for PyNUCTRAN's verifications.
#------------------------------------------------------------------------------------

# Prepare the transmutation matrix. In CRAM method, this matrix is the so-called
# matrix A. Recall that CRAM approximates the matrix exponential given in the 
# formula w(t) = exp(At) w0.
A = sim.prepare_transmutation_matrix()
w0_matrix = [np.float64('0.0') for i in range(len(isotopes))]
for key in w0.keys():
    w0_matrix[isotopes.index(key)] = np.float64(w0[key])
n0 = np.transpose(np.array(w0_matrix))
n_final_cram = cram.order48(A,n0,total_time)

# Prints the output of PyNUCTRAN solver and CRAM48, as well as their relative error.
print('%-5s   %-5s   %-21s   %-21s   %-21s' % ('ID', 'Name','Calculated','Reference (CRAM)', 'Rel. Error'))
for i in range(len(isotopes)):
    if n_final[i] > 1E-15:
        print('%i\t%s\t%+20.14e\t%+20.14e\t%+20.14e' % (i, isotopes[i],n_final[i],n_final_cram[i],\
                ((n_final[i])-n_final_cram[i])/n_final_cram[i]))
