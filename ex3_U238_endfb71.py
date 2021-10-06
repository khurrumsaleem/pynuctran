"""
  U-238 Depletion Problem [Using ENDFB71 nuclides library]
  
  Time step: 2.84E+4 (2 days).
  No. of sub-steps, t/dt: 1E+20.
   
"""

from pynuctran.solver import *
import numpy as np
from scipy.sparse.csr import csr_matrix

# Create a list of species involved in the depletion problem.
# It is possible to include all nuclides, but it is memory expensive.
# To include all nuclides replace the definition of isotopes with the
# following:
#
# isotopes = depletion_scheme.get_all_species_names('path\\to\\chains_endfb71.xml')
isotopes = [ 'U238','U239','Np239',
             'Pu239','Pu240','Pu241',
             'Pu242','Pu243','Am241',
             'Am243','Am244','Cm244'
           ]

# Initialize the PyNUCTRAN solver.
sim = solver(isotope_names = isotopes)

# Register the rates of neutron-induced reaction events.
# The pre-defined reaction ids are:
# [(n,gamma), (n,p), (n,a), (n,2n), (n,3n), (n,4n), fission]
rxn_rates = {
    'U238' : {'(n,gamma)': 1E-4},
    'Pu239': {'(n,gamma)': 1E-4},
    'Pu240': {'(n,gamma)': 1E-4},
    'Pu241': {'(n,gamma)': 1E-4},
    'Pu242': {'(n,gamma)': 1E-4},
    'Am243': {'(n,gamma)': 1E-4}
}

# Build the depletion scheme based on the nuclides data stored in chains_endfb71.xml.
# In this example, the XML data is located at E:\chain_endfb71.xml. You may want to
# replace this file location.
depletion_scheme.build_chain('E:\\chain_endfb71.xml',isotopes,sim,rxn_rates)

# Setup the initial concentration. The size of w0 must be equal to the number of 
# species defined in the depletion scheme.
w0 = [  1.00E+10, 1.00E+3, 0.00E+00, 
        0.00E+00, 0.00E+00, 0.00E+00,
        0.00E+00, 0.00E+00, 0.00E+00,
        0.00E+00, 0.00E+00, 0.00E+00
     ]
 
# Runs the calculation.
total_time = 8.64E+4
steps = int(1E15)
n_final = sim.solve(w0,total_time,steps)


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

n0 = np.transpose(np.array(w0))
n_final_cram = cram.order48(A,n0,total_time)

# Prints the output of PyNUCTRAN solver and CRAM48, as well as their relative error.
print('%21s---%21s---%21s' % ('-'*21,'-'*21,'-'*21))
print('%21s   %21s   %21s' % ('Calculated','Reference (CRAM)', 'Rel. Error'))
print('%21s---%21s---%21s' % ('-'*21,'-'*21,'-'*21))
for i in range(len(isotopes)):
    print('%+20.14e\t%+20.14e\t%+20.14e' % (n_final[i],n_final_cram[i],\
            (float(n_final[i])-n_final_cram[i])/n_final_cram[i]))
