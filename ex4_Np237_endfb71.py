from pynuctran.solver import *
import numpy as np
from scipy.sparse.csr import csr_matrix


isotopes = [
    
    'Np237',
    'Pa233',
    'U233',
    'Th229',
    'Ra225',
    'Ac225',
    'Fr221',
    'At217',
    'Bi213',
    'Po213',
    'Tl209',
    'Pb209',
    'Bi209'
]


# Initialize the PyNUCTRAN solver.
sim = solver(species_names = isotopes)

# Build the depletion scheme based on the nuclides data stored in chains_endfb71.xml.
# In this example, the XML data is located at E:\chain_endfb71.xml. You may want to
# replace this file location.
depletion_scheme.build_chains(sim, {}, 'E:\\chain_endfb71.xml')
 
w0 = [0.0 for _ in range(len(isotopes))]
w0[0] = 1.0

# Runs the calculation.
total_time = 3.1556926E+13
steps = int(1E20)
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
