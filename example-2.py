"""
  U-238 Depletion Problem with fission.
  
  U-238 Fission yields:
  Xe135 = 1.11541594E-04
  I135  = 0.013157
  
  U-238 Fission rate   : 1E-5 per secs.
  U-238 Absorption rate: 1E-4 per secs.
  Pu239 Absorption rate: 1E-5 per secs.
  
  Time step = 1E+5 seconds.
  No. of substeps = 1E+20.
  
"""

from pynuctran.solver import *
import numpy as np
from scipy.sparse.csr import csr_matrix

isotopes = [
    'U238',
    'U239',
    'Np239',
    'Pu239',
    'Pu240',
    'I135',
    'Xe135',
    'Cs135'
]


hl = [
    np.log(2) / 1.40902848e+17,
    np.log(2) / 1.4070E3,
    np.log(2) / 203558.4,
    np.log(2) / 7.609E11,
    np.log(2) / 2.0705E11,
    np.log(2) / 23688,
    np.log(2) / 32904,
    np.log(2) / 7.258E13
]

#---------------------------- PI DIST METHOD
sim = solver(isotope_names=isotopes)

sim.add_removal(0, hl[0], [-1])
sim.add_removal(1, hl[1], [ 2])
sim.add_removal(2, hl[2], [ 3])
sim.add_removal(3, hl[3], [-1])
sim.add_removal(4, hl[4], [-1])
sim.add_removal(5, hl[5], [ 6])
sim.add_removal(6, hl[6], [ 7])
sim.add_removal(7, hl[7], [-1])

sim.add_removal(0, 1E-4, [ 1])
sim.add_removal(3, 1E-5, [ 4])
sim.add_removal(0, 1E-5*(0.013157), [ 5])
sim.add_removal(0, 1E-5*(1.11541594E-04), [ 6])
sim.add_removal(0, 1E-5*(1-1.11541594E-04-0.013157), [-1])

w0 = [
    1E+12,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
]


total_time = 1E5
steps = int(1E10)
n_final = sim.solve(w0,total_time,steps)

#------------------ OBTAINING REFERENCE SOLUTION [CRAM48]----------------------------
# The CRAM solver can be obtained from MIT's CRPG codes repository:
# https://github.com/mit-crpg/opendeplete/blob/master/opendeplete/integrator/cram.py
#------------------------------------------------------------------------------------

A = sim.prepare_transmutation_matrix()
n0 = np.transpose(np.array(w0))
n_final_cram = CRAM48(A,n0,total_time)

# Prints the output of PyNUCTRAN solver and CRAM48, as well as their relative error.
for i in range(len(isotopes)):
    print('%+20.14e\t%+20.14e\t%+20.14e' % (n_final[i],n_final_cram[i],(float(n_final[i])-n_final_cram[i])/n_final_cram[i]))
