"""
  U-238 Depletion Problem
  
  The details of this problem can be found in this publication:
  http://dx.doi.org/10.1016/j.anucene.2016.09.004
  See pg. 268, Lago & Rahnema (2017), Annals of Nuclear Energy.
  
  Time step: 2.84E+4 (2 days).
  No. of sub-steps, t/dt: 1E+20.
   
"""

from pynuctran.solver import *
import numpy as np

isotopes = [ 'U-238','U-239','Np-239',
             'Pu-239','Pu-240','Pu-241',
             'Pu-242','Pu-243','Am241',
             'Am-243','Am244','Cm244'
           ]

hl = [   np.log(2)/1.4099935680E+17,  np.log(2)/1.4070E+03,  np.log(2)/2.0355840E+05, 
         np.log(2)/7.60853735110E+11,  np.log(2)/2.0704941360E+11, np.log(2)/4.509581040E+08,
         np.log(2)/1.178676360E+13,  np.log(2)/1.784160E+04,  np.log(2)/1.3651817760E+10, 
         np.log(2)/2.325795120E+11,  np.log(2)/3.6360E+04,  np.log(2)/5.715081360E+08
     ]

sim = solver(species_names = isotopes)
# Register decay events.
sim.add_removal(0, hl[0], [-1])
sim.add_removal(1, hl[1], [ 2])
sim.add_removal(2, hl[2], [ 3])
sim.add_removal(3, hl[3], [-1])
sim.add_removal(4, hl[4], [-1])
sim.add_removal(5, hl[5], [ 8])
sim.add_removal(6, hl[6], [-1])
sim.add_removal(7, hl[7], [ 9])
sim.add_removal(8, hl[8], [-1])
sim.add_removal(9, hl[9], [-1])
sim.add_removal(10,hl[10], [11])
sim.add_removal(11,hl[11], [4])

# Register absorption events.
sim.add_removal(0, 1E-4, [1])
sim.add_removal(3, 1E-4, [4])
sim.add_removal(4, 1E-4, [5])
sim.add_removal(5, 1E-4, [6])
sim.add_removal(6, 1E-4, [7])
sim.add_removal(9, 1E-4, [10])

# Setup the initial concentration.
w0 = [  1.00E+10, 1.00E+3, 0.00E+00, 
        0.00E+00, 0.00E+00, 0.00E+00,
        0.00E+00, 0.00E+00, 0.00E+00,
        0.00E+00, 0.00E+00, 0.00E+00
     ]

# Runs the calculation.
total_time = 8.64E+4
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