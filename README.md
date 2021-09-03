# PyNUCTRAN
## A Pseudo-Monte Carlo Simulation of Nuclear Transmutation.

### ```pip install pynuctran```

PyNUCTRAN is a Python library created by M. R. Omar for simulating various nuclear transmutations such as decays, fissions, and neutron absorptions. The code helps physicists avoid cumbersome numerical issues of solving the nuclide depletion equations (also known as the Bateman's equations). These issues include the stiffness of the Batemans equations due to the complex decay chain problems. To date, there are many numerical depletion solvers available such as CRAM, TTM etc. Unfortunately, these methods are complicated and requires sophisticated mathematical methods. It is also possible to simulates the actual transmutation processes using Monte Carlo method via iterations over a massive amount of nuclides. Alas, the simulation speed linearly increase with increasing accuracy, which makes it not practical. Thanks to the variance reduction technique in Monte Carlo method, the concept of isotope weight is applied in this algorithm to boost the CPU time and accuracy. PyNUCTRAN simulates the actual transmutation processes, however, it does not generate pseudo-random numbers nor any random sampling. By applying the standard Poisson statistics, isotope weights can be adjusted according the the pre-determined transmutation probabilities. Such a trick substantially improve the simulation time! The method is easy to understand, and it can be used by physicists from various mathematical background.

## Features

- Capable of simulating complex transmutation chains.
- PyNUCTRAN simulates the actual transmutations processes without losing the accuracy of the computation.
- Create a plot of various isotope concentrations.
- Helps nuclear physics students to understands transmutations processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

If you don't prefer dealing with complicated mathematical methods to solve Bateman's equations,  give PyNUCTRAN a try!

## Summary of the Stochastic Transmutation Method

A transmutation process involves the removal of a nuclide from a system. Then it leads to the creation of another daughter nuclide. For instance, the decay of U-238 into Th-234 involves removing U-238 from the system via alpha decay, which in fact mutating the U-238 nucleus into Th-234. In reality, such a transmutation process occur at a certain rate, [lambda]. PyNUCTRAN works by first accumulating the removal parameters from the user. The removal parameters include the rate, parent isotope and the daughter isotope(s). Practically, any number of removals can be added to PyNUCTRAN. Next, PyNUCTRAN prepares the removal probability (RP) table that stores the probability of various transmutation processes to occur. The RP-table boosts the computational speed by substantially reducing flops since the probabilities are pre-calculated before the simulation.

The simulation requires the division of time into a regular interval, dt, of N steps. Consider an isotope-i which consists of J(i) removal processes. The probability a removal of isotope-i from a system due to j-th removal process can be derived from Poisson statistics:

<img src="https://latex.codecogs.com/svg.image?p_{il}&space;=&space;\prod_{j=1}^{J_i}\left\{&space;\delta_{lj}&plus;(-1)^{\delta_{lj}}&space;e^{-\lambda_{j}\Delta&space;t}\right\}" title="p_{ij} = \prod_{l=1}^{J_i}\left\{ \delta_{lj}+(-1)^{\delta_{lj}} e^{-\lambda_{l}\Delta t}\right\}" />

The probability of isotope-i for not being removed from the system is given by

<img src="https://latex.codecogs.com/svg.image?p_{i0}&space;=&space;\prod_{j=1}^{J_i}&space;e^{-\lambda_{j}\Delta&space;t}" title="p_{i0} = \prod_{j=1}^{J_i} e^{-\lambda_{j}\Delta t}" />

At the beginning of the simulation, each isotope-i has an initial weight, w0(i), which corresponds to its initial concentration. During each time step, t, the weight of its daughter isotopes-k, w(k≠i) (due to all removals defined for isotope-i), is adjusted:

<img src="https://latex.codecogs.com/svg.image?w^{\(t&plus;1)}_k&space;=&space;w^{\(t)}_k&space;&plus;&space;w^{\(t)}_i&space;p_{ij}" title="w^{\(t+1)}_k = w^{\(t)}_k + w^{\(t)}_i p_{ij}" />

After processing all removal, the weight of isotope-i, w(i) is adjusted,

<img src="https://latex.codecogs.com/svg.image?w^{\(t&plus;1)}_i&space;=&space;w^{\(t)}_i&space;p_{i0}" title="w^{\(t+1)}_i = w^{\(t)}_i p_{i0}" />

These procedures repeat for N steps. Thus, the solution of the associated Bateman's equation can be viewed by plotting w versus the time steps. Note also that this algorithm does not require any random sampling. In essence, one may perceive this method as a Monte Carlo method because it simulates the actual transmutation events during each time step. Therefore, the method is a pseudo-Monte Carlo technique, and its solution is free from random errors. 

## Some Python Examples

_Installing PYNUCTRAN from PyPI repository._
```bash
pip install pynuctran
```
_Importing PYNUCTRAN library to your Python code._
```python
from pynuctran.solver import *
```

_Running the simulation._
```python
# Define isotope names.
iso_names = ['U235', 'U236','U237', 'Np237']
# Initialize the solver.
sim = solver(isotope_names=iso, time_interval=86.4, steps=1000)

# Add the removal processes for U235...
# U-235 decay: isotope_id=0 since it is the first iso_names element ...
# products=[-1] indicates the code should not monitor the daugther isotopes. The rate must be in per second.
sim.add_removal(isotope_id=0, rate=np.log(2)/2.2210238880E+16, products=[-1])
# U-235 neutron absorptions.
sim.add_removal(isotope_id=0, rate=1E-4, products=[1])

# Add the removal processes for U236, U237 and so on...
sim.add_removal(isotope_id=1, rate=np.log(2)/7.390789920E+14, [-1])
sim.add_removal(isotope_id=1, rate=1E-4, [2])
sim.add_removal(isotope_id=2, rate=np.log(2)/5.8320E+05, [3])
sim.add_removal(isotope_id=3, rate=np.log(2)/6.7659494310E+13, [-1])

# Prepare the RP-table (removal probability table)
sim.preprocess_p_table()

# Assign the initial weight of all isotopes. The length of w0 is equal to the number of isotopes being monitored.
w0 = [1.0, 0.0, 0.0, 0.0]

# Run the simulation.
final_w = sim.run(w0)

# Plot the concentrations of all isotopes.
sim.plot_concentrations(w=final_w, isotopes_to_plot=[0,1,2,3])
```

_Sample output_

<img src="https://user-images.githubusercontent.com/33319386/132009122-79e95a4e-0980-4185-af6e-eaabe323eedc.png" width="400">

Sample PYNUCTRAN output using ```solver.plot_concentrations()``` for Lago & Rahnema (2017) benchmark test #3. [doi: http://dx.doi.org/10.1016/j.anucene.2016.09.004] The simulation only takes 31ms after running 1000 time steps.


## License (MIT)

Permission is hereby granted,  free of charge,  to any person  obtaining  a copy of  PyNUCTRAN and associated documentation files (the "Code"), to deal in the Library without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Code, and to permit persons to whom the Library is furnished to do so,  subject to the following conditions:

The  above  copyright  notice  and  this permission notice  shall  be  included  in  all copies or substantial portions of the Software.

THE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  OR IMPLIED, INCLUDING  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT  SHALL  THE AUTHORS  OR COPYRIGHT  HOLDERS  BE LIABLE FOR ANY CLAIM,  DAMAGES OR  OTHER LIABILITY,  WHETHER  IN  AN  ACTION  OF CONTRACT,  TORT  OR OTHERWISE,  ARISING FROM,  OUT OF OR IN CONNECTION WITH THE CODE OR THE USE OR OTHER DEALINGS IN THE CODE.

## Library Documentation

Before reading this section, you must have the basic idea of object-oriented programming (OOP). This section uses the standard OOP terms to ensure effective explanation. Of course, thanks to the simplicity of the method implemented in PyNUCTRAN. Consequently, the library consists of a small number of classes and structures. You will discover that the library only has 600++ lines of Python code, which is a lot less than the state-of-the-art Monte Carlo transport codes!

#### ``` solver.solver(isotope_names: list, time_interval: float = 1.0, steps: int = 100, precision: float = 1E-20) ```
Initializes the solver. ```isotope_names``` is a list of isotope names involved in the simulation. An isotope ID corresponds to the index of ```isotope_names```. ```time_interval``` is the interval of between time steps (in seconds). ```steps``` is the total number of time steps to simulate and ``` precision``` is the tolerance value of isotope weights. Here, if an isotope weight falls below ```precision```, then PyNUCTRAN will stop updating the weight during the next time step (unless its weight increases to a value greater than ```precision``` at later times; for this case, PyNUCTRAN will re-monitor the isotope). 
#### ``` solver.add_removal(isotope_id: int, rate: float, products: list)```
Defines and adds a new removal process. ```isotope_id``` is the integer ID of the isotope subjected to the removal (the parent isotope). The ID corresponds to the index of ```isotope_names``` and  the isotope name is given by ```solver.isotope_names[isotope_id]```. ```rate``` is the rate of removal in /sec. For instance, ```rate``` is the decay constant of a decay process. ```products``` is a list of integers that corresponds to the IDs of the daughter isotopes. If the product is not known or not monitored, you must set ```products=[-1]```.
#### ``` solver.preprocess_p_table()```
Prepares the removal probabilities table.
#### ``` solver.advance_time_step(w: list)```
Computes the isotope weights of the current time steps. This method accepts ```w``` as the weights of the previous steps. Then, it computes the isotope weight of the current time steps and updates the values of ```w```.

#### ``` solver.run(w0: list) -> list```
Runs the simulation. ```w0``` is the initial isotope concentrations. The length of ```w0``` must equals to the total number of isotopes defined via ```solver.add_removal(...)```. Returns a 2D-list representing the isotope concentrations w versus time steps.


