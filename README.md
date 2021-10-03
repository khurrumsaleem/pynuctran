# PyNUCTRAN
## A Solver for Nuclear Transmutation Problems.

### ```pip install pynuctran```
<div align="justify">
PyNUCTRAN is a Python library created to simulate various nuclear transmutations such as decays, fissions, and neutron absorptions. The code helps physicists avoid cumbersome numerical issues of solving the nuclide depletion equations (also known as the Bateman's equations). These issues include the stiffness of the Bateman's equations due to the complex decay chain problems. The Bateman's equation of a depletion problem is given as follows:
</div>

\
<img src="https://latex.codecogs.com/svg.image?\frac{dn_{i}}{dt}=&space;\underset{\mathrm{production}}{\underbrace{\sum_{j=1}^{m}b_{ij}\lambda_{j}n_{j}&plus;\sum_{k=1}^{m}y_{ik}\Lambda_{ik}n_{k}}}&space;&space;-&space;&space;&space;\underset{\mathrm{removal}}{\underbrace{\left&space;(&space;\lambda_{i}&plus;\Lambda_{i}&space;\right)n_{i}}&space;" title="\frac{dn_{i}}{dt}= \underset{\mathrm{production}}{\underbrace{\sum_{j=1}^{m}b_{ij}\lambda_{j}n_{j}+\sum_{k=1}^{m}y_{ik}\Lambda_{ik}n_{k}}} - \underset{\mathrm{removal}}{\underbrace{\left ( \lambda_{i}+\Lambda_{i} \right)n_{i}} " />

<img src="https://latex.codecogs.com/svg.image?\inline&space;n_i" title="\inline n_i" /> is the atom density of isotope-*i*;\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\lambda_i" title="\inline \lambda_i" /> is the radioactive decay constant of isotope-*i* causing its removal from the system (the decay rate per second);\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\lambda_{ij}" title="\inline \lambda_{ij}" /> is the decay constant of isotope-*j* causing the production of isotope-*i* in the system; \
<img src="https://latex.codecogs.com/svg.image?\inline&space;b_{ij}" title="\inline b_{ij}" /> is the decay branching ratio from isotope-*j* into isotope-*i*;\
<img src="https://latex.codecogs.com/svg.image?\inline&space;y_{ij}" title="\inline y_{ij}" /> is the production yield from isotope-*k* into isotope-*i*; and,\
<img src="https://latex.codecogs.com/svg.image?\inline&space;\Lambda_{ij}" title="\inline \Lambda_{ij}" /> is the production rate of isotope-*i* due to the removal of isotope-*k* from the system (per second).

<div align="justify">
To date, there are many numerical depletion solvers available such as CRAM, TTM etc. Unfortunately, these methods are complicated and require sophisticated mathematical methods. It is also possible to simulates the actual transmutation processes using Monte Carlo method via iterations over a massive amount of nuclides. Alas, the simulation speed increases with increasing accuracy, which makes it not practical. Thanks to the variance reduction technique in Monte Carlo method, the concept of isotope weight is adopted. Here, the isotope weights can be adjusted and scaled using event probabilities, π, which will be described later. The method is easy to understand, and it can be used by physicists from various mathematical background.
</div>

## Features

- Capable of simulating complex transmutation chains.
- PyNUCTRAN simulates the actual transmutations processes without losing the accuracy of the computation.
- Create a plot of various isotope concentrations.
- Helps nuclear physics students to understands transmutation processes through simulations. They can create, study and design any depletion chains and simulate the transmutations.

If you don't prefer dealing with complicated mathematical methods to solve Bateman's equations,  give PyNUCTRAN a try!

## Summary of the Method

<div align="justify">
A transmutation process involves the removal of a nuclide from a system. Then it leads to the creation of another daughter nuclide. For instance, the decay of U-238 into Th-234 involves removing U-238 from the system via alpha decay, which in fact mutating the U-238 nucleus into Th-234. In reality, such a transmutation process occur at a certain rate, [lambda]. PyNUCTRAN works by first accumulating the removal parameters from the user. The removal parameters include the rate, parent isotope and the daughter isotope(s). Practically, any number of removals can be added to PyNUCTRAN. Next, PyNUCTRAN prepares the removal probability (RP) table that stores the probability of various transmutation processes to occur. The RP-table boosts the computational speed by substantially reducing flops since the probabilities are pre-calculated before the simulation.
</div>

<div align="justify">
The simulation requires the division of time into a regular interval, <i>dt</i>, of <i>N</i> steps. Consider an isotope-<i>i</i> which consists of <i>J</i>(<i>i</i>) removal processes. The probability a removal of isotope-<i>i</i> from a system due to <i>j</i>-th removal process can be derived from Poisson statistics:
</div>

\
<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{il}&space;=&space;\prod_{j=1}^{J_i}\left\{&space;\delta_{lj}&plus;(-1)^{\delta_{lj}}&space;e^{-\lambda_{j}\Delta&space;t}\right\}" title="\pi_{ij} = \prod_{l=1}^{J_i}\left\{ \delta_{lj}+(-1)^{\delta_{lj}} e^{-\lambda_{l}\Delta t}\right\}" />

The probability of isotope-*i* for not being removed from the system is given by

<img src="https://latex.codecogs.com/svg.image?\widetilde{\pi}_{i0}&space;=&space;\prod_{j=1}^{J_i}&space;e^{-\lambda_{j}\Delta&space;t}" title="\pi_{i0} = \prod_{j=1}^{J_i} e^{-\lambda_{j}\Delta t}" />

The normalized probability of removal-*l* to occur is given by (l=0 is for no-removal):

<img src="https://latex.codecogs.com/svg.latex?\pi_{il}&space;=&space;\frac{\widetilde{\pi}_{il}&space;}{\sum_{j=0}^{J_{i}}\widetilde{\pi}_{ij}}" title="P_{il} = \frac{f_{il} }{\sum_{j=0}^{J_{i}}f_{ij}}" />

At this point, we let <i>I</i> as the total number of species involved in the depletion problem, and we define the transfer matrix as 

<img src="https://latex.codecogs.com/svg.latex?\mathbf{A}&space;=&space;\begin{pmatrix}&space;\pi_{1\rightarrow&space;1}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;1}&space;\\&space;\vdots&space;&&space;\ddots&space;&&space;\vdots&space;\\&space;\pi_{1\rightarrow&space;I}&space;&&space;\cdots&space;&&space;\pi_{I\rightarrow&space;I}&space;\end{pmatrix}" title="\mathbf{A} = \begin{pmatrix} \pi_{1\rightarrow 1} & \cdots & \pi_{I\rightarrow 1} \\ \vdots & \ddots & \vdots \\ \pi_{1\rightarrow I} & \cdots & \pi_{I\rightarrow I} \end{pmatrix}" />

where π(k→i) is the transfer probability which is defined as

<img src="https://latex.codecogs.com/svg.latex?\pi_{k&space;\rightarrow&space;i}&space;=&space;\sum_{j\in&space;R}^{}&space;\pi_{kj}" title="\pi_{k \rightarrow i} = \sum_{j\in R}^{} \pi_{kj}" />

<div align="justify">
Here, R is a set of transmutation events that mutate species k into species i. Note that matrix A is a square matrix (<i>IxI</i>) with its columns as the parent species and the its rows as the daughter species. Now, let <b>w</b>(t) and <b>w</b>(0) be the column matrices representing the final and initial concentration of all species involved, respectively. Then, the final concentrations can be easily evaluated via the exponentiation of matrix A,
</div>
  
  
<img src="https://latex.codecogs.com/svg.latex?\mathbf{w}(t)&space;=&space;\mathbf{A}^{t/\Delta&space;t}&space;\mathbf{w}(0)" title="w(t) = \mathbf{A}^{t/\Delta t} w(0)" />

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
sim = solver(isotope_names=iso)

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


# Assign the initial weight of all isotopes. The length of w0 is equal to the number of isotopes being monitored.
w0 = [1.0, 0.0, 0.0, 0.0]

# Run the simulation.
final_w = sim.solve(w0, 1E+10, 10000000000)
```
_Sample output_

<img src="https://user-images.githubusercontent.com/33319386/132084062-98c99308-8e0e-49d9-a03d-1ce12b9a80a6.PNG" width="400"><img src="https://user-images.githubusercontent.com/33319386/132084073-abc31b3a-39dd-4491-920f-e27dd1d77d06.png" width="400">

Sample PYNUCTRAN output for Lago & Rahnema (2017) benchmark test #4. [doi: http://dx.doi.org/10.1016/j.anucene.2016.09.004] The simulation only takes 31ms after running 1000 time steps.


## License (MIT)

PyNUCTRAN is distributed under the MIT license.

## Library Documentation

<div align="justify">
Before reading this section, you must have the basic idea of object-oriented programming (OOP). This section uses the standard OOP terms to ensure effective explanation. Of course, thanks to the simplicity of the method implemented in PyNUCTRAN. Consequently, the library consists of a small number of classes and structures. You will discover that the algorithm requires <80 lines of Python code, which is a lot less than the state-of-the-art solvers!
</div>


#### ``` solver.solver(isotope_names: list) ```
Initializes the solver. ```isotope_names``` is a list of isotope names involved in the simulation. An isotope ID corresponds to the index of ```isotope_names```.  
#### ``` solver.add_removal(isotope_id: int, rate: float, products: list)```
Defines and adds a new removal process. ```isotope_id``` is the integer ID of the isotope subjected to the removal (the parent isotope). The ID corresponds to the index of ```isotope_names``` and  the isotope name is given by ```solver.isotope_names[isotope_id]```. ```rate``` is the rate of removal in /sec. For instance, ```rate``` is the decay constant of a decay process. ```products``` is a list of integers that corresponds to the IDs of the daughter isotopes. If the product is not known or not monitored, you must set ```products=[-1]```.
#### ``` solver.solve(w0: list, t: float, substeps: int) -> list```
Runs the simulation. ```w0``` is the initial isotope concentrations. The length of ```w0``` must equals to the total number of isotopes defined via ```solver.add_removal(...)```. Returns a column matrix representing the isotope concentrations w.


