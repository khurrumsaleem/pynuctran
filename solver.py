'''
Python Library for Stochastic Nuclear Transmutation Solver (PYNUCTRAN)

Initially developed, designed and proposed by Dr M. R. Omar for the purpose of simulating various
nuclear  transmutations such as  decays,  fissions as well as neutron absorptions. PYNUCTRAN  was
initially  developed  to  avoid cumbersome  numerical issue  of  solving  the  nuclide  depletion 
equations. To date, there  are many numerical  depletion  solvers available such as CRAM, TTM etc. 
Interestingly,  the  designed  stochastic  solver is  much more easier to code,  but it  consumes 
computational power. Future research involves parallelizing PYNUCTRAN to reduce CPU time.


Created on 28-08-21.

(c) M.R.Omar, School of Physics, Universiti Sains Malaysia, 11800 Penang, Malaysia.

Permission is hereby granted,  free of charge,  to any person  obtaining  a copy of  PYNUCTRAN and 
associated documentation files (the "Code"), to deal in the Library without restriction, including 
without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
sell copies of the Code, and to permit persons to whom the Library is furnished to do so,  subject
to the following conditions:

The  above  copyright  notice  and  this permission notice  shall  be  included  in  all copies or 
substantial portions of the Software.

THE Code IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS  OR IMPLIED, INCLUDING  BUT NOT 
LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT  SHALL  THE AUTHORS  OR COPYRIGHT  HOLDERS  BE LIABLE FOR ANY CLAIM,  DAMAGES OR  OTHER 
LIABILITY,  WHETHER  IN  AN  ACTION  OF CONTRACT,  TORT  OR OTHERWISE,  ARISING FROM,  OUT OF OR IN 
CONNECTION WITH THE CODE OR THE USE OR OTHER DEALINGS IN THE CODE.

'''


import numpy as np
import random
from enum import Enum
import matplotlib.pyplot as plt
import inspect as ins
from time import process_time, time
from progress.bar import Bar
from datetime import datetime

'''
PrintMode is an enumeration that indicates the printf() mode.

'''
class PrintMode(Enum):
    Normal = 0
    Warning = 1
    Error = 2
    FatalError = 3

'''
PRINTF is a subroutine that formats the style of the code console IO.
message: The message to print.
newline: (Optional) Set to True if the caret advances after the message print.
                    Set to False if the caret do not advance after the message print.
mode   : (Optional) Specifies the mode/style of the message print. The value must be derived from
                    PrintMode enum.

'''
def printf(message: str, newline: bool = True, mode: PrintMode = PrintMode.Normal, ofile = None):
    
    if newline == True:
        if mode == PrintMode.Normal:
            print('' + message)
            if not ofile is None:
                ofile.write('' + message + '\n')
        elif mode == PrintMode.Warning:
            print('[WARNING] @' + str(ins.stack()[1][3]) + ': ' + message)
            if not ofile is None:
                ofile.write('[WARNING] @' + str(ins.stack()[1][3]) + ': ' + message + '\n')
        elif mode == PrintMode.Error:
            print('[ERROR] @' + str(ins.stack()[1][3]) + ': ' + message)
            if not ofile is None:
                ofile.write('[ERROR] @' + str(ins.stack()[1][3]) + ': ' + message + '\n')
        elif mode == PrintMode.FatalError:
            print('[FATAL] @' + str(ins.stack()[1][3]) + ': ' + message)
            if not ofile is None:
                ofile.write('[FATAL] @' + str(ins.stack()[1][3]) + ': ' + message + '\n')
    else:
        if mode == PrintMode.Normal:
            print('__ ' + message, end='')
            if not ofile is None:
                ofile.write('' + message)
        elif mode == PrintMode.Warning:
            print('[WARNING] @' + str(ins.stack()[1][3]) + ': ' + message, end='')
            if not ofile is None:
                ofile.write('[WARNING] @' + str(ins.stack()[1][3]) + ': ' + message)
        elif mode == PrintMode.Error:
            print('[ERROR] @' + str(ins.stack()[1][3]) + ': ' + message, end='')
            if not ofile is None:
                ofile.write('[ERROR] @' + str(ins.stack()[1][3]) + ': ' + message)
        elif mode == PrintMode.FatalError:
            print('[FATAL] @' + str(ins.stack()[1][3]) + ': ' + message, end='')
            if not ofile is None:
                ofile.write('[FATAL] @' + str(ins.stack()[1][3]) + ': ' + message)

"""
    ============================================================================================
    THIS SECTION HANDLES THE DATA MANAGEMENT OF THE STOCHASTIC NUCLEAR TRANSMUTATION SOLVER.
    ============================================================================================
"""
'''
    ENUM RemovalType
    PURPOSE     To distinguish between the nuclide removal types.

'''
class RemovalType(Enum):
    Fission = 1
    Absorption = 2
    Decay = 3
    NoRemoval = 4

'''
    CLASS Removal
    PURPOSE     Defines a nuclide removal method. This is actually a data structure
                that stores removal parameters of a specific removal. For instance,
                A beta decay of isotope A is a removal method, where the decay removes
                isotope A from the system by transmutating A into other species.

    MEMBERS     __init__(..)            Class constructor.
                AddDaughters(..)        Add the products to the removal method. Only 
                                        a pair of isotopes is accepted. For removals that
                                        produce one product, set '' for the second product.
                SampleDaughters(..)     Randomly selects one of the daugter pair based
                                        on the specified product yield.

                Id (str)            The string identifier of the removal method.
                Lambda (float)      The rate of the removal (per second).
                Type (RemovalType)  The removal type, i.e. Decay, Fission or Absorption.
                Daughters (list)    List of product pair and its yield.

'''
class Removal:

    def __init__(self, Id: str, Lambda: float, Type: RemovalType):
        self.Id = Id
        self.Lambda = Lambda
        self.Type = Type
        self.Daughters = []


    def AddDaughters(self, Yield: float, NuclideIdA: str, NuclideIdB: str):
        # Checks whether there exists any daughter pair duplicates.
        for nuclide in self.Daughters:
            if NuclideIdA == nuclide[1] or NuclideIdB == nuclide[2]:
                printf('Duplicate nuclides detected. ' + str(nuclide), mode=PrintMode.Error)
                return

        NewPair = [Yield, NuclideIdA, NuclideIdB]
        self.Daughters.append(NewPair)    
        printf('Daughter nuclides ' + str(NewPair) + ' have been created.')


    def SampleDaughters(self):

        # Compute the total Yield.
        YTotal = 0.0
        for pair in self.Daughters:
            YTotal += pair[0]

        # Generate a random number in [0,1) and multiply with the total yield.
        rnd = random.random() * YTotal

        # Sample the nuclide(s).
        i = 0
        yield_cum = 0.0

        for pair in self.Daughters:
            if rnd > yield_cum and rnd <= (yield_cum + pair[0]):
                return pair[1], pair[2]
            else:
                yield_cum += pair[0]

'''
    CLASS       Isotope (Structure)
    PURPOSE     Isotope is a  data  structure  that store s various  information of an isotope.
    MEMBERS     __init__()      Class constructor.
                AddRemoval()    Add a removal method to the isotope.

                Id          A string identification of the isotope.
                Name        A string name of the isotope.
                Z           The proton number of the isotope.
                A           The nucleus number of the isotope
                Removals    A list of removal methods experienced by the isotope.
                

'''
class Isotope:
    
    def __init__(self, Id: str, Name: str, Z: float, A: float):
        self.Id = Id
        self.Name = Name.lower()
        self.Z = Z
        self.A = A
        self.Removals = []

    def AddRemoval(self, Method: Removal):
        self.Removals.append(Method)
        printf('Removal method ' + Method.Id + ' have been added to ' + self.Id + '.')

    
'''
    CLASS       Nuclide (Structure) 
    PURPOSE     Nuclide  is  a  data  structure  that  stores various information of a simulated
                nuclide.
    MEMBERS     __init__()      Class constructor.

                NId (int)       An integer identifier distinguishing between different nuclides.
                                This integer is auto-generated.
                Isotope         The isotope species of the nuclide.
                CreationStep    The time step in  which the nuclide is  created in  the computer 
                                memory.
                

'''
class Nuclide:

    CurrentNId = 0
    def __init__(self, NIsotope: Isotope, CreationStep: int):
        self.NId = Nuclide.CurrentNId
        Nuclide.CurrentNId += 1
        self.Isotope = NIsotope
        self.CreationStep = CreationStep
    

"""
    ============================================================================================
    THIS SECTION HANDLES THE PHYSICS OF THE STOCHASTIC NUCLEAR TRANSMUTATION SOLVER.
    ============================================================================================
"""

'''
    CLASS     Physics
    PURPOSE   Handles the necessary routines that helps to simulate the transmutation processes.

    MEMBERS     __init__()                  Class constructor.
                AddIsotope()                Adds a new isotope to the simulation.
                GenerateNuclides()          Generates and adds nuclides of the same species into 
                                            the source table.
                SelectRemoval()             Samples the  removal  method  to process  during the 
                                            simulation.
                DoesRemovalHappens()        Decides  whether a  removal  process  does  actually 
                                            occur.
                ProcessRemoval()            Simulates a removal event.
                UpdateStep()                Re-organizes   the   source   table,  preparing  the 
                                            simulation of the next time step.
                IterateOverNuclides()       Simulates the physics of each of the nuclides in the
                                            source table.
                Run()                       Begin the simulation.
                PlotConcentrations()        Plot the concentration of each isotope.

                Isotopes (dict(Isotopes))   A list of isotopes defined in the simulation.
                NuclidesTable               A dictionary consists of two arrays for storing the 
                                            source nuclides and the product nuclides.
                Id (str)                    Stores the string identification of the simulation.
                Dt (float)                  The time interval.
                Steps (int)                 The total number of steps to simulate.
                SourceTable (str)           A string specifying  the key of the current  source 
                                            table.
                ProductTable (str)          A string specifying the key of  the current product 
                                            table.
                IsotopesToMonitor           A list of  isotope  identifications  specifying the 
                                            isotopes should be tallied  during  the simulation.
                OutputSkip (int)            The number of steps  skipped  for console printing.
                                            Note: The DisplayProgressOnly property  must be set
                                            to False.
                IsotopeConcentrations       A dictionary containing the  concentrations of  all 
                                            isotopes over the entire simulated  time steps. The
                                            keys of the dictionary is equivalent to the isotope
                                            Id  defined  in  the  simulation. The  rows  of the 
                                            dictionary corresponds to the time step.
                                            
                                            [Isotope1]  [Isotope2]  [Isotope3] ...
                                            c1-1        c2-1        c3-1
                                            c1-2        c2-2        c3-2
                                            c1-3        c2-3        c3-3
                                            :           :           :
                                            :           :           :


                ProcessedCount              The total number of nuclides simulated during the 
                                            single simulation.
                InitialNuclidesCount        A dictionary storing the total number of initial 
                                            nuclides (when the time step = 0).
                IsotopeLifetime             A dictionary that stores the average lifetime of 
                                            all isotopes defined in the simulation.                           
                CountTable                  A dictionary that stores the CURRENT step's total
                                            number of nuclides.

'''
class Physics:

    def __init__(self, Id: str, TimeInterval: float, Steps: int):
        self.Isotopes = {}
        self.NuclidesTable = {'A': [], 'B': []}
        self.Id = Id
        self.Dt = TimeInterval
        self.Steps = Steps
        self.CurrentStep = 0
        self.SourceTable = 'A'
        self.ProductTable ='B'
        self.IsotopesToMonitor = []
        self.OutputSkip = 0
        self.DisplayProgressOnly = True
        self.IsotopeConcentrations = {}
        self.ProcessedCount = 0

        # Initialize the initial nuclides count table.
        self.InitialNuclidesCount = {}
        # Initialize the nuclides count table.
        self.CountTable = {}
        # Initialize the nuclides lifetime table.
        self.IsotopeLifetime = {}


    def AddIsotope(self, NewIsotope: Isotope):

        for item in self.Isotopes.keys():
            if self.Isotopes[item].Id == NewIsotope.Id:
                printf('Duplicate isotope is detected, ' + str(Isotope(item).Id) + '.', mode=PrintMode.Error)
                return
        self.Isotopes[NewIsotope.Id] = NewIsotope
        self.CountTable[NewIsotope.Id] = 0

        return

    """
        SUBROUTINE  GenerateNuclides
        
        PURPOSE     Unconditionally adds the desired nuclides with the specified isotope species
                    into the nuclide source table.
        PARAMETERS  [NuclideIsotope: Isotope] is the isotope species of the nuclides. 
                    [Count: int] is the number of nuclides to add.

        Last updated by M.R.Omar on 30/08/2021.
        
    """
    def GenerateNuclides(self, NuclideIsotope: Isotope, Count: int):
        for i in range(Count):
            NewNuclide = Nuclide(NuclideIsotope, -1)
            self.NuclidesTable[self.SourceTable].append(NewNuclide)
        printf(str(Count) + ' ' + str(NuclideIsotope.Id) +  ' nuclide(s) have been created.')
        
    """
        FUNCTION    SelectRemoval
        
        PURPOSE     To sample a nuclide removal method (Also selects the NoRemoval event).
        PARAMETERS  [SelectedIsotope: Isotope] is the isotope species of the simulated nuclide. 
        RETURN      [:Removal] Returns the selected removal method.

        Last updated by M.R.Omar on 30/08/2021.
        
    """
    def SelectRemoval(self, SelectedIsotope: Isotope) -> Removal:
        
        # Generate random number for selection.
        rnd = random.random()

        # Compute the probabilities for various events such as
        # (a) Nothing happens
        # (b) Removal-A, Removal-B and so on.
        Probs = [1.0]
        Exps  = []
        for rm in SelectedIsotope.Removals:
            Exps.append(np.exp(-rm.Lambda * self.Dt))
            Probs.append(1.0) # Ensure Probs size is equal to Exps size + 1.

        # Calculate the probability of nothing happens, and insert it in Probs.
        for i in range(1,len(Exps)+1):
            for j in range(len(Probs)):
                if j == i:
                    Probs[j] *= (1 - Exps[i-1])
                else:
                    Probs[j] *= Exps[i-1]

        CumProb = 0.0
        SumProbs = np.sum(Probs)
        for i in range(len(Probs)):
            if CumProb < (rnd*SumProbs) and (rnd*SumProbs) <= (Probs[i]+CumProb):
                if i == 0:
                    return None
                else:
                    return SelectedIsotope.Removals[i-1]
            CumProb += Probs[i]


    """
        FUNCTION    DoesRemovalHappens
        
        PURPOSE     Decides whether a removal process occurs.
        PARAMETERS  [rm: Removal] is the removal of the nuclide species. 
        RETURN      [:Bool] Returns True if the removal occurs, False otherwise.

        Last updated by M.R.Omar on 30/08/2021.
        
    """
    def DoesRemovalHappens(self, rm: Removal) -> bool:
        # Define the probability.
        Prob = 1 - np.exp(-rm.Lambda*self.Dt)
        rnd = random.random()
        if rnd > Prob:
            return True
        else:
            return False


    """
        SUBROUTINE  ProcessRemoval
        
        PURPOSE     Processes a removal event, i.e. kills the removed nuclide, add new
                    nuclides in the source table.
        PARAMETERS  [Type: RemovalType] is the type of the undergoing removal. 

        Last updated by M.R.Omar on 30/08/2021.
        
    """
    def ProcessRemoval(self, Type: RemovalType, NuclideToProcess: Nuclide):

        if Type == RemovalType.NoRemoval:
            self.NuclidesTable[self.ProductTable].append(NuclideToProcess)
            self.CountTable[NuclideToProcess.Isotope.Id] += 1
            return

        for rm in NuclideToProcess.Isotope.Removals:
            if rm.Type == Type:
                Products = rm.SampleDaughters()
                break
        if not Products[0] == '':
            DaughterNuclide = Nuclide(self.Isotopes[Products[0]], self.CurrentStep)
            self.NuclidesTable[self.ProductTable].append(DaughterNuclide)
            self.CountTable[Products[0]] += 1
        if not Products[1] == '' and Type == RemovalType.Fission:
            DaughterNuclide = Nuclide(self.Isotopes[Products[1]], self.CurrentStep)
            self.NuclidesTable[self.ProductTable].append(DaughterNuclide)
            self.CountTable[Products[1]] += 1
        return



    def UpdateStep(self):

        # Update the time step...
        self.CurrentStep += 1

        # Clear the source table...
        self.NuclidesTable[self.SourceTable] = []

        # Interchange the nuclides table...
        temp = self.SourceTable
        self.SourceTable = self.ProductTable
        self.ProductTable = temp

        # Reset the nuclides count table.
        self.CountTable = {}
        for item in self.Isotopes.keys():
            self.CountTable[item] = 0

        return

    '''

    '''
 
    def IterateOverNuclides(self, CurrentStep: int):

        # This is the loop header for iterating over the nuclides.
        for CurrentNuclide in self.NuclidesTable[self.SourceTable]:

            # Record the number of processed nuclides so far.
            self.ProcessedCount += 1

            # Compute the initial no of nuclides for all isotopes.
            if CurrentStep == 0:
                self.InitialNuclidesCount[CurrentNuclide.Isotope.Id] += 1

            # Begin selecting the face of the nuclide and process the fate.
            CurrentRemoval = self.SelectRemoval(CurrentNuclide.Isotope)
            if CurrentRemoval is None:
                self.ProcessRemoval(RemovalType.NoRemoval, CurrentNuclide)
            else:
                # Check whether the selected removal occur.
                if self.DoesRemovalHappens(CurrentRemoval) is True:
                    self.ProcessRemoval(CurrentRemoval.Type, CurrentNuclide)
                    # Tally the isotope lifetime.
                    self.IsotopeLifetime[CurrentNuclide.Isotope.Id].append(float(self.CurrentStep - CurrentNuclide.CreationStep))

    def ConcentrationToString(self) -> str:
        s = 'ISOTOPE CONCENTRATIONS vs TIME STEPS\n\n'
        s += str(self.IsotopesToMonitor) +'\n'
        for i in range(self.Steps+1):
            s += 'step = %-6i\t' % i
            for iso in self.IsotopesToMonitor:
                s += '%10g\t' % (self.IsotopeConcentrations[iso][i])
            s += '\n'
        return s

    def IsotopesToString(self) -> str:
        s = '** y = Product Yield [%], P1 = Removal Product 1, P2 = Removal Product 2.\n'
        s += '\n'
        for i in self.IsotopesToMonitor:
            iso = self.Isotopes[i]
            s += '%6s ________ [Name] %s\n' % (iso.Id, iso.Name)
            s += '          |____ [Z]    %-3i\n' % (iso.Z)
            s += '          |____ [A]    %-3i\n' % (iso.A)
            s += '          |____ [Removal Method(s)] %i type(s)\n' % (len(iso.Removals))
            s += '                       |\n'
            for i in range(len(iso.Removals)):
                rm = iso.Removals[i]
                s += '                       |____ %s\n' % (rm.Id)
                if len(iso.Removals) > 1 and i != (len(iso.Removals)-1):
                    s += '                       |      |____ [Type     ] %s\n' % (str(rm.Type))
                    s += '                       |      |____ [Lambda   ] %g /second\n' % (rm.Lambda)
                    s += '                       |      |____ [y, P1, P2] %s\n' % (str(rm.Daughters))
                    s += '                       |\n'
                else:
                    s += '                              |____ [Type     ] %s\n' % (str(rm.Type))
                    s += '                              |____ [Lambda   ] %g /second\n' % (rm.Lambda)
                    s += '                              |____ [y, P1, P2] %s\n' % (str(rm.Daughters))                    
            s += '\n\n'     
        return s


    def Run(self, Seed: int = 0):

        # Create a log file to print outputs.
        fname = str('out ' + self.Id + '.txt')
        out = open(file=fname, mode='w')
        out.write('PYNUCTRAN ' + str(datetime.now()) + '\n')
        out.write('Version 1.9.821\n')
        out.write('Simulation ID: ' + self.Id + '\n')
        out.write('\n\n')
        out.write(self.IsotopesToString() + '\n')

        # Set the random number seed.
        random.seed(Seed)

        # Record the start time.
        t_start = process_time()

        # Create plot points...
        self.IsotopeConcentrations = {}
        for i in self.IsotopesToMonitor:
            # Reserving the first row of all columns for the initial concentration @ 0s.
            self.IsotopeConcentrations[i] = np.array([0.0])
            # Resetting the initial nulides count to zero for all isotopes.
            # The purpose of computing the initial nuclides is to serve it as a 
            # normalization factor, i.e. N_normalized = N / N0, where N0 is the initial
            # concentration.
            self.InitialNuclidesCount[i] = 0
            # Initializing the isotope lifetime table.
            self.IsotopeLifetime[i] = []

        # Tallying the number of nuclides for each isotope species.
        for nuc in self.NuclidesTable[self.SourceTable]:
            self.IsotopeConcentrations[nuc.Isotope.Id][0] += 1

        # Loop header for progress bar reporting.
        with Bar('Running ' + self.Id, max=self.Steps, fill='░', suffix='%(percent)d%%') as bar:
            
            # Do not display progress bar if DisplayProgressOnly is False.
            if not self.DisplayProgressOnly:
                bar.clearln()
            
            # This is the loop header for iterating over time steps.
            for i in range(self.Steps):
                self.IterateOverNuclides(i)
            

                if len(self.IsotopesToMonitor) == 0:
                    printf('No isotope is being monitored. Please specify the IsotopesToMonitor list.', mode=PrintMode.FatalError, ofile=out)
                    return
                else:    
                    # Print the count monitoring for each isotope.
                    if i % self.OutputSkip == 0 and not self.DisplayProgressOnly:
                        # Print the current step.
                        printf('i = ' + str(i), newline=False, ofile=out)

                    for Isotope in self.IsotopesToMonitor:
                        if i % self.OutputSkip == 0 and not self.DisplayProgressOnly:
                            printf('  ' + Isotope + ' = ' + str(self.CountTable[Isotope]), newline=False, ofile=out)
                        self.IsotopeConcentrations[Isotope] = \
                            np.append(self.IsotopeConcentrations[Isotope], [self.CountTable[Isotope]])

                    if i % self.OutputSkip == 0 and not self.DisplayProgressOnly:
                        printf('', ofile=out)
                
                # Update the current step by preparing the nuclide source etc.
                self.UpdateStep()

                # Update the progress bar if DisplayProgressOnly=True.
                if self.DisplayProgressOnly:
                    bar.next()
            
            
        # Record the stop time.
        t_stop = process_time()

        # Print the computational stats.
        printf('\n\nSIMULATION STATISTICS\n', ofile=out)
        printf('Time elapsed           : ' + str(t_stop-t_start) + ' sec(s).', ofile=out)
        printf('Total steps            : %g.' % self.Steps, ofile=out)
        printf('Time interval [secs]   : %g.' % self.Dt, ofile=out)
        printf('RNG Seed               : %g.' % Seed, ofile=out)
        printf('Nuclides processed     : %g.' % self.ProcessedCount, ofile=out)
        printf('Initial concentrations : %-5s (%-g)' % (self.IsotopesToMonitor[0], self.InitialNuclidesCount[self.IsotopesToMonitor[0]]), ofile=out)
        if len(self.IsotopesToMonitor) > 1:
            for i in range(1, len(self.IsotopesToMonitor)):
                printf('                         %-5s (%-g)' % (self.IsotopesToMonitor[i], self.InitialNuclidesCount[self.IsotopesToMonitor[i]]), ofile=out)
        printf('Final concentrations   : %-5s (%-g)' % (self.IsotopesToMonitor[0], 
                self.IsotopeConcentrations[self.IsotopesToMonitor[0]][-1]), ofile=out)
        if len(self.IsotopesToMonitor) > 1:
            for i in range(1, len(self.IsotopesToMonitor)):
                printf('                         %-5s (%-g)' % (self.IsotopesToMonitor[i], 
                    self.IsotopeConcentrations[self.IsotopesToMonitor[i]][-1]), ofile=out)
        
        # Compute the average lifetime...
        mean_lifetime = np.Infinity
        if len(self.IsotopeLifetime[self.IsotopesToMonitor[0]]) != 0:
             mean_lifetime = np.average(self.IsotopeLifetime[self.IsotopesToMonitor[0]])
        printf('Average lifetime       : %-5s (%-g steps)' % \
            (self.IsotopesToMonitor[0], mean_lifetime), ofile=out)
        if len(self.IsotopesToMonitor) > 1:
            for i in range(1, len(self.IsotopesToMonitor)):
                mean_lifetime = np.Infinity
                if len(self.IsotopeLifetime[self.IsotopesToMonitor[i]]) != 0:
                    mean_lifetime = np.average(self.IsotopeLifetime[self.IsotopesToMonitor[i]])
                printf('                         %-5s (%-g steps)' % (self.IsotopesToMonitor[i], mean_lifetime), ofile=out)
        printf('%s simulation ended. [Total time elapsed: %-g secs.]' % (self.Id, process_time()-t_start), ofile=out)
        

        # Print the concentrations table into the output file and close the file.
        out.write('\n\n\n' + self.ConcentrationToString() +'\n')
        out.close()


    '''
        SUBROUTINE      PlotConcentrations

        PURPOSE         Plots the concentration curves over time steps.
        PARAMS          Ids: The list of isotope identification strings, specifying which isotope concentration to
                        include in the plot.
                        Range: The time step plot range. Default is from 0 to the specified maximum steps.
                        Normalize: A boolean specifying whether the concentrations should be normalized or not.
                        A normalized concentration is equals to [the nuclide count] / [total initial nuclides].
                        Color: A python dictionary specifying the curve colors that differentiate the isotopes. 
                        The key of the dictionary is equivalent to the isotope Id defined in the simulation.

    '''
    def PlotConcentrations(self, Ids: list = [], Range: tuple = (0, 0), Normalize: bool = False, Color: dict = {}):
        if Ids == []:
            printf('Plotting concentrations ' + str(self.IsotopesToMonitor) + '...')
        else:
            printf('Plotting concentrations ' + str(Ids) + '...')
        legend_strs = []
        for i in self.IsotopesToMonitor:
            if (i in Ids) or Ids == []:
                if Range == (0,0):
                    temp_array = np.array(self.IsotopeConcentrations[i])
                else:
                    temp_array = np.array(self.IsotopeConcentrations[i][Range[0]:Range[1]])
                total_initial_nuclides = 0
                for j in self.InitialNuclidesCount.keys():
                    total_initial_nuclides += self.InitialNuclidesCount[j]
                if Normalize:
                    if len(Color) == 0:
                        plt.plot(temp_array/total_initial_nuclides)
                    else:
                        if i in Color.keys():
                            plt.plot(temp_array/total_initial_nuclides, c = Color[i])
                        else:
                            printf('Plot color is not specified for isotope ' + i, mode=PrintMode.Warning)
                            plt.plot(temp_array/total_initial_nuclides)
                else:
                    if len(Color) == 0:
                        plt.plot(temp_array)
                    else:
                        if i in Color.keys():
                            plt.plot(temp_array, c = Color[i])
                        else:
                            printf('Plot color is not specified for isotope ' + i, mode=PrintMode.Warning)
                            plt.plot(temp_array)
                legend_strs.append(i)
        if Ids == []:
            plt.legend(self.IsotopesToMonitor)
        else:
            plt.legend(legend_strs)
        if Normalize:
            plt.ylabel('Isotope Fraction')
        else:
            plt.ylabel('Isotope Count')
        plt.xlabel('Time Step')
        plt.show()
