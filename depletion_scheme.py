import numpy as np
import decimal as dc
import xml.etree.ElementTree as ET
import time as tm
from pynuctran.solver import *
        
'''
    SECTION III: DEPLETION DATA PRE-PROCESSING ................................................ SEC. III
    
    *******************************************************************************************
    THIS SECTION ENABLES THE RETRIEVAL OF NUCLIDES DATA FROM ENDF.
    THE NUCLIDE DATA ARE STORED IN chain_endfb71.xml. 
    
    The XML file can be retrieved here:
    https://github.com/mit-crpg/opendeplete/blob/master/chains/chain_endfb71.xml
    *******************************************************************************************
'''
class depletion_scheme:

    '''
        Defines the depletion scheme in the code, based on ENDFB17 data. The nuclide data are
        stored in an xml file 'chains_endfb71.xml'. Here, the depletion chains are created
        based on the user specified reaction rates and species.

        Parameters:

        xml_data_location: A string specifying the location of chains_endfb71.xml on the disk.
        rxn_rates        : A 2D python dictionary containing the reaction rates of various
                           removal events. For example,

                           rxn_rates = {
                                'U238' : {'(n,gamma)': 1E-4, 'fission': 1E-5},
                                'Pu239': {'(n,gamma)': 1E-5},
                           }

    '''
    @staticmethod
    def build_chains(solver: solver, rxn_rates, xml_data_location: str = 'chain_endfb71.xml'):
        t0 = tm.process_time()

        species_names = solver.species_names
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
                        solver.add_removal(parent_id, decay_rate_adjusted, [solver.__no_product__])
                
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
                                removal_rate = dc.Decimal('%.15g' % rxn_rates[parent][removal.attrib['type']])
                                if daughter in species_names:
                                    daughter_id = species_names.index(daughter)
                                    solver.add_removal(parent_id, removal_rate, [daughter_id])
                                else:
                                    solver.add_removal(parent_id, removal_rate, [solver.__no_product__])

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
                                                    yields = [dc.Decimal(y) for y in param.text.split()]
                             
                                total_fission_rate = rxn_rates[parent]['fission']
                                yields_to_add = []
                                daughters_id_to_add = []
                                for product in products:
                                    if product in species_names:
                                        daughters_id_to_add.append(species_names.index(product))
                                        yields_to_add.append(yields[products.index(product)])
                                parent_id = species_names.index(species_name)
                                solver.add_removal(parent_id, total_fission_rate, daughters_id_to_add, yields_to_add)
  
                               
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

    @staticmethod
    def get_all_species_names_range(xml_data_location: str, AMin: int, AMax: int) -> list:
        tree = ET.parse(xml_data_location)
        root = tree.getroot()
        
        species_names = []
        for species in root:
            name = species.attrib['name']
            name = name.split('_')[0]
            x = ''
            for c in name:
                if c.isnumeric():
                    x += c
            A = int(x)
            if A >= AMin and A <= AMax:
                species_names.append(name)
        return species_names        
