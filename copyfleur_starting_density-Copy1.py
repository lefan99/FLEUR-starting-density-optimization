from aiida_fleur.data.fleurinp import FleurinpData
from aiida import load_profile
from aiida.orm import Dict, load_node , load_code
from aiida.engine import submit, run
from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation
from aiida.plugins import DataFactory, CalculationFactory
from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.data.fleurinpmodifier import FleurinpModifier
from itertools import product
import numpy as np


load_profile()

def gen_inp_itmax( structure , inpgen = load_code(1112), itmax = 1 ):
    
    FleurinpData = DataFactory('fleur.fleurinp')
    FleurCalculation = CalculationFactory('fleur.fleur')

    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine": 1}, 'queue_name': 'th1', 'withmpi': False}







    struct = structure


    inputs = FleurinputgenCalculation.get_builder()
    inputs.code = inpgen
    inputs.structure = struct
    inputs.metadata.options = options

    results , inpgen_node  = run.get_node(FleurinputgenCalculation , **inputs)


    print(results , inpgen_node)

    print('inpgen in progress pk is : ' , inpgen_node.pk)

    fleurinp = inpgen_node.outputs.fleurinpData

    print('fleurinp created from inpgen process, pk is: ' , fleurinp.pk)



    fleurmode = FleurinpModifier(fleurinp)
    dict_change = { 
        'itmax' : itmax }

    fleurmode.set_inpchanges(dict_change)
    inp_mod = fleurmode.freeze()
    print('PK of itmax:1 modified inp {}'.format(inp_mod.pk))
    
    return inp_mod


def starting_density_calc(inp , fleur = load_code(150) ,  step=0.5):
    
    inp_mod = inp
    
    constraints = { '(1s1/2)' : 2 ,
                   '(1d3/2)' : 4 ,
                   '(1d5/2)' : 6 ,
                   '(1p1/2)' : 2 ,
                   '(1p3/2)' : 4 ,
                   '(2s1/2)' : 2 ,
                   '(2d3/2)' : 4 ,
                   '(2d5/2)' : 6 ,
                   '(2p1/2)' : 2 ,
                   '(2p3/2)' : 4 ,
                   '(3s1/2)' : 2 ,
                   '(3d3/2)' : 4 ,
                   '(3d5/2)' : 6 ,
                   '(3p1/2)' : 2 ,
                   '(3p3/2)' : 4 ,
                   '(4s1/2)' : 2 ,
                   '(4d3/2)' : 4 ,
                   '(4d5/2)' : 6 ,
                   '(4p1/2)' : 2 ,
                   '(4p3/2)' : 4 ,
                   '(5s1/2)' : 2 ,
                   '(5d3/2)' : 4 ,
                   '(5d5/2)' : 6 ,
                   '(5p1/2)' : 2 ,
                   '(5p3/2)' : 4 ,
                   '(6s1/2)' : 2 ,
                   '(6d3/2)' : 4 ,
                   '(6d5/2)' : 6 ,
                   '(6p1/2)' : 2 ,
                   '(6p3/2)' : 4 ,
                   '(7s1/2)' : 2 ,
                   '(7d3/2)' : 4 ,
                   '(7d5/2)' : 6 ,
                   '(7p1/2)' : 2 ,
                   '(7p3/2)' : 4 ,
                   
                
                  }


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}



    number_valence_electrons = inp_mod.inp_dict['cell']['bzIntegration']['valenceElectrons']



    valence_orbs = inp_mod.inp_dict['atomSpecies'][0]['electronConfig']['valenceConfig']



    orbital_occupation = []

    occupation_combinations = {}

    for orbital in valence_orbs:
    
        if orbital in constraints:
        
            max_occupation = int(constraints[orbital] /2)
         
            occupation_combinations[orbital] = np.linspace(0 , max_occupation, round((max_occupation/step) + 1) ).tolist()
        
        else:
        
            print('error could not match valence electrons with saved list of constraint: ' , orbital)
        
        


        
        
    
    list_orb =[]

    for orbital in occupation_combinations:

        list_orb.append(occupation_combinations[orbital])






    all_comb = list(product(*list_orb))

    valid_lis = []

    for combination in all_comb:

        if sum(combination)*2 == number_valence_electrons:

            valid_lis.append(combination)


    print(valid_lis , 'number of possible combinations; ' , len(valid_lis))
    
    dict_out = {}

    for combination in valid_lis:
    
        index = 0
    
        change_list =[]
        fm = FleurinpModifier(inp_mod)

        for orbital in occupation_combinations:


            change_list.append({'state': orbital , 'spinup': combination[index] , 'spindown' : combination[index] })
            print(combination[index] , index , orbital)
        
            index += 1


        fm.set_species('all', { 'electronconfig': { 'stateoccupation': change_list}})
    
        change = fm.changes()
    
        mod_inp = fm.freeze()

    

        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
    
        print('started fleur calc')
    
        results , fleur_proc_node = run.get_node(FleurCalculation, **inputs)

        print('started fleur calc with pk:' , fleur_proc_node.pk)
        dict_out[str(combination)] = [fleur_proc_node.pk , fleur_proc_node.outputs.output_parameters.get_dict()['density_convergence']]
    
    


    return dict_out
