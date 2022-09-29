
from aiida_fleur.data.fleurinp import FleurinpData
from aiida import load_profile
from aiida.orm import Dict, load_node , load_code
from aiida.engine import submit, run
from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation
from aiida.plugins import DataFactory, CalculationFactory
from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.data.fleurinpmodifier import FleurinpModifier
import numpy as np



load_profile()

sample_number = 10

inp_mod = load_node(954)

constraints = { '(4s1/2)' : 2 ,
               '(3d3/2)' : 4 ,
               '(3d5/2)' : 6 }


options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
           'queue_name' : 'th1-2020-32',
           'withmpi' : True,
           'max_wallclock_seconds' : 600}
fleur = load_code(150)


number_valence_electrons = 11

step = 0.5


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





from itertools import product

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
    
    


print(dict_out)





















