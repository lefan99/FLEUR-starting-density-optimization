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
import json
import matplotlib.pyplot as plt
import numpy as np

load_profile()
####################################################################################3







def gen_inp_itmax( structure , inpgen = load_code(149), itmax = 1 ): # IN: structure data , number of wanted iterations OUT: input node of structre with itmax
    
    FleurinpData = DataFactory('fleur.fleurinp')
    FleurCalculation = CalculationFactory('fleur.fleur')

    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine": 1}, 'queue_name': 'th1-2020-32', 'withmpi': False}






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

####################################################################################3





#ONLY WORKS FOR UNARIES

def starting_density_calc(inp , fleur = load_code(150) ,  step=0.5): #IN: input node fleur: fleur code : step: grid size 0.01 means (0,01 , 0.02 , 0.03 etc)
    
    
    
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
                    '(1f5/2)' : 6 , 
                   '(1f7/2)' : 8 ,
                   '(2f5/2)' : 6 , 
                   '(2f7/2)' : 8 ,
                   '(3f5/2)' : 6 , 
                   '(3f7/2)' : 8 ,
                   '(4f5/2)' : 6 , 
                   '(4f7/2)' : 8 ,
                   '(5f5/2)' : 6 , 
                   '(5f7/2)' : 8 ,
                   '(6f5/2)' : 6 , 
                   '(6f7/2)' : 8 ,
                  
                
                  }

    number_atoms = len(inp_mod.inp_dict['atomGroups'][0]['labels'])  # how many atoms in unitcell
    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}



    number_valence_electrons = inp_mod.inp_dict['cell']['bzIntegration']['valenceElectrons']



    valence_orbs = inp_mod.inp_dict['atomSpecies'][0]['electronConfig']['valenceConfig'] #get all pre defined valence orbitals



    orbital_occupation = []

    occupation_combinations = {}

    for orbital in valence_orbs:                     #this part calculates all the possible combinations for the given gridsize and boundaries ( orbital size valence electrons
    
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

        if sum(combination)*2 * number_atoms  == number_valence_electrons:

            valid_lis.append(combination)


    print(valid_lis , 'number of possible combinations; ' , len(valid_lis))

    
    
    dict_out = {}
    bravais_matrix = inp.inp_dict['cell']['bulkLattice']['bravaisMatrix']
    atom_number = inp.inp_dict['atomSpecies'][0]['atomicNumber']

    
    
    
    for combination in valid_lis:
    
        index = 0
    
        change_list =[]
        fm = FleurinpModifier(inp_mod)
        
        combination_dict = {}

        for orbital in occupation_combinations:
            print(combination , index)


            change_list.append({'state': orbital , 'spinup': combination[index] , 'spindown' : combination[index] })
     #       print(combination[index] , index , orbital)
        
            index += 1
            
            
            
            #combination_dict[orbital] = combination[index]
            

        fm.set_species('all', { 'electronconfig': { 'stateoccupation': change_list}})
    
    #    change = fm.changes()
    
        mod_inp = fm.freeze()

    

        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
    

    
        fleur_proc_node = submit(FleurCalculation, **inputs)

        print('started fleur calc with pk:' , fleur_proc_node.pk)   
        dict_out[str(combination)] = fleur_proc_node
    
        fleur_proc_node.set_extra('configuration' , combination)
        fleur_proc_node.set_extra('atomNumber' , atom_number)
        fleur_proc_node.set_extra('Bravais' , bravais_matrix)

    return dict_out , bravais_matrix , atom_number , valence_orbs #dict out is a dictionary with orbital occupation as keys and the pk of the fleur calc as an item

####################################################################################3







def get_results( input_dict ): #after all the calculations are done use this function to retrieve the distances using the dict_out from the previous function
    
    results_dict = {}
    
    for orbital_combination in input_dict:

        results_dict[orbital_combination] = input_dict[orbital_combination].outputs.output_parameters.get_dict()['density_convergence']
        
    return results_dict








####################################################################################3




def num_possible_comb(inp , step = 0.1 , show_comb = False):   # calculate the number of possible calculations given an input and grid size 
    inp_mod = inp
    
    number_atoms = len(inp_mod.inp_dict['atomGroups'][0]['labels'])
    
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
                   '(1f5/2)' : 6 , 
                   '(1f7/2)' : 8 ,
                   '(2f5/2)' : 6 , 
                   '(2f7/2)' : 8 ,
                   '(3f5/2)' : 6 , 
                   '(3f7/2)' : 8 ,
                   '(4f5/2)' : 6 , 
                   '(4f7/2)' : 8 ,
                   '(5f5/2)' : 6 , 
                   '(5f7/2)' : 8 ,
                   '(6f5/2)' : 6 , 
                   '(6f7/2)' : 8 ,

                
                  }


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}


    number_valence_electrons = inp_mod.inp_dict['cell']['bzIntegration']['valenceElectrons']
    print(number_valence_electrons)


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
    #print(all_comb)
    valid_lis = []

    for combination in all_comb:
        #print(sum(combination)*2)


        if sum(combination)*2  *number_atoms == number_valence_electrons:
            
            
            valid_lis.append(combination)


    print('The number of all possible combinations is: ', len(valid_lis))
    
    if show_comb == True:
        
        return valid_lis
    


##################################################################################################

#doesnt work properly lol
def save_data( dif_dict , atom_number , bravais_matrix , valence_orbs, filename = 'test.json'):
    
    file = open(filename , 'w')
    
    file_data = json.load(file)
        
    atom_dict = {}
    
    bravais_dict = {}
    
    valence_dict  ={}

    valence_dict[str(valence_orbs)] = dif_dict

    bravais_dict[str(bravais_matrix)] = valence_dict
    
    atom_dict[str(atom_number)] = bravais_dict
    
    file_data.append(atom_dict)
    
    json.dump(file_data, file) 
    
    
########################################################################################################


# plot the data from the unaries using the saved json data for each orbital
def plot_unaries(results):
    
    for atom in results:
        
        results1 = results[atom]
        
        for bravais in results1:
            
            results2 = results1[bravais]
            
            for orbit in results2:
                
                orbits = eval(orbit) 
                
                number_orbits = len(orbits)
                
                results3 = results2[orbit]
                print(results3)
                
                
                for i in range(len(orbits)):
                    
                    f7 = [] 
                    distance = []

                    for occupation in results3:
    
                        f7.append(eval(occupation)[i])
                        distance.append(results3[occupation])
    
                    plt.plot(f7 , distance, 'ro')
                    plt.xlabel('orbital occupation')
                    plt.ylabel('distance')
                    plt.title('atom number: ' + str(atom) + 'orbital: ' + str(orbits[i]) )
                    plt.show()













#################################################################################################
#idk lol
def compare(results_dict , standard_dict ):
    
    return_dict ={}
    
    for atom in stadard_dict:
        
        standard_distance = standard_dict[atom]
        
        for bravais in results_dict[atom]: 
            
            for orbitals in results_dict[atom][bravais]:
                
                for occupations in results_dict[atom][bravais][orbitals]: 
                    
                    return_dict[occupations + atom]


##################################################################################################
#doesnt work 

def minimize_distance( list_inps , steps = 0.1 ):
    
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
                    '(1f5/2)' : 6 , 
                   '(1f7/2)' : 8 ,
                   '(2f5/2)' : 6 , 
                   '(2f7/2)' : 8 ,
                   '(3f5/2)' : 6 , 
                   '(3f7/2)' : 8 ,
                   '(4f5/2)' : 6 , 
                   '(4f7/2)' : 8 ,
                   '(5f5/2)' : 6 , 
                   '(5f7/2)' : 8 ,
                   '(6f5/2)' : 6 , 
                   '(6f7/2)' : 8 ,
                  
                
                  }


    sc_loop_running = True


    for inputs in list_inps:


        number_atoms = len(inputs.inp_dict['atomGroups'][0]['labels'])
        options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}



        number_valence_electrons = inputs.inp_dict['cell']['bzIntegration']['valenceElectrons']



        valence_orbs = inputs.inp_dict['atomSpecies'][0]['electronConfig']['valenceConfig']



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

            if sum(combination)*2 * number_atoms  == number_valence_electrons:

                valid_lis.append(combination)


        print(valid_lis , 'number of possible combinations; ' , len(valid_lis))

        
        
        dict_out = {}
        bravais_matrix = inp.inp_dict['cell']['bulkLattice']['bravaisMatrix']
        atom_number = inp.inp_dict['atomSpecies'][0]['atomicNumber']


        fleur_standard = inputs.inp_dict['atomSpecies'][0]['electronConfig']['stateOccupation']

        partial_occupation = [] 

        for dicts in fleur_standard:

            partial_occupation.append(dicts['state'])


        for orbitals in fleur_standard:

            if orbitals['spinUp'] != orbtials['spinDown']:

                print("error magnetic structure , spins do not match")




        fleur_standard_occupation = []

        for orbital in valence_orbs:

            if orbital in partial_occupation:

                for i in range(len(partial_occupation)):

                    if partial_occupation[i] == orbital:

                        fleur_standard_occupation.append(fleur_standard[i]['spinUp'])

            else: 

                fleur_standard_occupation.append(constraints[orbital])
                
        



        valid_lis2 = []

        for combination in valid_lis:

            valid_lis2.append(list(combination))

        if fleur_standard_occupation not in valid_lis2:

            valid_lis2.append(fleur_standard_occupation)

        valid_lis2.sort()

        for calc_index in range(len(valid_lis2)):

            if valid_lis2[calc_index] == fleur_standard_occupation:

                index = calc_index


                
        change_list =[]

        fm = FleurinpModifier(inputs)
        
        k = 0

        for orbital in valence_orbs:

        

            print(combination , index)


            change_list.append({'state': orbital , 'spinup': valid_lis2[index][k] , 'spindown' : valid_lis2[index][k] })
        
            k += 1
            
            
            

        fm.set_species('all', { 'electronconfig': { 'stateoccupation': change_list}})
    
    
        mod_inp = fm.freeze()

    

        inp = FleurCalculation.get_builder()
        inp.code = fleur
        inp.fleurinpdata = mod_inp
        inp.metadata.options = options 

        fleur_node = submit(FleurCalculation , **inp)

        fleur_node.set_extra('config' , valid_lis2[index])
        
        calc_nodes_list.append( fleur_node )
        bravais_list.append(bravais_matrix) 
        atom_number_list.append(atom_number)
        valence_orbs_list.append(valence_orbs)

    
    while sc_loop_running:

        if first:

            for process in range(len(calc_nodes_list)):

                if calc_nodes_list[process].is_finished:

                    distance = calc_nodes_list[process].outputs.output_parameters.get_dict()['density_convergence']

                    distance_list[process] = distance

                    combination_up = index_list[process] +  1
                    combination_down = index_list[process] -1

                    
                  #  for orbital in valence_orbs_list[process]:
#    list_running = []
#
#    proc_node_low = []
#    proc_node_high = []
#
#    occupation_low = []
#    occupation_high = []
#
#    proc_node_first = []
#    occupation_first = []
#
#    first = True
#
#
#    for inputs in list_inps:
#
#    list_finish.append(False)
#        
#        number_atoms = len(inputs.inp_dict['atomGroups'][0]['labels'])
#        
#        options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
#               'queue_name' : 'th1-2020-32',
#               'withmpi' : True,
#               'max_wallclock_seconds' : 600}
#
#
#
#        number_valence_electrons = inputs.inp_dict['cell']['bzIntegration']['valenceElectrons']
#    
#
#
#
#        valence_orbs = inputs.inp_dict['atomSpecies'][0]['electronConfig']['valenceConfig']
#
#        
#
#        state_occupation = inputs.inp_dict['atomSpecies'][0]['electronConfig']['stateOccupation']
#        
#        partial_occupation = []
#        
#        for dicts in state_occupation:
#            
#            
#            partial_occupation.append(dicts['state']) 
#
#
#        for orbitals in valence_orbs:
#            
#            if orbitals not in partial_occupation:
#                
#                state_occupation.append({'state': orbitals, 'spinUp': constraints[orbitals]/2 , 'spinDown': constraints[orbitals]/2})
#                                        
#        
#        fm = FleurinpModifier(inputs)
#        fm.set_species('all', { 'electronconfig': { 'stateoccupation': state_occupation}})
#        mod_inp = fm.freeze()
#                              
#                                        
#        
#
#        inputs = FleurCalculation.get_builder()
#        inputs.code = fleur
#        inputs.fleurinpdata = mod_inp
#        inputs.metadata.options  = options 
#        
#        
#
#    
#        fleur_proc_node = submit(FleurCalculation, **inputs)
#
#        print('started fleur calc with pk:' , fleur_proc_node.pk)
#        dict_out[str(combination)] = fleur_proc_node
#    
#        fleur_proc_node.set_extra('atomNumber' , atom_number)
#        fleur_proc_node.set_extra('Bravais' , bravais_matrix)
#        
#
#        proc_list.append(fleur_proc_node , state_occupation)
#
#
#    while sc_loop_finish:
#
#        if first:
#
#            for i in range(len(list_running)):
#
#                if list_running.finished[i]:
#
#                    distance = proc_node_first[i].outputs.output_parameters.get_dict()['density_convergence'] 
#
#                    distance_list.append(distance)
#
#                    step_up = state_occupation
#                    step_down = state_occupation
#
#                    for j in range(len(state_occupation)):
#                            
#                        if state_occupation[j]['spinUp'] <= constraints[state_occupation[j]['state']:
#                                                                            
#                            step_up[j]['spinUp'] += step
#                            step_up[j]['spinDown'] += step   
#                            break
#                                                                            
#                    fm = FleurinpModifier(list_inps[i])
#                    fm.set_species('all', { 'electronconfig': { 'stateoccupation': step_up}})
#                    mod_inp = fm.freeze()
#                    inputs = FleurCalculation.get_builder()
#                    inputs.code = fleur
#                    inputs.fleurinpdata = mod_inp
#                    inputs.metadata.options  = options                                    
#                                                                            
#                    fleur_proc_node = submit(FleurCalculation, **inputs)  
#                                                                            
#                    list_running_up[i] = fleur_proc_node                                                   
#                                                                            
#                    for k in range(len(state_occupation)):
#                            
#                        if state_occupation[k]['spinUp'] <= 0:
#                                                                            
#                            step_down[j]['spinUp'] -= step
#                            step_down[j]['spinDown'] -= step                           
#                            break
#                                                            
#                    fm = FleurinpModifier(list_inps[i])
#                    fm.set_species('all', { 'electronconfig': { 'stateoccupation': step_down}})
#                    mod_inp = fm.freeze()
#                    inputs = FleurCalculation.get_builder()
#                    inputs.code = fleur
#                    inputs.fleurinpdata = mod_inp
#                    inputs.metadata.options  = options                                    
#                                                                            
#                    fleur_proc_node = submit(FleurCalculation, **inputs)  
#                                                                            
#                    list_running_down[i] = fleur_proc_node                                                             
#                        
#    
#        for i in range(len(list_running_up)):
#
                        


#################################################################################################

import random

def p( l , per ): #permutation function for a given list or array l with a given permutation per
    
    return [l[i] for i in per]

def p_r(l , per):  # reverse the permutationfunction to get the orgianl array : p_r(p(list , per) , per) = list
    
    temp = list(range(len(per)))
    
    for i in range(len(temp)):
        temp[per[i]] = l[i]

    return temp

def random_grid( inp , samples  ,  fleur = load_code(150)): #this function uses a random grid approach for starting occupations instead of a grid of normalized size , randomly samples from possible combinations while respecting boundaries like valence electrons and orbital size
    
    con_file = open('constraints.json' , 'r')
    
    constraints_dict = json.load(con_file)
    
    number_atoms = len(inp.inp_dict['atomGroups'][0]['labels'])
    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}



    number_valence_electrons = inp.inp_dict['cell']['bzIntegration']['valenceElectrons']

    number_atoms = len(inp.inp_dict['atomGroups'][0]['labels'])

    valence_orbs = inp.inp_dict['atomSpecies'][0]['electronConfig']['valenceConfig']

    print(valence_orbs)
    
    constraints = []
    
    for orbs in valence_orbs:
    
        constraints.append(constraints_dict[orbs]/2)
        
    #print(constraints)
    
    
    configurations = [] 
    
    if 'stateOccupation' not in inp.inp_dict['atomSpecies'][0]['electronConfig']:
        samples = 1
    
    for i in range(samples): 
        
        sample = []
        
        rest = number_valence_electrons/ (2 * number_atoms)
        
        permutation = list(range(len(constraints)))
        
        random.shuffle(permutation)
        
        #get random permutation in order to shuffle the orbitals, the programm starts with a random orbital and picks a random occupation and then picks a random for the next one. without shuffling the occupation for the latter orbitals would be screwed due to boundaries
        
    
        #constr_p = constraints
        constr_p = p(constraints , permutation)
        
        print(permutation , constr_p)
        
        for j in range(len(constraints)): #this section makes sure that the occupations are picked with respect to the boundaries given by orbital size and valence electrons
            
            occ = 0
            
            if rest == 0:
                occ = 0
                
            elif j == len(constraints)-1:
                
                occ = rest
                        
            elif constr_p[j] < rest:
                
                if sum(constr_p[(j+1):]) < rest :
                    
                    occ = np.random.uniform( rest - sum(constr_p[(j+1):]) , constr_p[j])
                                                         
                          
                else:
            
                    occ = np.random.uniform(0 , constr_p[j])
    
    
            else:
            
                if sum(constr_p[(j+1):]) < rest :
                    
                    occ = np.random.uniform( rest - sum(constr_p[(j+1):]) , rest)
                                                      

                else:
            
                    occ = np.random.uniform(0 , rest)
    
    
            sample.append(occ)
            rest -= occ
            
     
            
        configurations.append(p_r(sample , permutation))
    
    
        
    
    bravais_matrix = inp.inp_dict['cell']['bulkLattice']['bravaisMatrix']
    atom_number = inp.inp_dict['atomSpecies'][0]['atomicNumber']    
    dict_out = {}
    for combination in configurations:
    
        index = 0
    
        change_list =[]
        fm = FleurinpModifier(inp)
        
        combination_dict = {}

        for orbital in valence_orbs:
            
            print(combination , index)


            change_list.append({'state': orbital , 'spinup': combination[index] , 'spindown' : combination[index] })
     #       print(combination[index] , index , orbital)
        
            index += 1
            
            
            
            #combination_dict[orbital] = combination[index]
            

        fm.set_species('all', { 'electronconfig': { 'stateoccupation': change_list}})
    
    #    change = fm.changes()
    
        mod_inp = fm.freeze()

    

        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
    

    
        fleur_proc_node = submit(FleurCalculation, **inputs)

        print('started fleur calc with pk:' , fleur_proc_node.pk)
        dict_out[str(combination)] = fleur_proc_node
    
        fleur_proc_node.set_extra('configuration' , combination)
        fleur_proc_node.set_extra('atomNumber' , atom_number)
        fleur_proc_node.set_extra('Bravais' , bravais_matrix)

    return dict_out , bravais_matrix , atom_number , valence_orbs    #dict_out is a dic with occupations as keys and the calc node of the corresponding calculation as an item
    
    
    






    
        

def oxides_charge1(inp_init ,fleur = load_code(150) , gridsize = 20 ):
    
    
    fleurmode = FleurinpModifier(inp_init)
    dict_change = { 
    'itmax' : 1 }

    fleurmode.set_inpchanges(dict_change)
    inp = fleurmode.freeze()
    dic = inp.inp_dict

    node_list = []
    


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}

    O_occ = 2
    X_occ = 1

    
    
    
    
    for species in dic['atomSpecies']:
        
        if species['name'] != 'Oxygen (O)':
            
            X_name = species['name']
            
            X_orb = species['electronConfig']['valenceConfig'][-1]
            
    for i in range(gridsize):
    
        fm = FleurinpModifier(inp)

        charge = 1/20
    
        changes = []
    

    
            
    
        #print(charge , Mg_occ , O_occ)
            
        changes.append({'state': '(2p1/2)' , 'spinup': O_occ/3 , 'spindown' : O_occ/3 })
        changes.append({'state': '(2p3/2)',
         'spinUp': O_occ * 2 / 3 ,
         'spinDown': O_occ * 2 / 3 })
     
        fm.set_species('Oxygen (O)', { 'electronconfig': { 'stateoccupation': changes}})
    
    
        changes = []
        changes.append({'state': X_orb, 'spinup': X_occ , 'spindown' : X_occ })
#    changes.append({'state': '(3p1/2)', 'spinUp': 0.0, 'spinDown': 0.0})
            
        fm.set_species(X_name, { 'electronconfig': { 'stateoccupation': changes}})
        
        
                   
        mod_inp = fm.freeze()
        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
        
        
        fleur_proc_node = submit(FleurCalculation, **inputs)
    
        print('started fleur calc with pk:' , fleur_proc_node.pk)
    
    
        node_list.append(fleur_proc_node) 
            
        O_occ += charge
        X_occ -= charge        

        
    return node_list


def oxides_charge2(inp_init ,fleur = load_code(150) , gridsize = 20 ):
    
    
    fleurmode = FleurinpModifier(inp_init)
    dict_change = { 
    'itmax' : 1 }

    fleurmode.set_inpchanges(dict_change)
    inp = fleurmode.freeze()
    dic = inp.inp_dict

    node_list = []
    


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}

    O_occ = 2
    X_occ = 0.5

    
    
    
    
    for species in dic['atomSpecies']:
        
        if species['name'] != 'Oxygen (O)':
            
            X_name = species['name']
            
            X_orb = species['electronConfig']['valenceConfig'][-1]
            
    for i in range(gridsize):
    
        fm = FleurinpModifier(inp)

        charge = 1/20
    
        changes = []
    

    
            
    
        #print(charge , Mg_occ , O_occ)
            
        changes.append({'state': '(2p1/2)' , 'spinup': O_occ/3 , 'spindown' : O_occ/3 })
        changes.append({'state': '(2p3/2)',
         'spinUp': O_occ * 2 / 3 ,
         'spinDown': O_occ * 2 / 3 })
     
        fm.set_species('Oxygen (O)', { 'electronconfig': { 'stateoccupation': changes}})
    
    
        changes = []
        changes.append({'state': X_orb, 'spinup': X_occ , 'spindown' : X_occ })
#    changes.append({'state': '(3p1/2)', 'spinUp': 0.0, 'spinDown': 0.0})
            
        fm.set_species(X_name, { 'electronconfig': { 'stateoccupation': changes}})
        
        
                   
        mod_inp = fm.freeze()
        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
        
        
        fleur_proc_node = submit(FleurCalculation, **inputs)
    
        print('started fleur calc with pk:' , fleur_proc_node.pk)
    
    
        node_list.append(fleur_proc_node) 
            
        O_occ += charge
        X_occ -= charge/2        

        
    return node_list



def oxides_charge3(inp_init ,fleur = load_code(150) , gridsize = 20 ):
    
    
    fleurmode = FleurinpModifier(inp_init)
    dict_change = { 
    'itmax' : 1 }

    fleurmode.set_inpchanges(dict_change)
    inp = fleurmode.freeze()
    dic = inp.inp_dict

    node_list = []
    


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}

    O_occ = 2
    X_occ = 0.5

    
    
    
    
    for species in dic['atomSpecies']:
        
        if species['name'] != 'Oxygen (O)':
            
            X_name = species['name']
            
            X_orb1 = species['electronConfig']['valenceConfig'][-1]
            
            X_orb2 = species['electronConfig']['valenceConfig'][-2]
            
            
            
    for i in range(gridsize):
    
        fm = FleurinpModifier(inp)

        charge = 1/20
    
        changes = []
    

    
            
    
        #print(charge , Mg_occ , O_occ)
            
        changes.append({'state': '(2p1/2)' , 'spinup': O_occ/3 , 'spindown' : O_occ/3 })
        changes.append({'state': '(2p3/2)',
         'spinUp': O_occ * 2 / 3 ,
         'spinDown': O_occ * 2 / 3 })
     
        fm.set_species('Oxygen (O)', { 'electronconfig': { 'stateoccupation': changes}})
    
    
        changes = []
        changes.append({'state': X_orb1, 'spinup': 0 , 'spindown' : 0 })
        changes.append({'state': X_orb2, 'spinup': X_occ , 'spindown' : X_occ })
#    changes.append({'state': '(3p1/2)', 'spinUp': 0.0, 'spinDown': 0.0})
            
        fm.set_species(X_name, { 'electronconfig': { 'stateoccupation': changes}})
        
        
                   
        mod_inp = fm.freeze()
        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
        
        
        fleur_proc_node = submit(FleurCalculation, **inputs)
    
        print('started fleur calc with pk:' , fleur_proc_node.pk)
    
    
        node_list.append(fleur_proc_node) 
            
        O_occ += charge
        X_occ -= charge/2        

        
    return node_list


def oxides_charge4(inp_init ,fleur = load_code(150) , gridsize = 20 ):
    
    
    fleurmode = FleurinpModifier(inp_init)
    dict_change = { 
    'itmax' : 1 }

    fleurmode.set_inpchanges(dict_change)
    inp = fleurmode.freeze()
    dic = inp.inp_dict

    node_list = []
    


    options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 2},
               'queue_name' : 'th1-2020-32',
               'withmpi' : True,
               'max_wallclock_seconds' : 600}

    O_occ = 2
    X_occ = 0.5

    
    
    
    
    for species in dic['atomSpecies']:
        
        if species['name'] != 'Oxygen (O)':
            
            X_name = species['name']
            
            X_orb1 = species['electronConfig']['valenceConfig'][-1]
            
            X_orb2 = species['electronConfig']['valenceConfig'][-2]
            
            
            
    for i in range(gridsize):
    
        fm = FleurinpModifier(inp)

        charge = 1/20
    
        changes = []
    

    
            
    
        #print(charge , Mg_occ , O_occ)
            
        changes.append({'state': '(2p1/2)' , 'spinup': O_occ/3 , 'spindown' : O_occ/3 })
        changes.append({'state': '(2p3/2)',
         'spinUp': O_occ * 2 / 3 ,
         'spinDown': O_occ * 2 / 3 })
     
        fm.set_species('Oxygen (O)', { 'electronconfig': { 'stateoccupation': changes}})
    
    
        changes = []
        changes.append({'state': X_orb1, 'spinup': 0 , 'spindown' : 0 })
        changes.append({'state': X_orb2, 'spinup': X_occ , 'spindown' : X_occ })
#    changes.append({'state': '(3p1/2)', 'spinUp': 0.0, 'spinDown': 0.0})
            
        fm.set_species(X_name, { 'electronconfig': { 'stateoccupation': changes}})
        
        
                   
        mod_inp = fm.freeze()
        inputs = FleurCalculation.get_builder()
        inputs.code = fleur
        inputs.fleurinpdata = mod_inp
        inputs.metadata.options  = options 
        
        
        fleur_proc_node = submit(FleurCalculation, **inputs)
    
        print('started fleur calc with pk:' , fleur_proc_node.pk)
    
    
        node_list.append(fleur_proc_node) 
            
        O_occ += charge
        X_occ -= charge/2        

        
    return node_list