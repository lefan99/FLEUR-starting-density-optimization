from aiida_fleur.data.fleurinp import FleurinpData
from aiida import load_profile
from aiida.orm import Dict, load_node , load_code
from aiida.engine import submit, run
from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation
from aiida.plugins import DataFactory, CalculationFactory
from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.data.fleurinpmodifier import FleurinpModifier


#this code uses a struct node and creates an inpgen with itmax = 1






load_profile() 

FleurinpData = DataFactory('fleur.fleurinp')
FleurCalculation = CalculationFactory('fleur.fleur')

options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine": 1}, 'queue_name': 'th1', 'withmpi': False}







inpgen = load_code(152)
struct = load_node(688)


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
    'itmax' : 1 }

fleurmode.set_inpchanges(dict_change)
inp_mod = fleurmode.freeze()
print('PK of itmax:1 modified inp {}'.format(inp_mod.pk))




















