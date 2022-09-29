from aiida_fleur.data.fleurinp import FleurinpData
from aiida import load_profile
load_profile();

from aiida.orm import Dict, load_node, load_code
from aiida.engine import submit
from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation

inpgen = load_code(152)

from aiida.plugins import DataFactory


FleurinpData = DataFactory('fleur.fleurinp')


from aiida_fleur.data.fleurinp import FleurinpData

options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine": 1}, 'queue_name': 'th1', 'withmpi': False}

inputs = FleurinputgenCalculation.get_builder()
inputs.code = inpgen
inputs.structure = struct
inputs.metadata.options = options

inpgen_process = submit(FleurinputgenCalculation, **inputs)
print('PK of inpgen job' , inpgen_process)

from aiida_fleur.calculation.fleur import FleurCalculation
from aiida.plugins import CalculationFactory
FleurCalculation = CalculationFactory('fleur.fleur')
fleur = load_code(150)



