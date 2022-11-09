from aiida import load_profile

load_profile()


from aiida.orm import Dict , load_code

from aiida.plugins import DataFactory

from aiida.orm import load_node
from aiida.engine import submit

from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation

from aiida.plugins import CalculationFactory

inpgen_code = load_code('inpgen_MaXR5_th1@iffslurm')

structure = load_node(535)

parameters = Dict(dict={
    'comp': {
        'kmax': 3.84
        },
    'kpt': {
        'div1': 2,
        'div2' : 2,
        'div3' : 2
        }})
options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine": 1}, 'queue_name': 'th1', 'withmpi': False}


inputs = FleurinputgenCalculation.get_builder()
inputs.parameters = parameters
inputs.code = inpgen_code
inputs.structure = structure
inputs.metadata.options = options  # note options not in inputs but inputs.metadat






inpgen_process = submit(FleurinputgenCalculation, **inputs)
print('The PK of submitted inpgen job is {}'.format(inpgen_process))





