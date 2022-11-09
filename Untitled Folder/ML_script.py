from aiida import load_profile

load_profile()


from aiida.orm import Dict , load_code

from aiida.plugins import DataFactory

from aiida.orm import load_node
from aiida.engine import submit

from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation

from aiida.plugins import CalculationFactory

inpgen_code = load_code('inpgen_MaXR5_th1@iffslurm')

bohr = 0.52917721092

a = bohr * 7.497

a1 = a * 0.7071068

cell = [[a1, 0.0 ,0.0 ], [0.0, a , 0.0 ] , [0.0 , 0.0 , a1]]

structure = StructureData(cell=cell)
structure.append_atom( position = (0. ,0.0 ,0.0), symbols = 'Fe')
structure.pbc = (True , True , False ) 
Fe_structure = structure 
Fe_structure.store() 
print('PK of MonoLayer Fe: {}'.format(Fe_structure.pk))


structure = StructureData(cell=cell)
structure.append_atom( position = (0. ,0.0 ,0.0), symbols = 'Ni')
structure.pbc = (True , True , False ) 
Ni_structure = structure 
Ni_structure.store() 
print('PK of MonoLayer Ni: {}'.format(Ni_structure.pk))






structure = StructureData(cell=cell)
structure.append_atom( position = (0. ,0.0 ,0.0), symbols = 'Co')
structure.pbc = (True , True , False ) 
Co_structure = structure 
Co_structure.store() 
print('PK of MonoLayer Co: {}'.format(Co_structure.pk))



structures = [Fe_structure , Ni_structure , Co_structure]


parameters = Dict(dict={
    'comp': {
        'kmax': 3.85
        },
    'kpt': {
        'div1': 4,
        'div2' : 4,
        'div3' : 1
        }})

# options

options = {'resources' : {"num_machines": 1, "num_mpiprocs_per_machine" : 1}, 'withmpi' : False}






for structure in structures:
    # assemble inputs in a single dictionary
    inputs = FleurinputgenCalculation.get_builder()
    inputs.parameters = parameters
    inputs.code = inpgen_code
    inputs.structure = structure
    inputs.metadata.options = options

    # submit
    inpgen_pk = submit(FleurinputgenCalculation, **inputs)
    print('The PK of submitted inpgen job is {} for {} structure'.format(res_pk, structure.formula))






