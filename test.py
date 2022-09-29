from aiida_fleur.data.fleurinp import FleurinpData
from aiida import load_profile
from aiida.orm import Dict, load_node , load_code
from aiida.engine import submit 
from aiida_fleur.calculation.fleurinputgen import FleurinputgenCalculation
from aiida.plugins import DataFactory, CalculationFactory
from aiida_fleur.calculation.fleur import FleurCalculation
from aiida_fleur.data.fleurinpmodifier import FleurinpModifier

load_profile()

FleurinpData = DataFactory('fleur.fleurinp')
FleurCalculation = CalculationFactory('fleur.fleur')

class FleurSDWorkchain(WorkChain):
 

    _default_options = {
        'optimize_resources': True,
        'resources': {
            'num_machines': 1,
            'num_mpiprocs_per_machine': 1
        },
        'max_wallclock_seconds': 6 * 60 * 60,
        'queue_name': '',
        'custom_scheduler_commands': '',
        'import_sys_environment': False,
        'environment_variables': {}
    }

    @classmethod
    def define(cls, spec):

        super().define(spec)
        spec.input('fleur', valid_type=Code, required=True)
        spec.input('inpgen', valid_type=Code, required=False)
        spec.input('wf_parameters', valid_type=Dict, required=False)
        spec.input('structure', valid_type=StructureData, required=False)
        spec.input('calc_parameters', valid_type=Dict, required=False)
        spec.input('fleurinp', valid_type=FleurinpData, required=False)
        spec.input('remote_data', valid_type=RemoteData, required=False)
        spec.input('options', valid_type=Dict, required=False)
        spec.input('settings', valid_type=Dict, required=False)
        spec.input('settings_inpgen', valid_type=Dict, required=False)
        spec.outline(cls.start, cls.inpgen , while_(cls.condition)(cls.change_fleurinp, cls.run_fleur , cls.get_result) , cls.return_results)

        spec.output('fleurinp', valid_type=FleurinpData)
        spec.output('distances', valid_type=Dict)
        spec.output('last_fleur_calc_output', valid_type=Dict)
        spec.output('fleur_calc_pk' , valid_type=list)
        spec.expose_outputs(FleurBaseWorkChain, namespace='last_calc')
        
    
    def start(self):
        
        self.report('began starting density workdflow')
        self.ctx.calc_cout = 0
        self.ctx.calcs = []
        
        options = self._default_options.copy()
        self.ctx.options = options
        
    def inpgen(self):
        
        structure =self.inputs.sructure
        self.ctx.formula = structure.get_formula
        label = 'SD workchain inpgen'
        inpgencode = self.inputs.inpgen
        
        options = {
            'max_wallclock_seconds': int(self.ctx.options.get('max_wallclock_seconds')),
            'resources': self.ctx.options.get('resources'),
            'queue_name': self.ctx.options.get('queue_name', '')
        }        
        
        

        inputs_build = get_inputs_inpgen(structure,
                                         inpgencode,
                                         options,
                                         label)
        
        self.report('INFO: run inpgen')
        future = self.submit(inputs_build)

        return ToContext(inpgen=future)

    def change_fleurinp(self):
        
        self.report('info: changing fleur inp')
        
        fleurin = self.ctx = self.ctx['inpgen'].outputs.fleurinpData
        
        fleurmode = FleurinpModifier(fleurin)
        itmax = 1
        











































