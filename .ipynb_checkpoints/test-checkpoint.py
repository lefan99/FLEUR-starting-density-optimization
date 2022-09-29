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

    def define(cls, spec):

        spec.in



    











































