#!/usr/bin/env python

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm import Code, load_node, DataFactory
from aiida.work import run
from aiida_kkr.tools.kkr_params import kkrparams
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from pprint import pprint
from scipy import array

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

def print_clean_inouts(node):
    from pprint import pprint
    inputs = node.get_inputs_dict()
    outputs = node.get_outputs_dict()
    for key in outputs.keys():
        try:
            int(key.split('_')[-1])
            has_id = True
        except:
            has_id = False
        if 'CALL' in key or 'CREATE' in key or has_id:
            outputs.pop(key)
    print 'inputs:'
    pprint(inputs)
    print '\noutputs:'
    pprint(outputs)

alat = 6.83 # in a_Bohr
abohr = 0.52917721067 # conversion factor to Angstroem units
# bravais vectors
bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])

a = 0.5*alat*abohr
Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

Cu.store()
print(Cu)

# here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
wfd = kkr_scf_wc.get_wf_defaults()
wfd['convergence_criterion'] = 10**-6
KKRscf_wf_parameters = ParameterData(dict=wfd)

# The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
VoroCode = Code.get_from_string('voronoi@slurmcontrol')
KKRCode = Code.get_from_string('KKRcode@slurmcontrol')

# Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a ParameterData object for the use in aiida
ParaNode = ParameterData(dict=kkrparams(LMAX=2, RMAX=7, GMAX=65, NSPIN=1, RCLUSTZ=1.9).get_dict())

label = 'KKR-scf for Cu bulk'
descr = 'KKR self-consistency workflow for Cu bulk'
run(kkr_scf_wc, structure=Cu, calc_parameters=ParaNode, voronoi=VoroCode, 
    kkr=KKRCode, wf_parameters=KKRscf_wf_parameters, _label=label, _description=descr)


