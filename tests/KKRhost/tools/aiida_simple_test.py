#!/usr/bin/env python

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm import Code, load_node, DataFactory
from aiida.work import run
from masci_tools.io.kkr_params import kkrparams
from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
from pprint import pprint
from numpy import array

ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

# make sure running the workflow exists after at most 5 minutes
import timeout_decorator
@timeout_decorator.timeout(300, use_signals=False)
def run_timeout(builder):
    from aiida.work.launch import run
    out = run(builder)
    return out

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
Cu.label = 'Cu_bulk_simple_test'

Cu.store()
print(Cu)

# here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
wfd = kkr_scf_wc.get_wf_defaults()

wfd['convergence_criterion'] = 10**-4
wfd['check_dos'] = False 
wfd['kkr_runmax'] = 5
wfd['nsteps'] = 50 
wfd['num_rerun'] = 2
wfd['natom_in_cls_min'] = 20

options = {'queue_name' : '', 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'use_mpi' : True, 'custom_scheduler_commands' : ''}
options = ParameterData(dict=options)


KKRscf_wf_parameters = ParameterData(dict=wfd)

# The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
VoroCode = Code.get_from_string('voronoi@slurmcontrol')
KKRCode = Code.get_from_string('KKRcode@slurmcontrol')
# for local tests:
#VoroCode = Code.get_from_string('voronoi@localhost')
#KKRCode = Code.get_from_string('KKRhost@localhost')

# Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a ParameterData object for the use in aiida
ParaNode = ParameterData(dict=kkrparams(LMAX=2, RMAX=7, GMAX=65, NSPIN=1, RCLUSTZ=1.9).get_dict())

label = 'KKR-scf for Cu bulk'
descr = 'KKR self-consistency workflow for Cu bulk'
try:
   # setup workflow process
   builder = kkr_scf_wc.get_builder()
   builder.calc_parameters = ParaNode
   builder.voronoi = VoroCode
   builder.kkr = KKRCode
   builder.structure = Cu
   builder.wf_parameters = KKRscf_wf_parameters
   builder.options = options
   builder.label = label
   builder.description = descr
   # now run calculation
   out = run_timeout(builder)
except:
   print 'some Error occured in run of kkr_scf_wc'


