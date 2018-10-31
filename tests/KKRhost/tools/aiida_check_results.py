#!/usr/bin/env python

# connect to db
from aiida import is_dbenv_loaded, load_dbenv
if not is_dbenv_loaded():
    load_dbenv()
# load aiida stuff
from aiida.orm import load_node
from aiida.orm  import QueryBuilder
from aiida.orm import DataFactory
from aiida.orm import WorkCalculation
# for nicer printing
from pprint import pprint

# get output dictionary quering in imported database (find structure first and then extract workflow results node)
StructureData = DataFactory('structure')
qb = QueryBuilder()
qb.append(StructureData, tag = 'struc')
qb.add_filter('struc', {'label': 'Cu_bulk_simple_test'})
structure_node = qb.get_results_dict()[0].get('struc').values()[0]
pk0 = structure_node.get_outputs(node_type=WorkCalculation)[0].pk
print 'found structure node:', structure_node
print 'Workflow pk=', pk0
# extract results node and get out dict
n = load_node(pk0)
for node in n.get_outputs():
    if node.label=='kkr_scf_wc_results':
        break
out = node.get_dict()
print '\n\noutput dictionary:\n-------------------------------------------------'
pprint(out)


# finally check some outputs
print '\n\ncheck values ...\n-------------------------------------------------'

print 'voronoi_step_success', out['voronoi_step_success']
assert out['voronoi_step_success']

print 'kkr_step_success', out['kkr_step_success']
assert out['kkr_step_success']

print 'successful', out['successful']
assert out['successful']

print 'error', out['errors'], out['errors'] == []
assert out['errors'] == []

print 'warning', out['warnings'], out['warnings'] == []
assert out['warnings'] == []

print 'convergence_reached', out['convergence_reached']
assert out['convergence_reached']

print 'convergence_value', out['convergence_value'], out['convergence_value'] < 10**-4
assert out['convergence_value'] < 10**-4

print 'charge_neutrality', abs(out['charge_neutrality']), abs(out['charge_neutrality']) < 5*10**-4
assert abs(out['charge_neutrality']) < 5*10**-4

print 'used_higher_accuracy', out['used_higher_accuracy']
assert out['used_higher_accuracy']

print '\ndone with checks\n'
