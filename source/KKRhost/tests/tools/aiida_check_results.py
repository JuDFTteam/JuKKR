#!/usr/bin/env python

from aiida import is_dbenv_loaded, load_dbenv
if not is_dbenv_loaded():
    load_dbenv()
from aiida.orm import load_node
from pprint import pprint

# load node of workflow
try:
   n = load_node(11)
   print '\noutputs of workflow\n-------------------------------------------------'
   pprint(n.get_outputs_dict())
   n = n.get_outputs()[-1]
except:
   n = load_node(10)
   print '\noutputs of workflow\n-------------------------------------------------'
   pprint(n.get_outputs_dict())
   n = n.get_outputs()[-1]

# get output dictionary
out = n.get_dict()
print '\n\noutput dictionary:\n-------------------------------------------------'
pprint(out)


Cu.label = 'Cu_bulk_simple_test'
from aiida.orm.QueryBuilder import QueryBuilder
from aiida.orm  import QueryBuilder
from aiida.orm import DataFactory
StructureData = DataFactory('structure')
qb = QueryBuilder()
qb.append(StructureData, tag = 'struc')
qb.add_filter('struc', {'label': 'Cu_bulk_simple_test'})
structure_node = qb.get_results_dict()[0].get('struc').values()[0]
from aiida.orm import WorkCalculation
pk0 = structure_node.get_outputs(node_type=WorkCalculation)[0].pk
print 'pk= pk0'

# finally check some output
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
