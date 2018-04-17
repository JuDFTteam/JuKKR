from aiida import is_dbenv_loaded, load_dbenv
if not is_dbenv_loaded():
    load_dbenv()
from aiida.orm import load_node

n = load_node(10)
n.get_outputs()
n = n.get_outputs()[-1]
out = n.get_dict()

print 'voronoi_step_success', out['voronoi_step_success']
assert out['voronoi_step_success']

print 'kkr_step_success', out['kkr_step_success']
assert out['kkr_step_success']

print 'successful', out['successful']
assert out['successful']

print 'error', out['errors']
assert out['errors'] == []

print 'warning', out['warnings']
assert out['warnings'] == []

print 'convergence_reached', out['convergence_reached']
assert out['convergence_reached']

print 'convergence_value', out['convergence_value']
assert out['convergence_value'] < 10**-4

print 'charge_neutrality', abs(out['charge_neutrality'])
assert abs(out['charge_neutrality']) < 5*10**-4

print 'used_higher_accuracy', out['used_higher_accuracy']
assert out['used_higher_accuracy']
