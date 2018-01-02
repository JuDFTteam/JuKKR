#!/usr/bin/env python

import pytest

def test_verify_results():
    import os
    from kkrparser_functions import *
    from pprint import pprint
   
    path0 = 'test_run/'
    print_all = False
   
    success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
   
    if not success or print_all:
        if not success:
            print 'WARNING parsing retured with error(s)'
            pprint(parser_msgs)
            print 
        print 'parsed out_dict:'
        pprint(out_dict)
   
    assert success
    assert out_dict['convergence_group']['rms'] < 10**-6
    assert abs(out_dict['convergence_group']['charge_neutrality']) < 10**-4

