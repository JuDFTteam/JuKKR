#!/usr/bin/env python

import pytest
import os, pprint
from kkrparser_functions import parse_kkr_outputfile

class Test_check_test_runs():
    """
    check results of the tests
    """

    def test_verify1_Au_bulk(self):
        path0 = './test_run1/'
        success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
	pprint.pprint(parser_msgs)
	pprint.pprint(out_dict)
        assert success
        assert out_dict['convergence_group']['rms'] < 10**-8
        assert abs(out_dict['convergence_group']['charge_neutrality']) <= 10**-6


    def test_verify2_Fe_slab(self):
        path0 = './test_run2/'
        success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
	pprint.pprint(parser_msgs)
	pprint.pprint(out_dict)
        assert success
        assert out_dict['convergence_group']['rms'] < 10**-8
        assert abs(out_dict['convergence_group']['charge_neutrality']) <= 10**-6


    def test_verify3_Si_lloyd(self):
        path0 = './test_run3/'
        success, parser_msgs, out_dict = parse_kkr_outputfile({}, path0+'out_kkr', path0+'output.0.txt', path0+'output.000.txt', path0+'out_timing.000.txt', path0+'out_potential', path0+'nonco_angle_out.dat')
	pprint.pprint(parser_msgs)
	pprint.pprint(out_dict)
        assert success
        assert out_dict['convergence_group']['rms'] < 10**-8
        assert abs(out_dict['convergence_group']['charge_neutrality']) <= 10**-6
