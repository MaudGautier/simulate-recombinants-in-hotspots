#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                           SCRIPT DESCRIPTION                          ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## SCRIPT SUMMARY:
# Writes a file containing parameter settings selected randomly among ranges 
# given as input.


## USAGE:
#  python $SRC/2_dBGC/7_Simulations_ABC/Scripts/write_parameters_file.py \
        #  --output $param_file \
        #  --tot_nb_simul 100000 \
        #  --method $method \
        #  --tot_nb_recomb $tot_nb_recomb \
        #  --nb_frag $nb_frag \
        #  --CO_asym $CO_asym \
        #  --NCO_asym $NCO_asym \
        #  --prop_CAST $prop_CAST \
        #  --name $name \
        #  --CO_l_min $CO_l_min \
        #  --CO_l_max $CO_l_max \
        #  --CO_sd_min $CO_sd_min \
        #  --CO_sd_max $CO_sd_max \
        #  --NCO_l_min $NCO_l_min \
        #  --NCO_l_max $NCO_l_max \
        #  --NCO_sd_min $NCO_sd_min \
        #  --NCO_sd_max $NCO_sd_max \
        #  --ratio_min $ratio_min \
        #  --ratio_max $ratio_max
        #  [-h,--help]


## FORMAT OF THE OUTPUT FILE
#name  method    nb_COs  CO_tl_mean  CO_tl_sd  CO_asym  nb_NCOs  NCO_tl_mean  NCO_tl_sd  NCO_asym  prob_CAST  nb_frag
#001   observed  10      300.0       10.0      0.5      10       50.0         5.0        0.5       0.5        10
#002   TEST      2       200.0       10.0      0.0      5        40.0         4.0        1         1          10
#003   TEST      2       200.0       10.0      0.0      5        40.0         4.0        1         1          10
#004   random    10000   566.0       277.0     0.5      10000    94.0         62.0       0.5       0.5        10


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              IMPORTATIONS                             ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

from __future__ import division
import argparse
import sys
from math import *

import random
import re


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                               FUNCTIONS                               ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def create_parser():
    """ Creates the parser for arguments. 
    """
    parser = argparse.ArgumentParser(description='Writes a file containing \
            the parameters of the simulations.')
    
    # Required arguments
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-o', '--output', dest='output_file', 
            metavar='output_file', type=str, required=True,
            help="""The name of the output file that will contain the list 
            of parameters.""")
    requiredNamed.add_argument('--CO_l_min', dest='CO_l_min', 
            metavar='CO_l_min', type=float, required=True,
            help="""The minimum value for CO length.""")
    requiredNamed.add_argument('--CO_l_max', dest='CO_l_max', 
            metavar='CO_l_max', type=float, required=True,
            help="""The maximum value for CO length.""")
    requiredNamed.add_argument('--NCO_l_min', dest='NCO_l_min', 
            metavar='NCO_l_min', type=float, required=True,
            help="""The minimum value for NCO length.""")
    requiredNamed.add_argument('--NCO_l_max', dest='NCO_l_max', 
            metavar='NCO_l_max', type=float, required=True,
            help="""The maximum value for NCO length.""")
    requiredNamed.add_argument('--CO_sd_min', dest='CO_sd_min', 
            metavar='CO_sd_min', type=float, required=True,
            help="""The minimum value for CO sd length.""")
    requiredNamed.add_argument('--CO_sd_max', dest='CO_sd_max', 
            metavar='CO_sd_max', type=float, required=True,
            help="""The maximum value for CO sd length.""")
    requiredNamed.add_argument('--NCO_sd_min', dest='NCO_sd_min', 
            metavar='NCO_sd_min', type=float, required=True,
            help="""The minimum value for NCO sd length.""")
    requiredNamed.add_argument('--NCO_sd_max', dest='NCO_sd_max', 
            metavar='NCO_sd_max', type=float, required=True,
            help="""The maximum value for NCO sd length.""")
    requiredNamed.add_argument('--ratio_min', dest='ratio_CO_NCO_min', 
            metavar='ratio_CO_NCO_min', type=float, required=True,
            help="""The maximum value for NCO sd length.""")
    requiredNamed.add_argument('--ratio_max', dest='ratio_CO_NCO_max', 
            metavar='ratio_CO_NCO_max', type=float, required=True,
            help="""The maximum value for NCO sd length.""")
    requiredNamed.add_argument('--tot_nb_simul', dest='tot_nb_simul', 
            metavar='tot_nb_simul', type=int, required=True,
            help="""The total number of simulations.""")
    
    # Optional arguments
    parser.add_argument('--method', dest = 'method', 
            metavar='method', type=str, 
            default='random',
            help = """The type of method to select hotspots.""")
    parser.add_argument('--tot_nb_recomb', dest = 'tot_nb_recomb', 
            metavar='tot_nb_recomb', type=int, 
            default=10000,
            help = """The total number of recombinants (try to get as close as possible).""")
    parser.add_argument('--nb_frag', dest = 'nb_frag', 
            metavar='nb_frag', type=int, 
            default=10,
            help = """The number of fragments.""")
    parser.add_argument('--CO_asym', dest = 'CO_asym', 
            metavar='CO_asym', type=float, 
            default=0.5,
            help = """The asymetry value for COs.""")
    parser.add_argument('--NCO_asym', dest = 'NCO_asym', 
            metavar='NCO_asym', type=float, 
            default=0.5,
            help = """The asymetry value for NCOs.""")
    parser.add_argument('--prop_CAST', dest = 'prop_CAST', 
            metavar='prop_CAST', type=float, 
            default=0.5,
            help = """The prop CAST value.""")
    parser.add_argument('--name', dest = 'name', 
            metavar='name', type=str, 
            default='count',
            help = """The name to give to the simulation.""")

    return parser.parse_args()


def count_sig_figs(answer):
    str_answer = str(answer).replace('.', '')
    m = re.match( r'[0-9](\d*[0-9])?', str_answer)
    return len(m.group())


def select_value_from_distribution(min_val, max_val, distribution):
    if distribution == "uniform":
        return round(random.uniform(min_val, max_val))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
####                              MAIN SCRIPT                              ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

## PROCEDURE FOLLOWED:
# 1. Extract each parameter. (start / stop + form)
# 2. For each parameter, select one value from the distribution.
# 3. Write to output file.
# 4. Re-do X times.

if __name__ == '__main__':
    
    # 0.   Get arguments
    args = create_parser()
    
    # 0.a. Files
    output_file = args.output_file
    
    # 0.b. Parameters
    sep = "\t"
    sep_name = "."
    # Required
    get_CO_l_min = args.CO_l_min
    get_CO_l_max = args.CO_l_max
    get_NCO_l_min = args.NCO_l_min
    get_NCO_l_max = args.NCO_l_max
    get_CO_sd_min = args.CO_sd_min
    get_CO_sd_max = args.CO_sd_max
    get_NCO_sd_min = args.NCO_sd_min
    get_NCO_sd_max = args.NCO_sd_max
    get_ratio_CO_NCO_min = args.ratio_CO_NCO_min
    get_ratio_CO_NCO_max = args.ratio_CO_NCO_max
    tot_nb_simul = args.tot_nb_simul
    # Optional
    get_method = args.method
    get_tot_nb_recomb = args.tot_nb_recomb
    get_nb_frag = args.nb_frag
    get_NCO_asym = args.NCO_asym
    get_CO_asym = args.CO_asym
    get_name = args.name
    get_prop_CAST = args.prop_CAST

    # 0.c. Remind user of their chosen parameters
    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "SCRIPT: write_parameters_file.py"
    print "PARAMETERS:\n"
    print "Output file:", output_file
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    
    # 1. Prepare the output file (header)
    with open(output_file, 'w') as filout:
        header = "#name\tmethod\t"
        header += "nb_COs\tCO_tl_mean\tCO_tl_sd\tCO_asym\t"
        header += "nb_NCOs\tNCO_tl_mean\tNCO_tl_sd\tNCO_asym\t"
        header += "prob_CAST\tnb_frag\n"
        filout.write(header)

    # 2. Process each fragment and add the required pieces of information
    with open(output_file, 'a') as filout:
        for nb_simul in xrange(tot_nb_simul):
            # Choose CO and NCO l and sd among distribution
            CO_l = select_value_from_distribution(get_CO_l_min, get_CO_l_max, 'uniform')
            NCO_l = select_value_from_distribution(get_NCO_l_min, get_NCO_l_max, 'uniform')
            CO_sd = select_value_from_distribution(get_CO_sd_min, get_CO_sd_max, 'uniform')
            NCO_sd = select_value_from_distribution(get_NCO_sd_min, get_NCO_sd_max, 'uniform')
            
            # Choose ratio
            power = random.uniform(get_ratio_CO_NCO_min, get_ratio_CO_NCO_max)
            ratio = 10**power
            ## RATIO = #CO / #NCO
            ##   TOT = #CO + #NCO  ==> #NCO = TOT - #CO   // #CO = TOT - #NCO
            ## ==> TOT - #NCO = #NCO * RATIO
            ## ==> #NCO = TOT / (1 + RATIO)
            ## ==> #CO = TOT - #NCO
            nb_NCOs = int(round(get_tot_nb_recomb / (ratio + 1)))
            nb_COs = get_tot_nb_recomb - nb_NCOs
            
            # Define name
            if get_name == "explicit":
                name = "CO_" + str(nb_COs) + sep_name + \
                        "m_" + str(CO_l) + sep_name + \
                        "sd_" + str(CO_sd) + sep_name + \
                        "asym_" + str(get_CO_asym) + sep_name + \
                        "NCO_" + str(nb_NCOs) + sep_name + \
                        "m_" + str(NCO_l) + sep_name + \
                        "sd_" + str(NCO_sd) + sep_name + \
                        "asym_" + str(get_NCO_asym)
            elif get_name == "count":
                nb_zeros = count_sig_figs(tot_nb_simul) - 1
                name = nb_zeros*"0" + str(nb_simul)

            # Create lineout
            lineout = name + sep + \
                    get_method + sep + \
                    str(nb_COs) + sep + \
                    str(CO_l) + sep + \
                    str(CO_sd) + sep + \
                    str(get_CO_asym) + sep + \
                    str(nb_NCOs) + sep + \
                    str(NCO_l) + sep + \
                    str(NCO_sd) + sep + \
                    str(get_NCO_asym) + sep + \
                    str(get_prop_CAST) + sep + \
                    str(get_nb_frag) + "\n"

            # Write line to output file
            filout.write(lineout)


