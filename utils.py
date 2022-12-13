#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022

@author: sshanker
"""

import importlib
from colorama import Fore, Style

def flag_validator(**kw):
    procede_after_after_flag_check = 1 # dfault to run the code
    
    # check openmm    
    if int(kw['minimize']) > 0:
        if importlib.find_loader('openmm') == None:
            print(f"{Fore.RED}OpenMM not found. Please install OpenMM or remove (or set value to 0) -nmin flag.{Style.RESET_ALL}")
            procede_after_after_flag_check = 0
        

    return procede_after_after_flag_check


def add_open_mm_flags(parser):
    parser.add_argument("-nmin", "--omm_nmin", type=int, default=0,
                       dest="minimize", help='Maximum number of poses for openmm minimization and energy calculations. Default is 0')
    
    parser.add_argument("-ra", "--omm_rearrange", type=int, default=1,
                       dest="omm_rearrange", help='Option for rearranging rescored poses. Default is true.')
    
    parser.add_argument("-omm_max_itr", "-omm_max_itr", type=int, default=5,
                       dest="omm_max_itr", help='Maximum steps for OpenMM minimization. Default is 5')
    
    parser.add_argument("-omm_env", "-omm_environment", type=int, default=0,
                       dest="omm_environment", help='0 for vaccum; 1 for implicit. Default is 0')
    
    
    return parser