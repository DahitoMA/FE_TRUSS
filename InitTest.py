#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# import Truss
# from Truss.ClampedBeam import *

###################### !!! TO MODIFY ACCORDING TO THE PROBLEM !!!
solver = 'mads' # 'cmaes', 'nsga2', 'mads'
" material variable type (integer or categorical) "
type_mat = 'Categ_type' # 'Int_type', 'Int_type2', 'Categ_type', 'Real_type' or 'Cycle_type'
" case "
case = 1 # 1 or 2
" Function(s) to minimize "
func = 'weight' # 'cost', 'weight', 'compliance', 'cost&weight', 'weight&compliance',
            # 'cost&compliance' or 'cost&weight&compliance'
nobj = func.count('&') + 1 # Number of objective functions
" Option to save data (plots and arrays) "
data_save = True # True or False
opt_pb = 'min ' + func
thick_lb = 0.01 # Thickness lower bound
thick_ub = 0.05 # Thickness upper bound

####################### Target values
target_per_mag = 10 # Number of target values per magnitude
target_precision = 1e-4

if '&' not in func: # a single objective --> target values
    if func == 'cost':
        maxv = 37 # sup of the max possible value
        if case == 1:
#            target_max = 9.5525
            target_min = 6.8889 # best min seen
        elif case == 2:
#            target_max = 37
            target_min = 4.2631
    elif func == 'weight':
        maxv = 72
        if case == 1:
#            target_max = 18.16
            target_min = 12.6887
        elif case == 2:
#            target_max = 72
            target_min = 7.2323
    elif func == 'compliance':
        if case == 1:
            maxv = 26825
#            target_max = 845
            target_min = 241.4213
        elif case == 2:
            maxv = 22223
#            target_max = 26825
            target_min = 199.9998
    