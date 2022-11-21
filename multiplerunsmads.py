#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
from InitTest import *
import Truss
from Truss.ClampedBeam import *
from io import StringIO
#import subprocess as sp

folder = os.path.dirname(__file__) + '/' # Working folder
folder = folder.replace("\\","/")
n = 2*nEls # Dimension of the problem

###################### !!! TO MODIFY ACCORDING TO THE PROBLEM !!!
it_max = 154#1000
nb_runs = 1#20
param_file = 'param.txt' # 'Fem_case'+str(case)+'/titan_mag_alu_steel/'+type_mat+'/'+func+'/'+solver+'/'+'param.txt'
bb_file = 'objectives.py' # 'Fem_case'+str(case)+'/titan_mag_alu_steel/'+type_mat+'/'+func+'/'+solver+'/'+'objectives.py'
neighbors_exe = 'neighbors_exe.py' # TODO
# min_poll_size = 'r1e-5' # variable range x 1e-5
min_mesh_size = '1e-11'
dir_type = 'ORTHO N+1 NEG'
nomad_type = 'ortho_nplus1_neg'
nb_dir = n+1 # number of poll directions
extended_poll_enabled = 'no'
max_bb_eval = 4000
savingpath = folder+'Fem_case'+str(case)+'/titan_mag_alu_steel/'+type_mat+'/'+func+'/'+solver+'/'
if data_save and not os.path.exists(savingpath):
    os.makedirs(savingpath)
######################

print('FE case '+str(case)+', '+opt_pb+' with '+solver+', '+str(type_mat)+ 
      ', ' + str(nomad_type) + ', '+str(nb_runs)+' runs')

lb = np.zeros(n)
ub = np.zeros(n)
lb[nEls:n] = thick_lb #0.01
ub[nEls:n] = thick_ub #0.05

" Problem definition "
if type_mat in ('Int_type', 'Int_type2'):
    bb_input_type = '(' + nEls*'I ' + nEls*'R ' + ')'
elif type_mat in ('Categ_type','Cycle_type'):
    bb_input_type = '(' + nEls*'I ' + nEls*'R ' + ')'
#    bb_input_type = '(' + nEls*'C ' + nEls*'R ' + ')'
#    extended_poll_enabled = 'yes' # Uncomment once neighbors_exe is fixed
elif type_mat == 'Real_type':
    bb_input_type = '* R'
    
if type_mat == 'Categ_type':
    lb[0:nEls] = 1
    ub[0:nEls] = 4
elif type_mat == 'Cycle_type':
    lb[0:nEls] = -np.inf
    ub[0:nEls] = np.inf
else:
    lb[0:nEls] = Mag
    ub[0:nEls] = Steel
    
" Create blackbox file "
with open(bb_file, 'w+') as in_file:
    in_file.write('#!/usr/bin/env python'
                  + '\n# -*- coding: utf-8 -*-'
                  + '\nimport sys'
                  + '\nfrom InitTest import *'
                  + '\nfrom Truss.ClampedBeam import ' + object_func.__name__ #object_func.func_name
                  + '\n\ndef obj_func(bb_args):'
                  + '\n\tfloat_args = []'
                  + '\n\tfor k in range(len(bb_args)):'
                  + '\n\t\tfloat_args.append(float(bb_args[k]))' 
                  + '\n\toutput = ' + object_func.__name__ + '(float_args)' #object_func.func_name + '(float_args)'
                  + '\n\treturn (output)'
                  + '\n\nif __name__ == \"__main__\":'
                  + '\n\tinput_args_txt = sys.argv[1]'
                  + '\n\targs = []'
                  + '\n\tfor line in open(input_args_txt,\"r\"):'
                  + '\n\t\tfor v in line.split(\" \"):'
                  + '\n\t\t\targs.append(v)'
                  + '\n\toutputs_bb = obj_func(args)'
                  + '\n\tbuffer = \"\\'+'n\".join([str(o) for o in outputs_bb]+[\" \"])'
                  + '\n\tsys.stdout.write(buffer)'
                  + '\n\tsys.exit(0)')
    
" Define target values "
maxpower = math.floor(math.log10(maxv-target_min)) # power of the max possible value
last = 10**maxpower
# range of magnitude between max possible value and target_min
mag_range = int(maxpower-math.floor(math.log10(target_precision)))
target_values = list()
for k in range(mag_range):
    nextv = last * 1e-1
    if k == 0 and maxpower != math.floor(math.log10(target_min)):
        values = np.logspace(maxpower,maxpower-1,num=target_per_mag+1)
        for j in range(len(values)-1):
            target_values.append(values[j])
    else: 
        values = np.linspace(last, nextv, num=target_per_mag+1)
        if k == mag_range-1:
            for j in range(len(values)):
                target_values.append(values[j])
        else:
            for j in range(len(values)-1):
                target_values.append(values[j])
    last = nextv
targets = copy.copy(target_values)
nb_of_targets = len(target_values)

mat = np.zeros((nb_dir*it_max+1,nb_runs)) # Matrix of best fvalues (or hypervolumes) per run
mat2 = np.zeros((it_max+1,nb_runs)) # Matrix of best fvalues per iteration
consv = np.zeros((nb_dir*it_max+1,nb_runs)) # Matrix of (maximum) constraint violations per run
consv2 = np.zeros((it_max+1,nb_runs)) # Matrix of (maximum) constraint violations per iteration
rates = np.zeros((it_max+1,nb_runs)) # Matrix of rates of reached targets

F1 = {}
# F2 = {}
# F3 = {}
NORMS = {}
COLOR = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 
         'olive', 'cyan', 'black', 'gold', 'yellowgreen', 'darkblue', 'crimson',
         'pink', 'deeppink', 'darkcyan', 'olivedrab', 'tan', 'orchid']

for it in range(nb_runs): # Launch nb_runs times the algorithm
    print('\nRun '+str(it+1)+'\n')
    statfile = 'Fem_case'+str(case)+'/titan_mag_alu_steel/'+type_mat+'/'+func+'/'+solver+'/'+'run'+str(it+1)+'_stats.txt'
    " Define target values "
    target_values = copy.copy(targets)
    target_nb = 0 # number of reached target values
    reached_targets = []
    " Define initial solution "
    x0 = np.zeros(n)
    for j in range(nEls):
        x0[j] = random.choice((Mag, Titan, Alu, Steel))
        x0[j+nEls] = random.uniform(thick_lb,thick_ub)
    " Create parameter file "        
    with open(param_file, 'w+') as in_file:
        in_file.write('X0 (')
        for k in range(n-1):
            in_file.write(str(x0[k]) + ' ')
        in_file.write(str(x0[n-1])+')')
        in_file.write('\nDIMENSION ' + str(n)
                      + '\nLOWER_BOUND (' + ' '.join(map(str, tuple(lb))) + ')'
                      + '\nUPPER_BOUND (' + ' '.join(map(str, tuple(ub))) + ')'
                      + '\nBB_EXE \"$python '+ bb_file +'\"'
                      + '\nBB_OUTPUT_TYPE OBJ CSTR'
                      + '\nMAX_BB_EVAL '+str(max_bb_eval)
#                      + '\nTMP_DIR /tmp'
                      + '\nDISPLAY_DEGREE NORMAL_DISPLAY'
                      + '\nDISPLAY_STATS BBE ( SOL ) OBJ'
                      + '\nDISPLAY_ALL_EVAL yes'
                      + '\nSTATS_FILE ' + statfile + ' BBE SOL BBO'
                      + '\nDIRECTION_TYPE ' + dir_type
                      # + '\nMIN_MESH_SIZE ' + min_mesh_size
                      + '\nBB_INPUT_TYPE ' + bb_input_type # R (real), B (binary), C (categorical), I (integer)
                      + '\nEXTENDED_POLL_ENABLED ' + extended_poll_enabled)
#                      + '\nNEIGHBORS_EXE ' + neighbors_exe) # TODO
                      # + '\nMIN_POLL_SIZE ' + min_poll_size)
                      # + '\nF_TARGET ' +  f_target # in () if multiobjectif pb
                      # + '\nLH_SEARCH p0 pi'
                      # + '\nCACHE_FILE cache.bin')
                      
#    run = sp.Popen(['nomad', param_file], stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
#    run.communicate()
    # os.system('C:/Users/mdahito/Desktop/nomad.3.9.1_Personal/bin/nomad.exe '+ param_file)
    os.system('nomad '+ param_file)
    " End of one run -- Store interesting quantities "
    try:
        mtx = np.loadtxt(savingpath+'run'+str(it+1)+'_stats.0.txt')
    except:
        with open(savingpath+'run'+str(it+1)+'_stats.0.txt') as in_file:
            lines = in_file.readlines()
            cplines = [] # will contain lines of the file bis
            for i in range(len(lines)-1):
                cplines.append(lines[i])
#                try:
#                    np.loadtxt(StringIO(lines[i]))
#                    cplines.append(lines[i])
#                except:
#                    pass                
        with open(savingpath+'run'+str(it+1)+'_stats.bis.txt', 'w+') as in_file:
            for i in range(len(cplines)):
                in_file.write(cplines[i])
        mtx = np.loadtxt(savingpath+'run'+str(it+1)+'_stats.bis.txt')
            
    mat[0,it] = mtx[0,-2]
    mat2[0,it] = mat[0,it]
    consv[0,it] = max(0,mtx[0,-1])
    consv2[0,it] = consv[0,it]
    left_lines = np.shape(mat)[0] - np.shape(mtx)[0] # left evaluation budget
    best_idx = 0 # index of the best point
    for l in range(1,np.shape(mtx)[0]):
        if max(0,mtx[l,-1]) < consv[l-1,it]: # smaller constraint violation
            best_idx = l
            consv[l,it] = max(0,mtx[l,-1])
            mat[l,it] = mtx[l,-2]
        elif max(0,mtx[l,-1]) == consv[l-1,it]: # same constraint violation
            if mtx[l,-2] < mat[l-1,it]:
                best_idx = l
                consv[l,it] = max(0,mtx[l,-1])
                mat[l,it] = mtx[l,-2]
            else:
                consv[l,it] = consv[l-1,it]
                mat[l,it] = mat[l-1,it]
        else:
            consv[l,it] = consv[l-1,it]
            mat[l,it] = mat[l-1,it]
            
    if left_lines > 0:
        for r in range(left_lines):
            consv[np.shape(mtx)[0]+r,it] = consv[np.shape(mtx)[0]+r-1,it]
            mat[np.shape(mtx)[0]+r,it] = mat[np.shape(mtx)[0]+r-1,it]
                    
    if consv[0,it] == 0: # initialize rates
        tmp = [target_values[i] for i in range(len(target_values)) if target_values[i] >= (mat[0,it]-target_min)]
        target_nb += len(tmp)
        rates[0,it] = target_nb/nb_of_targets
        if tmp != []:
            for el in range(len(tmp)):
                reached_targets.append(tmp[el])
                target_values.pop(0)
                
    for k in range(it_max):
        mat2[k+1,it] = mat[(k+1)*nb_runs,it]
        consv2[k+1,it] = consv[(k+1)*nb_runs,it]
            
    for i in range(1,np.shape(mat2)[0]):
        if consv2[i,it] == 0:
            tmp = [target_values[j] for j in range(len(target_values)) if target_values[j] >= (mat[i,it]-target_min)]
            target_nb += len(tmp)
            rates[i,it] = target_nb/nb_of_targets
            if tmp != []:
                for el in range(len(tmp)):
                    reached_targets.append(tmp[el])
                    target_values.pop(0)
    
    "Print mesh solution"
    plt.figure()
    FE.displaymaterials(mtx[best_idx,1:-2], Titan, Mag, Steel, Alu, connect, gcor, cycle=False)
    
    if data_save:
        plt.savefig(savingpath+'case'+str(case)+'_run'+str(it+1)
                    +'_mesh_solution')
        with open(savingpath+'case'+str(case)+'_nondominated_run'
                  +str(it+1)+'.txt', 'w+') as in_file:
            for k in range(n):
                in_file.write(str(mtx[best_idx,k+1]) + '\t')
            
" End of all runs "
" Plot best fvalues and constraint violations for all runs and save matrices 'mat' and 'consv' "
feval = np.linspace(1,mat.shape[0],num=mat2.shape[0])
plt.figure()
plt.title(solver+', '+str(opt_pb)+', case '+str(case)+ ', '+str(type_mat)+', best fvalues')
plt.xlabel('Function evaluations')
plt.ylabel('Best fvalues')
for k in range(nb_runs):
    plt.loglog(feval, mat2[:,k], COLOR[k])
if data_save:
    plt.savefig(savingpath+'best_fvalues_'+str(nb_runs)+'runs')

plt.figure()
plt.title(solver+', '+str(opt_pb)+', case '+str(case)+ ', '+str(type_mat)+', constraint violations')
plt.xlabel('Function evaluations')
plt.ylabel('Constraint violations')
for k in range(nb_runs):
    plt.loglog(feval, consv2[:,k], COLOR[k])
if data_save:
    plt.savefig(savingpath+'constr_viol_'+str(nb_runs)+'runs')
    np.savetxt(savingpath+'best_fvalues_'+str(nb_runs)+'runs.txt', mat2)
    np.savetxt(savingpath+'constr_viol_'+str(nb_runs)+'runs.txt', consv2)
    
plt.figure()
plt.title(solver+', '+str(opt_pb)+', case '+str(case)+', '+str(type_mat)+', ECDF')
plt.xlabel('Function evaluations')
plt.ylabel('Fraction of targets reached')
for k in range(nb_runs):
    plt.plot(feval, rates[:,k], color=COLOR[k])
if data_save:
    plt.savefig(savingpath+'ECDF_'+str(nb_runs)+'runs')
    np.savetxt(savingpath+'target_rates_'+str(nb_runs)+'runs.txt',rates)
    
