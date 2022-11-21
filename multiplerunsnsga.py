#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
from InitTest import *
import Truss
from Truss.ClampedBeam import *

folder = os.path.dirname(__file__) + '/' # Working folder
folder = folder.replace("\\","/")
n = 2*nEls # Dimension of the problem

###################### !!! TO MODIFY ACCORDING TO THE PROBLEM !!!
pop_size = 100#n
it_max = 260#1000
nb_runs = 20
savingpath = folder+'Fem_case'+str(case)+'/titan_mag_alu_steel/'+type_mat+'/'+func+'/'+solver+'/'
if data_save and not os.path.exists(savingpath):
    os.makedirs(savingpath)
######################

print('FE case '+str(case)+', '+opt_pb+' with '+solver+', '+str(type_mat)+ ', '+str(nb_runs)+' runs')
    
problem = Problem(n,nobj,1,object_func) # dim, nb obj, nb constraints, fitness

if type_mat == 'Int_type':
    problem.types[0:nEls] = Integer(Emag, Esteel)
elif type_mat == 'Int_type2':
    problem.types[0:nEls] = Integer(Emag2, Esteel2)
elif type_mat == 'Categ_type':
    problem.types[0:nEls] = Integer(1, 4)

problem.types[nEls:n] = Real(thick_lb,thick_ub)
problem.constraints[:] = '<=0'

algorithm = NSGAII(problem, pop_size, variator=CompoundOperator(SBX(), HUX(), PM(), BitFlip()))

" Define target values "
if nobj > 1:
    maxv = 37
    target_min = 6.8889
    
maxpower = math.floor(math.log10(maxv-target_min)) # power of the max possible value
last = 10**maxpower
# range of magnitude between max possible value value and target_min
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

mat = np.zeros((it_max+1,nb_runs)) # Matrix of best fvalues (or hypervolumes) per run
consv = np.zeros((it_max+1,nb_runs)) # Matrix of (maximum) constraint violations per run
rates = np.zeros((it_max+1, nb_runs)) # Matrix of rates of reached targets
" Dictionaries consisting of the fvalues on the pareto for each global run "
F1 = {}
F2 = {}
F3 = {}
NORMS = {}
COLOR = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 
         'olive', 'cyan', 'black', 'gold', 'yellowgreen', 'darkblue', 'crimson',
         'pink', 'deeppink', 'darkcyan', 'olivedrab', 'tan', 'orchid']

" To use cma-es "
#ub = np.zeros(n)
#ub[0:nEls] = int(Esteel)
#ub[nEls:n] = 0.05
#lb = np.zeros(n)
#lb[0:nEls] = int(Emag)
#lb[nEls:n] = 0.01
#opts = {'bounds': [lb,ub], 'CMA_stds': (ub-lb)/4, 'popsize': pop_size, 
#        'maxiter': it_max, 'tolfun': 1e-11}

for it in range(nb_runs): # Launch nb_runs times the algorithm
    " Define target values "
    target_values = copy.copy(targets)
    target_nb = 0 # number of reached target values
    " Using cma-es population "
#    algorithm.population = []
#    Thicknesses = np.zeros(nEls)
#    Es = np.zeros(nEls)
#    for k in range(nEls):
#        Es[k] = random.choice([Mag, Titan, Alu, Steel])
#        Thicknesses[k] = random.uniform(0.01,0.05)
#    x0 = np.concatenate((Es, Thicknesses))
#    sigma0 = 1
#    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)
#    cma_indivs = es.ask() # cma-es initial population
#    for r in range(pop_size):
#        indiv = platypus.Solution(problem)
#        cma_ind = cma_indivs[r]
#        MEncoding = [] # Binary encoding of Es
#        indiv.variables[nEls:] = cma_ind[nEls:]
#        for el in range(nEls):
#            indiv.variables[el] = problem.types[el].encode(min(Youngs, key=lambda v:abs(v - cma_ind[el])))
#        indiv.objectives[:] = object_func(cma_ind)[0]
#        U, F = FE.solveq(case,Fa,cma_ind[0:nEls],alphas,cma_ind[nEls:],gcor,nEls,nN,connect)
#        indiv.constraint_violation = FE.dispviolation(case, U, Uy3max, Ux5max)
#        algorithm.population.append(indiv)
         
    " Define initial population "
    algorithm.population = []
    for _ in range(pop_size):
        indiv = platypus.Solution(problem)
        M = np.zeros(nEls) # Materials
        RealM = np.zeros(nEls)
        MEncoding = [] # Binary encoding of Es
        Thicknesses = np.zeros(nEls)
        for k in range(nEls):
            M[k] = random.choice((Mag, Titan, Alu, Steel))
            if type_mat == 'Int_type2':
                RealM[k] = M[k] * 1e6
            elif type_mat == 'Categ_type':
                if M[k] == 1:
                    RealM[k] = Etitan
                elif M[k] == 2:
                      RealM[k] = Emag
                elif M[k] == 3:
                      RealM[k] = Esteel
                elif M[k] == 4:
                        RealM[k] = Ealu
            else:
                RealM[k] = M[k]
            MEncoding.append(problem.types[k].encode(M[k]))
            Thicknesses[k] = random.uniform(0.01,0.05)
        indiv.variables[0:nEls] = MEncoding
        indiv.variables[nEls:2*nEls] = Thicknesses
        indiv.objectives[:] = object_func(np.concatenate((M, Thicknesses)))[0]
        U, F = FE.solveq(case,Fa,RealM,alphas,Thicknesses,gcor,nEls,nN,connect)
        indiv.constraint_violation = FE.dispviolation(case, U, Uy3max, Ux5max)
        algorithm.population.append(indiv)
    " Info to plot convergence graphs "
    reached_targets = []
    pareto = platypus.nondominated(algorithm.population)
    if nobj == 1:
        best = pareto[0]
        fvalues, constr_violation = best.objectives[:], [best.constraint_violation]
        " number of reached target values "
        if constr_violation[-1] == 0:
            tmp = [target_values[i] for i in range(len(target_values)) if target_values[i] >= (fvalues[-1]-target_min)]
            target_nb += len(tmp)
            rates[0,it] = len(tmp)/nb_of_targets
            if tmp != []:
                for el in range(len(tmp)):
                    reached_targets.append(tmp[el])
                    target_values.pop(0)
        else:
            rates[0,it] = 0
    elif nobj == 2:
        constr_violation = [] # constraint violations of the current population
        max_constr_violation = [] # maximal constraint violation per iteration
        for x in algorithm.population:
            constr_violation.append(x.constraint_violation)
        max_constr_violation.append(max(constr_violation))
        if func == 'cost&weight':
            hyp = Hypervolume(minimum=[1,0.6] , maximum=[37,72])
        elif func == 'weight&compliance':
            hyp = Hypervolume(minimum=[0.6,199] , maximum=[72,26825])
        elif func == 'cost&compliance':
            hyp = Hypervolume(minimum=[1,199] , maximum=[37,26825])
        hypervol = [hyp.calculate(pareto)]
    elif nobj == 3:
        constr_violation = [] # constraint violations of the current population
        max_constr_violation = [] # maximal constraint violation per iteration
        for x in algorithm.population:
            constr_violation.append(x.constraint_violation)
        max_constr_violation.append(max(constr_violation))
        hyp = Hypervolume(minimum=[1,0.6,199], maximum=[37,72,26825])
        hypervol = [hyp.calculate(pareto)]
    " Run the algorithm "
    for k in range(it_max):
        algorithm.run(1)
        pareto = platypus.nondominated(algorithm.population)
        if nobj == 1:
            best = pareto[0]
            fvalues.append(best.objectives[0])
            constr_violation.append(best.constraint_violation)
            if constr_violation[-1] == 0:
                tmp = [target_values[i] for i in range(len(target_values)) if target_values[i] >= fvalues[-1]-target_min]
                target_nb += len(tmp)
                if tmp != []:
                    rates[k+1,it] = target_nb/nb_of_targets
                    for el in range(len(tmp)):
                        reached_targets.append(tmp[el])
                else:
                    rates[k+1,it] = rates[k,it]
                for j in range(len(tmp)):
                    target_values.pop(0)
            else:
                rates[k+1,it] = rates[k,it]
            if (k+1)%50 == 0:
                print('Run '+str(it+1)+', generation ' + str(k+1) + ':\n')
                for idx in range(nEls):
                    print(str(problem.types[0].decode(best.variables[idx])) + ', ')
                print(str(best.variables[nEls:2*nEls]))
        elif nobj == 2 or nobj == 3:
            hypervol.append(hyp.calculate(pareto))
            constr_violation = []
            for x in algorithm.population:
                constr_violation.append(x.constraint_violation)
            max_constr_violation.append(max(constr_violation))
            if (k+1)%50 == 0:
                print('Run '+str(it+1)+', generation '+str(k+1)+' hypervolume = '
                      +str(hypervol[k])+':\n')
    " End of one run -- Store interesting quantities "   
    " Draw pareto without duplicates "
    unique_pareto = list()
    for sublist in pareto:
        if sublist.variables not in (x.variables for x in unique_pareto):
            unique_pareto.append(sublist)
    for indiv in unique_pareto:
        for m in range(nEls):
            indiv.variables[m] = problem.types[0].decode(indiv.variables[m])
    " Print mesh solution "
    plt.figure()        
    for k in range(len(unique_pareto)):
        ind = unique_pareto[k]
        print('Solution ' + str(k+1))
        FE.displaymaterials(ind.variables, Titan, Mag, Steel, Alu, connect, gcor)
        if data_save:
            plt.savefig(savingpath+'case'+str(case)+'_run'+str(it+1)
                        +'_mesh_solution'+str(k+1))    
    if data_save:
        with open(savingpath+'case'+str(case)+'_nondominated_run'
                  +str(it+1)+'.txt', 'w+') as in_file:
            for k in range(len(unique_pareto)):
                in_file.write(str(unique_pareto[k]) + '\n')
            
    if nobj == 1:
        mat[:,it] = fvalues
        consv[:,it] = constr_violation
        
    elif nobj == 2 or nobj == 3:
        mat[:,it] = hypervol
        consv[:,it] = max_constr_violation
        f1 = []
        f2 = []
        f3 = []
        for sol in unique_pareto:
            f1.append(sol.objectives[0])
            f2.append(sol.objectives[1])
            if nobj == 3:
                f3.append(sol.objectives[2])
        if nobj == 2:              
            f = sorted(list((f1[i], f2[i]) for i in range(len(f1))), key=lambda v:v[0])    
            F1[it] = [f[i][0] for i in range(len(f))]
            F2[it] = [f[i][1] for i in range(len(f))]
        elif nobj == 3:
            f = sorted(list((f1[i], f2[i], f3[i]) for i in range(len(f1))), key=lambda v:v[0])
            F1[it] = [f[i][0] for i in range(len(f))]
            F2[it] = [f[i][1] for i in range(len(f))]
            F3[it] = [f[i][2] for i in range(len(f))]            
            NORMS[it]  = []
            for k in range(len(f1)):
                NORMS[it].append(np.linalg.norm([F1[it][k], 
                     F2[it][k], F3[it][k]]))
" End of all runs "
" Plot best fvalues and constraint violations for all runs sand save matrices 'mat' and 'consv' "
feval = np.linspace(pop_size,mat.shape[0]*pop_size,num=mat.shape[0])
if nobj == 1:
    plt.figure()
    plt.title(solver+', '+str(opt_pb)+', case '+str(case)+ ', '+str(type_mat)+', best fvalues')
    plt.xlabel('Function evaluations')
    plt.ylabel('Best fvalues')
    for k in range(nb_runs):
        plt.loglog(feval, mat[:,k], COLOR[k])
    if data_save:
        plt.savefig(savingpath+'best_fvalues_'+str(nb_runs)+'runs')
        
    plt.figure()
    plt.title(solver+', '+str(opt_pb)+', case '+str(case)+ ', '+str(type_mat)+', constraint violations')
    plt.xlabel('Function evaluations')
    plt.ylabel('Constraint violations')
    for k in range(nb_runs):
        plt.loglog(feval, consv[:,k], COLOR[k])
    if data_save:
        plt.savefig(savingpath+'constr_viol_'+str(nb_runs)+'runs')
        np.savetxt(savingpath+'best_fvalues_'+str(nb_runs)+'runs.txt', mat)
        np.savetxt(savingpath+'constr_viol_'+str(nb_runs)+'runs.txt', consv)
        
    plt.figure()
    plt.title(solver+', '+str(opt_pb)+', case '+str(case)+', '+str(type_mat)+', ECDF')
    plt.xlabel('Function evaluations')
    plt.ylabel('Fraction of targets reached')
    for k in range(nb_runs):
        plt.plot(feval, rates[:,k], color=COLOR[k])
    if data_save:
        plt.savefig(savingpath+'ECDF_'+str(nb_runs)+'runs')
        np.savetxt(savingpath+'target_rates_'+str(nb_runs)+'runs.txt',rates)
    
" Plot the pareto of each global run "
if nobj == 2  or nobj == 3:
    " Plot and save hypervolumes and maximal constraint violations "
    plt.figure()
    plt.title(solver+', '+str(opt_pb)+', case '+str(case)+', hypervolumes')
    plt.xlabel('Function evaluations')
    plt.ylabel('Hypervolumes')
    for k in range(nb_runs):
        plt.plot(feval, mat[:,k], COLOR[k])
    if data_save:
        plt.savefig(savingpath+'hypervolumes_'+str(nb_runs)+'runs')
    plt.figure()
    plt.title(solver+', '+str(opt_pb)+', case '+str(case)+', max constraint violations')
    plt.xlabel('Function evaluations')
    plt.ylabel('maximal constraint violations')
    for k in range(nb_runs):
        plt.loglog(feval, consv[:,k], COLOR[k])
    if data_save:
        plt.savefig(savingpath+'max_constr_viol_'+str(nb_runs)+'runs')
        np.savetxt(savingpath+'best_fvalues_'+str(nb_runs)+'runs.txt', mat)
        np.savetxt(savingpath+'max_constr_viol_'+str(nb_runs)+'runs.txt', consv)

    colors = ['#8ffe09','r','#0033ff','#003311','#993333','#21c36f','#c46210',
              '#ed3cca','#ffbf00','g','#000000']
    plt.figure()
    if nobj == 3:
        ax = plt.axes(projection='3d')
        ax.set_title(solver +', rounded offspring, '+str(opt_pb)+', case '+str(case))
        ax.set_xlabel('cost')
        ax.set_ylabel('weight')
        ax.set_zlabel('compliance')
    elif nobj == 2:
        plt.title(solver+', rounded offspring, ' + str(opt_pb) +', case '+ str(case))
        if func == 'cost&weight':
            plt.xlabel('cost')
            plt.ylabel('weight')
        elif func == 'weight&compliance':
            plt.xlabel('weight')
            plt.ylabel('compliance')
    for k in range(len(F1)):
        if nobj == 2:
            plt.plot(F1[k], F2[k], 'o', color=COLOR[k])
            plt.step(F1[k], F2[k], where='post', color=COLOR[k])
        elif nobj == 3:
            " Gradient of colours according to the distance to the origin "
            cNorm = matplotlib.colors.Normalize(vmin=min(NORMS[k]), vmax=max( NORMS[k]))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('jet'))
            scalarMap.set_array(NORMS[k])
            ax.scatter3D(F1[k], F2[k], F3[k], c=scalarMap.to_rgba(NORMS[k]))

    if data_save:
        plt.savefig(savingpath+'pareto_case'+str(case))
