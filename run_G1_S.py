from G1_S_v2 import *
from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants


# ***Generate ODEs and Plot***

declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
simulate_signal_degradation()
p16_p27_inhibition()
CDK2_activation()
Rb_E2F_activation()
DNA_damage_pathway()
     
# generate_equations(model, verbose=True)
#
# for monomers in model.monomers:
#     print monomers
# print
#   
# for i in range(len(model.parameters)):
#     print str(i) + ":"
#     print model.parameters[i]
# print
#  
# for initial_conditions in model.initial_conditions:
#     print initial_conditions
# print
# 
# for x in model.parameters_initial_conditions():
#     print x, ":", x.value
# print 
# print
# 
# for x in model.parameters_unused():
#     print x, ":", x.value
# print 
# 
# for i in range(len(model.parameters)):
# 	print str(i)+":", model.parameters[i], "=", model.parameters[i].value
# print
#
# for obs in model.observables:
#     print obs, ":", obs.species, ",", obs.coefficients
#     obs_string = ''
#     for i in range(len(obs.coefficients)):
#         if i > 0: obs_string += " + "
#         obs_string += "__s"+str(obs.species[i])
#         if obs.coefficients[i] > 1:
#              obs_string += "*"+str(obs.coefficients[i])
#     print obs_string
# print
# 
# for rules in model.rules:
#     print rules
# print
#   
# for i in range(len(model.parameters)):
#     print str(i)+":", model.parameters[i], model.parameters[i].value
# print
#
# for i in range(len(model.odes)):
#     print str(i) + ":", model.odes[i]
# print
# 
# for e in model.expressions:
#     print e, e.expand_expr()
# print
# 
# for obs in model.observables:
#     print obs, obs.species, obs.coefficients
# print
# 
# for i in range(len(model.species)):
#     print str(i)+":", model.species[i]
# print
# 
# quit()  
#
# from pysb.generator.bng import BngGenerator
# print BngGenerator(model).get_content()

def normalize_output(y):
	y['OBS_CycA']/=max(y['OBS_CycA']) if max(y['OBS_CycA']) > 0.0 else 1
	y['OBS_CycE']/=max(y['OBS_CycE']) if max(y['OBS_CycE']) > 0.0 else 1
	y['OBS_E2F']/=max(y['OBS_E2F']) if max(y['OBS_E2F']) > 0.0 else 1
	y['OBS_p53']/=max(y['OBS_p53']) if max(y['OBS_p53']) > 0.0 else 1
	y['OBS_p27']/= max(y['OBS_p27']) if max(y['OBS_p27']) > 0.0 else 1
	y['OBS_Int']/=max(y['OBS_Int']) if max(y['OBS_Int']) > 0.0 else 1
	y['OBS_Mdm2']/=max(y['OBS_Mdm2']) if max(y['OBS_Mdm2']) > 0.0 else 1
	y['OBS_p21']/=max(y['OBS_p21']) if max(y['OBS_p21']) > 0.0 else 1
	y['OBS_p16']/=max(y['OBS_p16']) if max(y['OBS_p16']) > 0.0 else 1
	y['OBS_Rb']/=max(y['OBS_Rb']) if max(y['OBS_Rb']) > 0.0 else 1
	y['OBS_CycD_CDK46']/=max(y['OBS_CycD_CDK46']) if max(y['OBS_CycD_CDK46']) > 0.0 else 1


t = linspace(0,3000,300)

vol = 1e-19
# ** Set No Damage **
# set_volume(vol)

Na_V = constants.N_A * vol

set_dna_damage(0.0) # For volume, * Na_V
y = odesolve(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
# set_dna_damage(0.0 * Na_V) # For volume, * Na_V
# y = run_ssa(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
normalize_output(y)

pl.figure()
for obs in ["OBS_p27", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.ylim(ymax=1.2)
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (No DNA Damage)", fontsize=22)
pl.savefig("G1_S Cell Cycle No DNA Damage1.png", format= "png")

pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.ylim(ymax=1.2)
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics (No DNA Damage)", fontsize=22)
pl.savefig("G1_S Cell Cycle No DNA Damage2.png", format= "png")

# ** Set Damage **

set_dna_damage(0.005) # For volume, * Na_V
y = odesolve(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
# set_dna_damage(0.005 * Na_V) # For volume, * Na_V
# y = run_ssa(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
normalize_output(y)

pl.figure()
for obs in ["OBS_p53", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
	pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.ylim(ymax=1.2)
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics with DNA Damage(0.005)", fontsize=22)
pl.savefig("G1_S Cell Cycle DNA Damage1.png", format= "png")

pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=re.match(r"OBS_(\w+)", obs).group(1), linewidth=3)
pl.legend(loc=0, prop={'size': 16})
pl.ylim(ymax=1.2)
pl.xlabel("Time (arbitrary units)", fontsize=22)
pl.ylabel("Protein Level", fontsize=22)
pl.xticks(fontsize=18)
pl.yticks(fontsize=18)
pl.title("Protein Dynamics with DNA Damage(0.005)", fontsize=22)
pl.savefig("G1_S Cell Cycle DNA Damage2.png", format= "png")

pl.show()

"""
# ***Generate ODEs and Plot***

declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
simulate_signal_degradation()
p16_p27_inhibition()
CDK2_activation()
Rb_E2F_activation()
DNA_damage_pathway()

set_volume(1.0e-19)
     
generate_equations(model, verbose=True)

# for i in range(len(model.species)):
#     print i, model.species[i]
# print
# 
# for monomers in model.monomers:
#     print monomers
# print
#  
# for initial_conditions in model.initial_conditions:
#     print initial_conditions
# print
# 
# for x in model.parameters_initial_conditions():
#     print x, ":", x.value
# print 
# print
# 
# for x in model.parameters_unused():
#     print x, ":", x.value
# print 
# 
# for x in model.parameters_rules():
#     print x
# print
# 
# for obs in model.observables:
#     print obs, ":", obs.species, ",", obs.coefficients
#     obs_string = ''
#     for i in range(len(obs.coefficients)):
#         if i > 0: obs_string += " + "
#         obs_string += "__s"+str(obs.species[i])
#         if obs.coefficients[i] > 1:
#              obs_string += "*"+str(obs.coefficients[i])
#     print obs_string
# print    
#    
# for rules in model.rules:
#     print rules
# print
#   
# for i in range(len(model.parameters)):
#     print str(i)+":", model.parameters[i], model.parameters[i].value
# print
# 
# for e in model.expressions:
#     print e, e.expand_expr()
# print
# 
# for obs in model.observables:
#     print obs, obs.species, obs.coefficients
# print
#
# quit()  
#  
# from pysb.generator.bng import BngGenerator
# print BngGenerator(model).get_content()

def normalize_output(y):
    y['OBS_CycA']/=max(y['OBS_CycA']) if max(y['OBS_CycA']) > 0.0 else 1
    y['OBS_CycE']/=max(y['OBS_CycE']) if max(y['OBS_CycE']) > 0.0 else 1
    y['OBS_E2F']/=max(y['OBS_E2F']) if max(y['OBS_E2F']) > 0.0 else 1
    y['OBS_p53']/=max(y['OBS_p53']) if max(y['OBS_p53']) > 0.0 else 1
    y['OBS_p27']/= max(y['OBS_p27']) if max(y['OBS_p27']) > 0.0 else 1
    y['OBS_Int']/=max(y['OBS_Int']) if max(y['OBS_Int']) > 0.0 else 1
    y['OBS_Mdm2']/=max(y['OBS_Mdm2']) if max(y['OBS_Mdm2']) > 0.0 else 1
    y['OBS_p21']/=max(y['OBS_p21']) if max(y['OBS_p21']) > 0.0 else 1
    y['OBS_p16']/=max(y['OBS_p16']) if max(y['OBS_p16']) > 0.0 else 1
    y['OBS_Rb']/=max(y['OBS_Rb']) if max(y['OBS_Rb']) > 0.0 else 1
    y['OBS_CycD_CDK46']/=max(y['OBS_CycD_CDK46']) if max(y['OBS_CycD_CDK46']) > 0.0 else 1

t = linspace(0,3000,300)

## ** Set No DNA Damage **

set_dna_damage(0.0 * Na_V) # For volume, * Na_V
y = run_ssa(model,t,verbose=True)
normalize_output(y)
 
pl.figure()
for obs in ["OBS_p27", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Stochastic Protein Dynamics with No DNA Damage")
pl.savefig("Stochastic G1_S Cell Cycle No DNA Damage1.png", format= "png")

pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Stochastic Protein Dynamics with No DNA Damage")
pl.savefig("Stochastic G1_S Cell Cycle No DNA Damage2.png", format= "png")

## ** Set DNA Damage **

set_dna_damage(0.005 * Na_V) # For volume, * Na_V
y = run_ssa(model,t,verbose=True)
normalize_output(y)

pl.figure()
for obs in ["OBS_p53", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Stochastic Protein Dynamics with DNA Damage(0.005)")
pl.savefig("Stochastic G1_S Cell Cycle DNA Damage1.png", format= "png")

pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Stochastic Protein Dynamics with DNA Damage(0.005)")
pl.savefig("Stochastic G1_S Cell Cycle DNA Damage2.png", format= "png")

pl.show()
"""