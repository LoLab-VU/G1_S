from fan_modules_mod import *
from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.core import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify

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


print model.rules
# quit()
# declare_monomers()
# declare_parameters()
# declare_initial_conditions()
# declare_observables()
# declare_functions()
# declare_rules()

     
generate_equations(model, verbose=True)
for i in range(len(model.species)):
	print i, model.species[i]
print
# quit()
  
for monomers in model.monomers:
    print monomers
print
  
for i in range(len(model.parameters)):
    print str(i) + ":"
    print model.parameters[i]
 
 
for initial_conditions in model.initial_conditions:
    print initial_conditions
print

for x in model.parameters_initial_conditions():
    print x, ":", x.value
print 
print

for x in model.parameters_unused():
    print x, ":", x.value
print 

for x in model.parameters_rules():
    print x

# for obs in model.observables:
#     print obs, ":", obs.species, ",", obs.coefficients
#     obs_string = ''
#     for i in range(len(obs.coefficients)):
#         if i > 0: obs_string += " + "
#         obs_string += "__s"+str(obs.species[i])
#         if obs.coefficients[i] > 1:
#              obs_string += "*"+str(obs.coefficients[i])
#     print obs_string
#   
# for rules in model.rules:
#     print rules
# print
#  
for i in range(len(model.species)):
    print str(i)+":", model.species[i]
print
 
for i in range(len(model.odes)):
    print str(i)+":", model.odes[i]
print

for i in range(len(model.parameters)):
    print str(i)+":", model.parameters[i], model.parameters[i].value
# quit()  
 
from pysb.generator.bng import BngGenerator
print BngGenerator(model).get_content()

def normalize_output(y):
	y['OBS_CycA']/=max(y['OBS_CycA']) if max(y['OBS_CycA']) > 0.0 else y['OBS_CycA']
	y['OBS_CycE']/=max(y['OBS_CycE']) if max(y['OBS_CycE']) > 0.0 else y['OBS_CycE']
	y['OBS_E2F']/=max(y['OBS_E2F']) if max(y['OBS_E2F']) > 0.0 else y['OBS_E2F']
	y['OBS_p53']/=max(y['OBS_p53']) if max(y['OBS_p53']) > 0.0 else y['OBS_p53']
	y['OBS_p27']/= max(y['OBS_p27']) if max(y['OBS_p27']) > 0.0 else y['OBS_p27']
	y['OBS_Int']/=max(y['OBS_Int']) if max(y['OBS_Int']) > 0.0 else y['OBS_Int']
	y['OBS_Mdm2']/=max(y['OBS_Mdm2']) if max(y['OBS_Mdm2']) > 0.0 else y['OBS_Mdm2']
	y['OBS_p21']/=max(y['OBS_p21']) if max(y['OBS_p21']) > 0.0 else y['OBS_p21']
	y['OBS_p16']/=max(y['OBS_p16']) if max(y['OBS_p16']) > 0.0 else y['OBS_p16']
	y['OBS_Rb']/=max(y['OBS_Rb']) if max(y['OBS_Rb']) > 0.0 else y['OBS_Rb']
	y['OBS_CycD_CDK46']/=max(y['OBS_CycD_CDK46']) if max(y['OBS_CycD_CDK46']) > 0.0 else y['OBS_CycD_CDK46']

  
t = linspace(0,3000,300)



# Observable("OBS_p27", p27(b=None))
# 	Observable("OBS_p53", p53(b=None))
# 	Observable("OBS_E2F", E2F(b=None))
# 	Observable("OBS_CycE", CycE(c=None))
# 	Observable("OBS_CycA", CycA(c=None))
# 	Observable("OBS_Int", I())
# 	Observable("OBS_Mdm2", Mdm2(i=None))
# 	Observable("OBS_p21", p21(b=None))
# 	Observable("OBS_p16", p16(b=None))
# 	Observable("OBS_Rb", Rb(b=None, Y='u'))
# 	Observable("OBS_CycD_CDK46", CycD(c=1) % CDK46(Y='p',b=None,c=1))
# 	Observable("signal", Signal())
# 	Observable("signal_damp", SignalDamp())
#####
# for monomers in model.monomers:
# 	print monomers
# print
#   
# for parameters in model.parameters:
# 	print parameters
# print
#   
# for initial_conditions in model.initial_conditions:
# 	print initial_conditions
# print
#   
# for observables in model.observables:
# 	print observables
# print
  
# for rules in model.rules:
# 	print rules
# print
#   
# for i in range(len(model.species)):
#     print str(i) + ":"
#     print model.species[i]
#    
# for i in range(len(model.species)):
#     print str(i) + ":"
#     print model.odes[i]
# quit()
#  
#   
set_dna_damage(0.0)
y = odesolve(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
normalize_output(y)
 
pl.figure()
for obs in ["OBS_p27", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper right')
pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper right')
#####
set_dna_damage(0.005)
y = odesolve(model,t,verbose=True, rtol = 1e-15, atol = 1e-15)
normalize_output(y)

pl.figure()
for obs in ["OBS_p53", "OBS_E2F", "OBS_CycE", "OBS_CycA"]:
	pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper left')

pl.figure()
for obs in ["OBS_p53", "OBS_Mdm2", "OBS_Int", "OBS_p21"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper right')

pl.show()    

# import fan_modules as m
# import pylab as pl
# from pysb.integrate import odesolve
# from sympy import sympify
 
# Set to False to plot I, Mdm2, p21, p53
# cyclins = True
 
#################################
## Force PySB to generate ODEs ##
#################################
# t = pl.linspace(0,3000, num = 300)
# y = odesolve(m.model,t)
 
# ###############
# ## Overrides ##
# ###############
# _source = "s30"
# #m.model.odes[24] = sympify("k60*s32 + k61*s28 - k62*s24 - k73*s24*s25*s28 + k73*s24*s25*s31 - k74*s24*s25")
# m.model.odes[25] = sympify("k63*" + _source + " + k66*s27**9/(k65**9 + s27**9) - k64*s25") #!!!
# m.model.odes[27] = sympify("(k70*s24*s28)/(1 + k71*s24*s25) - k67*s27") #!!!
# m.model.odes[28] = sympify("-k72*s28") #!!!
# #m.model.odes[18] = sympify("k40*s32 + k41/(1 + k42*s23) - k43*s18 - k44*s18*s5")
# #m.model.odes[23] = sympify("-k45*s21*s23 + k55*s22 + k56*s32 - k57*s23 + k58/(1 + k59*s18)")
# print m.model.odes[25]
# print m.model.odes[27]
# print m.model.odes[28]
# print
# print m.model.odes[24]
 
 
# ###########
# ## Solve ##
# ###########
# t = pl.linspace(0,3000, num = 300)
# print t
# print len(t)
# print 
#  
# y = odesolve(m.model,t)
#  

 	
 	
 
 
### Checking ODEs ###
 
 
 
#### NOTE: NORMALIZING SO WE CAN COMPARE TO IWAMOTO ET AL.

#  
# ##########
# ## Plot ##
# ##########
# # pl.ion()
# pl.figure()
#  
# pl.plot(y['OBSp53'], label="p53")
# if cyclins:
# 	pl.plot(y['OBSCycA'], label="CycA")
# 	pl.plot(y['OBSCycE'], label="CycE")
# 	pl.plot(y['OBSE2F'], label="E2F")
# 	pl.plot(y['OBSp27'], label="p27")
# else: 
# 	pl.plot(y['OBSI'], label="I")
# 	pl.plot(y['OBSMdm'], label="Mdm2")
# 	pl.plot(y['OBSp21'], label="p21")
# pl.legend(loc='lower left')
# pl.show()

