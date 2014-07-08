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
from scipy import constants

def set_volume(vol):
    global Na_V
    Na_V = constants.N_A * vol              #1/[]
        
    # ***2nd Order Reactions*** (divide by Avogadro's Number * Volume)
    model.parameters['k3'].value /= Na_V
    model.parameters['k20'].value /= Na_V
    model.parameters['k24'].value /= Na_V
    model.parameters['k30'].value /= Na_V
    model.parameters['k13'].value /= Na_V
    model.parameters['k44'].value /= Na_V
    model.parameters['k35'].value /= Na_V
    model.parameters['k36'].value /= Na_V
    model.parameters['k7'].value /= Na_V
    model.parameters['k11'].value /= Na_V
    model.parameters['k17'].value /= Na_V
    model.parameters['k22'].value /= Na_V
    model.parameters['k28'].value /= Na_V
    model.parameters['k46'].value /= Na_V
    model.parameters['k47'].value /= Na_V
    model.parameters['k48'].value /= Na_V
    model.parameters['k49'].value /= Na_V
    model.parameters['k50'].value /= Na_V
    model.parameters['k54'].value /= Na_V
    model.parameters['k45'].value /= Na_V
    model.parameters['k18'].value /= Na_V
    model.parameters['k26'].value /= Na_V
    model.parameters['k32'].value /= Na_V   
        
    # ***1st Order Reactions*** - LEAVE ALONE!
#     model.parameters['k72'].value
#     model.parameters['k4'].value
#     model.parameters['k2'].value
#     model.parameters['k21'].value
#     model.parameters['k25'].value
#     model.parameters['k31'].value
#     model.parameters['k6'].value
#     model.parameters['k10'].value
#     model.parameters['k69'].value
#     model.parameters['k8'].value
#     model.parameters['k12'].value
#     model.parameters['k68'].value
#     model.parameters['k5'].value
#     model.parameters['k9'].value
#     model.parameters['k14'].value
#     model.parameters['k15'].value
#     model.parameters['k16'].value
#     model.parameters['k23'].value
#     model.parameters['k29'].value
#     model.parameters['k53'].value
#     model.parameters['k57'].value
#     model.parameters['k51'].value
#     model.parameters['k55'].value
#     model.parameters['k39'].value
#     model.parameters['k62'].value
#     model.parameters['k64'].value
#     model.parameters['k67'].value
#     model.parameters['k19'].value
#     model.parameters['k27'].value
#     model.parameters['k33'].value
#     model.parameters['k38'].value
#     model.parameters['k61'].value
        
    # ***0th Order Reactions*** (multiply by Avogadro's Number * Volume)
    model.parameters['k1'].value *= Na_V
    model.parameters['k34'].value *= Na_V
    model.parameters['k40'].value *= Na_V
    model.parameters['k52'].value *= Na_V
    model.parameters['k56'].value *= Na_V
    model.parameters['k37'].value *= Na_V
    model.parameters['k60'].value *= Na_V
    model.parameters['k63'].value *= Na_V
     
        
    # ***In Functions *** (act differently)
    model.parameters['k41'].value *= (Na_V * Na_V)
    model.parameters['k44'].value /= Na_V
    model.parameters['k58'].value *= (Na_V * Na_V)
    model.parameters['k65'].value *= Na_V
    model.parameters['k66'].value *= Na_V
    model.parameters['k71'].value /= Na_V
    model.parameters['k73'].value /= (Na_V * Na_V)
    model.parameters['k74'].value /= Na_V
    model.parameters['k75'].value /= Na_V
    model.parameters['One'].value *= Na_V
    
    # ** Should not be changed **
#     model.parameters['k42'].value /= Na_V
#     model.parameters['k43'].value /= Na_V
#     model.parameters['k70'].value 
#     model.parameters['k59'].value /= Na_V
        
    # ***Initial Condition Parameters*** (multiply by Avogadro's Number * Volume)
    model.parameters['Y0_0'].value *= Na_V
    model.parameters['Y1_0'].value *= Na_V
    model.parameters['Y2_0'].value *= Na_V
    model.parameters['Y3_0'].value *= Na_V
    model.parameters['Y4_0'].value *= Na_V
    model.parameters['Y5_0'].value *= Na_V
    model.parameters['Y6_0'].value *= Na_V
    model.parameters['Y7_0'].value *= Na_V
    model.parameters['Y8_0'].value *= Na_V
    model.parameters['Y9_0'].value *= Na_V
    model.parameters['Y10_0'].value *= Na_V
    model.parameters['Y11_0'].value *= Na_V
    model.parameters['Y12_0'].value *= Na_V
    model.parameters['Y13_0'].value *= Na_V
    model.parameters['Y14_0'].value *= Na_V
    model.parameters['Y15_0'].value *= Na_V
    model.parameters['Y16_0'].value *= Na_V
    model.parameters['Y17_0'].value *= Na_V
    model.parameters['Y18_0'].value *= Na_V
    model.parameters['Y19_0'].value *= Na_V
    model.parameters['Y20_0'].value *= Na_V
    model.parameters['Y21_0'].value *= Na_V
    model.parameters['Y22_0'].value *= Na_V
    model.parameters['Y23_0'].value *= Na_V
    model.parameters['Y24_0'].value *= Na_V
    model.parameters['Y25_0'].value *= Na_V
    model.parameters['Y26_0'].value *= Na_V
    model.parameters['Y27_0'].value *= Na_V
    model.parameters['DDS_0'].value *= Na_V

###
# Kinetic parameters of proposed model
#     Parameter("k1", 5.00e-03)     #
#     Parameter("k2", 5.00e-04)     #
#     Parameter("k3", 5.00e-03)     # 
#     Parameter("k4", 2.50e-03)     # 
#     Parameter("k5", 7.50e-02)     # 
#     Parameter("k6", 2.50e-03)     #
#     Parameter("k7", 1.25e-03)     #
#     Parameter("k8", 2.50e-04)     #
#     Parameter("k9", 8.00e-04)     #
#     Parameter("k10", 5.00e-04)    #
#     Parameter("k11", 1.00e-03)    #
#     Parameter("k12", 2.00e-04)    #
#     Parameter("k13", 5.00e-04)    #
#     Parameter("k14", 5.00e-04)    #
#     Parameter("k15", 5.00e-04)    #
#     Parameter("k16", 5.00e-04)    #
#     Parameter("k17", 2.00e-03)    #
#     Parameter("k18", 5.00e-04)    #
#     Parameter("k19", 5.00e-03)    #
#     Parameter("k20", 5.00e-04)    #
#     Parameter("k21", 5.00e-05)    #
#     Parameter("k22", 2.50e-02)    #
#     Parameter("k23", 1.75e-03)    #
#     Parameter("k24", 2.25e-02)    # 
#     Parameter("k25", 1.75e-04)    #
#     Parameter("k26", 2.25e-02)    #
#     Parameter("k27", 1.75e-04)    #
#     Parameter("k28", 1.90e-02)    #
#     Parameter("k29", 5.00e-04)    #
#     Parameter("k30", 2.50e-03)    # 
#     Parameter("k31", 1.75e-04)    #
#     Parameter("k32", 2.50e-03)    #
#     Parameter("k33", 1.75e-04)    #
#     Parameter("k34", 5.00e-08)    #
#     Parameter("k35", 1.00e-02)    # 
#     Parameter("k36", 1.50e-03)    # 
#     Parameter("k37", 5.00e-05)    #
#     Parameter("k38", 1.00e-02)    #
#     Parameter("k39", 5.00e-03)    #
#     Parameter("k40", 2.00e-03)    #
#     Parameter("k41", 5.00e-05)    #
#     Parameter("k42", 1.00e-04)    #
#     Parameter("k43", 5.00e-04)    #
#     Parameter("k44", 5.00e-04)    #
#     Parameter("k45", 5.00e-05)    #
#     Parameter("k46", 2.50e-03)    #
#     Parameter("k47", 2.50e-03)    #
#     Parameter("k48", 2.50e-03)    #
#     Parameter("k49", 4.00e-02)    #
#     Parameter("k50", 2.50e-03)    #
#     Parameter("k51", 5.00e-08)    #
#     Parameter("k52", 5.00e-07)    #
#     Parameter("k53", 5.00e-05)    #
#     Parameter("k54", 1.00e-02)    #
#     Parameter("k55", 5.00e-08)    #
#     Parameter("k56", 5.00e-05)    #
#     Parameter("k57", 5.00e-03)    #
#     Parameter("k58", 5.00e-05)    #
#     Parameter("k59", 5.00e-04)    #
#     Parameter("k60", 1.00e-04)    #
#     Parameter("k61", 1.50e+00)    #
#     Parameter("k62", 1.00e-03)    #
#     Parameter("k63", 9.40e-04)    #
#     Parameter("k64", 2.00e-02)    #
#     Parameter("k65", 9.50e+00)    #
#     Parameter("k66", 1.00e+01)    #
#     Parameter("k67", 5.00e-03)    #
#     Parameter("k68", 5.00e-02)    #
#     Parameter("k69", 8.00e-04)    #
#     Parameter("k70", 6.00e+00)    #
#     Parameter("k71", 4.00e-03)    #
#     Parameter("k72", 1.00e-08)    #
#     Parameter("k73", 7.72e-01)    #
#     Parameter("k74", 5.56e-02)    #
#     Parameter("k75", 2.00e-02)    #

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

set_volume(1.0e-20)


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
for i in range(len(model.parameters)):
	print str(i)+":", model.parameters[i], "=", model.parameters[i].value
print
#   
# for initial_conditions in model.initial_conditions:
# 	print initial_conditions
# print
#   

   
# for rules in model.rules:
# 	print rules
# print
#   
# for i in range(len(model.species)):
#     print str(i) + ":"
#     print model.species[i]
#    
for i in range(len(model.odes)):
    print str(i) + ":", model.odes[i]
print

for e in model.expressions:
    print e, e.expand_expr()
print

for obs in model.observables:
    print obs, obs.species, obs.coefficients
print
# quit()
#  
#   

set_dna_damage(0.0 * Na_V)
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
set_dna_damage(0.005 * Na_V)
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

