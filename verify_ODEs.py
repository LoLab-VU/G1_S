"""
This is a utility script that compares ODEs generated for G1/S models to the 
ODEs present in Table B.1 in Iwamoto et al. 

Usage:

run run_fan_model.py	# Run model once to generate ODEs, or force it with
						#	pysb.bng.generate_equations(model)
run verify_ODEs.py

IMPORTANT: "_source" needs to be set to the variable corresponding to the 
species "__source()". If other species are added it may be something other than
"s30".
"""

import fan_modules as m
import pylab as pl
from sympy import sympify

# Stuff for colors
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        
# Hard-coded ODEs
_source = "s30"
ode = []
ode.append(sympify("k1*" + _source + "+ k3*s5 - (k2*s0 + k4*s0*s3)"))
ode.append(sympify("k5*s21 + k8*s6 - (k6*s1 + k7*s1*s4)"))
ode.append(sympify("k9*s26 + k12*s8 - (k10*s2 + k11*s2*s4)"))
ode.append(sympify("k3*s5 + k13*s5 - (k4*s0*s3)"))
ode.append(sympify("k8*s6 + k12*s8 + k14*s9 + k15*s8 + k16*s6 + k17*s7**2 - (k7*s1*s4 + k11*s2*s4)"))
ode.append(sympify("k4*s0*s3 + k19*s15 + k21*s11 - (k3*s5 + k13*s5 + k18*s5*s14 + k20*s5*s10 + k44*s5*s18)"))
ode.append(sympify("k7*s1*s4 + k23*s7- (k8*s6 + k22*s6*s7 + k16*s6)"))
ode.append(sympify("k22*s6*s7 + k25*s12 + k27*s16 - (k23*s7 + k24*s7*s10 + k26*s7*s14 + k17*s7**2)"))
ode.append(sympify("k11*s2*s4 + k29*s9 - (k12*s8 + k28*s8*s9 + k15*s8)"))
ode.append(sympify("k28*s8*s9 + k31*s13 + k33*s17 - (k29*s9 + k30*s9*s10 + k32*s9*s14 + k14*s9)"))
ode.append(sympify("k34*" + _source + " + k21*s11 + k25*s12 + k31*s13 - (k35*s7*s10 + k36*s9*s10 + k20*s5*s10 + k24*s7*s10 + k30*s9*s10)"))
ode.append(sympify("k20*s5*s10 - (k21*s11)"))
ode.append(sympify("k24*s7*s10 - (k25*s12)"))
ode.append(sympify("k30*s9*s10 - (k31*s13)"))
ode.append(sympify("k37*" + _source + " + k38*s24 + k19*s15 + k27*s16 + k33*s17 - (k39*s14 + k18*s5*s14 + k26*s7*s14 + k32*s9*s14)"))
ode.append(sympify("k18*s5*s14 - (k19*s15)"))
ode.append(sympify("k26*s7*s14 - (k27*s16)"))
ode.append(sympify("k32*s9*s14 - (k33*s17)"))
ode.append(sympify("k40*" + _source + "+ k41/(1+k42*s23)-(k43*s18+k44*s5*s18)"))
ode.append(sympify("k45*s21*s23 - (k46*s5*s19 + k47*s11*s19 + k48*s15*s19)"))
ode.append(sympify("k46*s5*s19 + k47*s11*s19 + k48*s15*s19 - (k49*s7*s20 + k50*s9*s20)"))
ode.append(sympify("k49*s7*s20 + k50*s9*s20 + k51*s21 + k52*" + _source + "- (k45*s21*s23 + k53*s21 + k54*s9*s21)"))
ode.append(sympify("k49*s7*s20 + k50*s9*s20 - (k55*s22)"))
ode.append(sympify("k56*" + _source + "+ k58/(1+k59*s18)+k55*s22 - (k57*s23 + k45*s21*s23)"))
ode.append(sympify("k60*" + _source + "+ k61*s28 - (degradations24*s25 + k62*s24)"))
ode.append(sympify("k63*" + _source + "+ k66*s27^9/(k65^9 + s27^9) - (k64*s25)"))
ode.append(sympify("k68*s21 - (k69*s26)"))
ode.append(sympify("k70*s24*s28/(1+k71*s24*s25) - (k67*s27)"))

# Check each ODE
for i in range(28):
	if ode[i] != m.model.odes[i]:
		print bcolors.WARNING + "Warning: ODE " + str(i) + " does not match:" + bcolors.ENDC
		print "Original:  " + ode[i]
		print "Generated: " + m.model.odes[i]
	
	
