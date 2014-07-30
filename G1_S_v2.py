"""An implementation of the model from:

Mathematical modeling and sensitivity analysis of G1/S phase in the cell cycle
including the DNA-damage signal transduction pathway.  Kazunari Iwamoto, 
Yoshihiko Tashima, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
BioSystems. 2008 May 23;94(2008):109-117.
doi:10.1016/j.biosystems.2008.05.016

http://www.sciencedirect.com/science/article/pii/S0303264708001330

Implemented by: Corey Hayford (from David Fan)
"""

from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants

# Model()

# def set_dna_damage(damage):
# 	model.parameters['DDS_0'].value = damage
	
def declare_monomers():
	"""Declare the monomers in the Iwamoto model
	'phos' is the phosphorylation state
	'b' is the binding site between p## and CDKs
	'c' is the binding site between cyclins and CDKs
	'i' is the binding site for inhibitors"""
	
# 	Monomer("Signal")
# 	Monomer("SignalDamp")
	# **Cyclins**
# 	Monomer("CycD", [        'c'])
# 	Monomer("CycE", [        'c'])
# 	Monomer("CycA", [        'c'])
# 	
# 	# **Cyclin-dependent kinases**
# 	Monomer("CDK4_6",['phos','b','c'], {'phos':['u','p']})
# 	Monomer("CDK2",	['phos','b','c'], {'phos':['u','p']})
# 	
# 	# **Regulatory proteins**
# 	Monomer("p27", 	[    'b'])
# 	Monomer("p21", 	[    'b'])
# 	Monomer("p16", 	[    'b'])
# 	Monomer("p53",	[    'b'])
# 	Monomer("Rb", 	['phos','b'], {'phos':['u','pp','pppp']})
# 	Monomer("E2F", 	[    'b'])
# 	Monomer("Mdm2",	['i'])
# 	Monomer("X")
# 	Monomer("I")

def declare_parameters():
	Parameter("Y0_0", 3.00e-2)		# 	CycD
	Parameter("Y1_0", 1.00e-3)		# 	CycE
	Parameter("Y2_0", 4.00e-5)		# 	CycA
	Parameter("Y3_0", 5.00e+00)		# 	CDK4/6
	Parameter("Y4_0", 1.35e+01)		# 	CDK2
	Parameter("Y5_0", 2.00e+00)		# 	CycD/CDK4/6
	Parameter("Y6_0", 1.00e-03)		# 	CycE/CDK2
	Parameter("Y7_0", 1.00e-03)		# 	CycE/CDK2-P
	Parameter("Y8_0", 4.00e-04)		# 	CycA/CDK2
	Parameter("Y9_0", 1.00e-04)		# 	CycA/CDK2-P
	Parameter("Y10_0", 1.00e+01)	# 	p27
	Parameter("Y11_0", 1.00e-03)	# 	p27/CycD/CDK4/6 (-P)
	Parameter("Y12_0", 1.00e+00)	# 	p27/CycE/CDK2 (-P)
	Parameter("Y13_0", 1.00e-04)	# 	p27/CycA/CDK2 (-P)
	Parameter("Y14_0", 0.00e+00)	# 	p21
	Parameter("Y15_0", 0.00e+00)	# 	p21/CycD/CDK4/6 (-P)
	Parameter("Y16_0", 0.00e+00)	# 	p21/CycE/CDK2 (-P)
	Parameter("Y17_0", 0.00e+00)	# 	p21/CycA/CDK2 (-P)
	Parameter("Y18_0", 1.00e-03)	# 	p16
	Parameter("Y19_0", 1.95e+00)	# 	Rb/E2F
	Parameter("Y20_0", 1.00e-03)	# 	Rb-PP/E2F
	Parameter("Y21_0", 0.00e+00)	# 	E2F
	Parameter("Y22_0", 1.00e-02)	# 	Rb-PPPP
	Parameter("Y23_0", 5.00e-02)	# 	Rb
	Parameter("Y24_0", 2.65e-02)	# 	p53
	Parameter("Y25_0", 2.35e-04)	# 	Mdm2
	Parameter("Y26_0", 1.00e-04)	# 	X
	Parameter("Y27_0", 0.00e+00)	# 	I
# 	Parameter("DDS_0")				# 	DDS
	
	# Kinetic parameters of proposed model
	Parameter("G1_S_k1", 5.00e-03)	#
	Parameter("G1_S_k2", 5.00e-04)	#
	Parameter("G1_S_k3", 5.00e-03)	# 
	Parameter("G1_S_k4", 2.50e-03)	# 
	Parameter("G1_S_k5", 7.50e-02)	# 
	Parameter("G1_S_k6", 2.50e-03)	#
	Parameter("G1_S_k7", 1.25e-03)	#
	Parameter("G1_S_k8", 2.50e-04)	#
	Parameter("G1_S_k9", 8.00e-04)	#
	Parameter("G1_S_k10", 5.00e-04)	#
	Parameter("G1_S_k11", 1.00e-03)	#
	Parameter("G1_S_k12", 2.00e-04)	#
	Parameter("G1_S_k13", 5.00e-04)	#
	Parameter("G1_S_k14", 5.00e-04)	#
	Parameter("G1_S_k15", 5.00e-04)	#
	Parameter("G1_S_k16", 5.00e-04)	#
	Parameter("G1_S_k17", 2.00e-03)	#
	Parameter("G1_S_k18", 5.00e-04)	#
	Parameter("G1_S_k19", 5.00e-03)	#
	Parameter("G1_S_k20", 5.00e-04)	#
	Parameter("G1_S_k21", 5.00e-05)	#
	Parameter("G1_S_k22", 2.50e-02)	#
	Parameter("G1_S_k23", 1.75e-03)	#
	Parameter("G1_S_k24", 2.25e-02)	# 
	Parameter("G1_S_k25", 1.75e-04)	#
	Parameter("G1_S_k26", 2.25e-02)	#
	Parameter("G1_S_k27", 1.75e-04)	#
	Parameter("G1_S_k28", 1.90e-02)	#
	Parameter("G1_S_k29", 5.00e-04)	#
	Parameter("G1_S_k30", 2.50e-03)	# 
	Parameter("G1_S_k31", 1.75e-04)	#
	Parameter("G1_S_k32", 2.50e-03)	#
	Parameter("G1_S_k33", 1.75e-04)	#
	Parameter("G1_S_k34", 5.00e-08)	#
	Parameter("G1_S_k35", 1.00e-02)	# 
	Parameter("G1_S_k36", 1.50e-03)	# 
	Parameter("G1_S_k37", 5.00e-05)	#
	Parameter("G1_S_k38", 1.00e-02)	#
	Parameter("G1_S_k39", 5.00e-03)	#
	Parameter("G1_S_k40", 2.00e-03)	#
	Parameter("G1_S_k41", 5.00e-05)	#
	Parameter("G1_S_k42", 1.00e-04)	#
	Parameter("G1_S_k43", 5.00e-04)	#
	Parameter("G1_S_k44", 5.00e-04)	#
	Parameter("G1_S_k45", 5.00e-05)	#
	Parameter("G1_S_k46", 2.50e-03)	#
	Parameter("G1_S_k47", 2.50e-03)	#
	Parameter("G1_S_k48", 2.50e-03)	#
	Parameter("G1_S_k49", 4.00e-02)	#
	Parameter("G1_S_k50", 2.50e-03)	#
	Parameter("G1_S_k51", 5.00e-08)	#
	Parameter("G1_S_k52", 5.00e-07)	#
	Parameter("G1_S_k53", 5.00e-05)	#
	Parameter("G1_S_k54", 1.00e-02)	#
	Parameter("G1_S_k55", 5.00e-08)	#
	Parameter("G1_S_k56", 5.00e-05)	#
	Parameter("G1_S_k57", 5.00e-03)	#
	Parameter("G1_S_k58", 5.00e-05)	#
	Parameter("G1_S_k59", 5.00e-04)	#
	Parameter("G1_S_k60", 1.00e-04)	#
	Parameter("G1_S_k61", 1.50e+00)	#
	Parameter("G1_S_k62", 1.00e-03)	#
	Parameter("G1_S_k63", 9.40e-04)	#
	Parameter("G1_S_k64", 2.00e-02)	#
	Parameter("G1_S_k65", 9.50e+00)	#
	Parameter("G1_S_k66", 1.00e+01)	#
	Parameter("G1_S_k67", 5.00e-03)	#
	Parameter("G1_S_k68", 5.00e-02)	#
	Parameter("G1_S_k69", 8.00e-04)	#
	Parameter("G1_S_k70", 6.00e+00)	#
	Parameter("G1_S_k71", 4.00e-03)	#
	Parameter("G1_S_k72", 1.00e-08)	#
	Parameter("G1_S_k73", 7.72e-01)	#
	Parameter("G1_S_k74", 5.56e-02)	#
	Parameter("G1_S_k75", 2.00e-02)	#
	Parameter("One", 1)			#
	
	alias_model_components()
	
def declare_initial_conditions():
	x=1
# 	Initial(CycD(c=None), Y0_0)	#0
# 	Initial(CycE(c=None), Y1_0)	#1
# 	Initial(CycA(c=None), Y2_0)	#2
# 	Initial(CDK4_6(phos='u',b=None,c=None), Y3_0)	#3
# 	Initial(CDK2(phos='u',b=None,c=None), Y4_0)	#4
# 	Initial(CycD(c=1) % CDK4_6(phos='p',b=None,c=1), Y5_0)	#5
# 	Initial(CycE(c=1) % CDK2(phos='u',b=None,c=1), Y6_0)	#6
# 	Initial(CycE(c=1) % CDK2(phos='p',b=None,c=1), Y7_0)	#7
# 	Initial(CycA(c=1) % CDK2(phos='u',b=None,c=1), Y8_0)	#8
# 	Initial(CycA(c=1) % CDK2(phos='p',b=None,c=1), Y9_0)	#9
# 	Initial(p27(b=None), Y10_0)	#10
# 	Initial(CycD(c=2) % CDK4_6(phos='p',b=1,c=2) % p27(b=1), Y11_0)	#11
# 	Initial(CycE(c=2) % CDK2(phos='p',b=1,c=2) % p27(b=1), Y12_0)	#12
# 	Initial(CycA(c=2) % CDK2(phos='p',b=1,c=2) % p27(b=1), Y13_0)	#13
# 	Initial(p21(b=None), Y14_0)	#14
# 	Initial(CycD(c=2) % CDK4_6(phos='p',b=1,c=2) % p21(b=1), Y15_0)	#15
# 	Initial(CycE(c=2) % CDK2(phos='p',b=1,c=2) % p21(b=1), Y16_0)	#16
# 	Initial(CycA(c=2) % CDK2(phos='p',b=1,c=2) % p21(b=1), Y17_0)	#17
# 	Initial(p16(b=None), Y18_0)	#18
# 	Initial(Rb(b=1,phos='u') % E2F(b=1), Y19_0)	#19
# 	Initial(Rb(b=1,phos='pp') % E2F(b=1), Y20_0)	#20
# 	Initial(E2F(b=None), Y21_0)	#21
# 	Initial(Rb(b=None,phos='pppp'), Y22_0)	#22
# 	Initial(Rb(b=None,phos='u'), Y23_0)	#23
# 	Initial(p53(b=None), Y24_0)	#24
# 	Initial(Mdm2(b=None), Y25_0)	#25
# 	Initial(X(), Y26_0)	#26
# 	Initial(I(), Y27_0)	#27
# 	Initial(Signal(), DDS_0) #28
# 	Initial(SignalDamp(), DDS_0) #29



# Figure 6 (CycE, CycA, p27, E2F with no DNA-damage)
def declare_observables():
	x=1
# 	Observable("OBS_p27", p27(b=None))
# 	Observable("OBS_p53", p53(b=None))
# 	Observable("OBS_E2F", E2F(b=None))
# 	Observable("OBS_CycE", CycE(c=None))
# 	Observable("OBS_CycA", CycA(c=None))
# 	Observable("OBS_Int", I())
# 	Observable("OBS_Mdm2", Mdm2(b=None))
# 	Observable("OBS_p21", p21(b=None))
# 	Observable("OBS_p16", p16(b=None))
# 	Observable("OBS_Rb", Rb(b=None, phos='u'))
# 	Observable("OBS_CycD_CDK4_6", CycD(c=1) % CDK4_6(phos='p',b=None,c=1))
# 	Observable("signal", Signal())
# 	Observable("signal_damp", SignalDamp())

###############
## Functions ##
## Implemented by Leonard Harris and Corey Hayford ##

def declare_functions():
	x=1
# 	Expression("create_p16", sympify("G1_S_k41/((One + G1_S_k42*OBS_Rb) - (G1_S_k43*OBS_p16 + G1_S_k44*OBS_p16*OBS_CycD_CDK4_6))"))
# 	Expression("create_Rb", sympify("G1_S_k58/(One + G1_S_k59*OBS_p16)"))
# 	Expression("create_Mdm2", sympify("G1_S_k66*OBS_Int**9/(G1_S_k65**9 + OBS_Int**9)"))
# 	Expression("create_Int", sympify("(G1_S_k70*OBS_p53*signal)/(One + G1_S_k71*OBS_p53*OBS_Mdm2)"))
# 	Expression("sig_deg", sympify("G1_S_k74 - G1_S_k73*(signal-signal_damp)"))
# 	Expression("kdamp_DDS0", sympify("G1_S_k75*DDS_0"))

###########
## Rules ##
###########

def simulate_signal_degradation():
	x=1
# 	##### DNA Damage Signal ####
# 	Rule('Signal_Degrade', Signal() >> None, G1_S_k72)
# 	Rule('Signal_Damp', SignalDamp() >> None, kdamp_DDS0)
	
# Figure 1 in Iwamoto et al. 
# Reactions 1-14

def p16_p27_inhibition():
	# Synthesis and Degradation
	synthesize_degrade_table([[CycD(c=None),	G1_S_k1,		G1_S_k2],	#1,2
							  [p27(b=None),		G1_S_k34,	None],	#10
							  [p16(b=None),		G1_S_k40,	None]])	#11,13
	# Inhibition by p27
	bind_table_complex([[										p27(b=None)	],
						[CycD(c=1) % CDK4_6(phos='p',b=None,c=1),	(G1_S_k20, G1_S_k21)	],	#5
						[CycE(c=1) % CDK2( phos='p',b=None,c=1),	(G1_S_k24, G1_S_k25)	],	#7
						[CycA(c=1) % CDK2( phos='p',b=None,c=1),	(G1_S_k30, G1_S_k31)	]],	#9
						'b', 'b')
						
	#TODO: Clean up below
	equilibrate(CycD(c=None) + CDK4_6(phos='u',b=None,c=None), CycD(c=1) % CDK4_6(phos='p',b=None,c=1), [G1_S_k3, G1_S_k4])			#3	Activation of CDK4/6 by CycD
	Rule("R4", CycD(c=1) % CDK4_6(phos='p',b=None,c=1) >> CDK4_6(phos='u',b=None,c=None), G1_S_k13)	#4	Degradation of bound CycD
	
	Rule("Degrade_p16_from_CycD_CDK4_6", p16(b=None) + CycD(c=1) % CDK4_6(phos='p',b=None,c=1) >>  p16(b=None), G1_S_k44)   		#14 Degradation of p16/CycD/CDK4/6 complex
	catalyze_one_step(CycE(c=1) % CDK2(phos='p',b=None,c=1), p27(b=None), None, G1_S_k35)		#6	Degradation of p27 by CycE/CDK2 complex
# 	catalyze_one_step(CycA(c=1) % CDK2(phos='p',b=None,c=1), p27(b=None), None, G1_S_k36)		#8	Degradation of p27 by CycA/CDK2 complex
	Rule("Degrade_p27_from_CycA_CDK2P", p27(b=None) + CycA(c=1) % CDK2(phos='p',b=None,c=1) >> CycA(c=1) % CDK2(phos='p',b=None,c=1), G1_S_k36)
# 	catalyze_one_step(Rb(b=None,phos='u'), p16(b=None), None, G1_S_k41)							#12 Inhibition of p16 by pRb
	Rule("Create_p16", None >> p16(b=None), create_p16)
	

# Figure 2 in Iwamoto et al.
# Reactions 15-30
def CDK2_activation():
	# Synthesis and Degradation
	synthesize_degrade_table([[CycE(c=None),	None,	G1_S_k6],	#16
							  [CycA(c=None),	None,	G1_S_k10],	#19
							  [X(), 			None,	G1_S_k69]])	#30
	# Activation of CDK2
	bind_table_complex([[				CDK2(phos='u',b=None,c=None)	],	
						[CycE(c=None),	(G1_S_k7, 	G1_S_k8)					],	#17
						[CycA(c=None),	(G1_S_k11, 	G1_S_k12)				]],	#20
						'c','c')
	#TODO: Clean up below
	catalyze_one_step(E2F(b=None), None, X(), 			G1_S_k68)	#29 Activation of "X" synthesis by E2F
	catalyze_one_step(E2F(b=None), None, CycE(c=None), 	G1_S_k5)		#15 Activation of CycE synthesis by E2F
	catalyze_one_step(X(), 		   None, CycA(c=None), 	G1_S_k9)		#18 Activation of CycA synthesis by "X"
	Rule("R21", CycA(c=1) % CDK2(phos='p',b=None,c=1) >> CDK2(phos='u',b=None,c=None), G1_S_k14)			#21	Degradation of bound CycA
	Rule("R22", CycA(c=1) % CDK2(phos='u',b=None,c=1) >> CDK2(phos='u',b=None,c=None), G1_S_k15)			#22	Degradation of bound CycA
	Rule("R23", CycE(c=1) % CDK2(phos='u',b=None,c=1) >> CDK2(phos='u',b=None,c=None), G1_S_k16)			#23	Degradation of bound CycE
	Rule("R26", CycE(c=1) % CDK2(phos='p',b=None,c=1) >> CycE(c=1) % CDK2(phos='u',b=None,c=1), G1_S_k23)	#26	Dephosphorylation of CycE/CDK2 complex
	Rule("R28", CycA(c=1) % CDK2(phos='p',b=None,c=1) >> CycA(c=1) % CDK2(phos='u',b=None,c=1), G1_S_k29)	#28	Dephosphorylation of CycA/CDK2 complex
	catalyze_one_step(CycE(c=1) % CDK2(phos='p',b=None,c=1), CycE(c=1) % CDK2(phos='p',b=None,c=1), CDK2(phos='u',b=None,c=None), G1_S_k17)			#24 Self-inactivation of CycE/CDK2
	catalyze_one_step(CycE(c=1) % CDK2(phos='p',b=None,c=1), CycE(c=1) % CDK2(phos='u',b=None,c=1), CycE(c=1) % CDK2(phos='p',b=None,c=1), G1_S_k22)	#25 Self-phosphorylation of CycE/CDK2
	catalyze_one_step(CycA(c=1) % CDK2(phos='p',b=None,c=1), CycA(c=1) % CDK2(phos='u',b=None,c=1), CycA(c=1) % CDK2(phos='p',b=None,c=1), G1_S_k28)	#27 Self-phosphorylation of CycA/CDK2

# Figure 3 in Iwamoto et al.
# Reactions 31-44
def Rb_E2F_activation():
	# Synthesis and Degradation
	synthesize_degrade_table([[E2F(b=None),			G1_S_k52, G1_S_k53],	#38,39
							  [Rb(phos='u',b=None),	G1_S_k56, G1_S_k57]])	#42,43

	#32	Phosphorylation of Rb/E2F complex by CycD/CDK4/6 complex	TODO: change to "catalyze_state"
	catalyze_one_step(CycD(c=1) % CDK4_6(phos='p',b=None,c=1), Rb(phos='u',b=2) % E2F(b=2), Rb(phos='pp',b=2) % E2F(b=2), G1_S_k46)
	#33 Phosphorylation of Rb/E2F complex by p27/CycD/CDK4/6 complex	TODO: change to "catalyze_state"
	catalyze_one_step(CycD(c=1) % CDK4_6(phos='p',b=2,c=1) % p27(b=2), Rb(phos='u',b=2) % E2F(b=2), Rb(phos='pp',b=2) % E2F(b=2), G1_S_k47)
	#34 Phosphorylation of Rb/E2F complex by p21/CycD/CDK4/6 complex	TODO: change to "catalyze_state"
	catalyze_one_step(CycD(c=1) % CDK4_6(phos='p',b=2,c=1) % p21(b=2), Rb(phos='u',b=2) % E2F(b=2), Rb(phos='pp',b=2) % E2F(b=2), G1_S_k48)
	#35 Cleavage of Rb/E2F complex by CycE/CDK2 complex
	#catalyze_one_step(CycE(c=1) % CDK2(phos='p',b=None,c=1), Rb(phos='pp',b=2) % E2F(b=2), Rb(phos='pppp',b=None) + E2F(b=None), G1_S_k49)
	Rule("R35", CycE(c=1) % CDK2(phos='p',b=None,c=1) + Rb(phos='pp',b=2) % E2F(b=2) >> Rb(phos='pppp',b=None) + E2F(b=None) + CycE(c=1) % CDK2(phos='p',b=None,c=1), G1_S_k49)	#35
	#36 Cleavage of Rb/E2F complex by CycA/CDK2 complex
	#catalyze_one_step(CycA(c=1) % CDK2(phos='p',b=None,c=1), Rb(phos='pp',b=2) % E2F(b=2), Rb(phos='pppp',b=None) + E2F(b=None), G1_S_k50)
	Rule("R36", CycA(c=1) % CDK2(phos='p',b=None,c=1) + Rb(phos='pp',b=2) % E2F(b=2) >> Rb(phos='pppp',b=None) + E2F(b=None) + CycA(c=1) % CDK2(phos='p',b=None,c=1), G1_S_k50)	#36
	#37 Self-activation of E2F
	catalyze_one_step(E2F(b=None), None, E2F(b=None), G1_S_k51)
	#40 Degradation of E2F by CycA/CDK2 complex
	catalyze_one_step(CycA(c=1) % CDK2(phos='p',b=None,c=1), E2F(b=None), None, G1_S_k54)
	#44 Inhibition of Rb by p16
# 	catalyze_one_step(p16(b=None), Rb(phos='u',b=None), None, G1_S_k58)
	Rule("R31", E2F(b=None) + Rb(b=None,phos='u') >> E2F(b=2) % Rb(b=2,phos='u'), G1_S_k45)	#31
	Rule("R41", Rb(b=None, phos='pppp') >> Rb(b=None,phos='u'), G1_S_k55)	#41
	Rule("Create_Rb", None >> Rb(b=None, phos='u'), create_Rb)
	
# Figure 4 in Iwamoto et al.
# Reactions 45-58
def DNA_damage_pathway():
	# Synthesis and Degradation
	synthesize_degrade_table([[p21(b=None),		G1_S_k37, G1_S_k39],		#48,50
							  [p53(b=None),		G1_S_k60, G1_S_k62],		#51,53
							  [Mdm2(b=None),	G1_S_k63, G1_S_k64],		#54,55
							  [I(),				None, G1_S_k67]])	#58,57
							  
	# Inhibition by p21
	bind_table_complex([[										p21(b=None)	],
						[CycD(c=1) % CDK4_6(phos='p',b=None,c=1),	(G1_S_k18, G1_S_k19)	],	#45
						[CycE(c=1) % CDK2( phos='p',b=None,c=1),	(G1_S_k26, G1_S_k27)	],	#46
						[CycA(c=1) % CDK2( phos='p',b=None,c=1),	(G1_S_k32, G1_S_k33)	]],	#47
						'b','b')

	#49 Activation of p21 synthesis by p53
	catalyze_one_step(p53(b=None), None, p21(b=None), G1_S_k38)
	#52 Activation of p53 synthesis by DNA Damage Signal
	catalyze_one_step(Signal(), None, p53(b=None), G1_S_k61)
	#56 Activation of Mdm2 synthesis by "I"
# 	catalyze_one_step(I(), None, Mdm2(b=None), G1_S_k66)
	#59 Degradation of p53 by Mdm2
# 	catalyze_one_step(Mdm2(b=None), p53(b=None), None, G1_S_k74)
	Rule("Create_Mdm2", None >> Mdm2(b=None), create_Mdm2)
	Rule("Create_Int", None >> I(), create_Int)
	Rule('p53_Create_Mdm2', p53(b=None) + Mdm2(b=None) >> Mdm2(b=None), sig_deg)	
	
def set_volume(vol):
    Na_V = constants.N_A * vol              #1/[]
         
    # ***2nd Order Reactions*** (divide by Avogadro's Number * Volume)
    model.parameters['G1_S_k3'].value /= Na_V
    model.parameters['G1_S_k20'].value /= Na_V
    model.parameters['G1_S_k24'].value /= Na_V
    model.parameters['G1_S_k30'].value /= Na_V
    model.parameters['G1_S_k13'].value /= Na_V
    model.parameters['G1_S_k44'].value /= Na_V
    model.parameters['G1_S_k35'].value /= Na_V
    model.parameters['G1_S_k36'].value /= Na_V
    model.parameters['G1_S_k7'].value /= Na_V
    model.parameters['G1_S_k11'].value /= Na_V
    model.parameters['G1_S_k17'].value /= Na_V
    model.parameters['G1_S_k22'].value /= Na_V
    model.parameters['G1_S_k28'].value /= Na_V
    model.parameters['G1_S_k46'].value /= Na_V
    model.parameters['G1_S_k47'].value /= Na_V
    model.parameters['G1_S_k48'].value /= Na_V
    model.parameters['G1_S_k49'].value /= Na_V
    model.parameters['G1_S_k50'].value /= Na_V
    model.parameters['G1_S_k54'].value /= Na_V
    model.parameters['G1_S_k45'].value /= Na_V
    model.parameters['G1_S_k18'].value /= Na_V
    model.parameters['G1_S_k26'].value /= Na_V
    model.parameters['G1_S_k32'].value /= Na_V   
         
    # ***1st Order Reactions*** - LEAVE ALONE!
#     model.parameters['G1_S_k72'].value
#     model.parameters['G1_S_k4'].value
#     model.parameters['G1_S_k2'].value
#     model.parameters['G1_S_k21'].value
#     model.parameters['G1_S_k25'].value
#     model.parameters['G1_S_k31'].value
#     model.parameters['G1_S_k6'].value
#     model.parameters['G1_S_k10'].value
#     model.parameters['G1_S_k69'].value
#     model.parameters['G1_S_k8'].value
#     model.parameters['G1_S_k12'].value
#     model.parameters['G1_S_k68'].value
#     model.parameters['G1_S_k5'].value
#     model.parameters['G1_S_k9'].value
#     model.parameters['G1_S_k14'].value
#     model.parameters['G1_S_k15'].value
#     model.parameters['G1_S_k16'].value
#     model.parameters['G1_S_k23'].value
#     model.parameters['G1_S_k29'].value
#     model.parameters['G1_S_k53'].value
#     model.parameters['G1_S_k57'].value
#     model.parameters['G1_S_k51'].value
#     model.parameters['G1_S_k55'].value
#     model.parameters['G1_S_k39'].value
#     model.parameters['G1_S_k62'].value
#     model.parameters['G1_S_k64'].value
#     model.parameters['G1_S_k67'].value
#     model.parameters['G1_S_k19'].value
#     model.parameters['G1_S_k27'].value
#     model.parameters['G1_S_k33'].value
#     model.parameters['G1_S_k38'].value
#     model.parameters['G1_S_k61'].value
         
    # ***0th Order Reactions*** (multiply by Avogadro's Number * Volume)
    model.parameters['G1_S_k1'].value *= Na_V
    model.parameters['G1_S_k34'].value *= Na_V
    model.parameters['G1_S_k40'].value *= Na_V
    model.parameters['G1_S_k52'].value *= Na_V
    model.parameters['G1_S_k56'].value *= Na_V
    model.parameters['G1_S_k37'].value *= Na_V
    model.parameters['G1_S_k60'].value *= Na_V
    model.parameters['G1_S_k63'].value *= Na_V
      
         
    # ***In Functions *** (act differently)
    model.parameters['G1_S_k41'].value *= (Na_V * Na_V)
    model.parameters['G1_S_k44'].value /= Na_V
    model.parameters['G1_S_k58'].value *= (Na_V * Na_V)
    model.parameters['G1_S_k65'].value *= Na_V
    model.parameters['G1_S_k66'].value *= Na_V
    model.parameters['G1_S_k71'].value /= Na_V
    model.parameters['G1_S_k73'].value /= (Na_V * Na_V)
    model.parameters['G1_S_k74'].value /= Na_V
    model.parameters['G1_S_k75'].value /= Na_V
    model.parameters['One'].value *= Na_V
     
    # ** In functions - should not be changed **
#     model.parameters['G1_S_k42'].value
#     model.parameters['G1_S_k43'].value
#     model.parameters['G1_S_k70'].value 
#     model.parameters['G1_S_k59'].value
         
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