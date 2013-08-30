import fan_modules as m
import pylab as pl
from pysb.integrate import odesolve
from sympy import sympify

# Set to False to plot I, Mdm2, p21, p53
cyclins = True

#################################
## Force PySB to generate ODEs ##
#################################
t = pl.linspace(0,3000, num = 300)
y = odesolve(m.model,t)

###############
## Overrides ##
###############
_source = "s30"
#m.model.odes[24] = sympify("k60*s32 + k61*s28 - k62*s24 - k73*s24*s25*s28 + k73*s24*s25*s31 - k74*s24*s25")
m.model.odes[25] = sympify("k63*" + _source + " + k66*s27**9/(k65**9 + s27**9) - k64*s25") #!!!
m.model.odes[27] = sympify("(k70*s24*s28)/(1 + k71*s24*s25) - k67*s27") #!!!
m.model.odes[28] = sympify("-k72*s28") #!!!
#m.model.odes[18] = sympify("k40*s32 + k41/(1 + k42*s23) - k43*s18 - k44*s18*s5")
#m.model.odes[23] = sympify("-k45*s21*s23 + k55*s22 + k56*s32 - k57*s23 + k58/(1 + k59*s18)")

###########
## Solve ##
###########
t = pl.linspace(0,3000, num = 300)
y = odesolve(m.model,t)

##### NOTE: NORMALIZING SO WE CAN COMPARE TO IWAMOTO ET AL.
y['OBSCycA']=y['OBSCycA']/max(y['OBSCycA'])
y['OBSCycE']=y['OBSCycE']/max(y['OBSCycE'])
y['OBSE2F']=y['OBSE2F']/max(y['OBSE2F'])
y['OBSp53']=y['OBSp53']/max(y['OBSp53'])
y['OBSp27']=y['OBSp27']/max(y['OBSp27'])
y['OBSI']=y['OBSI']/max(y['OBSI'])
y['OBSMdm']=y['OBSMdm']/max(y['OBSMdm'])
y['OBSp21']=y['OBSp21']/max(y['OBSp21'])

##########
## Plot ##
##########
pl.ion()
pl.figure()

pl.plot(y['OBSp53'], label="p53")
if cyclins:
	pl.plot(y['OBSCycA'], label="CycA")
	pl.plot(y['OBSCycE'], label="CycE")
	pl.plot(y['OBSE2F'], label="E2F")
	pl.plot(y['OBSp27'], label="p27")
else: 
	pl.plot(y['OBSI'], label="I")
	pl.plot(y['OBSMdm'], label="Mdm2")
	pl.plot(y['OBSp21'], label="p21")
#pl.legend()

