"""An stochastic implementation of the model from:

Mathematical modeling and sensitivity analysis of G1/S phase in the cell cycle
including the DNA-damage signal transduction pathway.  Kazunari Iwamoto, 
Yoshihiko Tashima, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
BioSystems. 2008 May 23;94(2008):109-117.
doi:10.1016/j.biosystems.2008.05.016

http://www.sciencedirect.com/science/article/pii/S0303264708001330

Implemented by: Corey Hayford (from David Fan)
"""

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
     
    # ** In functions - should not be changed **
#     model.parameters['k42'].value
#     model.parameters['k43'].value
#     model.parameters['k70'].value 
#     model.parameters['k59'].value
         
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