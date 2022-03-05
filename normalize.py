from parameters import *
import logging
import math

DEBYE_LENGTH = math.sqrt(epsilon0*epsilon*k*T/q/q/ni)

logging.info("PAR: The Debye Length is %.3f [cm]" % DEBYE_LENGTH)
# normalize parameter
## for location
logging.info("NORMALIZATION: The normalization constant:")
LOC_NORMAL           = DEBYE_LENGTH
logging.info("NORMALIZATION: Locatoin X:                %.5e [cm]" % LOC_NORMAL)
## for time
TIME_NORMAL          =   DEBYE_LENGTH*DEBYE_LENGTH/D_0
# TIME_NORMAL          =   4.483e-4
logging.info("NORMALIZATION: Time t:                    %.5e [s]" % TIME_NORMAL)
## for potential
PHI_NORMAL            =   Vt
logging.info("NORMALIZATION: ElectroStatic potential :  %.5e [V]" % PHI_NORMAL)
## for applied voltage
VA_NORMAL           =   Vt
logging.info("NORMALIZATION: Applied Voltage Va :       %.5e [V]" % VA_NORMAL)
## for electric field
E_NORMAL            =   Vt/DEBYE_LENGTH
# E_NORMAL            =   7.295
logging.info("NORMALIZATION: ElectricField E :          %.5e [V/cm]" % E_NORMAL)
## for carrier density
CARRIER_NORMAL      =   ni
logging.info("NORMALIZATION: Carrier Density n & p :    %.5e [cm^-3]" % CARRIER_NORMAL)
## for doping
DOPPING_NORMAL      =   ni
logging.info("NORMALIZATION: Dopping Concentration:     %.5e [cm^-3]" % DOPPING_NORMAL)
## for current density
CURRENT_DEN_NORMAL  =   q*D_0*ni/DEBYE_LENGTH
# CURRENT_DEN_NORMAL  =   1.76e-8
logging.info("NORMALIZATION: Current Density:           %.5e [A/cm^2]" % CURRENT_DEN_NORMAL)
DIFF_DEN_NORMAL     =   D_0
logging.info("NORMALIZATION: diffusion coefficient:     %.5e [cm^2/s]" % DIFF_DEN_NORMAL)
MOBILITY_NROMAL     =   D_0/Vt
# MOBILITY_NROMAL     =   1
logging.info("NORMALIZATION: mobility:                  %.5e [cm^2/V/s]" % MOBILITY_NROMAL)
