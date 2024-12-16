import numpy as np
from amuse.units import units

#--------------------------------------------------------------------------------------------------------------------
# TRES general settings
REPORT_USER_WARNINGS = True
REPORT_DEBUG = False
REPORT_DT = False 
REPORT_SN_EVOLUTION = False
REPORT_TRIPLE_EVOLUTION = False 
MAKE_PLOTS = False

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = False
REPORT_MASS_TRANSFER_STABILITY = False

GET_GYRATION_RADIUS_FROM_STELLAR_CODE = False
GET_AMC_FROM_STELLAR_CODE = False
no_stellar_evolution = False
USE_MESA_AS_STELLAR_CODE = False

#--------------------------------------------------------------------------------------------------------------------
#TRES constants
time_step_factor_stable_mt = 0.01 #1% mass loss during mass transfer
# lowering this to 0.005 makes the code twice as slow
time_step_factor_ecc = 0.01
#Rl_fraction = 0.8
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
#maximum_radius_change_factor = 0.005
error_dr = 0.05 #0.01
minimum_time_step = 1.e-9 |units.Myr
if USE_MESA_AS_STELLAR_CODE:
    minimum_time_step = 1.e-7 |units.Myr

max_mass = 100 |units.MSun
min_mass = 0.08 |units.MSun # for primary stars

maximum_time_step_factor = 100.
maximum_time_step_factor_after_stable_mt = 5. 
time_step_factor_find_RLOF = 0.5
#Rl_fraction = 0.9#1.0-10.*error_dr # ratio or star radius over Roche lobe at which time step is decreased
                              # radius grows maximally by error_dr
time_step_factor_kozai = 0.025 # 0.2*0.1, 0.2-> for error in kozai timescale, 0.1-> 10 steps per cycle
kozai_type_factor = 10.
maximum_time_step = np.inf|units.Myr

kanonical_neutron_star_mass = 1.4|units.MSun
fall_back_mass = 41 |units.MSun

#--------------------------------------------------------------------------------------------------------------------
#TPS general settings
REPORT_TPS = False
REPORT_USER_WARNINGS_TPS = False
EXCLUDE_SSO = True #in order to not simulate systems with exoplanet or brown dwarf secondaries and tertiaries

#--------------------------------------------------------------------------------------------------------------------
#TPS constants
precision = 1.e-10 
absolute_max_mass = 100 |units.MSun
#  for secondaries and tertiaries
if EXCLUDE_SSO:
    absolute_min_mass = 0.0075|units.MSun 
else:
    absolute_min_mass = 0.2|units.MJupiter  

#--------------------------------------------------------------------------------------------------------------------
