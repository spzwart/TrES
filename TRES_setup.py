from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles
from amuse.units import units
from seculartriple_TPS.interface import SecularTriple

from TRES_options import max_mass, absolute_min_mass, USE_MESA_AS_STELLAR_CODE
from MESA_setup import setup_mesa_stars

from interactions import *
from tidal_friction_constant import *

import numpy as np


#-------
#to initialize the triple object
def make_stars(inner_primary_mass, inner_secondary_mass, outer_mass):
    stars = Particles(3)
    stars.is_star = True
    stars.is_donor = False
    
    if inner_primary_mass < inner_secondary_mass:
        spare = inner_primary_mass
        inner_primary_mass = inner_secondary_mass
        inner_secondary_mass = spare     

    stars[0].mass = inner_primary_mass
    stars[1].mass = inner_secondary_mass
    stars[2].mass = outer_mass

    stars[0].initial_mass = inner_primary_mass
    stars[1].initial_mass = inner_secondary_mass
    stars[2].initial_mass = outer_mass

    return stars 
     
def make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
        inner_eccentricity, outer_eccentricity,
        inner_argument_of_pericenter, outer_argument_of_pericenter,
        inner_longitude_of_ascending_node, outer_longitude_of_ascending_node):     

    bins = Particles(2)
    bins.is_star = False
    bins.is_mt_stable = True
    bins.part_dt_mt = 1.
    bins.bin_type = bin_type['unknown'] #Unknown

    bins[0].child1 = stars[0]
    bins[0].child2 = stars[1]
    bins[0].child1.parent = bins[0]
    bins[0].child2.parent = bins[0]

    bins[0].semimajor_axis = inner_semimajor_axis
    bins[0].eccentricity = inner_eccentricity
    bins[0].argument_of_pericenter = inner_argument_of_pericenter
    bins[0].longitude_of_ascending_node = inner_longitude_of_ascending_node
    
    bins[0].mass_transfer_rate = 0.0 | units.MSun/units.yr
    bins[0].accretion_efficiency_mass_transfer = 1.0
    bins[0].accretion_efficiency_wind_child1_to_child2 = 0.0
    bins[0].accretion_efficiency_wind_child2_to_child1 = 0.0

    bins[1].child1 = stars[2]
    bins[1].child2 = bins[0]
    bins[1].child1.parent = bins[1]
    bins[1].child2.parent = bins[1]
    
    bins[1].semimajor_axis = outer_semimajor_axis
    bins[1].eccentricity = outer_eccentricity
    bins[1].argument_of_pericenter = outer_argument_of_pericenter                
    bins[1].longitude_of_ascending_node = outer_longitude_of_ascending_node
    
    bins[1].mass_transfer_rate = 0.0 | units.MSun/units.yr        
    bins[1].accretion_efficiency_mass_transfer = 1.0
    bins[1].accretion_efficiency_wind_child1_to_child2 = 0.0
    bins[1].accretion_efficiency_wind_child2_to_child1 = 0.0

    # binary evolutionary settings
    bins[0].specific_AM_loss_mass_transfer = 2.5 
    bins[1].specific_AM_loss_mass_transfer = 2.5

    return bins
    
def test_initial_parameters(inner_primary_mass, inner_secondary_mass, outer_mass,
        inner_semimajor_axis, outer_semimajor_axis,
        inner_eccentricity, outer_eccentricity,
        relative_inclination,
        inner_argument_of_pericenter, outer_argument_of_pericenter,
        inner_longitude_of_ascending_node):

    if max(inner_primary_mass, outer_mass) > max_mass:  
        print('error: masses not in allowed range')
        print('m1=',inner_primary_mass, 'm2=',inner_secondary_mass, 'm3=',outer_mass)
        print('should be below:', max_mass)
        print('max_mass settable in TRES_options.py')
        return False, 0,0
               
    if min(inner_secondary_mass, outer_mass) <= absolute_min_mass:  
        print('error: masses not in allowed range')
        print('m1=',inner_primary_mass, 'm2=',inner_secondary_mass, 'm3=',outer_mass)
        print('should be at least above:', absolute_min_mass)
        print('absolute_min_mass settable in TRES_options.py')
        print('substellar objects can be included through EXCLUDE_SSO in TRES_options.py')
        return False, 0,0

    if inner_semimajor_axis >= outer_semimajor_axis:
        print('error input parameters, should be:')
        print('inner_semimajor_axis < outer_semimajor_axis' )
        return False, 0,0
    if (inner_semimajor_axis < 0.|units.RSun):
        print('error: inner separation not in allowed range')
        return False, 0,0
    if (outer_semimajor_axis < 0.|units.RSun):
        print('error: outer separation not in allowed range')
        return False, 0,0

    if (inner_eccentricity < 0.) or (inner_eccentricity > 1.):
        print('error: inner eccentricity not in allowed range')
        return False, 0,0
    if (outer_eccentricity < 0.) or (outer_eccentricity > 1.):
        print('error: outer eccentricity not in allowed range')
        return False, 0,0
    if (inner_eccentricity < minimum_eccentricity):
        inner_eccentricity = minimum_eccentricity
    if (outer_eccentricity < minimum_eccentricity):
        outer_eccentricity = minimum_eccentricity

    if (relative_inclination < 0.) or (relative_inclination > np.pi):
        print('error: relative inclination not in allowed range')
        return False, 0,0

    if (inner_argument_of_pericenter < -1.*np.pi) or (inner_argument_of_pericenter > np.pi):
        print('error: inner argument of pericenter not in allowed range')
        return False, 0,0
    if (outer_argument_of_pericenter < -1.*np.pi) or (outer_argument_of_pericenter > np.pi):
        print('error: outer argument of pericenter not in allowed range')
        return False, 0,0

    if (inner_longitude_of_ascending_node < -1.*np.pi) or (inner_longitude_of_ascending_node > np.pi):
        print('error: inner longitude of ascending node not in allowed range')
        return False, 0,0
        
    return True, inner_eccentricity, outer_eccentricity 
    
def make_particle_sets(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node):
            
        correct_params, inner_eccentricity, outer_eccentricity = test_initial_parameters(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis, inner_eccentricity, outer_eccentricity,
            relative_inclination, inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node)  
            
        outer_longitude_of_ascending_node = inner_longitude_of_ascending_node - np.pi
                        
        stars = make_stars(inner_primary_mass, inner_secondary_mass, outer_mass)
        bins = make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)

        return stars, bins, correct_params
    
#-------
#setup community codes

def setup_stellar_code(stellar_code, stars):
    stellar_code.particles.add_particles(stars)
    if USE_MESA_AS_STELLAR_CODE:
        setup_mesa_stars(stellar_code)
    return stellar_code

                    
def setup_secular_code(triple, secular_code, stop_at_semisecular_regime):
    triple_set = triple.as_set()
    triple_time = triple_set.time
    secular_code.triples.add_particles(triple_set)
    secular_code.parameters.verbose = False
#    secular_code.parameters.verbose = True
    
    #needed for initialisation in some circumstances 
    secular_code.model_time = triple_time
    
    secular_code.parameters.equations_of_motion_specification = 0
    secular_code.parameters.roche_radius_specification = 0
    #0: eccentric eggleton, 1: sepinsky, 2: classical circular eggleton
    secular_code.parameters.stability_limit_specification = 0
    #for stars 0, 5-6, for exoplanets 1-4
    #0: mardling & aarseth 2001, 1:petrovich et al. 2015 simple, 2:petrovich et al. 2015 
    #3: holman et al. 98 s-type, 4: holman et al. 98 p-type,  
    #5: vynatheya+ 22
    #6: tory+ 22

    secular_code.parameters.ignore_tertiary = False

    secular_code.parameters.include_quadrupole_terms = True
    secular_code.parameters.include_octupole_terms = True        
    secular_code.parameters.include_inner_wind_terms = True
    secular_code.parameters.include_outer_wind_terms = True
    secular_code.parameters.include_inner_RLOF_terms = True
    secular_code.parameters.include_outer_RLOF_terms = True
    secular_code.parameters.include_magnetic_braking_terms = False # not tested

    secular_code.parameters.include_inner_tidal_terms = True
    secular_code.parameters.include_outer_tidal_terms = True
    
    secular_code.parameters.include_1PN_inner_terms = True
    secular_code.parameters.include_1PN_outer_terms = True
    secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
    secular_code.parameters.include_25PN_inner_terms = True
    secular_code.parameters.include_25PN_outer_terms = True

    secular_code.parameters.check_for_dynamical_stability = True
    secular_code.parameters.check_for_dynamical_stability_at_initialisation = True

    secular_code.parameters.check_for_semisecular_regime = stop_at_semisecular_regime
    secular_code.parameters.check_for_semisecular_regime_at_initialisation = stop_at_semisecular_regime
    
    secular_code.parameters.check_for_inner_collision = True
    secular_code.parameters.check_for_outer_collision = True

    secular_code.parameters.check_for_inner_RLOF = True 
    secular_code.parameters.check_for_outer_RLOF = True 
    
    secular_code.parameters.include_spin_radius_mass_coupling_terms_star1 = True
    secular_code.parameters.include_spin_radius_mass_coupling_terms_star2 = True
    secular_code.parameters.include_spin_radius_mass_coupling_terms_star3 = True
    
     # accuracy of secular code
#        secular_code.parameters.input_precision = 1.0e-10#1.0e-5
#        secular_code.parameters.relative_tolerance = 1.0e-10
#        secular_code.parameters.threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero = 1.0e-12
    secular_code.parameters.threshold_value_of_spin_angular_frequency_for_setting_spin_angular_frequency_dot_moment_of_inertia_plus_wind_changes_zero = 1.0e-7|units.Myr**-1

    secular_code.parameters.include_linear_mass_change = True #needed for Jspin conservation
    secular_code.parameters.include_linear_radius_change = True #needed for Jspin conservation

#    channel_from_secular = secular_code.triples.new_channel_to(triple_set)
#    channel_to_secular = triple_set.new_channel_to(secular_code.triples)

    return secular_code 
    
#-------
