# to do 
# min teken in mean anomaly

## Triple:      Triple evolution
##              computes the evolution of a given triple
##              given any initial conditions (M, m, l, A, a, E, e, i, G, g, O, o, T, z).
 
import sys

import numpy as np
# from interactions import *
# from tidal_friction_constant import *
from amuse.units import units
from amuse.support.console import set_printing_strategy

from amuse.community.seba.interface import SeBa
from seculartriple_TPS.interface import SecularTriple


from triple_class import Triple_Class
from TRES_plotting import plot_data_container, plot_function
from TRES_setup import make_particle_sets, setup_stellar_code
from TRES_options import REPORT_DEBUG, \
                         REPORT_TRIPLE_EVOLUTION, \
                         MAKE_PLOTS, \
                         REPORT_USER_WARNINGS, \
                         USE_MESA_AS_STELLAR_CODE
from interactions import corotating_spin_angular_frequency_binary, \
                        lang_spin_angular_frequency, \
                        break_up_angular_frequency, \
                        criticial_angular_frequency_CHE

if USE_MESA_AS_STELLAR_CODE:
    from amuse.community.mesa.interface import Mesa

def initialize_triple_class(stars, bins, correct_params,
                            stellar_code, secular_code, relative_inclination = 80.0*np.pi/180.0,
                            metallicity = 0.02, tend = 5.0 |units.Myr, tinit = 0.0|units.Myr, 
                            number = 0, maximum_radius_change_factor = 0.005,
                            stop_at_mass_transfer = True, stop_at_init_mass_transfer = True, stop_at_outer_mass_transfer = True,
                            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
                            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False, which_common_envelope = 2,
                            stop_at_no_CHE = False, include_CHE = False, 
                            stop_at_merger = True, stop_at_disintegrated = True, stop_at_inner_collision = True, stop_at_outer_collision = True, 
                            stop_at_dynamical_instability = True, stop_at_semisecular_regime = False, 
                            stop_at_SN = False, SN_kick_distr = 2, impulse_kick_for_black_holes = True, fallback_kick_for_black_holes = True, 
                            stop_at_CPU_time = False, max_CPU_time = 3600.0, file_name = "TRES.hdf", file_type = "hdf5", dir_plots = ""):

    triple = Triple_Class(stars, bins, correct_params, stellar_code, secular_code,
                          relative_inclination, tend, tinit,
                          number, maximum_radius_change_factor,  
                          stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
                          stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
                          stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
                          stop_at_no_CHE, include_CHE, 
                          stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
                          stop_at_dynamical_instability, stop_at_semisecular_regime, 
                          stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes, fallback_kick_for_black_holes,
                          stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots)
    triple.stellar_code.parameters.metallicity = metallicity

    return triple


#-----
#for running TRES.py from other routines
def main(inner_primary_mass = 1.3|units.MSun, inner_secondary_mass = 0.5|units.MSun, outer_mass = 0.5|units.MSun,
            inner_semimajor_axis = 1.0 |units.AU, outer_semimajor_axis = 100.0 |units.AU,
            inner_eccentricity = 0.1, outer_eccentricity= 0.5,
            relative_inclination = 80.0*np.pi/180.0,
            inner_argument_of_pericenter = 0.1, outer_argument_of_pericenter = 0.5,
            inner_longitude_of_ascending_node = 0.0,
            metallicity = 0.02, tend = 5.0 |units.Myr, tinit = 0.0|units.Myr, 
            number = 0, maximum_radius_change_factor = 0.005,
            stop_at_mass_transfer = True, stop_at_init_mass_transfer = True, stop_at_outer_mass_transfer = True,
            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False, which_common_envelope = 2,
            stop_at_no_CHE = False, include_CHE = False, 
            stop_at_merger = True, stop_at_disintegrated = True, stop_at_inner_collision = True, stop_at_outer_collision = True, 
            stop_at_dynamical_instability = True, stop_at_semisecular_regime = False, 
            stop_at_SN = False, SN_kick_distr = 2, impulse_kick_for_black_holes = True, fallback_kick_for_black_holes = True, 
            stop_at_CPU_time = False, max_CPU_time = 3600.0, file_name = "TRES.hdf", file_type = "hdf5", dir_plots = "",
            stellar_code=None, secular_code=None):

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    relative_inclination = float(relative_inclination)
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)

    stars, bins, correct_params = make_particle_sets(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node)

    clean_up_stellar_code = False
    clean_up_secular_code = False
    if stellar_code is None:
        if USE_MESA_AS_STELLAR_CODE:
            stellar_code = Mesa(redirection='none')
        else:
            stellar_code = SeBa()
    #    stellar_code = SeBa(redirection='none')
    #    stellar_code = SeBa(redirection='file', redirect_file='output_SeBa_TRES.txt')
        clean_up_stellar_code = True

    stellar_code.parameters.metallicity = metallicity
    if secular_code is None:
        secular_code = SecularTriple()
    #    secular_code = SecularTriple(redirection='none')
    #    secular_code = SecularTriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')
        clean_up_secular_code = True

    triple_class_object = Triple_Class(stars, bins, correct_params, stellar_code, secular_code,
            relative_inclination, tend, tinit,
            number, maximum_radius_change_factor,  
            stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
            stop_at_no_CHE, include_CHE, 
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime, 
            stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes, fallback_kick_for_black_holes,
            stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots)


    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given triple are incorrect.')
        return triple_class_object # no codes initialized yet
    elif stop_at_semisecular_regime == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is in the semisecular regime at initialization.')
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is dynamically unstable at initialization.')
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given triple at initialization.')
    elif stop_at_no_CHE == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()
            

    stellar_code.particles.remove_particles(stars)
    triple_set = triple_class_object.triple.as_set()    
    secular_code.triples.remove_particles(triple_set)            
    del stars, bins, triple_set

    if clean_up_stellar_code:
        triple_class_object.stellar_code.stop()
        print('cleaning se')
    if clean_up_secular_code:
        triple_class_object.secular_code.stop()
        print('cleaning sec')
    
    return triple_class_object


def main_developer(stars, bins, correct_params, stellar_code, secular_code,
            relative_inclination = 80.0*np.pi/180.0, 
            metallicity = 0.02, tend = 5.0 |units.Myr, tinit = 0.0 |units.Myr, 
            number = 0, maximum_radius_change_factor = 0.005,
            stop_at_mass_transfer = True, stop_at_init_mass_transfer = True, stop_at_outer_mass_transfer = True,
            stop_at_stable_mass_transfer = True, stop_at_eccentric_stable_mass_transfer = True,
            stop_at_unstable_mass_transfer = False, stop_at_eccentric_unstable_mass_transfer = False, which_common_envelope = 2, 
            stop_at_no_CHE = False, include_CHE = False, 
            stop_at_merger = True, stop_at_disintegrated = True, stop_at_inner_collision = True, stop_at_outer_collision = True, 
            stop_at_dynamical_instability = True, stop_at_semisecular_regime = False, 
            stop_at_SN = False, SN_kick_distr = 2, impulse_kick_for_black_holes = True, fallback_kick_for_black_holes = True,
            stop_at_CPU_time = False, max_CPU_time = 3600.0, file_name = "TRES.hdf", file_type = "hdf5", dir_plots = ""):


    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    bins.eccentricity[0] = float(bins.eccentricity[0])
    bins.eccentricity[1] = float(bins.eccentricity[1])
    bins.argument_of_pericenter[0] = float(bins.argument_of_pericenter[0])
    bins.argument_of_pericenter[1] = float(bins.argument_of_pericenter[1])
    bins.longitude_of_ascending_node[0] = float(bins.longitude_of_ascending_node[0])
    bins.longitude_of_ascending_node[1] = float(bins.longitude_of_ascending_node[1])    
    relative_inclination = float(relative_inclination)

    triple_class_object = Triple_Class(stars, bins, correct_params, stellar_code, secular_code,
                                       relative_inclination, tend, tinit,
                                       number, maximum_radius_change_factor,  
                                       stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
                                       stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
                                       stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
                                       stop_at_no_CHE, include_CHE,
                                       stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
                                       stop_at_dynamical_instability, stop_at_semisecular_regime, 
                                       stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes, fallback_kick_for_black_holes,
                                       stop_at_CPU_time, max_CPU_time, file_name, file_type, dir_plots)


    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given triple are incorrect.')
        return triple_class_object # no codes initialized yet
    elif stop_at_semisecular_regime == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is in the semisecular regime at initialization.')
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is dynamically unstable at initialization.')
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given triple at initialization.')
    elif stop_at_no_CHE == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        triple_class_object.evolve_model(tend)
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, dir_plots)
            triple_class_object.print_stellar_system()
   
    return triple_class_object



#-----

#-----
#for running triple.py from the commandline
def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-M", "--M1", unit=units.MSun, 
                      dest="inner_primary_mass", type="float", default = 1.3|units.MSun,
                      help="inner primary mass [%default]")
    parser.add_option("-m", "--M2", unit=units.MSun, 
                      dest="inner_secondary_mass", type="float", default = 0.5|units.MSun,
                      help="inner secondary mass [%default]")
    parser.add_option("-l", "--M3", unit=units.MSun, 
                      dest="outer_mass", type="float", default = 0.5|units.MSun,
                      help="outer mass [%default]")

    parser.add_option("-A", "--Ain",  unit=units.RSun,
                      dest="inner_semimajor_axis", type="float", 
                      default = 200.0 |units.RSun,
                      help="inner semi major axis [%default]")
    parser.add_option("-a", "--Aout",unit=units.RSun,
                      dest="outer_semimajor_axis", type="float", 
                      default = 20000.0 |units.RSun,
                      help="outer semi major axis [%default]")
    parser.add_option("-E", "--Ein",
                      dest="inner_eccentricity", type="float", default = 0.1,
                      help="inner eccentricity [%default]")
    parser.add_option("-e", "--Eout",
                      dest="outer_eccentricity", type="float", default = 0.5,
                      help="outer eccentricity [%default]")
    parser.add_option("-i","-I",
                      dest="relative_inclination", type="float", default = 80.0*np.pi/180.0,
                      help="relative inclination [rad] [%default]")
    parser.add_option("-G", "--Gin",
                      dest="inner_argument_of_pericenter", type="float", default = 0.1,
                      help="inner argument of pericenter [rad] [%default]")
    parser.add_option("-g","--Gout",
                      dest="outer_argument_of_pericenter", type="float", default = 0.5,
                      help="outer argument of pericenter [rad] [%default]")
    parser.add_option("-O", "--Oin",
                      dest="inner_longitude_of_ascending_node", type="float", default = 0.0,
                      help="inner longitude of ascending node [rad] [%default]")
##             outer longitude of ascending nodes = inner - pi               
#    parser.add_option("-o",
#                      dest="outer_longitude_of_ascending_node", type="float", default = 0.0,
#                      help="outer longitude of ascending node [rad] [%default]")

    parser.add_option("-z", "-Z", dest="metallicity", type="float", default = 0.02,
                      help="metallicity [%default] %unit")
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 5.0 |units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("--initial_time", unit=units.Myr, 
                      dest="tinit", type="float", default = 0.0 |units.Myr,
                      help="initial time [%default] %unit")
    parser.add_option("-N", dest="number", type="int", default = 0,
                      help="number ID of system [%default]")
    parser.add_option("-r", dest="maximum_radius_change_factor", type="float", default = 0.01,
                      help="maximum_radius_change_factor [%default] %unit")

#    parser.add_option("--tidal", dest="tidal_terms", action="store_false", default = True, 
#                      help="tidal terms included [%default] %unit")

    parser.add_option("--no_stop_at_mass_transfer", dest="stop_at_mass_transfer", action="store_false", default = True,
                      help="stop at mass transfer [%default] %unit")
    parser.add_option("--no_stop_at_init_mass_transfer", dest="stop_at_init_mass_transfer", action="store_false", default = True,
                      help="stop if initially mass transfer[%default] %unit")
    parser.add_option("--no_stop_at_outer_mass_transfer", dest="stop_at_outer_mass_transfer", action="store_false", default = True,
                      help="stop at triple mass transfer [%default] %unit")
                      
#   if stop_at_mass_transfer is False, the following 4 stopping conditions can be used to further specify.
#   if stop_at_mass_transfer is True, the following 4 are ignored.
    parser.add_option("--stop_at_stable_mass_transfer", dest="stop_at_stable_mass_transfer", action="store_true", default = False,
                      help="stop at stable mass transfer [%default] %unit")
    parser.add_option("--stop_at_eccentric_stable_mass_transfer", dest="stop_at_eccentric_stable_mass_transfer", action="store_true",                                               
                    default = False, help="stop at eccentric stable mass transfer [%default] %unit")
    #unstable mass transfer leads to common-envelope evolution
    parser.add_option("--stop_at_unstable_mass_transfer", dest="stop_at_unstable_mass_transfer", action="store_true", 
                    default = False, help="stop at unstable mass transfer [%default] %unit")
    parser.add_option("--stop_at_eccentric_unstable_mass_transfer", dest="stop_at_eccentric_unstable_mass_transfer", 
                    action="store_true", default = False, help="stop at eccentric unstable mass transfer [%default] %unit")
    #0  alpha-ce + alpha-dce
    #1  gamma-ce + alpha-dce
    #2  seba style; combination of gamma-ce, alpha-ce & alpha-dce
    parser.add_option("--CE", dest="which_common_envelope",  type="int", default = 2,
                      help="which common envelope modeling [%default]")                      

    parser.add_option("--stop_at_no_CHE", dest="stop_at_no_CHE", 
                    action="store_true", default = False, help="stop if no chemically homogeneous evolution [%default] %unit")
    parser.add_option("--include_CHE", dest="include_CHE", 
                    action="store_true", default = False, help="include chemically homogeneous evolution in the stellar evolution [%default] %unit")

    parser.add_option("--no_stop_at_merger", dest="stop_at_merger", action="store_false", default = True, 
                      help="stop at merger [%default] %unit")
    parser.add_option("--no_stop_at_disintegrated", dest="stop_at_disintegrated", action="store_false", default = True,
                      help="stop at disintegrated [%default] %unit")
    parser.add_option("--no_stop_at_inner_collision", dest="stop_at_inner_collision", action="store_false",default = True,
                      help="stop at collision in inner binary[%default] %unit")
    parser.add_option("--no_stop_at_outer_collision", dest="stop_at_outer_collision", action="store_false",default = True,
                      help="stop at collision in outer binary[%default] %unit")
    parser.add_option("--no_stop_at_dynamical_instability", dest="stop_at_dynamical_instability", action="store_false", default = True,
                      help="stop at dynamical instability [%default] %unit")
    parser.add_option("--stop_at_semisecular_regime", dest="stop_at_semisecular_regime", action="store_true", default = False,
                      help="stop at semisecular regime [%default] %unit")

    parser.add_option("--stop_at_SN", dest="stop_at_SN", action="store_true", default = False,
                      help="stop at supernova [%default] %unit")
    #0  No kick 
    #1  Hobbs, Lorimer, Lyne & Kramer 2005, 360, 974  
    #2  Arzoumanian ea 2002, 568, 289
    #3  Hansen & Phinney 1997, 291, 569
    #4  Paczynski 1990, 348, 485
    #5  Verbunt, Igoshev & Cator, 2017, 608, 57
    parser.add_option("--SN_kick_distr", dest="SN_kick_distr",  type="int", default = 5,
                      help="which supernova kick distribution [%default]")                      
    parser.add_option("--no_impulse_kick_for_black_holes", dest="impulse_kick_for_black_holes",  action="store_false", default = True,
                      help="do not rescale the BH SN kick by mass -> impulse kick [%default]")                      
    parser.add_option("--no_fallback_kick_for_black_holes", dest="fallback_kick_for_black_holes",  action="store_false", default = True,
                      help="do not rescale the BH SN kick with fallback  [%default]")                      
                      
                      
                      
    parser.add_option("--stop_at_CPU_time", dest="stop_at_CPU_time", action="store_true", default = False,
                      help="stop at CPU time [%default] %unit")
    parser.add_option("--max_CPU_time", dest="max_CPU_time", type="float", default = 3600.0,
                      help="max CPU time [%default] %unit")
                      
                      
    parser.add_option("-f", dest="file_name", type ="string", default = "TRES.hdf",#"TRES.txt"
                      help="file name[%default]")
    parser.add_option("-F", dest="file_type", type ="string", default = "hdf5",#"txt"
                      help="file type[%default]")
    parser.add_option("--dir_plots", dest="dir_plots", type ="string", default = "",#"txt"
                      help="directory for plots for debugging mode [%default]")
                                           
    options, args = parser.parse_args()
    return options.__dict__


if __name__ == '__main__':
    opt = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    stars, bins, correct_params = make_particle_sets(opt["inner_primary_mass"], opt["inner_secondary_mass"], opt["outer_mass"],
            opt["inner_semimajor_axis"], opt["outer_semimajor_axis"],
            opt["inner_eccentricity"], opt["outer_eccentricity"],
            opt["relative_inclination"],
            opt["inner_argument_of_pericenter"], opt["outer_argument_of_pericenter"],
            opt["inner_longitude_of_ascending_node"])

    if USE_MESA_AS_STELLAR_CODE:
        stellar_code = Mesa(redirection='none')
    else:
        stellar_code = SeBa()
#    stellar_code = SeBa(redirection='none')
#    stellar_code = SeBa(redirection='file', redirect_file='output_SeBa_TRES.txt')
    stellar_code.parameters.metallicity = opt["metallicity"]
    secular_code = SecularTriple()
#    secular_code = SecularTriple(redirection='none')
#    secular_code = SecularTriple(redirection='file', redirect_file='output_SecularTriple_TRES.txt')

    triple_class_object = Triple_Class(stars, bins, correct_params, stellar_code, secular_code,
            opt["relative_inclination"], opt["tend"], opt["tinit"],
            opt["number"], opt["maximum_radius_change_factor"],  
            opt["stop_at_mass_transfer"], opt["stop_at_init_mass_transfer"], opt["stop_at_outer_mass_transfer"],
            opt["stop_at_stable_mass_transfer"], opt["stop_at_eccentric_stable_mass_transfer"],
            opt["stop_at_unstable_mass_transfer"], opt["stop_at_eccentric_unstable_mass_transfer"], opt["which_common_envelope"],
            opt["stop_at_no_CHE"], opt["include_CHE"],
            opt["stop_at_merger"], opt["stop_at_disintegrated"], opt["stop_at_inner_collision"], opt["stop_at_outer_collision"], 
            opt["stop_at_dynamical_instability"], opt["stop_at_semisecular_regime"], 
            opt["stop_at_SN"], opt["SN_kick_distr"], opt["impulse_kick_for_black_holes"], opt["fallback_kick_for_black_holes"],
            opt["stop_at_CPU_time"], opt["max_CPU_time"], opt["file_name"], opt["file_type"], opt["dir_plots"])  

    if triple_class_object.correct_params == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The parameters of the given triple are incorrect.' )   
        # no codes initialized yet
        sys.exit('Choose a different system. The parameters of the given triple are incorrect.')
    elif opt['stop_at_semisecular_regime'] == True and triple_class_object.semisecular_regime_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is in the semisecular regime at initialization.')
    elif triple_class_object.dynamical_instability_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. The given triple is dynamically unstable at initialization.')
    elif triple_class_object.mass_transfer_at_initialisation == True:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. There is mass transfer in the given triple at initialization.')
    elif opt["stop_at_no_CHE"] == True and triple_class_object.CHE_at_initialisation == False:
        if REPORT_USER_WARNINGS:
            print('Choose a different system. No chemically homogeneous evolution at initialization')
    else:    
        triple_class_object.evolve_model(opt["tend"])
        if REPORT_DEBUG or MAKE_PLOTS:
            plot_function(triple_class_object, opt['dir_plots'])
            triple_class_object.print_stellar_system()



        if REPORT_TRIPLE_EVOLUTION:
            print('Simulation has finished succesfully')
            
    print('\nYou have used the TRES triple evolution code. Literature reference:')
    print('** Toonen, Hamers & Portegies Zwart 2016, ComAC, 3, 6T:')
    print('... "The evolution of hierarchical triple star-systems" ')
            
    triple_class_object.stellar_code.stop()
    triple_class_object.secular_code.stop()

