    
from interactions import *
from tidal_friction_constant import *

import sys
import time
from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.io import write_set_to_file
from amuse.units import quantities
from scipy.stats import maxwell
from scipy import optimize
# from math import isnan
import numpy as np

# from tres_options import REPORT_USER_WARNINGS, \
#                          GET_GYRATION_RADIUS_FROM_STELLAR_CODE, \
#                          GET_AMC_FROM_STELLAR_CODE, \
#                          REPORT_TRIPLE_EVOLUTION, \
#                          REPORT_DEBUG, \
#                          MAKE_PLOTS, \
#                          REPORT_DT, \
#                          kozai_type_factor, \
#                          maximum_time_step_factor, \
#                          minimum_time_step

from TRES_options import *
from TRES_setup import setup_secular_code, setup_stellar_code
from TRES_plotting import plot_data_container


class Triple_Class:
    #-------
    #setup stellar system
    def __init__(self, stars, bins, correct_params,
            stellar_code, secular_code, relative_inclination, 
            tend, tinit, number, maximum_radius_change_factor, 
            stop_at_mass_transfer, stop_at_init_mass_transfer, stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, which_common_envelope,
            stop_at_no_CHE, include_CHE, 
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime,  
            stop_at_SN, SN_kick_distr, impulse_kick_for_black_holes,fallback_kick_for_black_holes,
            stop_at_CPU_time, max_CPU_time,
            file_name, file_type, dir_plots):
        
        self.correct_params = correct_params
        if correct_params == False:
            self.triple = Particles(1) 
            return
            
        self.tend = tend 
        self.tinit = tinit
        self.previous_time = 0.0|units.yr
        self.previous_dt = 0.0|units.yr
        self.instantaneous_evolution = False 
        if tinit > 0|units.Myr:
            self.instantaneous_evolution = True # no secular    
        
        self.maximum_radius_change_factor = maximum_radius_change_factor
        self.fixed_timestep = -1|units.Myr

        self.file_name = file_name
        self.file_type = file_type
        self.which_common_envelope = which_common_envelope
        self.include_CHE = include_CHE      
        if REPORT_USER_WARNINGS  and self.include_CHE:
            print("Note: For CHE evolution to be included, it also needs to be switched on manually in SeBa.") 
            print("\t This can be done by setting include_CHE=True in SeBa's sstar/starclass/constants.C.")
            print("\t AMUSE developer mode is required to access SeBa files.")
        self.SN_kick_distr = SN_kick_distr
        self.impulse_kick_for_black_holes = impulse_kick_for_black_holes
        self.fallback_kick_for_black_holes = fallback_kick_for_black_holes
        self.max_CPU_time = max_CPU_time

        self.triple = bins[1]        
        self.triple.time = 0.0|units.yr
        self.triple.relative_inclination = relative_inclination 
        self.triple.is_star = False #maybe not necessary?
        self.triple.dynamical_instability = False 
        self.triple.number = number 
        self.triple.error_flag_secular = 0
        self.triple.CPU_time = 0.0
        self.triple.child2.delta_e_in = 0.0
        self.triple.child2.max_delta_e_in = 0.0
                        
        self.dynamical_instability_at_initialisation = False
        self.semisecular_regime_at_initialisation = False
        self.mass_transfer_at_initialisation = False
        self.CHE_at_initialisation = False

        self.set_stopping_conditions(stop_at_mass_transfer, stop_at_init_mass_transfer,stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, stop_at_no_CHE,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime,  stop_at_SN, stop_at_CPU_time)

        self.initialize_stellar(stellar_code, stars)
        self.initialize_secular(secular_code, stop_at_semisecular_regime, 
                                              stop_at_dynamical_instability)

        self.check_RLOF() 
        if self.has_tertiary_donor() and (self.stop_at_outer_mass_transfer or self.stop_at_mass_transfer or self.stop_at_init_mass_transfer): 
            self.mass_transfer_at_initialisation = True
            self.triple.bin_type = bin_type['rlof']
            return
        
        self.check_OLOF()
        if self.stop_at_no_CHE and (not self.check_CHE()):
            self.CHE_at_initialization = False
            return
                    
        if (self.has_donor() or self.has_OLOF_donor()) and (self.stop_at_mass_transfer or self.stop_at_init_mass_transfer):            
            self.mass_transfer_at_initialisation = True
            #assuming object is triple as is triple constructor
            if self.is_binary(self.triple.child1):
                bin = self.triple.child1
            else:
                bin = self.triple.child2

            if bin.child1.is_OLOF_donor or bin.child2.is_OLOF_donor:    
                bin.bin_type = bin_type['olof'] 
                return             
            elif bin.child1.is_donor and bin.child2.is_donor:
                bin.bin_type = bin_type['contact'] 
                return
            else:
                bin.bin_type = bin_type['rlof']
                return

            return

        self.triple.kozai_type = self.get_kozai_type()
        self.update_stellar_parameters() 
        self.update_time_derivative_of_radius()
        self.update_previous_stellar_parameters()

    def initialize_stellar(self, stellar_code, stars):
        self.stellar_code = setup_stellar_code(stellar_code, stars)
        self.channel_from_stellar = stellar_code.particles.new_channel_to(stars)
        self.channel_to_stellar = stars.new_channel_to(stellar_code.particles)
        self.copy_from_stellar()        
        self.initial_angular_frequency() 

    def initialize_secular(self, secular_code, 
                            stop_at_semisecular_regime,
                            stop_at_dynamical_instability):

        triple_set = self.triple.as_set()
        self.secular_code = setup_secular_code(self.triple, secular_code, stop_at_semisecular_regime)      
        self.channel_from_secular = self.secular_code.triples.new_channel_to(triple_set)
        self.channel_to_secular = triple_set.new_channel_to(self.secular_code.triples)
        self.channel_to_secular.copy()

        self.secular_code.check_for_dynamical_stability()
        if stop_at_dynamical_instability == True and self.secular_code.triples[0].dynamical_instability == True:
            self.dynamical_instability_at_initialisation = True
            self.triple.dynamical_instability = True
            self.set_bintype_to_dynamical_instability()
            return 

        self.secular_code.check_for_semisecular_regime()
        if stop_at_semisecular_regime == True and self.secular_code.triples[0].semisecular_regime == True:
            self.semisecular_regime_at_initialisation = True
            self.triple.semisecular_regime = True        
            return
   
        
    def set_stopping_conditions(self, stop_at_mass_transfer,stop_at_init_mass_transfer,stop_at_outer_mass_transfer,
            stop_at_stable_mass_transfer, stop_at_eccentric_stable_mass_transfer,
            stop_at_unstable_mass_transfer, stop_at_eccentric_unstable_mass_transfer, stop_at_no_CHE,
            stop_at_merger, stop_at_disintegrated, stop_at_inner_collision, stop_at_outer_collision, 
            stop_at_dynamical_instability, stop_at_semisecular_regime, stop_at_SN, stop_at_CPU_time):

        if stop_at_disintegrated == False:
            sys.exit('stop_at_disintegrated = False not possible yet. After the disintegration of the triple, further evolution can be done with stellar code directly. ')
        if stop_at_outer_mass_transfer == False:
            sys.exit('stop_at_outer_mass_transfer = False not possible yet. Methodology is as of yet non-existent.' )
        if stop_at_outer_collision == False:
            sys.exit('stop_at_outer_collision = False not possible. Non-hierarchical triples can not be simulated using the secular equations as used in TRES. Further evolution should be done by other means, e.g. one of the N-body codes implemented in AMUSE.' )
        if stop_at_dynamical_instability == False:
            sys.exit('stop_at_dynamical_instability = False not possible. Unstable triples can not be simulated using the secular equations as used in TRES. Further evolution should be done by other means, e.g. one of the N-body codes implemented in AMUSE.') 

                            
        self.stop_at_mass_transfer = stop_at_mass_transfer            
        self.stop_at_init_mass_transfer = stop_at_init_mass_transfer
        self.stop_at_outer_mass_transfer = stop_at_outer_mass_transfer            

        self.stop_at_stable_mass_transfer =  stop_at_stable_mass_transfer
        self.stop_at_eccentric_stable_mass_transfer = stop_at_eccentric_stable_mass_transfer
        self.stop_at_unstable_mass_transfer = stop_at_unstable_mass_transfer
        self.stop_at_eccentric_unstable_mass_transfer = stop_at_eccentric_unstable_mass_transfer
        self.stop_at_no_CHE = stop_at_no_CHE
        
        self.stop_at_merger = stop_at_merger            
        self.stop_at_disintegrated = stop_at_disintegrated            
        self.stop_at_inner_collision = stop_at_inner_collision            
        self.stop_at_outer_collision = stop_at_outer_collision            

        self.stop_at_dynamical_instability = stop_at_dynamical_instability            
        self.stop_at_semisecular_regime = stop_at_semisecular_regime
        self.stop_at_SN = stop_at_SN
        self.stop_at_CPU_time = stop_at_CPU_time
    
            
    #-------

    def copy_from_stellar(self):
#        self.channel_from_stellar.copy()        
        self.channel_from_stellar.copy_attributes(["age", "mass", "core_mass", "radius", "core_radius", "convective_envelope_radius",  "convective_envelope_mass",  "stellar_type", "luminosity", "wind_mass_loss_rate",  "temperature"]) 
        if GET_GYRATION_RADIUS_FROM_STELLAR_CODE:
            self.channel_from_stellar.copy_attributes(["gyration_radius"]) 
        if GET_AMC_FROM_STELLAR_CODE:
            self.channel_from_stellar.copy_attributes(["apsidal_motion_constant"]) 
        

    #-------        
        
    def initial_angular_frequency(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.previous_time = self.triple.time
        if stellar_system.is_star:
            if stellar_system.stellar_type in stellar_types_planetary_objects:  
                stellar_system.spin_angular_frequency = 0.125 * break_up_angular_frequency(stellar_system)  
            else:
                stellar_system.spin_angular_frequency = lang_spin_angular_frequency(stellar_system)
                if self.include_CHE: #sets initial spin to corotation
                    stellar_system.spin_angular_frequency = corotating_spin_angular_frequency_binary(stellar_system.parent.semimajor_axis, self.get_mass(stellar_system.parent.child1), self.get_mass(stellar_system.parent.child2))
                    stellar_system.rotation_period = (2*np.pi/stellar_system.spin_angular_frequency)

        else:
            self.initial_angular_frequency(stellar_system.child1)        
            self.initial_angular_frequency(stellar_system.child2)

        if self.include_CHE: #only needed when including CHE
            self.channel_to_stellar.copy_attributes(['rotation_period']) 
                               
    #-------

    #-------
    def refresh_memory(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.triple.time = self.previous_time
        if stellar_system.is_star:
            stellar_system.mass = stellar_system.previous_mass 
            stellar_system.radius = stellar_system.previous_radius           
            stellar_system.stellar_type = stellar_system.previous_stellar_type
            stellar_system.moment_of_inertia_of_star = stellar_system.previous_moment_of_inertia_of_star
            stellar_system.time_derivative_of_radius = stellar_system.previous_time_derivative_of_radius
            stellar_system.spin_angular_frequency = stellar_system.previous_spin_angular_frequency
            
        else:
            self.refresh_memory(stellar_system.child1)        
            self.refresh_memory(stellar_system.child2)
            stellar_system.mass = stellar_system.previous_mass
            stellar_system.bin_type = stellar_system.previous_bin_type
            
            if self.is_triple(stellar_system):
                stellar_system.kozai_type = stellar_system.previous_kozai_type

    def update_previous_stellar_parameters(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.previous_time = self.triple.time
        if stellar_system.is_star:
            stellar_system.previous_mass = self.get_mass(stellar_system)      
            stellar_system.previous_radius = stellar_system.radius
            stellar_system.previous_stellar_type = stellar_system.stellar_type
            stellar_system.previous_moment_of_inertia_of_star = stellar_system.moment_of_inertia_of_star
            stellar_system.previous_spin_angular_frequency = stellar_system.spin_angular_frequency 
            if self.triple.time == quantities.zero: #initialization
               stellar_system.previous_time_derivative_of_radius = 0.0 | units.RSun/units.yr
            else:
               stellar_system.previous_time_derivative_of_radius = stellar_system.time_derivative_of_radius
        else:
            self.update_previous_stellar_parameters(stellar_system.child1)        
            self.update_previous_stellar_parameters(stellar_system.child2)
            stellar_system.previous_mass = self.get_mass(stellar_system) 
            stellar_system.previous_bin_type = stellar_system.bin_type
            
            if self.is_triple(stellar_system):
                stellar_system.previous_kozai_type = stellar_system.kozai_type
                
    #-------

    #-------
    def update_time_derivative_of_radius(self, stellar_system = None):
        #update time_derivative_of_radius for effect of wind on spin
        #radius change due to stellar evolution, not mass transfer
        if stellar_system == None:
            stellar_system = self.triple
                
        time_step = self.triple.time - self.previous_time

        if self.triple.time == quantities.zero:
            #initialization
            self.triple.child2.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child2.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
        else:     
            if stellar_system.is_star:
                stellar_system.time_derivative_of_radius = (stellar_system.radius - stellar_system.previous_radius)/time_step
            else:
                self.update_time_derivative_of_radius(stellar_system.child1)        
                self.update_time_derivative_of_radius(stellar_system.child2)
    #-------
    
    #-------
    def update_time_derivative_of_moment_of_inertia(self, stellar_system = None):
        #update time_derivative_of_radius for effect of changing Jspin
        if stellar_system == None:
            stellar_system = self.triple
                
        time_step = self.triple.time - self.previous_time

        if self.triple.time == quantities.zero:
            #initialization
            self.triple.child2.child1.time_derivative_of_moment_of_inertia = 0.0 | units.MSun*units.RSun**2/units.yr
            self.triple.child2.child2.time_derivative_of_moment_of_inertia = 0.0 | units.MSun*units.RSun**2/units.yr
            self.triple.child1.time_derivative_of_moment_of_inertia = 0.0 | units.RSun/units.yr
        else:     
            if stellar_system.is_star:
                stellar_system.time_derivative_of_moment_of_inertia = (stellar_system.moment_of_inertia_of_star - stellar_system.previous_moment_of_inertia_of_star)/time_step
            else:
                self.update_time_derivative_of_moment_of_inertia(stellar_system.child1)        
                self.update_time_derivative_of_moment_of_inertia(stellar_system.child2)
    #-------

    #-------
    def update_stellar_parameters(self, stellar_system = None):
        # for the convective envelope mass:
        # the prescription of Hurley, Pols & Tout 2000 is implemented in SeBa, however note that the prescription in BSE is different
        # for the convective envelope radius:
        # the prescription of Hurley, Tout & Pols 2002 is implemented in SeBa, however note that the prescription in BSE is different
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            
            if not GET_GYRATION_RADIUS_FROM_STELLAR_CODE:
                stellar_system.gyration_radius = set_gyration_radius(stellar_system.stellar_type, stellar_system.mass)
            if not GET_AMC_FROM_STELLAR_CODE:
                stellar_system.apsidal_motion_constant = self.apsidal_motion_constant(stellar_system) 
            if stellar_system.core_radius > stellar_system.radius:
                #can happen very late on the agb before WD formation
                stellar_system.core_radius = stellar_system.radius                
            stellar_system.moment_of_inertia_of_star = self.moment_of_inertia(stellar_system)

            if stellar_system.convective_envelope_radius < 0|units.RSun:
                sys.exit('convective_envelope_radius < 0')
            if stellar_system.convective_envelope_radius == 0|units.RSun:
                stellar_system.convective_envelope_mass = 1.e-10 |units.MSun    
                stellar_system.convective_envelope_radius = 1.e-10 |units.RSun  
                                 
            #When the GW inspiral time is shorter than the inner orbit, the numerical solver crashes
            #This is only possible for BH & NS as other stars would fill their RL earlier
            #To avoid this articially increase stellar radii of BH/NS in secular code 
            #Does not affect any other processes
            if stellar_system.stellar_type in stellar_types_SN_remnants:
                stellar_system.radius = stellar_system.radius*10                
                 
        else:
            self.update_stellar_parameters(stellar_system.child1)        
            self.update_stellar_parameters(stellar_system.child2)
            if self.is_triple(stellar_system):
                stellar_system.kozai_type = self.get_kozai_type()           
    #-------

    #-------
    # useful functions general
    
    #whether or not a stellar system consists of just two stars
    def is_binary(self, stellar_system=None):
        if stellar_system == None:
            stellar_system = self.triple    
    
        if not stellar_system.is_star and stellar_system.child1.is_star and stellar_system.child2.is_star:
            return True
        else:
            return False

    def is_triple(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if not stellar_system.is_star:
            if stellar_system.child1.is_star and self.is_binary(stellar_system.child2):
                return True
            elif stellar_system.child2.is_star and self.is_binary(stellar_system.child1):
                return True

        return False

    def has_donor(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            if stellar_system.is_donor:
                return True
        else:
            if self.has_donor(stellar_system.child1) or self.has_donor(stellar_system.child2):
                return True                        
            
        return False            


    def has_OLOF_donor(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
 
        if stellar_system.is_star:
            if stellar_system.is_OLOF_donor:
                return True
        else:
            if self.has_OLOF_donor(stellar_system.child1) or self.has_OLOF_donor(stellar_system.child2):
                return True                                               

        return False   

    def has_contact_system(self, stellar_system = None):
       if stellar_system == None:
           stellar_system = self.triple

       if stellar_system.is_star:
           return False
       elif self.is_binary(stellar_system):
           if stellar_system.child1.is_donor and stellar_system.child2.is_donor:
               return True
       else:
           if self.has_contact_system(stellar_system.child1):
               return True
           if self.has_contact_system(stellar_system.child2):
               return True

       return False

# if a merger is currently taking place, not if a merger has happened in the past
    def has_merger(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_merger(stellar_system.child1):
                return True
            if self.has_merger(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['merger']:  
                return True    

        return False            

# if a disruption is currently taking place, not if a disruption has happened in the past
    def has_disintegrated(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_disintegrated(stellar_system.child1):
                return True
            if self.has_disintegrated(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['disintegrated']:  
                return True    
            
        return False            



# if a dynamical instability is currently taking place, not if an instability has happened in the past
    def has_dynamical_instability(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_dynamical_instability(stellar_system.child1):
                return True
            if self.has_dynamical_instability(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['dyn_inst']:  
                return True    
            
        return False            

    def set_bintype_to_dynamical_instability(self, stellar_system = None): 
        if self.triple.dynamical_instability: 
            if stellar_system == None:
                stellar_system = self.triple
                
            if stellar_system.is_star == False:
                stellar_system.bin_type = bin_type['dyn_inst']       
                self.set_bintype_to_dynamical_instability(stellar_system.child1)
                self.set_bintype_to_dynamical_instability(stellar_system.child2)           


#doesn't work well, as it uses bin_types that are set later -> use has_tertiary_donor
# if a mass transfer in the outer binary of the triple is currently taking place, not if a mass transfer has happened in the past
#    def has_outer_mass_transfer(self, stellar_system = None): 
#        if stellar_system == None:
#            stellar_system = self.triple
#            
#        if stellar_system.is_star:
#            return False
#        elif self.is_binary(stellar_system):
#            return False
#        else:
#            if self.has_outer_mass_transfer(stellar_system.child1):
#                return True
#            if self.has_outer_mass_transfer(stellar_system.child2):
#                return True
#            if stellar_system.bin_type != bin_type['unknown'] and stellar_system.bin_type != bin_type['detached']:  
#                return True    
#            
#        return False            


# if a mass transfer in the outer binary of the triple is currently taking place, not if a mass transfer has happened in the past
    def has_tertiary_donor(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        elif self.is_binary(stellar_system):
            return False
        else:
            if self.has_tertiary_donor(stellar_system.child1): 
                return True
            if self.has_tertiary_donor(stellar_system.child2): 
                return True
            if stellar_system.child1.is_star and stellar_system.child1.is_donor and not stellar_system.child2.is_star:
                return True 
            if stellar_system.child2.is_star and stellar_system.child2.is_donor and not stellar_system.child1.is_star:
                return True             
        return False            
       
    def contains_SN_remnant(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type in stellar_types_SN_remnants:
                return True
        else:
            if self.contains_SN_remnant(stellar_system.child1) or self.contains_SN_remnant(stellar_system.child2):
                return True                        
            
        return False            
    

    def has_stellar_type_changed_into_SN_remnant(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                return True
        else:
            if self.has_stellar_type_changed_into_SN_remnant(stellar_system.child1) or self.has_stellar_type_changed_into_SN_remnant(stellar_system.child2):
                return True                        
            
        return False            


    def has_stellar_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type:
                return True
        else:
            if self.has_stellar_type_changed(stellar_system.child1) or self.has_stellar_type_changed(stellar_system.child2):
                return True                        
            
        return False            
    
    

    def has_kozai_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if self.is_triple(stellar_system):
            if stellar_system.kozai_type != stellar_system.previous_kozai_type:
                return True
        else: #binaries and single stars do not have a kozai timescale
            return False            

#obsolete?
#    def is_system_stable(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#            
#        if stellar_system.is_star:
#            return True
#        elif self.is_binary(stellar_system):
#            return stellar_system.is_mt_stable
#        else:
#            if stellar_system.is_mt_stable and self.is_system_stable(stellar_system.child1) and self.is_system_stable(stellar_system.child2):
#                return True                        
#            
#        return False                    

    def get_mass(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_system.mass
        else:
            M1 = self.get_mass(stellar_system.child1)        
            M2 = self.get_mass(stellar_system.child2)
            return M1 + M2 

    def get_size(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_system.radius
        else:
            return stellar_system.semimajor_axis
 

    #-------

    #-------
    # useful functions general
            
    def orbital_period(self, bs):
        if not bs.is_star:
            Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / self.get_mass(bs))
            return Porb
        else:
            sys.exit('orbital_period: single star does not have a period')

    def orbital_angular_momentum(self, bs):
        if not bs.is_star:
            M = self.get_mass(bs.child1)
            m = self.get_mass(bs.child2)
            a = bs.semimajor_axis
            e = bs.eccentricity
            J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))
        
            if REPORT_BINARY_EVOLUTION:
                print('Jorb:', M, m, a, e, J)
                
            return J
        else:
            sys.exit('orbital_angular_momentum: single star does not have an orbit')
    
    def spin_angular_momentum(self, ss):
        if ss.is_star:
            return ss.moment_of_inertia_of_star * ss.spin_angular_frequency
        else:
            sys.exit('spin_angular_momentum: structure stellar system unknown')        

            
            
    def apsidal_motion_constant(self, star):
    
        if star.stellar_type in [13]|units.stellar_type: #ns
            #based on Brooke & Olle 1955, for n=1 polytrope
            return 0.260
    
        elif star.stellar_type in [13, 14]|units.stellar_type: #bh
            # Hamers et al. 2013
            return 0.
        elif star.stellar_type in [1,7,10,11,12]|units.stellar_type:#ms, he-ms, wd
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144            

        elif star.stellar_type in [0,2,3,4,5,6,8,9,17]|units.stellar_type:#low-mass ms, hg, gb, cheb, agb, he-g, pre-ms
            #based on Brooke & Olle 1955, for n=3 polytrope
#            return 0.143 
            #based on Claret & Gimenez 1992, 96, 225 the value should be smaller, try:
            return 0.05
        elif star.stellar_type in [18]|units.stellar_type:#planet
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144            
        elif star.stellar_type in [19]|units.stellar_type:#bd
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144
        else:
            print(star.stellar_type)
            sys.exit('apsidal motion constant: stellar_type unknown')
            

    def moment_of_inertia(self, star):
        if star.is_star:

            if GET_GYRATION_RADIUS_FROM_STELLAR_CODE:
                I = star.gyration_radius**2 * (star.mass)*star.radius**2                     
            else: 
                k2 = 0.1
                k3 = 0.21
                
                if star.stellar_type in stellar_types_remnants:
                    I = k3*(star.mass)*star.radius**2 
                else:            
                    I = k2*(star.mass - star.core_mass)*star.radius**2 + k3*star.core_mass*star.core_radius**2
                        
            return I                   
        else:
            sys.exit('moment_of_inertia: structure stellar system unknown')        





    def octupole_parameter(self):
        if self.is_triple():
           if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
           else: 
                star = self.triple.child2
                bin = self.triple.child1

           return (self.get_mass(bin.child1)-self.get_mass(bin.child2))/self.get_mass(bin) * bin.semimajor_axis/self.triple.semimajor_axis * self.triple.eccentricity/(1-self.triple.eccentricity**2)           

        else:
            print('Octupole parameter needs triple system')
            return np.nan   




    def kozai_timescale(self):
        if self.is_triple():
           alpha_kozai = 1.
           if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
           else: 
                star = self.triple.child2
                bin = self.triple.child1

            
           P_in = self.orbital_period(bin) #period inner binary 
           P_out = self.orbital_period(self.triple)#period outer binary 
           return alpha_kozai * P_out**2 / P_in * (self.get_mass(self.triple) / self.get_mass(star)) * (1-self.triple.eccentricity**2)**1.5       

        else:
            print('Kozai timescale needs triple system')
            return np.nan   
            
    def get_kozai_type(self):
        if self.is_triple():
            if self.secular_code.parameters.ignore_tertiary == True:
                return False

            t_kozai = self.kozai_timescale()
            if t_kozai >kozai_type_factor * self.tend:
                return False
            
            t_ev = self.get_min_stellar_evolution_timescale_of_system()
            if t_kozai < kozai_type_factor * t_ev:
                return True
            else:
                return False
        else:
           sys.exit('Kozai type needs triple system')
  

    def get_min_stellar_evolution_timescale_of_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_evolution_timescale(stellar_system)
        else:
            t1 = self.get_min_stellar_evolution_timescale_of_system(stellar_system.child1)        
            t2 = self.get_min_stellar_evolution_timescale_of_system(stellar_system.child2)
            return max(t1, t2)


    def check_RLOF(self):
        if self.triple.is_star:
            return
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2)
                print('Stellar radii:', self.triple.child1.radius, self.triple.child2.radius)
            
            self.triple.child1.is_donor = False
            self.triple.child2.is_donor = False
            
            if self.triple.child1.radius >= Rl1 - (1.0 * small_numerical_error|units.RSun):
                self.triple.child1.is_donor = True               
            if self.triple.child2.radius >= Rl2 - (1.0 * small_numerical_error|units.RSun):
                self.triple.child2.is_donor = True                             
 
        elif self.is_triple() and self.secular_code.parameters.ignore_tertiary == True:
            # for disrupted binary
            if self.triple.child1.is_star:
                bin = self.triple.child2
            else:
                bin = self.triple.child1
        
            Rl1 = roche_radius(bin, bin.child1, self)
            Rl2 = roche_radius(bin, bin.child2, self)
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2)
                print('Stellar radii:', bin.child1.radius, bin.child2.radius)
            
            bin.child1.is_donor = False
            bin.child2.is_donor = False
            
            if bin.child1.radius >= Rl1 - (1.0 * small_numerical_error|units.RSun):
                bin.child1.is_donor = True               
            if bin.child2.radius >= Rl2 - (1.0 * small_numerical_error|units.RSun):
                bin.child2.is_donor = True                             
 
 
        elif self.is_triple():
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1

            #assumping secular code always returns inner binary first
            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
    
            if REPORT_TRIPLE_EVOLUTION:
                print('Roche lobe radii:', Rl1, Rl2, Rl3)
                print('Stellar radii:', bin.child1.radius, bin.child2.radius, star.radius)
                print('binary Roche lobe radii:', roche_radius(bin, bin.child1, self), roche_radius(bin, bin.child2, self), roche_radius(self.triple, star, self))
                print('eccentric binary Roche lobe radii:', roche_radius(bin, bin.child1, self)* (1-bin.eccentricity), roche_radius(bin, bin.child2, self)* (1-bin.eccentricity), roche_radius(self.triple, star, self)*(1-self.triple.eccentricity))
                print('Masses:', bin.child1.mass, bin.child2.mass, star.mass)
                print('Semi:', bin.semimajor_axis, self.triple.semimajor_axis)
                print('Ecc:', bin.eccentricity, self.triple.eccentricity)
                print('Stellar type:', bin.child1.stellar_type, bin.child2.stellar_type, star.stellar_type)
                print('Spin:', bin.child1.spin_angular_frequency, bin.child2.spin_angular_frequency, star.spin_angular_frequency)
                

            bin.child1.is_donor = False
            bin.child2.is_donor = False
            star.is_donor = False
                
            if bin.child1.radius >= Rl1 - (1.0 * small_numerical_error|units.RSun):
                bin.child1.is_donor = True
            if bin.child2.radius >= Rl2 - (1.0 * small_numerical_error|units.RSun):
                bin.child2.is_donor = True
            if star.radius >= Rl3 - (1.0 * small_numerical_error|units.RSun):
                star.is_donor = True

                
        else:
            sys.exit('check_RLOF: structure stellar system unknown')


    def check_OLOF(self):
        if self.triple.is_star:
            return
        elif self.is_binary():
            Rl2_1 = L2_radius(self, self.child1)
            Rl2_2 = L2_radius(self, self.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print('L2 lobe radii:', Rl2_1, Rl2_2)
                print('Stellar radii:', self.triple.child1.radius, self.triple.child2.radius)
            
            self.triple.child1.is_OLOF_donor = False
            self.triple.child2.is_OLOF_donor = False
            
            if self.triple.child1.radius >= Rl2_1 - (1.0 * small_numerical_error|units.RSun):
                self.triple.child1.is_OLOF_donor = True               
            if self.triple.child2.radius >= Rl2_2 - (1.0 * small_numerical_error|units.RSun):
                self.triple.child2.is_OLOF_donor = True                             
 
        elif self.is_triple():
            # for disrupted binary
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1
        
            Rl2_1 = L2_radius(bin, bin.child1, self)
            Rl2_2 = L2_radius(bin, bin.child2, self)
            
            bin.child1.is_OLOF_donor = False
            bin.child2.is_OLOF_donor = False
            
            if bin.child1.radius >= Rl2_1 - (1.0 * small_numerical_error|units.RSun):
                bin.child1.is_OLOF_donor = True               
            if bin.child2.radius >= Rl2_2 - (1.0 * small_numerical_error|units.RSun):
                bin.child2.is_OLOF_donor = True                             
    
            if REPORT_TRIPLE_EVOLUTION:
                print('L2 lobe radii:', Rl2_1, Rl2_2)
                print('Stellar radii:', bin.child1.radius, bin.child2.radius)
                print('Masses:', bin.child1.mass, bin.child2.mass, star.mass)
                print('Semi:', bin.semimajor_axis, self.triple.semimajor_axis)
                print('Ecc:', bin.eccentricity, self.triple.eccentricity)
                print('Stellar type:', bin.child1.stellar_type, bin.child2.stellar_type, star.stellar_type)
                print('Spin:', bin.child1.spin_angular_frequency, bin.child2.spin_angular_frequency, star.spin_angular_frequency)
               
        else:
            print('check_OLOF: structure stellar system unknown')
            sys.exit('check_OLOF: structure stellar system unknown')     


    def check_CHE(self): #future option: potentially use: che_flag in SeBa
        #returns true when one or both of the inner binary components are chemically homogeneously evolving
        #if !include_CHE, then return False
        #problematic for quadruples - what if one binary is CHE, and other is not
        
        metallicity = self.stellar_code.parameters.metallicity
        if self.triple.is_star:
            return False
        elif self.is_binary():  
            bin = self.triple
            if self.include_CHE and ((bin.child1.spin_angular_frequency >= criticial_angular_frequency_CHE(bin.child1.mass, metallicity) and bin.child1.stellar_type <= 1|units.stellar_type) or 
                    (bin.child2.spin_angular_frequency >= criticial_angular_frequency_CHE(bin.child2.mass, metallicity) and bin.child2.stellar_type <= 1|units.stellar_type)):
                 return True
            else: 
                return False
        
        elif self.is_triple():
            if self.triple.child1.is_star:
                bin = self.triple.child2
            else:
                bin = self.triple.child1

            if self.include_CHE and ((bin.child1.spin_angular_frequency >= criticial_angular_frequency_CHE(bin.child1.mass, metallicity) and bin.child1.stellar_type <= 1|units.stellar_type) or 
                    (bin.child2.spin_angular_frequency >= criticial_angular_frequency_CHE(bin.child2.mass, metallicity) and bin.child2.stellar_type <= 1|units.stellar_type)):
                 return True
            else: 
                return False
        else:
            print('check_CHE: structure stellar system unknown')
            sys.exit('check_CHE: structure stellar system unknown')



   
                     
#            
#    def determine_partial_timestep_stable_mass_transfer(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#
#        if stellar_system.is_star:
#            return np.inf |units.Myr 
#        else:
#            dt1 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child1)        
#            dt2 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child2)
#            dt =  stellar_system.part_dt_mt
#            return min(dt, min(dt1, dt2))
                   
    #-------

    #-------
    # useful functions for printing
    def print_star(self, star):
        if star.is_star:
            print('star:')
            print(star.age, )
            print(star.stellar_type, )
            print(star.mass, )
            print(star.radius, )
            print(star.core_mass, )
            print(star.core_radius,)
            print(star.convective_envelope_mass,)
            print(star.convective_envelope_radius,)
            print(star.luminosity,)
            print(star.temperature,)
            print(star.wind_mass_loss_rate,)
            print(star.spin_angular_frequency,)
            print(star.is_donor)
            print('\t')
        else:
            sys.exit('print_star needs a star')
    
    
    def print_binary(self, binary):
        if not binary.is_star:
            print(self.get_mass(binary), )
            print(binary.semimajor_axis,) 
            print(binary.eccentricity, )
            print(binary.argument_of_pericenter, )
            print(binary.longitude_of_ascending_node,)
            print(binary.mass_transfer_rate,)
#            print(binary.mass_transfer_timescale,)
            print(binary.accretion_efficiency_mass_transfer,)
            print(binary.accretion_efficiency_wind_child1_to_child2,)
            print(binary.accretion_efficiency_wind_child2_to_child1,)
            print(binary.specific_AM_loss_mass_transfer,)
            print(binary.is_mt_stable)
            print('\t')
        else:
            sys.exit('print_binary needs a binary')        

    
    
    def print_stellar_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            print(stellar_system.number,)
            print(stellar_system.relative_inclination)

        if stellar_system.is_star:
            self.print_star(stellar_system)
        else:
            print('binary star: ')
            self.print_binary(stellar_system)
            self.print_stellar_system(stellar_system.child1)
            self.print_stellar_system(stellar_system.child2)
        print('\t')
    #-------
    #-------
    #don't change this unless you know what you're doing
    def remove_parents(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        parents = []
        if stellar_system.is_star:
            try:
                p = stellar_system.parent
                stellar_system.parent = 0
                return p
            except AttributeError: #when there is no parent
                return 0
        else:
            parents.append(self.remove_parents(stellar_system.child1))
            parents.append(self.remove_parents(stellar_system.child2))

            p = stellar_system.parent
#            except AttributeError: #when there is no parent=
            if p != None:
                stellar_system.parent = 0
                parents.append(p)                                    

            return parents
            
    #don't change this unless you know what you're doing
    def set_parents(self, parents, stellar_system=None):            
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
                stellar_system.parent = parents
        else:
            self.set_parents(parents[0], stellar_system.child1)
            self.set_parents(parents[1], stellar_system.child2)
            if len(parents) == 3:
                stellar_system.parent = parents[2]
            elif len(parents) != 2:
                sys.exit('set_parents: structure stellar system unknown') 

 
    def save_snapshot(self):
        file_name = self.file_name

        snapshot_triple = self.triple.copy()
        #make nice recursive loop
        #find star & bin particles sets 
        if snapshot_triple.is_star:
            star_particle_set = snapshot_triple.particles_set
            bin_particle_set = None
        else:
            bin_particle_set = snapshot_triple.particles_set
            
            if snapshot_triple.child1.is_star:
                star_particle_set = snapshot_triple.child1.particles_set
            else:
                sys.exit('save_snapshot: structure stellar system unknown')        

        #deleting redundant parameters of stellar particle set
        del star_particle_set.apsidal_motion_constant
        del star_particle_set.convective_envelope_mass
        del star_particle_set.convective_envelope_radius
        del star_particle_set.core_radius
        del star_particle_set.gyration_radius
        del star_particle_set.initial_mass
        del star_particle_set.moment_of_inertia_of_star
        del star_particle_set.previous_moment_of_inertia_of_star
        del star_particle_set.previous_radius
        del star_particle_set.previous_spin_angular_frequency
        del star_particle_set.previous_stellar_type
        del star_particle_set.previous_time_derivative_of_radius
        del star_particle_set.time_derivative_of_radius
        del star_particle_set.wind_mass_loss_rate

        #deleting redundant parameters of binary particle set 
        del bin_particle_set.accretion_efficiency_mass_transfer
        del bin_particle_set.accretion_efficiency_wind_child1_to_child2
        del bin_particle_set.accretion_efficiency_wind_child2_to_child1
        del bin_particle_set.mass_transfer_rate
        del bin_particle_set.part_dt_mt
        del bin_particle_set.previous_bin_type
        del bin_particle_set.previous_kozai_type
        del bin_particle_set.previous_mass
        del bin_particle_set.specific_AM_loss_mass_transfer


    #double parameters
#            del snapshot_triple.child2.particles_set.CPU_time
#            del snapshot_triple.child2.particles_set.time
#            del snapshot_triple.child2.particles_set.dynamical_instability
#            del snapshot_triple.child2.particles_set.error_flag_secular
#            del snapshot_triple.child2.particles_set.delta_e_in
#            del snapshot_triple.child2.particles_set.max_delta_e_in
#            del snapshot_triple.child2.particles_set.number

    #some minor parameters are not saved:
#        self.instantaneous_evolution = False 
#        self.tend = tend #...
#        self.triple.time = 0.0|units.yr
#        self.previous_time = 0.0|units.yr
#        self.file_name
#        self.file_type


        if self.file_type == 'txt':
            print(self.file_name,self.file_type)
            bin_particle_set.parent = 0
            star_particle_set.parent = 0
            write_set_to_file(snapshot_triple.as_set(), self.file_name, self.file_type) 
        else:
            write_set_to_file(snapshot_triple.as_set(), self.file_name, self.file_type, version='2.0', append_to_file=True)

        del snapshot_triple
            
        self.triple.child2.max_delta_e_in = 0


    #-------
        
    #-------
    # determining time steps
    def determine_time_step_wind(self, stellar_system = None):
    #note: returned value can be inf when the wind_mass_loss_rate = 0
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr

            if stellar_system.wind_mass_loss_rate * -1. > quantities.zero:
                dt = maximum_wind_mass_loss_factor * (stellar_system.mass - stellar_system.core_mass)/stellar_system.wind_mass_loss_rate*-1.  
            if REPORT_DT:
                print("Dt_wind_star = ", dt)
            return dt 
        else:
            dt1 = self.determine_time_step_wind(stellar_system.child1)        
            dt2 = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_DT:
                print("Dt_wind_binary = ", dt1, dt2)
            return min(dt1, dt2) 
    


    def determine_time_step_radius_change(self, stellar_system = None):
        #note: returned value can be inf when the change in radius <= 0
        #radius is only necessary for tides

        if not self.secular_code.parameters.include_inner_tidal_terms and not self.secular_code.parameters.include_outer_tidal_terms:
            return np.inf |units.Myr

        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            if stellar_system.time_derivative_of_radius == quantities.zero:
                dt = np.inf |units.Myr
            else:     
                if stellar_system.previous_time_derivative_of_radius == quantities.zero:
                    growth_factor = 0.1
                elif stellar_system.previous_time_derivative_of_radius * stellar_system.time_derivative_of_radius < quantities.zero:
                    growth_factor = 0.01
                else:
                    growth_factor = stellar_system.previous_time_derivative_of_radius/stellar_system.time_derivative_of_radius 
                    if growth_factor > 1:
                        growth_factor = 1./growth_factor
                
                if stellar_system.stellar_type > 1|units.stellar_type and stellar_system.time_derivative_of_radius < stellar_system.previous_time_derivative_of_radius: #not a MS star and Rdot < Rdot_prev
                    growth_factor = 1. 

                dt = abs(growth_factor * self.maximum_radius_change_factor*stellar_system.radius / stellar_system.time_derivative_of_radius)


            if REPORT_DT:
                print("Dt_radius_change_star = ", dt)
            return dt 
        else:
            dt1 = self.determine_time_step_radius_change(stellar_system.child1)        
            dt2 = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_DT:
                print("Dt_radius_change_binary = ", dt1, dt2)
            return min(dt1, dt2) 

    def determine_time_step_kozai(self):
    #note: returned value can be inf when the system is a binary or single star            
        if self.is_triple():
            dt = self.kozai_timescale()*time_step_factor_kozai
            if REPORT_DT:
                print("Dt_kozai = ", dt)
        else:
            dt = np.inf |units.Myr
        
        P_out = self.orbital_period(self.triple)
        P_in = self.orbital_period(self.triple.child2)        
        return dt

 
    def determine_time_step_stable_mt(self, stellar_system = None):
    #note: returned value can be inf when the mass_transfer_rate = 0 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.is_donor:
                dt = abs(time_step_factor_stable_mt*stellar_system.mass/stellar_system.parent.mass_transfer_rate)
            if REPORT_DT:
                print("Dt_mt_star = ", dt, time_step_factor_stable_mt, stellar_system.mass, stellar_system.parent.mass_transfer_rate )
            return dt
        else:
            dt1 = self.determine_time_step_stable_mt(stellar_system.child1)        
            dt2 = self.determine_time_step_stable_mt(stellar_system.child2)
            if REPORT_DT:
                print("Dt_mt_binary = ", dt1, dt2)
            return min(dt1, dt2) 
         
         
    #based on secular code      
    def e_dot_tides(self, star,m_comp, semi, eccentricity):
        e_p2 = eccentricity**2
        l = np.sqrt(1.0 - e_p2)        
        f_tides3 = 1.0 + e_p2*(15./4. + e_p2*(15./8. + e_p2*5./64.))
        f_tides4 = 1.0 + e_p2*(3./2 + e_p2/8.)
    
        R_div_a = star.radius/semi        
        spin = star.spin_angular_frequency
        k_div_T_tides = tidal_friction_constant(star.stellar_type, star.mass, m_comp, semi, star.radius, star.convective_envelope_mass, star.convective_envelope_radius, star.luminosity, spin, star.gyration_radius, star.apsidal_motion_constant)
        n = corotating_spin_angular_frequency_binary(semi, star.mass, m_comp) # mean orbital angular speed

        e_dot = -27.0*(1.0+m_comp/star.mass)*(m_comp/star.mass) * R_div_a**6 * k_div_T_tides * R_div_a**2 *eccentricity* (l**-13.0) * (f_tides3 - 11.0/18.0*(l**3)*f_tides4*(spin/n))
#        print(star.mass, m_comp, semi, eccentricity, k_div_T_tides, e_dot)
        return e_dot
         
         
         
    def determine_time_step_tides(self, stellar_system = None):
#        print("determine_time_step_tides")
        if stellar_system == None:
            stellar_system = self.triple


        if self.triple.is_star:
            return False, 0        
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if self.triple.child1.radius >= Rl1 or self.triple.child2.radius >= Rl2:
                return abs(time_step_factor_stable_mt*min(self.triple.child1.mass, self.triple.child2.mass)/self.triple.mass_transfer_rate)
            else: 
                de_dt_child1 = abs(self.e_dot_tides(self.triple.child1, self.triple.child2.mass, self.triple.semimajor_axis, self.triple.eccentricity))
                de_dt_child2 = abs(self.e_dot_tides(self.triple.child2, self.triple.child1.mass, self.triple.semimajor_axis, self.triple.eccentricity))
                return time_step_factor_ecc * self.triple.eccentricity / max(de_dt_child1, de_dt_child2)
        elif self.is_triple():
        
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1

            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
            if star.radius >= Rl3:
                sys.exit('R>Rl3')
            elif bin.child1.radius >= Rl1 or bin.child2.radius >= Rl2:
                return abs(time_step_factor_stable_mt*min(bin.child1.mass, bin.child2.mass)/self.triple.mass_transfer_rate)
            else:
                de_dt_star = abs(self.e_dot_tides(star, self.get_mass(bin), self.triple.semimajor_axis, self.triple.eccentricity))
                de_dt_bin_child1 = abs(self.e_dot_tides(bin.child1, bin.child2.mass, bin.semimajor_axis, bin.eccentricity))
                de_dt_bin_child2 = abs(self.e_dot_tides(bin.child2, bin.child1.mass, bin.semimajor_axis, bin.eccentricity))
                dt_star = time_step_factor_ecc * self.triple.eccentricity / de_dt_star
                dt_bin = time_step_factor_ecc * bin.eccentricity / max(de_dt_bin_child1, de_dt_bin_child2)
                return min(dt_star, dt_bin)

        else: 
            sys.exit('determine_time_step_tides: structure stellar system unknown')        
   

    
         
         
#    def close_to_RLOF(self):
#
#        if self.triple.is_star:
#            return False, 0        
#        elif self.is_binary():
#            Rl1 = roche_radius(self, self.child1)
#            Rl2 = roche_radius(self, self.child2)
#            if self.triple.child1.radius >= Rl1 or self.triple.child2.radius >= Rl2:
#                print('close to rlof already rlof')
#                exit(0)
#            elif self.triple.child1.radius >= Rl_fraction*Rl1 and self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = maximum(self.triple.child1.radius/Rl1, self.triple.child2.radius/Rl2)
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.child2.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*min(self.triple.child1.mass, self.triple.child2.mass)/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof = self.triple.child1.radius/Rl1
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*self.triple.child1.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = self.triple.child2.radius/Rl2
#                print('dt_close_to_mt', time_step_factor_stable_mt, self.triple.child2.mass, self.triple.mass_transfer_rate)
#                dt = abs(time_step_factor_stable_mt*self.triple.child2.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt              
#            else:
#                return False, 0, 0    
#        elif self.is_triple():
#            if self.triple.child1.is_star:
#                star = self.triple.child1
#                bin = self.triple.child2
#            else:
#                star = self.triple.child2
#                bin = self.triple.child1
#
#            #assumping secular code always returns inner binary first
#            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
#            if bin.child1.radius >= Rl1 or bin.child2.radius >= Rl2 or star.radius >= Rl3:
#                print('close to rlof already rlof')
#                exit(0)
#
#            dt=[]
#            ratio_rad_rlof=[]
#            if star.radius >= Rl_fraction*Rl3:
#                ratio_rad_rlof.append(star.radius / Rl3)
#                print('dt_close_to_mt', time_step_factor_stable_mt, star.mass, star.parent.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*star.mass/star.parent.mass_transfer_rate))
#            if bin.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof.append(bin.child1.radius / Rl1)
#                print('dt_close_to_mt', time_step_factor_stable_mt, bin.child1.mass, bin.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*bin.child1.mass/bin.mass_transfer_rate))
#            if bin.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof.append(bin.child2.radius / Rl2)
#                print('dt_close_to_mt', time_step_factor_stable_mt, bin.child2.mass, bin.mass_transfer_rate)
#                dt.append(abs(time_step_factor_stable_mt*bin.child2.mass/bin.mass_transfer_rate))
#                    
#            if len(dt) > 0:
#                print('close to rlof', ratio_rad_rlof, dt)
#                return True, max(ratio_rad_rlof), min(dt)
#            else:
#                return False, 0, 0
#        else:
#            print('close_to_RLOF: structure stellar system unknown')        
#            exit(2)    
    
        
    def determine_time_step(self):         
        if REPORT_DT:
            print("Dt = ", self.stellar_code.particles.time_step, self.tend, self.previous_time)

        #maximum time_step            
        time_step_max = self.tend - self.triple.time      
        
        #make sure we hit the timestamp of tinit
        if self.triple.time < self.tinit:
            time_step_max = min(time_step_max, self.tinit - self.triple.time)
                               
        # time_step of stellar evolution
        time_step_stellar_code = self.stellar_code.particles.time_step

        # small time_step after type change
        # done automatically by SeBa when due to wind mass losses, but not during RLOF
        if self.has_stellar_type_changed():
            time_step_stellar_code.append(minimum_time_step)                    
                    
        # small time_step during heavy wind mass losses
        time_step_wind =  self.determine_time_step_wind()

        # small time_step during phases of fast growth (note: not yet during fast shrinkage)
        time_step_radius_change = self.determine_time_step_radius_change()             
        
#        #during stable mass transfer   
#        if self.has_donor():
##            print(time_step, self.determine_time_step_stable_mt())
#            time_step = min(time_step, self.determine_time_step_stable_mt())
#            if REPORT_DT or REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
#                print('donor time:',  self.determine_time_step_stable_mt())

        #tides - ADD ALSO DURING MT
        time_step_tides = np.inf |units.Myr 
	#interesting alternative, slows down code 
#        if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:    
#             time_step_tides = self.determine_time_step_tides()  	
                
        if REPORT_DT or REPORT_DEBUG:
            print('time:', self.triple.time, time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step_tides)



        if not self.has_donor():
            if REPORT_DT:
                print('no RLOF' )
            time_step = min(time_step_radius_change, min(time_step_wind, min( min(time_step_stellar_code), time_step_max)))
            if time_step_tides < time_step:
                Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
                if REPORT_DT:
                    print(time_step, time_step_tides)
                    print(Rl1,Rl2,Rl3)
                    print(self.triple.child2.child1.radius, self.triple.child2.child2.radius, self.triple.child1.radius )
                time_step = time_step_tides

            if REPORT_DT or REPORT_DEBUG:
                print('min increase', time_step, maximum_time_step_factor*self.previous_dt)   
                
            #maybe not relevant with MESA                                 
            time_step = min(time_step, maximum_time_step_factor*self.previous_dt) 


        elif self.has_donor() and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached':
            #find beginning of RLOF carefully
            time_step = time_step_factor_find_RLOF*self.previous_dt
            if REPORT_DT:
                print('find rlof', time_step_factor_find_RLOF,self.previous_dt)

            #resetting is_donor
            self.check_RLOF()

        else:
            #during stable mass transfer   
            time_step = min(time_step_wind, min( min(time_step_stellar_code), time_step_max))
            time_step = min(time_step, self.determine_time_step_stable_mt())
            if REPORT_DT or REPORT_DEBUG:
                print('donor time:',  self.determine_time_step_stable_mt())

            if REPORT_DT or REPORT_DEBUG:
                print('min increase', time_step, maximum_time_step_factor*self.previous_dt)                                
            time_step = min(time_step, maximum_time_step_factor*self.previous_dt)  


        if self.triple.child2.bin_type == 'detached' and self.triple.child2.previous_bin_type == 'stable_mass_transfer':
            corotation_spin = corotating_spin_angular_frequency_binary(self.triple.child2.semimajor_axis, self.triple.child2.child1.mass, self.triple.child2.child2.mass)
            self.triple.child2.child1.spin_angular_frequency = corotation_spin
            self.triple.child2.child2.spin_angular_frequency = corotation_spin
            
#            if REPORT_DT or REPORT_DEBUG:
#                print('prev timestep', time_step, previous_time_step)
            previous_time_step = self.triple.time - self.previous_time
            time_step = min(time_step, maximum_time_step_factor_after_stable_mt*previous_time_step)  



        if self.triple.time == quantities.zero:
            #initialization (e.g. time_derivative_of_radius)
            P_out = self.orbital_period(self.triple) #period outer binary 
            # do not take 0.1*P_in -> resonance -> large error
            time_step = min(min(P_out, time_step), 1.|units.yr)


        time_step = max(time_step, minimum_time_step)  
        if self.triple.time >= self.tinit:
            time_step = min(time_step, maximum_time_step)  





#        else:
            #during run-up towards mass transfer 
#            close_to_RLOF_bool, ratio_rad_rlof, time_step_close_to_mt = self.close_to_RLOF()     
#            if close_to_RLOF_bool:
#                if ratio_rad_rlof >= 1.0:
#                    print('ratio_rad_rlof >1?')
#                    exit(0)
#                else:
#                    print('close to rlof', time_step_close_to_mt, time_step, ratio_rad_rlof)
#                    exit(0)
#                    time_step = min(time_step, time_step_close_to_mt)


      # in case secular code finds RLOF 
        if self.fixed_timestep > 0.|units.yr:
            if REPORT_DT or REPORT_DEBUG:
                print('time:', self.fixed_timestep)
                print('donor time:',  self.determine_time_step_stable_mt())

            t_donor = self.determine_time_step_stable_mt()*time_step_factor_stable_mt # extra small for safety
            t_donor_lim = max(minimum_time_step, min(time_step, t_donor))#
            # although fixed_timestep < time_step, fixed_timestep can be > time_step_stable_mt

            if t_donor == np.inf|units.Myr:
                #bh or ns donors
                time_step = self.fixed_timestep
            elif self.fixed_timestep < t_donor_lim:
                time_step = t_donor_lim
                self.secular_code.parameters.check_for_inner_RLOF = False
            else:
                time_step = self.fixed_timestep
                           
            self.fixed_timestep = -1|units.Myr
            return time_step # as previous timestep was effectively zero -> min increase is zero           

        return time_step
    #-------

    #-------
    #SN functions
    #based on SeBa (Portegies Zwart et al. 1996, Toonen et al. 2012)
    
    def anomaly_converter(self, e_a, ecc, m_a):
    #    mean_anomaly = eccentric_anomaly - ecc*sin(eccentric_anomaly)
        return e_a - ecc * np.sin(e_a) - m_a

    #Paczynski 1990, 348, 485
    def kick_velocity_paczynski(self):
        return self.random_paczynski_velocity(270|units.km/units.s)

    def kick_velocity_SeBa_std(self):
        return self.random_paczynski_velocity(600|units.km/units.s)

    #Paczynski 1990, 348, 485
    #taken from SeBa (Portegies Zwart et al. 1996, Toonen et al. 2012)
    def random_paczynski_velocity(self, v_disp):
    #Velocity distribution used by Paczynski, B., 1990, ApJ 348, 485.
    #with a dispersion of 270 km/s.
    #Phinney uses the same distribution with a dispersion of 600 km/s.
    #The distribution:
    #P(u)du = \frac{4}{\pi} \cdot \frac{du}{(1+u^2)^2},
    #u=v/v_d,
    #v_d is the dispersion velocity = 270 or 600 km/s respectively.

        max_distr = 4./np.pi
        v_max = 4
        prob = max_distr*np.random.uniform(0., 1.)
        velocity = v_max*np.random.uniform(0., 1.)
        dist_value = max_distr/(1+velocity**2.)**2.
    
        while(prob > dist_value):
            prob = max_distr*np.random.uniform(0., 1.)
            velocity = v_max*np.random.uniform(0., 1.)
            dist_value = max_distr/(1+velocity**2.)**2.
        
        vector = self.random_direction()
        return velocity*v_disp*vector

    #Hansen & Phinney 1997, 291, 569
    # V_3d = 300, V_1d = 190 km/s
    def kick_velocity_hansen(self):
        vector = self.random_direction()
        return vector*maxwell.rvs(loc=0, scale = 190)|units.km/units.s
    
    #Hobbs, Lorimer, Lyne & Kramer, 2005, 360, 974
    # V_3d = 400, V_1d = 265 km/s
    def kick_velocity_hobbs(self):
        vector = self.random_direction()
        return vector*maxwell.rvs(loc=0, scale = 265)|units.km/units.s
     
     #Arzoumanian ea 2002, 568, 289
    def kick_velocity_arzoumanian(self):
        vector = self.random_direction()
        if np.random.uniform(0,1) < 0.4:
            return vector*maxwell.rvs(loc=0, scale = 90)|units.km/units.s
        else:
            return vector*maxwell.rvs(loc=0, scale = 500)|units.km/units.s


     # Verbunt, Igoshev & Cator, 2017, 608, 57 - combination of two maxwellians
    def kick_velocity_verbunt(self):
        vector = self.random_direction()
        if np.random.uniform(0,1) < 0.42:
            return vector*maxwell.rvs(loc=0, scale = 75)|units.km/units.s
        else:
            return vector*maxwell.rvs(loc=0, scale = 316)|units.km/units.s


    
    def random_direction(self):
        theta_kick = np.arccos(np.random.uniform(-1, 1))
        phi_kick   = np.random.uniform(0, 2*np.pi)
        x_kick = np.sin(theta_kick)*np.cos(phi_kick)
        y_kick = np.sin(theta_kick)*np.sin(phi_kick)
        z_kick = np.cos(theta_kick)
        return np.array([x_kick, y_kick, z_kick])
	    

    def get_SN_kick(self, star):
        v_kick = [0.,0.,0.]|units.kms
        if star.stellar_type != star.previous_stellar_type and star.stellar_type in stellar_types_SN_remnants:
            if self.SN_kick_distr == 0: 
                v_kick =  [0.,0.,0.]|units.kms
            elif self.SN_kick_distr in [1]: # Hobbs, Lorimer, Lyne & Kramer, 2005, 360, 974
                v_kick = self.kick_velocity_hobbs()
            elif self.SN_kick_distr in [2]: #Arzoumanian ea 2002, 568, 289
                v_kick = self.kick_velocity_arzoumanian()
            elif self.SN_kick_distr in [3]: #Hansen & Phinney 1997, 291, 569
                v_kick = self.kick_velocity_hansen()
            elif self.SN_kick_distr in [4]: # Paczynski 1990, 348, 485
                v_kick = self.kick_velocity_paczynski()
            elif self.SN_kick_distr in [5]: # Verbunt, Igoshev & Cator, 2017, 608, 57
                v_kick = self.kick_velocity_verbunt()
            else:
                print("SN kick distribution not specified, assuming no kick")
                v_kick =  [0.,0.,0.]|units.kms


            print(self.impulse_kick_for_black_holes, self.fallback_kick_for_black_holes, v_kick)
            #reduce kick of BH by conservation of momentum
            if self.impulse_kick_for_black_holes and star.stellar_type == 14|units.stellar_type:
                v_kick *= (kanonical_neutron_star_mass / star.mass)
#                print(star.mass, kanonical_neutron_star_mass)
            if self.fallback_kick_for_black_holes and star.stellar_type == 14|units.stellar_type:
#                self.channel_from_stellar.copy_attributes(["fallback"]) 
#                v_kick *= (1-star.fallback)

                star_in_stellar_code = star.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]                
                fallback = star_in_stellar_code.get_fallback() 
                v_kick *= (1-fallback)
#                print(fallback)                    

        return v_kick

    def save_mean_anomalies_at_SN(self, inner_mean_anomaly, outer_mean_anomaly, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                stellar_system.inner_mean_anomaly = inner_mean_anomaly
                stellar_system.outer_mean_anomaly = outer_mean_anomaly                
        else:
            self.save_mean_anomalies_at_SN(inner_mean_anomaly, outer_mean_anomaly,stellar_system.child1)        
            self.save_mean_anomalies_at_SN(inner_mean_anomaly, outer_mean_anomaly,stellar_system.child2)

    def adjust_spin_after_supernova(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants:
                stellar_system.spin_angular_frequency *= stellar_system.previous_moment_of_inertia_of_star/stellar_system.moment_of_inertia_of_star
        else:
            self.adjust_spin_after_supernova(stellar_system.child1)        
            self.adjust_spin_after_supernova(stellar_system.child2)
 
        
    def adjust_system_after_supernova_kick(self):
        if self.triple.is_star:
            #only velocity kick needs to be taken into account here
            #not followed currently
            return False
        elif not self.is_triple():
            sys.exit('SN only implemented in triple')
                    
        #SN in triple
        if self.triple.child1.is_star:
            star = self.triple.child1
            bin = self.triple.child2
        else:
            star = self.triple.child2
            bin = self.triple.child1

        #needs BH mass    
        vel_kick1 = self.get_SN_kick(bin.child1)
        vel_kick2 = self.get_SN_kick(bin.child2)
        vel_kick3 = self.get_SN_kick(star)
    
        #convention in secular code: dm > 0
        dm1 = bin.child1.previous_mass - bin.child1.mass
        dm2 = bin.child2.previous_mass - bin.child2.mass
        dm3 = star.previous_mass - star.mass
 
        #secular code requires the mass of the stars to be the mass prior to the SN
        bin.child1.mass = bin.child1.previous_mass 
        bin.child2.mass = bin.child2.previous_mass
        star.mass = star.previous_mass 
        
        #reset parameters to before SN 
        #to get correct values for the snapshot
        bin_child1_proper_stellar_type = bin.child1.stellar_type
        bin_child2_proper_stellar_type = bin.child2.stellar_type
        star_proper_stellar_type = star.stellar_type
        bin.child1.stellar_type = bin.child1.previous_stellar_type
        bin.child2.stellar_type = bin.child2.previous_stellar_type
        star.stellar_type = star.previous_stellar_type                

        bin_child1_proper_radius = bin.child1.radius
        bin_child2_proper_radius = bin.child2.radius
        star_proper_radius = star.radius
        bin.child1.radius = bin.child1.previous_radius
        bin.child2.radius = bin.child2.previous_radius
        star.radius = star.previous_radius                
 
        self.save_snapshot()                    
 
         #reset orbital parameters tertiary in case of binary or ignore_tertiary
        a_out_init = self.triple.semimajor_axis
        e_out_init = self.triple.eccentricity              

        # determine stellar anomaly
        # mean anomaly increases uniformly from 0 to 2\pi radians during each orbit
        inner_ecc = bin.eccentricity

        inner_mean_anomaly = np.random.uniform(0., 2.*np.pi)  
        inner_eccentric_anomaly = optimize.brentq(self.anomaly_converter, 0., 2.*np.pi, args=(inner_ecc, inner_mean_anomaly))
        inner_true_anomaly = 2.* np.arctan2(np.sqrt(1-inner_ecc) * np.cos(inner_eccentric_anomaly/2.), np.sqrt(1+inner_ecc) * np.sin(inner_eccentric_anomaly/2.))
        # min teken verschil met cos_true_anomaly = (np.cos(eccentric_anomaly)-ecc) / (1-ecc*np.cos(eccentric_anomaly))

        outer_ecc = self.triple.eccentricity
        outer_mean_anomaly = np.random.uniform(0., 2.*np.pi)  
        outer_eccentric_anomaly = optimize.brentq(self.anomaly_converter, 0., 2.*np.pi, args=(outer_ecc, outer_mean_anomaly))
        outer_true_anomaly = 2.* np.arctan2(np.sqrt(1-outer_ecc) * np.cos(outer_eccentric_anomaly/2.), np.sqrt(1+outer_ecc) * np.sin(outer_eccentric_anomaly/2.))

        if REPORT_SN_EVOLUTION:
            print('before SN:', self.triple.number) 
            print('\n\ntime:', self.triple.time, self.previous_time, self.triple.time-self.previous_time)
            print('eccentricity:', inner_ecc, outer_ecc)
            print('semi-major axis:', bin.semimajor_axis, self.triple.semimajor_axis)
            print('inner anomaly:', inner_mean_anomaly, inner_eccentric_anomaly, inner_true_anomaly)
            print('outer anomaly:', outer_mean_anomaly, outer_eccentric_anomaly, outer_true_anomaly)
            print('vel_kicks:', vel_kick1, vel_kick2, vel_kick3  )
            print('dm:', dm1, dm2, dm3)        

        self.channel_to_secular.copy()
        V1, V2, V3 = self.secular_code.compute_effect_of_SN_on_triple(vel_kick1, vel_kick2, vel_kick3, dm1, dm2, dm3, inner_true_anomaly, outer_true_anomaly)
        self.channel_from_secular.copy() 
        
        #as secular code required the mass of the stars to be the mass prior to the SN
        #reset mass to after SN
        bin.child1.mass = bin.child1.previous_mass - dm1
        bin.child2.mass = bin.child2.previous_mass - dm2
        star.mass = star.previous_mass - dm3

        #reset the stellar types to after SN 
        bin.child1.stellar_type = bin_child1_proper_stellar_type
        bin.child2.stellar_type = bin_child2_proper_stellar_type
        star.stellar_type =  star_proper_stellar_type

        bin.child1.radius = bin_child1_proper_radius
        bin.child2.radius = bin_child2_proper_radius
        star.radius =  star_proper_radius

        if self.secular_code.parameters.ignore_tertiary:
            self.triple.semimajor_axis = a_out_init
            self.triple.eccentricity = e_out_init

        if REPORT_SN_EVOLUTION:
            print('after SN')
            print('eccentricity:', bin.eccentricity, self.triple.eccentricity)
            print('semi-major axis:', bin.semimajor_axis, self.triple.semimajor_axis)

        if bin.eccentricity >= 1.0 or bin.eccentricity < 0.0 or bin.semimajor_axis <=0.0|units.RSun or np.isnan(bin.semimajor_axis.value_in(units.RSun)):
            if REPORT_SN_EVOLUTION:
                print("Inner orbit dissociated by SN at time = ",self.triple.time)                 
            bin.bin_type = bin_type['disintegrated']   
            return False
        elif self.triple.eccentricity >= 1.0 or self.triple.eccentricity < 0.0 or self.triple.semimajor_axis <=0.0|units.RSun or np.isnan(self.triple.semimajor_axis.value_in(units.RSun)):
            if REPORT_SN_EVOLUTION:
                print("Outer orbit dissociated by SN at time = ",self.triple.time) 
            self.triple.bin_type = bin_type['disintegrated']   

	    # When the outer orbit has disintegrated, change its orbital parameters such
	    # that it has no influence on the inner orbit for the remainder of the simulation
            self.triple.semimajor_axis = 1e100|units.RSun
            self.triple.eccentricity = 0

            self.secular_code.parameters.ignore_tertiary = True
            self.secular_code.parameters.check_for_dynamical_stability = False
            self.secular_code.parameters.check_for_outer_collision = False
            self.secular_code.parameters.check_for_outer_RLOF = False

            # The evolution could be followed further by uncommenting the next line
            # For now we stop the simulation if we do not have a bound triple system anymore
            return False 

        self.adjust_spin_after_supernova() 

        #bh spin frequency has no meaning. make sure it doesn't affect the evolution
        bin.child1.previous_moment_of_inertia_of_star = bin.child1.moment_of_inertia_of_star
        bin.child2.previous_moment_of_inertia_of_star = bin.child2.moment_of_inertia_of_star
        star.previous_moment_of_inertia_of_star =  star.moment_of_inertia_of_star
                
        if bin.child1.stellar_type in stellar_types_SN_remnants:
            self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star1 = False
        if bin.child2.stellar_type in stellar_types_SN_remnants:
            self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star2 = False
        if star.stellar_type in stellar_types_SN_remnants:
            self.secular_code.parameters.include_spin_radius_mass_coupling_terms_star3 = False

        if self.stop_at_SN:
            if REPORT_SN_EVOLUTION:
                print("Supernova at time = ",self.triple.time )
                print(self.triple.child2.child1.mass, self.triple.child2.child1.stellar_type)
                print(self.triple.child2.child2.mass, self.triple.child2.child2.stellar_type)
                print(self.triple.child1.mass, self.triple.child1.stellar_type)
            return False   
            
        #check for rlof -> if rlof than stop?
#        self.check_RLOF()

        return True
        
    #-------

    #-------
    #evolution
 
    # uses instantaneous eccentricity (from global timestep) to calculate an orbit averaged stellar flux, 
    # that can evaporise part of the planetary envelope. 
    # would be best to do in secular code instead, such that eccentricity variations due to KL are taken into account properly
    def planetary_mass_evaporation(self, dt, stellar_system = None):        
        if stellar_system == None:
           stellar_system = self.triple
        
        if stellar_system.is_star:
            return
        elif self.is_binary(stellar_system):
            if stellar_system.child1.stellar_type in stellar_types_planetary_objects and stellar_system.child2.stellar_type not in stellar_types_planetary_objects:
                donor = stellar_system.child1
                dm = mass_lost_due_to_evaporation_in_binary(stellar_system, dt, donor, stellar_system.child2, self) 
                donor.previous_mass = donor.mass
                donor_in_stellar_code = donor.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
                donor_in_stellar_code.change_mass(-1*dm, 0.|units.yr)
                print("mass lost child1:", dm)
            elif stellar_system.child2.stellar_type in stellar_types_planetary_objects and stellar_system.child1.stellar_type not in stellar_types_planetary_objects:
                donor = stellar_system.child2
                dm = mass_lost_due_to_evaporation_in_binary(stellar_system, dt, donor, stellar_system.child1, self) 
                donor.previous_mass = donor.mass
                donor_in_stellar_code = donor.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
                donor_in_stellar_code.change_mass(-1*dm, 0.|units.yr)
                print("mass lost child2:", dm)
        elif self.is_triple(stellar_system):
            if stellar_system.child1.is_star: #child1 is the tertiary
                # effect on tertiary
                if stellar_system.child1.stellar_type in stellar_types_planetary_objects: #evaporation if tertiary is a planet
                    donor = stellar_system.child1
                    dm = mass_lost_due_to_evaporation_tertiary(stellar_system, dt, donor, stellar_system.child2, self) 
                    donor.previous_mass = donor.mass
                    donor_in_stellar_code = donor.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
                    donor_in_stellar_code.change_mass(-1*dm, 0.|units.yr)
                    print("mass lost child1:", dm)
                # check the inner binary: for now we only take into account evaporation due the closest stellar companion 
                self.planetary_mass_evaporation(dt, stellar_system.child2) 
            elif self.child2.is_star: #child2 is the tertiary
                # effect on tertiary
                if self.child2.stellar_type in stellar_types_planetary_objects: #evaporation if tertiary is a planet
                    donor = stellar_system.child2
                    dm = mass_lost_due_to_evaporation_tertiary(stellar_system, dt, donor, stellar_system.child1, self) 
                    donor.previous_mass = donor.mass
                    donor_in_stellar_code = donor.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
                    donor_in_stellar_code.change_mass(-1*dm, 0.|units.yr)
                    print("mass lost child2:", dm)
                # check the inner binary: for now we only take into account evaporation due the closest stellar companion 
                self.planetary_mass_evaporation(dt, stellar_system.child1) 
            else:             
                sys.exit('planetary_mass_evaporation: structure stellar system unknown')
    
        else:         
            sys.exit('planetary_mass_evaporation: structure stellar system unknown')
                 
    
    
    def resolve_stellar_interaction(self, stellar_system = None):
    # the most inner binary is calculated first, and then move outwards

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print('single stellar evolution')
            sys.exit('for now no single stellar evolution - exiting program')
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print('\n perform stellar interaction: binary')
#            stellar_system = perform_stellar_interaction(stellar_system, self)
            stopping_condition = perform_stellar_interaction(stellar_system, self)
            return stopping_condition #stellar interaction
        else:
            if REPORT_TRIPLE_EVOLUTION:
                print('\n perform stellar interaction')
      
            stopping_condition = perform_stellar_interaction(stellar_system, self)            
            if not stopping_condition: #stellar interaction
                return False                                     

            if stopping_condition>0:
                #in case of outer mass transfer, skip inner evolution, is taken care of in outer mass transfer
                if not stellar_system.child1.is_star: #child1 is a multiple
                    stopping_condition = self.resolve_stellar_interaction(stellar_system.child1)  
                    if not stopping_condition: #stellar interaction
                        return False                                     
                elif not stellar_system.child2.is_star: #child2 is a multiple
                    stopping_condition = self.resolve_stellar_interaction(stellar_system.child2)        
                    if not stopping_condition: #stellar interaction
                        return False                                     
                else:
                    sys.exit('resolve_stellar_interaction: structure stellar system unknown, both children are binaries')

            return True                  
            
                                                                
    def determine_mass_transfer_timescale(self, stellar_system = None):

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print('single stellar evolution')
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print('\n determine_mass_transfer_timescale: binary - double star')
            mass_transfer_stability(stellar_system, self)
            if REPORT_TRIPLE_EVOLUTION:
                print('mt rate double star:', stellar_system.mass_transfer_rate)
        else:
            self.determine_mass_transfer_timescale(stellar_system.child1)        
            self.determine_mass_transfer_timescale(stellar_system.child2)        

            if REPORT_TRIPLE_EVOLUTION:
                print('\n determine_mass_transfer_timescale: binary')
            mass_transfer_stability(stellar_system, self)            
            if REPORT_TRIPLE_EVOLUTION:
                print('mt rate binary:', stellar_system.mass_transfer_rate)




    def solve_for_partial_time_step_stable_mass_transfer(self):
        full_dt = self.triple.time - self.previous_time
        self.secular_code.evolve_model(self.previous_time + full_dt * self.triple.child2.part_dt_mt)

        self.channel_from_secular.copy()
        self.check_RLOF() 
        if self.has_donor():
            print('After partial timestep the system should be detached...')
            print(self.has_donor())
            print(self.triple.child1.mass, self.triple.child2.child1.mass, self.triple.child2.child2.mass)
            return False

        self.secular_code.model_time = self.triple.time # not redundant!
        self.triple.child2.part_dt_mt = 1.                    
        return True                      
                    
    
    def safety_check_time_step(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        successfull_dr = True
        
        if stellar_system.is_star:
            dm = (stellar_system.previous_mass - stellar_system.mass)/stellar_system.mass
#            if REPORT_DEBUG:
#                print('wind mass loss rate:', stellar_system.wind_mass_loss_rate, stellar_system.mass,stellar_system.core_mass,)
#                print('relative wind mass losses:', dm)
                        
            if (abs(dm) > error_dm) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type ) and not (stellar_system.stellar_type in stellar_types_SN_remnants):
                successfull_dr = False 
                print('should this statement be here? does this give problems during mt?')

            if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:
                dr = (stellar_system.radius - stellar_system.previous_radius)/stellar_system.radius
    
#                if REPORT_DEBUG:    
#                    print('change in radius over time:',  stellar_system.time_derivative_of_radius,)
#                    print('relative change in radius:', dr)
                                     

                if (abs(dr) > error_dr) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_dr) and  not stellar_system.parent.bin_type == 'stable_mass_transfer':
                    #during stable mass transfer we reset the system to synchronization, 
                    #so tides and spin change due to moment of inertia change are not important, 
                    #therefore the radius change is not important

                    successfull_dr = False                    
            return successfull_dr, 1, stellar_system

        else:
            successfull_dr1, nr_dr1, star_dr1 = self.safety_check_time_step(stellar_system.child1)        
            successfull_dr2, nr_dr2, star_dr2 = self.safety_check_time_step(stellar_system.child2)
            if successfull_dr1==False and successfull_dr2==False:
                if nr_dr1 > 1:
                    return False, 3, np.array([star_dr1[0], star_dr1[1], star_dr2])
                elif nr_dr2 > 1:
                    return False, 3, np.array([star_dr1, star_dr2[0], star_dr2[1]])
                else:
                    return False, 2, np.array([star_dr1, star_dr2])
            elif successfull_dr1==False: 
                return False, nr_dr1, star_dr1
            elif successfull_dr2==False: 
                return False, nr_dr2, star_dr2
            else:
                return True, 0, 0    
            
          
    def recall_memory_one_step_stellar(self, nr_unsuccessfull, star_unsuccessfull):          
        if REPORT_DEBUG:
            print('recall_memory_one_step_stellar')

        self.update_time_derivative_of_radius()
        dt = self.triple.time-self.previous_time
        #crashes for M=1.4Msun: MS hook at stellar type change	
        #dt_new = max(minimum_time_step, min(self.determine_time_step(), 0.9*dt))
        dt_new = max(minimum_time_step, 0.5*dt)
        self.stellar_code.particles.recall_memory_one_step()
        self.triple.time += (dt_new - dt)                     
        
        self.stellar_code.evolve_model(self.triple.time)
        
        self.copy_from_stellar()
        self.update_stellar_parameters()          
                
        if nr_unsuccessfull > 1:
            for i_rec in range(nr_unsuccessfull):
                triple_number = self.triple.number
                if dt_new <= minimum_time_step:
                    triple_number += 0.5
                dr = (star_unsuccessfull[i_rec].radius - star_unsuccessfull[i_rec].previous_radius)/star_unsuccessfull[i_rec].radius        
        else:
            dr = (star_unsuccessfull.radius - star_unsuccessfull.previous_radius)/star_unsuccessfull.radius        
            triple_number = self.triple.number
            if dt_new <= minimum_time_step:
                triple_number += 0.5

        if dt_new <= minimum_time_step:
            return True, 0, 0    

        return self.safety_check_time_step()                        
        
    #when the stellar code finds RLOF
    def rewind_to_begin_of_rlof_stellar(self, dt):          
        if REPORT_DEBUG:
            print('rewind to begin of rlof stellar')
        dt_new = dt
        dt_min = max(minimum_time_step, self.determine_time_step_stable_mt())
        while self.has_donor() == True and dt_new > dt_min:
            dt_old = self.triple.time - self.previous_time
            dt_new = max(minimum_time_step, 0.5*dt_old)

            self.stellar_code.particles.recall_memory_one_step()
            self.triple.time += (dt_new - dt_old)   
            if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
                print('\n\ntime:', self.triple.time, self.has_donor())
            self.stellar_code.evolve_model(self.triple.time)
            self.copy_from_stellar()
            self.update_stellar_parameters()              
            self.check_RLOF()                    
                    

    #when the secular code finds RLOF
    def rewind_to_begin_of_rlof_secular(self):          
        if REPORT_DEBUG:
            print('rewind to begin of rlof secular')

        #find beginning of RLOF
        self.fixed_timestep = self.secular_code.model_time - self.previous_time
        self.triple.time = self.previous_time
        self.secular_code.model_time = self.previous_time      
                         
        if self.fixed_timestep < 0.|units.yr:
            print(self.triple.time, self.secular_code.model_time, self.fixed_timestep)      
            sys.exit('fixed_timestep < 0: should not be possible')           
            
        #rewind system
        self.stellar_code.particles.recall_memory_one_step()
        self.refresh_memory() 
                                           
            
    def check_stopping_conditions_stellar_interaction(self):             
        if self.stop_at_merger and self.has_merger():
            if REPORT_TRIPLE_EVOLUTION:
                print('Merger at time = ',self.triple.time )                             
            return False    
        if self.stop_at_disintegrated and self.has_disintegrated():
            if REPORT_TRIPLE_EVOLUTION:
                print('Disintegration of system at time = ',self.triple.time )             
            return False
        return True            
            
            
#    def check_stopping_conditions_stellar(self):  
#        if self.stop_at_outer_mass_transfer and self.has_tertiary_donor():
#            if self.secular_code.model_time < self.triple.time:
#                self.triple.time = self.secular_code.model_time
#
#            if REPORT_TRIPLE_EVOLUTION:
#                print('Mass transfer in outer binary of triple at time = ',self.triple.time)
#            self.triple.bin_type = bin_type['rlof']
#            
#            self.determine_mass_transfer_timescale() # to set the stability #obsolete?
#            return False                                   
#
#        elif self.stop_at_mass_transfer and self.has_donor():
#            if self.secular_code.model_time < self.triple.time:
#                self.triple.time = self.secular_code.model_time
#        
#            if REPORT_TRIPLE_EVOLUTION:
#                print('Mass transfer at time = ',self.triple.time)
#            if self.is_binary(self.triple.child2):
#                self.triple.child2.bin_type = bin_type['rlof'] 
#            elif self.is_binary(self.triple.child1):
#                self.triple.child1.bin_type = bin_type['rlof']    
#            else:
#                print('currently not implemented')
#                exit(-1)    
#                
#            self.determine_mass_transfer_timescale() # to set the stability #obsolete? 
#            return False
#            
#        return True



    def check_stopping_conditions_stellar(self, stellar_system = None):  
        #before check stopping conditions_stellar, always make sure the stability labels are up to date
        #by running self.determine_mass_transfer_timescale()

        if stellar_system == None:
            stellar_system = self.triple

        if self.is_triple(stellar_system):
            if self.stop_at_outer_mass_transfer and self.has_tertiary_donor():
                if self.secular_code.model_time < self.triple.time:
                    self.triple.time = self.secular_code.model_time
                    
                if REPORT_TRIPLE_EVOLUTION:
                    print('Mass transfer in outer binary of triple at time = ',self.triple.time)

                if stellar_system.is_mt_stable:
                    stellar_system.bin_type = bin_type['stable_mass_transfer']
                else:
                    stellar_system.bin_type = bin_type['common_envelope']

                return False                                   
            else:
                if not self.check_stopping_conditions_stellar(stellar_system.child1):
                    return False
                if not self.check_stopping_conditions_stellar(stellar_system.child2):
                    return False
                return True
        elif stellar_system.is_star:
            return True
        elif self.is_binary(stellar_system): 
            if (self.has_donor(stellar_system) or self.has_OLOF_donor(stellar_system)) and (self.stop_at_mass_transfer or 
                (self.stop_at_stable_mass_transfer and stellar_system.is_mt_stable) or
                (self.stop_at_unstable_mass_transfer and not stellar_system.is_mt_stable) or
                (self.stop_at_eccentric_stable_mass_transfer and stellar_system.is_mt_stable and stellar_system.eccentricity > minimum_eccentricity*5.) or
                (self.stop_at_eccentric_unstable_mass_transfer and not stellar_system.is_mt_stable and stellar_system.eccentricity > minimum_eccentricity*5.)):

                if REPORT_TRIPLE_EVOLUTION:
                    print('Mass transfer in inner binary at time = ',self.triple.time)
                    print(self.stop_at_mass_transfer,self.stop_at_stable_mass_transfer, self.stop_at_unstable_mass_transfer, self.stop_at_eccentric_stable_mass_transfer, self.stop_at_eccentric_unstable_mass_transfer, stellar_system.is_mt_stable)                       

                if stellar_system.is_mt_stable:
                    stellar_system.bin_type = bin_type['stable_mass_transfer']
                else:
                    stellar_system.bin_type = bin_type['common_envelope']

                return False    
            else:
                return True
            
        return True
        
            
    #when these processes are implemented, also the bintypes need to be set in those functions
    def check_stopping_conditions(self):
        if self.check_stopping_conditions_stellar()==False:
            return False

        if self.stop_at_dynamical_instability and self.triple.dynamical_instability == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time  
            self.set_bintype_to_dynamical_instability()            
            if REPORT_TRIPLE_EVOLUTION:
                print("Dynamical instability at time = ",self.triple.time)
            return False
        if self.stop_at_inner_collision and self.triple.inner_collision == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.child2.bin_type = bin_type['collision']
            if REPORT_TRIPLE_EVOLUTION:
                print("Inner collision at time = ",self.triple.time)
            return False
        if self.stop_at_outer_collision and self.triple.outer_collision == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.bin_type = bin_type['collision']
            if REPORT_TRIPLE_EVOLUTION:
                print("Outer collision at time = ",self.triple.time)
            return False
        if self.stop_at_semisecular_regime and self.triple.semisecular_regime == True:
            if self.secular_code.model_time < self.triple.time:
                self.triple.time = self.secular_code.model_time
            self.triple.bin_type = bin_type['semisecular']
            if REPORT_TRIPLE_EVOLUTION:
                print("Semisecular regime at time = ",self.triple.time)
            return False
        if self.triple.time - self.secular_code.model_time > numerical_error|units.Myr:
            print('triple time > sec time: should not be possible', self.triple.time, self.secular_code.model_time)
            print(self.has_donor(), self.secular_code.triples[0].error_flag_secular, self.secular_code.triples[0].dynamical_instability)
            return False
            
        return True


    def check_spin_angular_frequency(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.spin_angular_frequency < 0.0|1./units.Myr: 
                if REPORT_DEBUG:
                  print("Error in spin at time = ",self.triple.time)
               

                return False
            return True                
        else:
            if self.check_spin_angular_frequency(stellar_system.child1) == False or self.check_spin_angular_frequency(stellar_system.child2)==False:
                return False
        return True

    def check_error_flag_secular(self):    
        if self.triple.error_flag_secular < 0:
            if REPORT_TRIPLE_EVOLUTION:
                print("Error in secular code at time = ",self.triple.time)
                print(self.triple.error_flag_secular)
            return False
        return True
        
    def calculate_maximum_change_eccentricity(self):
        if self.triple.child2.delta_e_in > self.triple.child2.max_delta_e_in:
            self.triple.child2.max_delta_e_in = self.triple.child2.delta_e_in
            
            
#-------------------------    

    def evolve_model(self, tend):
        self.tend = tend
        CPU_start_time = time.time()
    
        if REPORT_DEBUG or MAKE_PLOTS:
            # for plotting data
            times_array = quantities.AdaptingVectorQuantity() 
            a_in_array = quantities.AdaptingVectorQuantity()
            e_in_array = []
            g_in_array = []
            o_in_array = []
            a_out_array = quantities.AdaptingVectorQuantity()
            e_out_array = []
            g_out_array = []
            o_out_array = []
            i_relative_array = []
            m1_array = quantities.AdaptingVectorQuantity()
            m2_array = quantities.AdaptingVectorQuantity()
            m3_array = quantities.AdaptingVectorQuantity()
            spin1_array = []
            spin2_array = []
            spin3_array = []
            r1_array = []
            r2_array = []
            r3_array = []
            moi1_array = []
            moi2_array = []
            moi3_array = []
            RL1_array = []
            RL2_array = []
            RL3_array = []
            delta_e_in_array = []
        
            times_array.append(self.triple.time)
            e_in_array.append(self.triple.child2.eccentricity)
            a_in_array.append(self.triple.child2.semimajor_axis)
            g_in_array.append(self.triple.child2.argument_of_pericenter) 
            o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
            e_out_array.append(self.triple.eccentricity)
            a_out_array.append(self.triple.semimajor_axis)
            g_out_array.append(self.triple.argument_of_pericenter) 
            o_out_array.append(self.triple.longitude_of_ascending_node)               
            i_relative_array.append(self.triple.relative_inclination)        
            m1_array.append(self.triple.child2.child1.mass)
            m2_array.append(self.triple.child2.child2.mass)
            m3_array.append(self.triple.child1.mass)
            spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
            spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
            spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
            r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
            r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
            r3_array.append(self.triple.child1.radius.value_in(units.RSun))
            moi1_array.append(self.triple.child2.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            moi2_array.append(self.triple.child2.child2.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            moi3_array.append(self.triple.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
            RL1, RL2, RL3 = self.secular_code.give_roche_radii(self.triple)
            RL1_array.append(RL1.value_in(units.RSun))
            RL2_array.append(RL2.value_in(units.RSun))
            RL3_array.append(RL3.value_in(units.RSun))
            delta_e_in_array.append(self.triple.child2.delta_e_in)

        include_inner_RLOF_term_initial = self.secular_code.parameters.include_inner_RLOF_terms
        include_outer_RLOF_term_initial = self.secular_code.parameters.include_outer_RLOF_terms
        
        if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
            print('kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.tend)

        self.determine_mass_transfer_timescale()
        self.save_snapshot()          
        while self.triple.time<self.tend: 
            if REPORT_TRIPLE_EVOLUTION or REPORT_DEBUG:
                print('\n\n kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.octupole_parameter())

            #resetting and testing some parameters     
                         
	        #if the maximum CPU time has been exceeded for the system, stop the evolution and continue with the next system
            CPU_end_time = time.time()
            self.triple.CPU_time = CPU_end_time - CPU_start_time
            if (self.stop_at_CPU_time == True) and float(CPU_end_time - CPU_start_time) > self.max_CPU_time:
                print('stopping conditions maximum CPU time')
                print("CPU time: ", self.triple.CPU_time)
                break
                
            #turning RLOF terms back on if they were on in the first place
            #they might have been turned off for contact systems with is_mt_stable in perform_mass_equalisation_for_contact
            if include_inner_RLOF_term_initial == True:
                self.secular_code.parameters.include_inner_RLOF_terms = True 
            if include_outer_RLOF_term_initial == True:
                self.secular_code.parameters.include_outer_RLOF_terms = True 
                

            #setting time parameters
            if self.has_stellar_type_changed() or self.has_kozai_type_changed():
                self.save_snapshot()  

            dt = self.determine_time_step()  
            if not no_stellar_evolution: 
                self.update_previous_stellar_parameters()
                if not USE_MESA_AS_STELLAR_CODE: # for now, refresh_memory and recall_memory_one_step not available in MESA interface
                    self.stellar_code.particles.refresh_memory()
            self.triple.time += dt   
            self.previous_dt = dt         
            if REPORT_DEBUG or REPORT_DT:
                print('\t time:', self.triple.time, self.previous_time, self.previous_dt)
                
            if self.triple.time <self.tinit:
                self.instantaneous_evolution = True
                         
                            
            #do stellar evolution 
            if not no_stellar_evolution: 
                if REPORT_DEBUG:
                    print('Stellar evolution')

                if self.include_CHE:#only needed when including CHE
                    self.channel_to_stellar.copy_attributes(['rotation_period'])
                
                self.planetary_mass_evaporation(dt)                    
                self.stellar_code.evolve_model(self.triple.time)
                self.copy_from_stellar()
                self.update_stellar_parameters() 
                                
                successfull_step, nr_unsuccessfull, star_unsuccessfull = self.safety_check_time_step() 
                while successfull_step == False:
                    if not USE_MESA_AS_STELLAR_CODE: # for now, refresh_memory and recall_memory_one_step not available in MESA interface
                        successfull_step, nr_unsuccessfull, star_unsuccessfull = self.recall_memory_one_step_stellar(nr_unsuccessfull, star_unsuccessfull)
                    else:
                        successfull_step = True
                # if SN has taken place
                if self.has_stellar_type_changed_into_SN_remnant():
                    if REPORT_TRIPLE_EVOLUTION:
                        print("Supernova at time = ", self.triple.time, dt)
                    if self.adjust_system_after_supernova_kick()==False:
                        break  
                    self.instantaneous_evolution = True #skip secular evolution

                self.check_RLOF()                                       
                self.check_OLOF()                                       

                if (self.has_donor() or self.has_OLOF_donor()) and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached' and dt > minimum_time_step:
#                    self.rewind_to_begin_of_rlof_stellar(dt) 
#                    print('RLOF:', self.triple.child2.child1.is_donor, self.triple.bin_type , self.triple.child2.bin_type )

                    if not USE_MESA_AS_STELLAR_CODE: # for now, refresh_memory and recall_memory_one_step not available in MESA interface
                        self.stellar_code.particles.recall_memory_one_step()
                    self.copy_from_stellar()
                    self.update_stellar_parameters()   
                    if self.include_CHE:#only needed when including CHE
                        self.channel_to_stellar.copy_attributes(['rotation_period']) #for CHE                           
                    self.refresh_memory()                     
                    #note that 'previous' parameters cannot be reset 
                    #resetting is_donor in determine_time_step                                    
                    continue

                                         

            #needed for nucleair timescale 
            self.update_time_derivative_of_radius() 
          
            #do stellar interaction
            if REPORT_DEBUG:
                print('Stellar interaction')

            self.determine_mass_transfer_timescale()
            if self.check_stopping_conditions_stellar()==False:
                print('stopping conditions stellar')
                break  
            if not self.resolve_stellar_interaction():
                print('stopping conditions stellar interaction')
                break            
            self.update_time_derivative_of_radius() # includes radius change from wind and ce, not stable mt
            self.update_time_derivative_of_moment_of_inertia() # includes mass and radius change from wind and mass transfer            
            if self.check_stopping_conditions_stellar_interaction()==False: 
                print('stopping conditions stellar interaction 2')
                break
            
            
                                                           
            # do secular evolution
            if self.instantaneous_evolution == False:
                if REPORT_DEBUG:
                    print('Secular evolution')

                #needed for refreshing memory in case secular finds RLOF    
                previous_semimajor_axis_in = self.triple.child2.semimajor_axis
                previous_eccentricity_in = self.triple.child2.eccentricity
                previous_argument_of_pericenter_in = self.triple.child2.argument_of_pericenter
                previous_longitude_of_ascending_node_in = self.triple.child2.longitude_of_ascending_node
                previous_spin1 = self.triple.child2.child1.spin_angular_frequency
                previous_spin2 = self.triple.child2.child2.spin_angular_frequency
    

                self.channel_to_secular.copy()
                
                # if mass transfer should only take place for a fraction of the timestep
                # e.g. at the end of mass transfer when the envelope is thin
                if self.triple.child2.part_dt_mt < 1: # inner binary, see function determine_partial_time_step_stable_mass_transfer
                    print('skipping')
                    if self.solve_for_partial_time_step_stable_mass_transfer() == False:
                        break
                else:
                    self.secular_code.evolve_model(self.triple.time)
                    
                if REPORT_DEBUG:
                    print('Secular evolution finished')
                
                # to differentiate between semi-detached and contact 
                self.check_RLOF()
                self.check_OLOF()
		
                if self.triple.time - self.secular_code.model_time < -1*numerical_error|units.Myr and self.secular_code.triples[0].error_flag_secular >= 0:
                    print('triple time < sec time: should not be possible', self.triple.time, self.secular_code.model_time)
                    print(self.has_donor(), self.secular_code.triples[0].error_flag_secular)
                    break
                elif (self.has_donor() or self.has_OLOF_donor()) and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached' and self.secular_code.model_time < self.triple.time-max(minimum_time_step, 0.01*self.determine_time_step_stable_mt()):
                    #factor 0.01 times time_step_stable_mt as used mass_transfer_timescale from previous timestep
                    
                    self.determine_mass_transfer_timescale()            
                    if self.check_stopping_conditions_stellar()==False:
                        print('stopping conditions stellar 2')                    
                        break                   

                    if not USE_MESA_AS_STELLAR_CODE: # for now, refresh_memory and recall_memory_one_step not available in MESA interface
                        self.rewind_to_begin_of_rlof_secular()
                    self.triple.child2.semimajor_axis = previous_semimajor_axis_in
                    self.triple.child2.eccentricity = previous_eccentricity_in
                    self.triple.child2.previous_argument_of_pericenter_in = previous_argument_of_pericenter_in
                    self.triple.child2.previous_longitude_of_ascending_node_in = previous_longitude_of_ascending_node_in
                    self.triple.child2.child1.spin_angular_frequency = previous_spin1
                    self.triple.child2.child2.spin_angular_frequency = previous_spin2
                    continue
                elif self.has_donor() and self.triple.bin_type == 'detached' and self.triple.child2.bin_type == 'detached':
                    #time difference secular code & tres small enough to continue mass transfer 
                    self.secular_code.model_time = self.triple.time
                                  
                self.channel_from_secular.copy() 



                if self.check_error_flag_secular()==False:
                    break
                if self.check_spin_angular_frequency()==False:
                    break
                self.determine_mass_transfer_timescale()
                if self.check_stopping_conditions()==False:
                    print('stopping conditions')                
                    break                
                if not self.stop_at_inner_collision and self.triple.inner_collision == True:
                    perform_inner_collision(self)
                                        
            else:
                if REPORT_DEBUG:
                    print('skip secular')
                self.secular_code.model_time = self.triple.time
                self.instantaneous_evolution = False
                
            self.calculate_maximum_change_eccentricity()


            if REPORT_DEBUG or MAKE_PLOTS:
                # for plotting data
                times_array.append(self.triple.time)
                e_in_array.append(self.triple.child2.eccentricity)
                a_in_array.append(self.triple.child2.semimajor_axis)
                g_in_array.append(self.triple.child2.argument_of_pericenter) 
                o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
                e_out_array.append(self.triple.eccentricity)
                a_out_array.append(self.triple.semimajor_axis)
                g_out_array.append(self.triple.argument_of_pericenter) 
                o_out_array.append(self.triple.longitude_of_ascending_node)               
                i_relative_array.append(self.triple.relative_inclination)        
                m1_array.append(self.triple.child2.child1.mass)
                m2_array.append(self.triple.child2.child2.mass)
                m3_array.append(self.triple.child1.mass)
                spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
                spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
                spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
                r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
                r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
                r3_array.append(self.triple.child1.radius.value_in(units.RSun))
                moi1_array.append(self.triple.child2.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
                moi2_array.append(self.triple.child2.child2.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
                moi3_array.append(self.triple.child1.moment_of_inertia_of_star.value_in(units.RSun**2*units.MSun))
                RL1, RL2, RL3 = self.secular_code.give_roche_radii(self.triple)
                RL1_array.append(RL1.value_in(units.RSun))
                RL2_array.append(RL2.value_in(units.RSun))
                RL3_array.append(RL3.value_in(units.RSun))
                delta_e_in_array.append(self.triple.child2.delta_e_in)
                                    
        self.save_snapshot()        
            
            
        if (REPORT_DEBUG or MAKE_PLOTS) and self.triple.time >= self.tinit:
            # for plotting data
            e_in_array = np.array(e_in_array)
            g_in_array = np.array(g_in_array)
            o_in_array = np.array(o_in_array)
            e_out_array = np.array(e_out_array)
            g_out_array = np.array(g_out_array)
            o_out_array = np.array(o_out_array)
            i_relative_array = np.array(i_relative_array)
            spin1_array = np.array(spin1_array)
            spin2_array = np.array(spin2_array)
            spin3_array = np.array(spin3_array)
            r1_array = np.array(r1_array)
            r2_array = np.array(r2_array)
            r3_array = np.array(r3_array)
            moi1_array = np.array(moi1_array)
            moi2_array = np.array(moi2_array)
            moi3_array = np.array(moi3_array)
            RL1, RL2, RL3 = self.secular_code.give_roche_radii(self.triple)
            RL1_array = np.array(RL1_array)
            RL2_array = np.array(RL2_array)
            RL3_array = np.array(RL3_array)
    
            self.plot_data = plot_data_container()
            self.plot_data.times_array = times_array
            self.plot_data.a_in_array = a_in_array
            self.plot_data.e_in_array = e_in_array
            self.plot_data.g_in_array = g_in_array
            self.plot_data.o_in_array = o_in_array        
            self.plot_data.a_out_array = a_out_array
            self.plot_data.e_out_array = e_out_array
            self.plot_data.g_out_array = g_out_array
            self.plot_data.o_out_array = o_out_array        
            self.plot_data.i_relative_array = i_relative_array
            self.plot_data.m1_array = m1_array
            self.plot_data.m2_array = m2_array
            self.plot_data.m3_array = m3_array
            self.plot_data.spin1_array = spin1_array
            self.plot_data.spin2_array = spin2_array
            self.plot_data.spin3_array = spin3_array
            self.plot_data.r1_array = r1_array
            self.plot_data.r2_array = r2_array
            self.plot_data.r3_array = r3_array
            self.plot_data.moi1_array = moi1_array
            self.plot_data.moi2_array = moi2_array
            self.plot_data.moi3_array = moi3_array
            self.plot_data.RL1_array = RL1_array
            self.plot_data.RL2_array = RL2_array
            self.plot_data.RL3_array = RL3_array
            self.plot_data.delta_e_in_array = delta_e_in_array
            
        
    #-------
