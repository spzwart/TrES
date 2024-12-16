##########################################################################################
# A script to set a few MESA controls for using TRES with MESA.                          #
#                                                                                        #
# When used with its default setting, MESA does not account for several processes (mass  #
# loss, overshoot, semi-convection, ...)                                                 #
#                                                                                        #
# Here we provided a set of MESA controls appropriate for massive stars.                 #
#                                                                                        #
# The choice of input physics is made in order to be the closest to the SeBa defaults.   #   
# More details can be found in Sciarini et al. 2025 (in prep.)                           #
#                                                                                        #
# The list of MESA controls can be found in                                              #
# https://github.com/MESAHub/mesa/blob/r15140/star/defaults/controls.defaults            #
##########################################################################################


def setup_mesa_stars(stellar_code):   
    
    for i in range(len(stellar_code.particles)):
          
        # Opacity
        stellar_code.particles[i].set_kap('use_Type2_opacities',True)
        stellar_code.particles[i].set_kap('Zbase', 0.02)

            
        # End phase of simulation (here end of central C burning)
        stellar_code.particles[i].set_control('xa_central_lower_limit_species(1)','c12')
        stellar_code.particles[i].set_control('xa_central_lower_limit(1)',1e-12)
            
        # Winds
        # The 'Dutch' prescription combines Vink, de Jager and Nugis & Lamers
        stellar_code.particles[i].set_control('hot_wind_scheme','Dutch')
        
        # RGB winds
        stellar_code.particles[i].set_control('cool_wind_rgb_scheme','Dutch')
      
        # AGB winds    
        stellar_code.particles[i].set_control('cool_wind_agb_scheme','Dutch')
        
        # Default factor in SeBa
        stellar_code.particles[i].set_control('Dutch_scaling_factor',0.333333)

            
            
        # Overshoot  
        # Here a step overshoot is chosen. This is not the only option
        # 'exponential' overshoot is also available in MESA

        # Overshoot during central H burning 
        stellar_code.particles[i].set_control("overshoot_scheme(1)", 'step')
        stellar_code.particles[i].set_control("overshoot_zone_type(1)", 'burn_H')
        stellar_code.particles[i].set_control("overshoot_zone_loc(1)", 'core')
        stellar_code.particles[i].set_control("overshoot_bdy_loc(1)", 'top')


        # alpha_ov calibrated so that maximum radius in the MS matches that of SeBa 
        # for the 50Msun star, which corresponds to the upper limit of Pols+98 grid.
        # This value of a_ov also allows a pretty good match for all the models with M <= 50Msun.        
        stellar_code.particles[i].set_control("overshoot_f(1)",0.31)# alpha_ov/f_ov for step/exponential overshot
        stellar_code.particles[i].set_control("overshoot_f0(1)", 0.03)# insert a value strictly smaller than 10% of alpha_ov 
        
        # Overshoot during central He burning
        stellar_code.particles[i].set_control("overshoot_scheme(2)", 'step')
        stellar_code.particles[i].set_control("overshoot_zone_type(2)", 'burn_He')
        stellar_code.particles[i].set_control("overshoot_zone_loc(2)", 'core')
        stellar_code.particles[i].set_control("overshoot_bdy_loc(2)", 'top')

        stellar_code.particles[i].set_control("overshoot_f(2)",0.31)# alpha_ov/f_ov for step/exponential overshot
        stellar_code.particles[i].set_control("overshoot_f0(2)", 0.03)# insert a value strictly smaller than 10% of alpha_ov 

            
        # Semiconvection
        stellar_code.particles[i].set_control('alpha_semiconvection',0)
            
        # Semiconvection only works with Ledoux (Schwarzschild by default)
        stellar_code.particles[i].set_control('use_Ledoux_criterion', False)
            
            
        # Mixing length
        stellar_code.particles[i].set_control('mixing_length_alpha', 2.0)
            
        # Nuclear network
        stellar_code.particles[i].set_control('default_net_name','mesa_49.net')
        	
		#Time and space resolution
        stellar_code.particles[i].set_control('mesh_delta_coeff', 0.5)
        stellar_code.particles[i].set_control('time_delta_coeff', 0.4)
        stellar_code.particles[i].set_control('dH_hard_limit', 0.001)
        stellar_code.particles[i].set_control('varcontrol_target', 0.01)
            
        
        # MLT ++ ('okay_to_reduce_gradT_excess' = True)
        # We recommend to use MLT++ for massive stars (>= 10 Msun)
        # For very massive stars (>= 50 Msun), the stars usually will contract a bit at the beginning of the simulation
        # during model convergence when MLT++ is used, which may lead to TRES imposing a timestep reduction (because of too high radius_change).
        # A possible solution may be to prevent the timestep reduction in this specific case.
        stellar_code.particles[i].set_control('okay_to_reduce_gradT_excess', True)
        
        # Convergence
        stellar_code.particles[i].set_control('Pextra_factor', 2.0)

