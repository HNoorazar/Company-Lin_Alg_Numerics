import numpy as np
import machinery.core as my_core
import machinery.config as my_cfg
import machinery.surface as my_surf
import machinery.boundary as my_bnd
import machinery.definitions as my_def
import machinery.state as my_stat

#
# process command line
# This just will read configuration/input file 
# which defines the problem and parameters associated with it.
#
cmdline =  my_cfg.CmdLineArguments()

prob_parameters = 'wave.cfg'

config = my_cfg.staticParameters()
config.readFromFile(prob_parameters)
ibound = 0

state = my_stat.GridState(config.init_x_vel_scalar,
                          config.init_y_vel_scalar,
                          config.init_press_scalar,
                          ibound=0,
                          current_time=0 )
state = my_bnd.set_boun_cond(state, config)
state = my_bnd.set_specific_conditions(state, config)

""" --------------------
     t i m e    l o o p 
    --------------------
"""

t = 0; cycle=0;
while (t < config.final_t):
    config.delta_t = compute_time_step(config,
                                       x_grid_vel, y_grid_vel, 
                                       Reynolds)
    if ((problem != "drop") or (problem != "dam") or 
    	(problem != "molding") or (problem != "wave")):
        
        

    t += config.delta_t
    cycle += 1




