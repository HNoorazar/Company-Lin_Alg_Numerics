import numpy as np
import machinery.core as my_core
import machinery.config as my_cfg
import machinery.surface as my_surf
import machinery.boundary as my_bnd
import machinery.definitions as my_def


config = my_cfg.staticParameters()
config.readFromFile('staticParameters.cfg')
config.printOut()

#  Read initial values from file "infile" 
# I have no idea how read_bin fucntion is supposed to work!
Lines 63 to 81 of main.c is skipped here.



state = ?

state = my_bnd.set_boun_cond(state, config)
state = my_bnd.set_specific_conditions(state, config)


"""     --------------------
         t i m e    l o o p 
        --------------------
"""

for time in "What the hell is that now?"
    config.delta_t = compute_time_step(config,
                                       x_grid_vel, y_grid_vel, 
                                       Reynolds)