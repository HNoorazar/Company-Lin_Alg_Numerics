"""
Company. 
Very First code.
"""

import numpy as np
"""
1D discrete derivatives.
"""
##########################
########################## Forward difference derivative
##########################
def 1D_forward_diff_deriv(grids, func):
    # "grids" is a numpy array of grid points. It is a column vector.
    # "func" is the function we want to 
    # find use its derivatives at grid points given by "grids".
    no_derivs = np.size(grids, 1) - 1
    for_derivatives = np.zeros((1, no_derivs))

    for grid_counter in range(0, no_derivs):
        difference  = func(grids[1, grid_counter+1]) - func(grids[1,grid_counter])
        step_size = grids[1, grid_counter+1] - grids[1, grid_counter]
        for_derivatives[0, grid_counter] = difference / step_size
    return for_derivatives

##########################
########################## Backward difference derivative
##########################
def 1D_backward_diff_deriv(grids, func):
    # "grids" is a numpy array of grid points. It is a column vector.
    # "func" is the function we want to 
    # find use its derivatives at grid points given by "grids".
    no_derivs = np.size(grids, 1) - 1
    back_derivatives = np.zeros((1, no_derivs))

    for grid_counter in range(1, no_derivs+1):
        difference  = func(grids[1, grid_counter]) - func(grids[1, grid_counter-1])
        step_size = grids[1, grid_counter] - grids[1, grid_counter-1]
        back_derivatives[0, grid_counter] = difference / step_size
    return back_derivatives[0,1:]
        
        
##########################        
########################## centered difference derivative
##########################
