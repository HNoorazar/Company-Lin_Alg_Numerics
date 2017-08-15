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
    # "grids" is a numpy array of grid points. It is a column vector of size n-by-1.
    # "func" is the function we want to 
    # find use its derivatives at grid points given by "grids".
    no_derivs = np.size(grids, 1) - 1
    for_derivatives = np.zeros((1, no_derivs))
    step_sizes = grids[0, 0:-1] - grids[0, 1:]
    for deriv_counter in range(0, no_derivs):
        difference  = func(grids[1, deriv_counter+1]) - func(grids[1, deriv_counter])
        for_derivatives[0, deriv_counter] = difference / step_sizes[deriv_counter]
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
    step_sizes = grids[0, 0:-1] - grids[0, 1:]
    for deriv_counter in range(1, no_derivs+1):
        difference  = func(grids[1, deriv_counter]) - func(grids[1, deriv_counter-1])
        back_derivatives[0, deriv_counter-1] = difference / step_sizes[deriv_counter-1]
    return back_derivatives

##########################        
########################## centered difference derivative
##########################
def 1D_central_diff_deriv(grids, func):
    no_derivs = np.size(grids, 1) - 2
    step_sizes = grids[0, 0:-1] - grids[0, 1:]
    central_derivatives = np.zeros((1, no_derivs))
    for deriv_counter in range(1, no_derivs+1):
        difference  = func(grids[1, deriv_counter+1]) - func(grids[1, deriv_counter-1])
        current_step_size = step_sizes[deriv_counter] + step_sizes[deriv_counter-1]
        central_derivatives[0, deriv_counter] = difference / current_step_size
    return central_derivatives

##########################        
########################## Second derivative - centered difference
##########################
def second_centre_deriv(grids, func):
    # here I assume all step_sizes are of the same size.
    no_derivs = np.size(grids, 1) - 2
    step_size = np.abs(grids[0, 1] - grids[0, 1])
    half_size = step_size/2
    second_derivatives = np.zeros((1, no_derivs))
    for deriv_counter in range(1, no_derivs+1):
        numerator = func(grids(deriv_counter-1)) - 2 * func(grids(deriv_counter)) + func(grids(deriv_counter+1))
        second_derivatives[0, deriv_counter-1] = numerator / (step_size**2)
    return second_derivatives










































