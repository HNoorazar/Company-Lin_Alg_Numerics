"""
Company. 
Very First code.
Functions would be written here
and then separated as needed into
different files, modules, etc.
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

"""
Problem 2 of page 43 of the book.
Initialize the velocities and pressures
at interior grid points:
"""
def initialize_vel_press(imax, jmax, initial_x_vel, initial_y_vel, initial_pressure):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         initial_x_vel is initial velocity in x-direction
    #         initial_y_vel is initial velocity in y-direction
    #         initial_pressure is initial initial pressure
    # output: velocities in x- and y- directions and pressure on all
    #         interior grid points.
    x_velocities = initial_x_vel * np.ones(imax+2, jmax+2)
    y_velocities = initial_y_vel * np.ones(imax+2, jmax+2)
    pressures = initial_pressure * np.ones(imax+2, jmax+2)
    return x_velocities, y_velocities, pressures
    
"""
Problem 3 of page 43 of the book.
The stepsize delt for the 
next time step is calculated 
according to (3.50). 
In case of neg- ative tau 
the stepsize read in READ-PARAMETER 
is to be use
"""
def compute_time_step(imax, jmax, delx, dely, x_velocities, y_velocities, Reynolds, safety_tau):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         delx is step sizes in x-direction
    #         dely is step sizes in y-direction
    #         x_velocities are velocities at interior grids in x-direction
    #         y_velocities are velocities at interior grids in y-direction
    #         Reynolds is Reynolds number.
    #         safety_tau is safety factor for time step size control tau.
    
    # output: The stepsize delt for the next time step is calculated.
    first_element = np.reciprocal(1./(delx ** 2) + 1./(dely ** 2)) * (Reynolds / 2.)
    second_element = delx / np.abs(np.max(x_velocities))
    third_element = dely / np.abs(np.max(y_velocities)) 
    delt = safety_tau * np.minimum(first_element, second_element, third_element)
    
    
    


































