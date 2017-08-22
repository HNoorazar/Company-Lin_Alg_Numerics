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
def initialize_grid_state(config, initial_x_vel, initial_y_vel, initial_pressure):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         initial_x_vel is initial velocity in x-direction
    #         initial_y_vel is initial velocity in y-direction
    #         initial_pressure is initial initial pressure
    # output: velocities in x- and y- directions and pressure on all
    #         interior grid points.
    init_grid_x_vel = initial_x_vel * np.ones(config.imax-1, config.jmax-1)
    init_grid_y_vel = initial_y_vel * np.ones(config.imax-1, config.jmax-1)
    init_grid_press = initial_pressure * np.ones(config.imax-1, config.jmax-1)
    return init_grid_x_vel, init_grid_y_vel, init_grid_press
    
"""
Problem 3 of page 43 of the book.
The stepsize delt for the 
next time step is calculated 
according to (3.50). 
In case of neg- ative tau 
the stepsize read in READ-PARAMETER 
is to be use
"""
def compute_time_step(config,
                      x_velocities, y_velocities, 
                      Reynolds):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         delx is step sizes in x-direction
    #         dely is step sizes in y-direction
    #         x_velocities are velocities at interior grids in x-direction
    #         y_velocities are velocities at interior grids in y-direction
    #         Reynolds is Reynolds number.
    #         safety_tau is safety factor for time step size control tau.
    
    # output: The stepsize delta_t for the next time step is calculated.
    if config.safety_tau >= 0:
        first_element = np.reciprocal(1./(config.delta_x ** 2) + 1./(config.delta_y ** 2)) * (Reynolds / 2.)
        second_element = config.delta_x / np.abs(np.max(x_velocities))
        third_element = config.delta_y / np.abs(np.max(y_velocities)) 
        delta_t = config.safety_tau * np.minimum(first_element, second_element, third_element)
    else:
        delta_t = config.deta_t
    return delta_t    
    
    

"""
-----------------------------------------------------------------
 Setting the boundary conditions at the boundary strip.          
 The flags wW,wE,wN, and wS can have the values:                 
 1 = slip               2 = no-slip                              
 3 = outflow            4 = periodic                             
 Moreover, no-slip conditions are set at internal obstacle cells 
 by default.                                                     
 For temperature, adiabatic boundary conditions are set.         
-----------------------------------------------------------------
"""
# I have to look into this and vectorize it, if possible!
# this is just copied from C++ code!
def set_boun_cond(state, config, temp, flag):
    # Left and right boundary
    for jj in xrange(0, config.jmax+2):
        if config.wW == 1:
            state.x_velocities[0, jj] = 0.
            state.y_velocities[0, jj] = state.y_velocities[1, jj]
        elif config.wW == 2:
            state.x_velocities[0, jj] = 0.
            state.y_velocities[0, jj] = -state.y_velocities[1, jj]
        elif config.wW == 3:
            state.x_velocities[0, jj] = state.x_velocities[1, jj]
            state.y_velocities[0, jj] = state.y_velocities[1, jj]
        elif config.wW == 4:
            state.x_velocities[0, jj] = state.x_velocities[config.imax - 1, jj]
            state.y_velocities[0, jj] = state.y_velocities[config.imax - 1, jj]
            state.y_velocities[1, jj] = state.y_velocities[config.imax, jj]            
            state.pressures[1, jj] = state.pressures[imax, jj]
        temp[0, jj] = temp[1, jj]
        
        if config.wE == 1:
            state.x_velocities[config.imax, jj] = 0.
            state.y_velocities[config.imax+1, jj] = state.y_velocities[config.imax, jj]
        elif config.wE == 2:
            state.x_velocities[config.imax, jj] = 0.
            state.y_velocities[config.imax+1, jj] = -state.y_velocities[config.imax, jj]
        elif config.wE == 3:
            state.x_velocities[config.imax, jj] = state.x_velocities[config.imax-1, jj]
            state.y_velocities[config.imax+1, jj] = state.y_velocities[config.imax, jj]
        elif config.wE == 4:
            state.x_velocities[config.imax, jj] = state.x_velocities[1, jj]
            state.y_velocities[config.imax+1, jj] = state.y_velocities[2, jj]
        temp[config.imax+1, jj] = temp[config.imax, jj]
    # Northern and Southern boundary conditions
    for ii in xrange(0, config.imax+2):
       if config.wN == 1:
           state.y_velocities[ii, config.jmax] = 0.
           state.x_velocities[ii, config.jmax+1] = state.x_velocities[ii, config.jmax]
        if config.wN == 2:
           state.y_velocities[ii, config.jmax] = 0.
           state.x_velocities[ii, config.jmax+1] = -state.x_velocities[ii, config.jmax]
        if config.wN == 3:
            state.y_velocities[ii, config.jmax] = state.y_velocities[ii, config.jmax-1]
            state.x_velocities[ii, config.jmax+1] = state.x_velocities[ii, config.jmax]
        if config.wN == 4:
            state.y_velocities[ii, config.jmax] = state.y_velocities[ii, 1]
            state.x_velocities[ii, config.jmax+1] = state.x_velocities[ii, 2]
        temp[ii,0] = temp[ii,1]
        
        if config.wS == 1:
            state.y_velocities[ii, 0] = 0.
            state.x_velocities[ii, 0] = state.x_velocities[ii, 1]
        elif config.wS == 2:
            state.y_velocities[ii, 0] = 0.
            state.x_velocities[ii, 0] = -state.x_velocities[ii, 1]
        elif config.wS == 3:
            state.y_velocities[ii, 0] = state.y_velocities[ii, 1]
            state.x_velocities[ii, 0] = state.x_velocities[ii, 1]
        elif config.wS == 4:
            state.y_velocities[ii, 0] = state.y_velocities[ii, config.jmax-1]
            state.x_velocities[ii, 0] = state.x_velocities[ii, config.jmax-1]
            state.x_velocities[ii, 1] = state.x_velocities[ii, config.jmax]
            state.pressures[ii, 1] = state.pressures[ii, jmax]
        temp[ii, jmax+1] = temp[ii, jmax]
#  /* setting the boundary values at inner obstacle cells */
#  /*                  (only no-slip)                     */
#  /*-----------------------------------------------------*/
    for ii in xrange(0, imax+1):
        for jj in xrange(0, jmax+1):
            if flag[ii, jj] and 
            """ What is 0x000f """
            
            

"""
/*-----------------------------------------------------------------*/
/* Setting specific boundary conditions, depending on "problem"    */
/*-----------------------------------------------------------------*/
Again, this is not consistent with the book.
Problem 5 has five input, the C++ code has nine!
"""
def set_specific_conditions(state, config, temp, problem)







    