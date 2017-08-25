"""
Company. 
Very First code.
Functions would be written here
and then separated as needed into
different files, modules, etc.
"""

import numpy as np
from math import log, sqrt
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
        numerator = func(grids(deriv_counter-1)) - 2 * func(grids(deriv_counter)) + 
                                                     func(grids(deriv_counter+1))
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
    init_grid_x_vel = config.init_x_vel_scalar * np.ones(config.imax-1, config.jmax-1)
    init_grid_y_vel = config.init_y_vel_scalar * np.ones(config.imax-1, config.jmax-1)
    init_grid_press = config.init_press_scalar * np.ones(config.imax-1, config.jmax-1)
    if config.problem != 'backstep':
        init_grid_x_vel[:, 0: 1+int(config.jmax/2)] = 0.
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
                      x_grid_vel, y_grid_vel, 
                      Reynolds):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         delx is step sizes in x-direction
    #         dely is step sizes in y-direction
    #         x_grid_vel are velocities at interior grids in x-direction
    #         y_grid_vel are velocities at interior grids in y-direction
    #         Reynolds is Reynolds number.
    #         safety_tau is safety factor for time step size control tau.
    
    # output: The stepsize delta_t for the next time step is calculated.
    if config.safety_tau >= 0:
        first_element = np.reciprocal(1./(config.delta_x ** 2) + 
                                               1./(config.delta_y ** 2)) * (Reynolds / 2.)
        second_element = config.delta_x / np.abs(np.max(x_grid_vel))
        third_element = config.delta_y / np.abs(np.max(y_grid_vel)) 
        delta_t = config.safety_tau * 
                                  np.minimum(first_element, second_element, third_element)
    else:
        delta_t = config.deta_t
    return delta_t    
    
"""    
/*----------------------------------------------------------------*/
/* Computation of tentative velocity field (F,G)                  */
/*----------------------------------------------------------------*/
Problem 6 of the book. 
Computation of F and G according to (3.36) and (3.37). 
At the boundary the formulas (3.42) must be applied.
"""
# I might be able to vectorize this! look @ it closely!!!
def compute_FG(config, state, flag, temp, F, G):
    for ii in range(1, config.imax):
        for jj in range(1, config.jmax+1):
            if (( ((flag[ii, jj]   & C_F) and (flag[ii, jj]   < C_E)) and
               ((   flag[ii+1, jj] & C_F) and (flag[ii+1, jj] < C_E)) ):
                DU2DX = (
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii+1, jj]) +
	                     config.gamma * np.abs(state.x_grid_vel[ii, jj] + 
	                     state.x_grid_vel[ii+1, jj]) * 
	                     (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii+1, jj]) -
	                     (state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii, jj]) * 
	                     (state.x_grid_vel[i-1][j] + state.x_grid_vel[ii, jj]) -
                         config.gamma * np.abs(state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii, jj]) * 
                         (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii, jj])
                         )
                         /(4.0 * config.delta_x);
                
                DUVDY = (
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj+1]) +
                         config.gamma * np.abs(state.y_grid_vel[ii, jj] + 
                         state.y_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[i][j] - state.x_grid_vel[i][j+1]) -
	                     (state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii+1, jj-1]) * 
	                     (state.x_grid_vel[ii, jj-1] + state.x_grid_vel[ii, jj])-
	                     config.gamma * np.abs(state.y_grid_vel[i][j-1] + state.y_grid_vel[ii+1, jj-1]) * 
	                     (state.x_grid_vel[ii, jj-1] - state.x_grid_vel[ii, jj])
	                     )	                     
                         /(4.0 * config.delta_y)
                
                LAPLU = (state.x_grid_vel[ii+1, jj] - 2.0 * state.x_grid_vel[ii, jj] + state.x_grid_vel[ii-1, jj]) / (config.delta_x ** 2) +  
	                    (state.x_grid_vel[ii, jj+1] - 2.0 * state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj-1]) / (config.delta_y ** 2)
                  
                
                F[ii, jj] = state.x_grid_vel[ii, jj] + 
                            config.delta_t * (LAPLU / config.Ray_no - DU2DX - DUVDY + config.GX) -
                            config.delta_t * config.beta * config.GX * (temp[ii, jj] + 
                            temp[ii+1, jj])/2
            else:
                F[ii, jj] = state.x_grid_vel[ii, jj]
                
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax):
            # /* only if both adjacent cells are fluid cells */
            if( ((flag[ii, jj]   & C_F) and (flag[ii, jj]   < C_E)) and
                ((flag[ii, jj+1] & C_F) and (flag[ii, jj+1] < C_E)) ):
                
                DUVDX = (
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj+1]) * 
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii+1, jj]) +
                          config.gamma * np.abs(state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj+1]) *
                         (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii+1, jj]) -
                         (state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii-1, jj+1]) * 
                         (state.y_grid_vel[ii-1, jj] + state.y_grid_vel[ii, jj])-
                          config.gamma * np.abs(state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii-1, jj+1]) * 
                          (state.y_grid_vel[ii-1, jj]-state.y_grid_vel[ii, jj])
                        )
                         / (4.0 * config.delta_x)
                
                DV2DY = (
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii, jj+1]) *
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii, jj+1]) +
                          gamma * np.abs(state.y_grid_vel[ii, jj] + state.y_grid_vel[ii, jj+1]) * 
                          (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj+1]) -
                         (state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii, jj]) * 
                         (state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii, jj]) -
                          gamma * np.abs(state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii, jj]) * 
                          (state.y_grid_vel[ii, jj-1] - state.y_grid_vel[ii, jj])
                        )
                        /(4.0 * config.delta_y)

                LAPLV = (state.y_grid_vel[ii+1, jj] - 2.0 * state.y_grid_vel[ii, jj] + 
                         state.y_grid_vel[ii-1 ,jj]) / (config.delta_x ** 2) +
                        (state.y_grid_vel[ii, jj+1] - 2.0 * state.y_grid_vel[ii, jj] + 
                          state.y_grid_vel[ii, jj-1]) / (config.delta_y ** 2)

                G[ii, jj] = state.y_grid_vel[ii, jj] + config.delta_t * ( LAPLV / config.Ray_no - DUVDX - DV2DY + config.GY)
                            -config.delta_t * config.beta * config.GY * (temp[ii, jj] + temp[ii, jj+1]) / 2
            else:
                G[ii, jj] = state.y_grid_vel[ii, jj];
 """
  /* F und G at external boundary */
  /*------------------------------*/ 
 """
    for ii in xrange(1, config.jmax+1):
        F[0, jj] = state.x_grid_vel[0, jj]
        F[config.imax, jj] = state.x_grid_vel[imax, jj]
        
    for ii in xrange(1, config.imax+1):
        G[ii, 0] = state.y_grid_vel[i][0]
        G[ii, config.jmax] = state.y_grid_vel[ii, config.jmax]
    return F, G

"""
/*-------------------------------------------------------------*/
/* Computation of the right hand side of the pressure equation */
/*-------------------------------------------------------------*/
Problem 7 of the book.
Computation ofthe right-hand side of the pressure equation (3.38).
"""
def compute_RHS(config, F, G, RHS, flag):
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if ((flag[ii, jj] & C_F) and (flag[ii, jj] < 0x0100))  
        # /* only for fluid and non-surface cells */
                RHS[ii, jj] = (
                               (F[ii, jj] - F[ii-1, jj]) / config.delta_x + 
                               (G[ii, jj] - G[ii, jj-1]) / config.delta_y
                              )/config.delta_t
    return RHS

"""
/*-------------------------------------------------------------*/
/* Computation of the right hand side of the pressure equation */
/*-------------------------------------------------------------*/

Problem 8 of the book:
SOR iteration for the pressure Poisson 
equation according to (3.44). 
The iteration is terminated once 
the residual norm res drops below 
the tolerance limit eps 
(absolute or relative, multiplied 
by the norm of the initial pressure) 
or once the maximal number of 
iterations itermax is reached. 
Upon completion, the number of 
steps taken is returned and the 
current residual norm is stored in res.
If the pressure boundary values 
are treated using the second method, 
then the boundary values must be 
set according to (3.48) prior to 
each iteration step.
"""
config.iteration_max
def POISSON(config, state, RHS, flag, press_residual, ifull):
    # config.res_norm = press_residual \* If press_residual is a parameter
    # provided by user we can use config. but it that will be changing during
    # dynamics then it has to be taken out of config. I still do not know.
    # Have to read the book! */
    p0 = 0.0
    rdx2 = 1./(config.delta_x ** 2)
    rdy2 = 1./(config.delta_y ** 2)
    beta_2 = -config.relax_param /(2.0*(rdx2+rdy2))
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if (flag[ii, jj] & C_F):
                p0 += (state.pressures[ii,jj] ** 2)
    p0 = sqrt(p0/ifull)
    if p0 < 0.0001:
        p0 = 1.0
    """
    /* SOR-iteration */
    /*---------------*/
    """
    for (iter=1;iter<=itermax;iter++):
        if (p_bound == 1):
        """
        /* modify the equation at the boundary */
        /*-------------------------------------*/
        /* relaxation for fluid cells */
        /*----------------------------*/
        """
            for ii in xrange(1, config.imax+1):
                for jj in xrange(1, config.jmax+1):
                    # /* five point star for interior fluid cells */
                    if (flag[ii, jj] == 0x001f):
                        state.pressures[ii, jj] = (1. - config.relax_param) * state.pressures[ii, jj] - 
                                                  beta_2 * 
                                                  (
                                                   (state.pressures[ii+1, jj] + 
                                                    state.pressures[ii-1, jj]) * rdx2 +
                                                   (state.pressures[ii, jj+1] + 
                                                    state.pressures[ii, jj-1]) * rdy2 - 
                                                    RHS[ii, jj]
                                                  )
                    # /* modified star near boundary */
                    elif ((flag[ii, jj] & C_F) and (flag[ii, jj] < 0x0100)):
                        beta_mod = -config.relax_param/((eps_E + eps_W) * rdx2 + (eps_N + eps_S) * rdy2);
                        state.pressures[ii, jj] = (1. - config.relax_param) * state.pressures[ii, jj] -
                                   beta_mod * ( (eps_E * state.pressures[ii+1, jj] + 
                                   eps_W * state.pressures[ii-1, jj]) * rdx2 +
                                   (eps_N * state.pressures[ii, jj+1] + 
                                   eps_S * state.pressures[ii, jj-1]) * rdy2 -
                                   RHS[i][j])
            """
            /* computation of residual */
            /*-------------------------*/
            """     







