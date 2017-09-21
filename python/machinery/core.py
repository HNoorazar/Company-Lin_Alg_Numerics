"""
Company. 
Very First code.
Functions would be written here
and then separated as needed into
different files, modules, etc.
"""

import numpy as np
from math import log, sqrt
import definitions as defn
import IO as IO
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
This is the same as "INIT_UVP" in init.py of C++ code.
/*---------------------------------------------------------------*/
/* Setting the initial values for U,V,P, and TEMP                */
/*---------------------------------------------------------------*/
"""
def initialize_grid_state(config, state):
    # input:  imax is number of interior cells in x-direction
    #         jmax is number of interior cells in y-direction
    #         initial_x_vel_scalar is initial velocity in x-direction
    #         initial_y_vel_scalar is initial velocity in y-direction
    #         initial_press_scalar is initial initial pressure
    # output: velocities in x- and y- directions, pressure and temperature on all
    #         interior grid points.
    state.x_grid_vel = config.init_x_vel_scalar * np.ones(config.imax+2, config.jmax+2)
    state.y_grid_vel = config.init_y_vel_scalar * np.ones(config.imax+2, config.jmax+2)
    state.pressures = config.init_press_scalar * np.ones(config.imax+2, config.jmax+2)
    state.temp  = config.initial_temp_scalar * np.ones(config.imax+2, config.jmax+2)
    # /* Set U=0.0 in the lower half for the flow past a backward facing step */
    # /*----------------------------------------------------------------------*/
    if config.problem == 'backstep':
        state.x_grid_vel[:, 0: 1+int(config.jmax/2)] = 0.
    return state
"""
\*----------------------------------------------------------------------*/
\* Initializing the integer array FLAG, dependent of the problem type   */
\*----------------------------------------------------------------------*/
"""
def initialize_flag(config, state):
    state.flag = np.zeros((config.imax+2, config.jmax+2))
    """"
    /* boundary strip to C_B */
    /*-----------------------*/
    """"
    state.flag[:, 0] = defn.C_B
    state.flag[:, config.jmax+1] = defn.C_B
    state.flag[0, 1:config.jmax+1] = defn.C_B
    state.flag[config.imax+1, 1:config.jmax+1] = defn.C_B
    """
    /* all inner cells fluid cells */
    /*-----------------------------*/
    """
    state.flag[1:config.imax+1, 1:config.jmax+1] = defn.C_F
    """
        \* problem dependent obstacle cells in the interior */
        \*--------------------------------------------------*/
    """
    # why these are weird conditions?
    if config.problem == 'fluidtrap':
        low = 1 + (9 * config.imax / 22) 
        up  = (13* config.imax / 22) 
        state.flag[low : up+1, 1 : 4*config.jmax/11 + 1 ] = defn.C_B
        state.flag[low : up+1, (8*config.jmax/11)+1 : config.jmax+1 ] = defn.C_B

    if config.problem == 'plate':
        """
        /* flow past an inclined plate */
        /* lower and upper bound of the plate */
        """
        low = 2 * config.jmax / 5;
        up  = 3 * config.jmax / 5;
        state.flag[low, low]   = defn.C_B
        state.flag[low, low+1] = defn.C_B
        state.flag[up, up-1]   = defn.C_B
        state.flag[up, up]     = defn.C_B
        # might be doable by reshape and some tricks!:
        # Toeplitz and/or band stuff.
        for ii in range(low+1, up):
            for jj in range(ii-1, ii+2):
                state.flag[ii, jj] = defn.C_B

    if (config.problem == 'backstep') or (config.problem == 'wave'):
              # \* flow past a backward facing step */
        state.flag[1:config.jmax+1, 1:jmax/2+1] = defn.C_B
    
    if config.problem == 'circle':
       # \* flow past a cylinder/circle */
        mx = 20.0/41.0 * config.jmax * config.delta_y
        my = mx
        rad1 = mx / 4.
        
        xMatrix = np.broadcast_to(np.arange(1, config.imax+1), (config.jmax, config.imax))
        xMatrix = np.transpose(xMatrix)
        xMatrix = np.array(xMatrix)
        yMatrix = np.broadcast_to(np.arange(1, config.jmax+1), (config.imax, config.jmax))
        yMatrix = np.array(yMatrix)
        xMatrix = (xMatrix - 0.5) * config.delta_x
        yMatrix = (yMatrix - 0.5) * config.delta_y
        
        A = ( ((xMatrix - mx) ** 2 + (yMatrix - my) ** 2) <= rad1 )
        state.flag[A == True] = defn.C_B
        
    if config.problem == 'molding':
        # \* circular obstacle */
        mx = config.jmax * config.delta_y/2
        my = config.jmax * config.delta_y/2
        rad1=config.jmax * config.delta_y/6
        xMatrix = np.broadcast_to(np.arange(1, config.imax+1), (config.jmax, config.imax))
        xMatrix = np.transpose(xMatrix)
        xMatrix = np.array(xMatrix)
        yMatrix = np.broadcast_to(np.arange(1, config.jmax+1), (config.imax, config.jmax))
        yMatrix = np.array(yMatrix)
        xMatrix = (xMatrix - 0.5) * config.delta_x
        yMatrix = (yMatrix - 0.5) * config.delta_y
        
        A = ( ((xMatrix - mx) ** 2 + (yMatrix - my) ** 2) <= (rad1**2) )
        state.flag[A == True] = defn.C_B
    # /* Printing the geometry of the fluid domain */
    print "nGeometry of the fluid domain:\n\n"
    print ""
    print ""
    for jj in sorted(range(0, config.jmax+1), reverse=True):
        for ii in range(0,config.imax+1):
            if not (flag[ii, jj] & defn.C_F):
                print "**"
            else:
                print "  "
        print "\n"
    print "\n"
    print "\n"
                  # \* flags for boundary cells */
    state.ibound = 0
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if not (flag[ii, jj] & defn.C_F):
                state.ibound += 1
            flag[ii, jj] += ((flag[ii-1, jj] & defn.C_F) * defn.B_W + 
                             (flag[ii+1, jj] & defn.C_F) * defn.B_O +
                             (flag[ii, jj-1] & defn.C_F) * defn.B_S + 
                             (flag[ii, jj+1] & defn.C_F) * defn.B_N)/defn.C_F
            if (flag[ii, jj] == 0x0003) or (flag[ii, jj] == 0x0007) or
               (flag[ii, jj] == 0x000b) or (flag[ii, jj] == 0x000c) or
               (flag[ii, jj] == 0x000d) or (flag[ii, jj] == 0x000e) or
               (flag[ii, jj] == 0x000f):
                print('Illegal obstacle cell {}{}'.format(ii, jj))
                break # the C++ code has exit(0) here.  
                      # I'm still not sure what that is supposed to do!
    return state


def write_state(state, filename):
    system_state = {'x_grid_vel': state.x_grid_vel,
                    'y_grid_vel': state.y_grid_vel,
                    'pressures': state.pressures,
                    'temp': state.temp,
                    'flag': state.flag
                    }
    try:
        IO.save_matrix(filename, system_state)
    except:
        print("ERROR: could not write matrix file "+filename)                    

"""
Problem 3 of page 43 of the book.
The stepsize delt for the 
next time step is calculated 
according to (3.50). 
In case of negative tau 
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
    
#                                  UVP.C 
"""
/*---------------------------------------------------------------*/
/* Computation of new temperature                                */
/*---------------------------------------------------------------*/
"""
def compute_temp(config, state):
    T2 = np.zeros((config.imax+2, config.jmax+2))
    indelx2 = 1. / (config.delta_x ** 2)
    indely2 = 1. / (config.delta_y ** 2)
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if ( (state.flag[ii, jj] & defn.C_F) and (state.flag[ii, jj] < defn.C_E) ):
                LAPLT = (
                         state.temp[ii+1, jj] - 
                         2.0 * state.temp[ii, jj] + 
                         state.temp[ii-1, jj]
                         ) * indelx2 +
		                (
		                 state.temp[ii, jj+1] - 
		                 2.0 * state.temp[ii, jj] + state.temp[ii, jj-1]
		                 ) * indely2
		        DUTDX = ((
                          state.x_grid_vel[ii, jj] * 0.5 * 
                          (state.temp[ii, jj]+state.temp[ii+1, jj]) -
                          state.x_grid_vel[ii-1, jj] * 0.5 * 
                          (state.temp[ii-1, jj] + state.temp[ii, jj])
				          ) +
		                  config.gamma * 
		                   ( np.abs(state.x_grid_vel[ii, jj]) * 0.5 * 
		                   (state.temp[ii, jj] - state.temp[ii+1, jj]) -
				             np.abs(state.x_grid_vel[ii-1, jj]) * 0.5 * 
				             (state.temp[ii-1, jj] - state.temp[ii, jj])
				           )) / config.delta_x;
				DVTDY = (( state.y_grid_vel[ii, jj] * 0.5 * 
				           ( state.temp[ii, jj] + state.temp[ii, jj+1]) -
                           state.y_grid_vel[ii, jj-1] * 0.5 * 
                           ( state.temp[ii, jj-1] + state.temp[ii, jj])) +
                           gamma * (np.abs(state.y_grid_vel[ii, jj]) * 0.5 * 
                                               (state.temp[ii, jj] - state.temp[ii, jj+1]) -
                           np.abs(state.y_grid_vel[ii, jj-1]) * 0.5 * 
                               ( state.temp[ii, jj-1] - state.temp[ii, jj]))
                           ) / config.delta_y;
                T2[ii, jj] = state.temp[ii, jj] + 
                             config.delta_t * 
                             ( LAPLT / config.Ray_no / config.Pr - DUTDX - DVTDY);
    """
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if ( (state.flag[ii, jj] & defn.C_F) and (state.flag[ii, jj] < defn.C_E) ):
                state.temp[ii, jj] = T2[ii, jj]
    """
    # This is vectorized version of the code above. (last 4 lines)
    A =  state.flag[1:config.imax+1, 1:config.jmax+1] & defn.C_F
    B = (state.flag[1:config.imax+1, 1:config.jmax+1] < defn.C_E) * 1
    C = A * B
    D = np.zeros((config.imax+2, config.jmax+2))
    D[1:config.imax+1, 1:config.jmax+1] = C
    state.temp[D == True] = T2[D == True]
    del T2, A, B, C, D # free up memory
    return state

"""
/*----------------------------------------------------------------*/
/* Computation of tentative velocity field (F,G)                  */
/*----------------------------------------------------------------*/
Problem 6 of the book. 
Computation of F and G according to (3.36) and (3.37). 
At the boundary the formulas (3.42) must be applied.
"""
# I might be able to vectorize this! look @ it closely!!!
def compute_FG(config, state, F, G):
    for ii in range(1, config.imax):
        for jj in range(1, config.jmax+1):
            if ( ((state.flag[ii, jj]   & C_F) and (state.flag[ii, jj]   < C_E)) and
                 ((state.flag[ii+1, jj] & C_F) and (state.flag[ii+1, jj] < C_E)) 
               ):
                DU2DX = (
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii+1, jj]) +
	                     config.gamma * np.abs(state.x_grid_vel[ii, jj] + 
	                     state.x_grid_vel[ii+1, jj]) * 
	                     (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii+1, jj]) -
	                     (state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii, jj]) * 
	                     (state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii, jj]) -
                         config.gamma * np.abs(state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii, jj]) * 
                         (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii, jj])
                         )
                         /(4.0 * config.delta_x);
                
                DUVDY = (
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj+1]) +
                         config.gamma * np.abs(state.y_grid_vel[ii, jj] + 
                         state.y_grid_vel[ii+1, jj]) * 
                         (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii, jj+1]) -
	                     (state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii+1, jj-1]) * 
	                     (state.x_grid_vel[ii, jj-1] + state.x_grid_vel[ii, jj])-
	                     config.gamma * np.abs(state.y_grid_vel[ii, jj-1] + state.y_grid_vel[ii+1, jj-1]) * 
	                     (state.x_grid_vel[ii, jj-1] - state.x_grid_vel[ii, jj])
	                     )	                     
                         /(4.0 * config.delta_y)
                
                LAPLU = (state.x_grid_vel[ii+1, jj] - 2.0 * state.x_grid_vel[ii, jj] + state.x_grid_vel[ii-1, jj]) / (config.delta_x ** 2) +  
	                    (state.x_grid_vel[ii, jj+1] - 2.0 * state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj-1]) / (config.delta_y ** 2)
                  
                
                F[ii, jj] = state.x_grid_vel[ii, jj] + 
                            config.delta_t * (LAPLU / config.Ray_no - DU2DX - DUVDY + config.GX) -
                            config.delta_t * config.beta * config.GX * (state.temp[ii, jj] + 
                            state.temp[ii+1, jj])/2
            else:
                F[ii, jj] = state.x_grid_vel[ii, jj]
                
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax):
            # /* only if both adjacent cells are fluid cells */
            if( ((state.flag[ii, jj]   & C_F) and (state.flag[ii, jj]   < C_E)) and
                ((state.flag[ii, jj+1] & C_F) and (state.flag[ii, jj+1] < C_E)) ):
                
                DUVDX = (
                         (state.x_grid_vel[ii, jj] + state.x_grid_vel[ii, jj+1]) * 
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii+1, jj]) +
                          config.gamma * np.abs(state.x_grid_vel[ii, jj] + 
                          state.x_grid_vel[ii, jj+1]) *
                         (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii+1, jj]) -
                         (state.x_grid_vel[ii-1, jj] + state.x_grid_vel[ii-1, jj+1]) * 
                         (state.y_grid_vel[ii-1, jj] + state.y_grid_vel[ii, jj])-
                          config.gamma * np.abs(state.x_grid_vel[ii-1, jj] + 
                          state.x_grid_vel[ii-1, jj+1]) * 
                          (state.y_grid_vel[ii-1, jj]-state.y_grid_vel[ii, jj])
                        )
                         / (4.0 * config.delta_x)
                
                DV2DY = (
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii, jj+1]) *
                         (state.y_grid_vel[ii, jj] + state.y_grid_vel[ii, jj+1]) +
                          gamma * np.abs(state.y_grid_vel[ii, jj] + 
                          state.y_grid_vel[ii, jj+1]) * 
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
                    
                G[ii, jj] = state.y_grid_vel[ii, jj] + config.delta_t * 
                            ( LAPLV / config.Ray_no - DUVDX - DV2DY + config.GY)
                            -config.delta_t * config.beta * config.GY * 
                            (state.temp[ii, jj] + state.temp[ii, jj+1]) / 2
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
        G[ii, 0] = state.y_grid_vel[ii, 0]
        G[ii, config.jmax] = state.y_grid_vel[ii, config.jmax]
    return F, G

"""
/*-------------------------------------------------------------*/
/* Computation of the right hand side of the pressure equation */
/*-------------------------------------------------------------*/
Problem 7 of the book.
Computation ofthe right-hand side of the pressure equation (3.38).
"""
def compute_RHS(config, state, F, G, RHS):
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if ((state.flag[ii, jj] & C_F) and (state.flag[ii, jj] < 0x0100))  
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
def POISSON(config, state, RHS, press_residual, ifull):
    # config.res_norm = press_residual \* If press_residual is a parameter
    # provided by user we can use config. but it that will be changing during
    # dynamics then it has to be taken out of config. I still do not know.
    # Have to read the book! */
    p0 = 0.
    rdx2 = 1./(config.delta_x ** 2)
    rdy2 = 1./(config.delta_y ** 2)
    beta_2 = -config.relax_param / (2.0 * (rdx2 + rdy2))
    # Can be vectorized:
    for ii in xrange(1, config.imax+1):
        for jj in xrange(1, config.jmax+1):
            if (state.flag[ii, jj] & C_F):
                p0 += (state.pressures[ii,jj] ** 2)
    p0 = sqrt(p0/ifull)
    if p0 < 0.0001:
        p0 = 1.
    """
    /* SOR-iteration */
    /*---------------*/
    """
    for iter in xrange(1, config.iteration_max+1)
        if (config.p_bound == 1):
        """
        /* modify the equation at the boundary */
        /*-------------------------------------*/
        /* relaxation for fluid cells */
        /*----------------------------*/
        """
            for ii in xrange(1, config.imax+1):
                for jj in xrange(1, config.jmax+1):
                    # /* five point star for interior fluid cells */
                    if (state.flag[ii, jj] == 0x001f):
                        state.pressures[ii, jj] = (1.-config.relax_param)*state.pressures[ii, jj] - 
                                                  beta_2 * 
                                                  (
                                                   (state.pressures[ii+1, jj] + 
                                                    state.pressures[ii-1, jj]) * rdx2 +
                                                   (state.pressures[ii, jj+1] + 
                                                    state.pressures[ii, jj-1]) * rdy2 - 
                                                    RHS[ii, jj]
                                                  )
                    # /* modified star near boundary */
                    elif ((state.flag[ii, jj] & C_F) and (state.flag[ii, jj] < 0x0100)):
                        beta_mod = -config.relax_param/((defn.eps_E + defn.eps_W) * rdx2 + 
                                    (defn.eps_N + defn.eps_S) * rdy2)
                        
                        state.pressures[ii, jj] = (1. - config.relax_param) * state.pressures[ii, jj] -
                                                  beta_mod * ( (eps_E * state.pressures[ii+1, jj] + 
                                                  defn.eps_W * state.pressures[ii-1, jj]) * rdx2 +
                                                  (defn.eps_N * state.pressures[ii, jj+1] + 
                                                  eps_S * state.pressures[ii, jj-1]) * rdy2 -
                                                  RHS[ii, jj])
                
                # \* computation of residual */
                res = 0.
                """
                # The following would be a vectorized attemp of
                # version of the function Possion, lines 197-208, in upv.c
                A = (state.flag[1:config.imax+1, 1:config.jmax+1] & defn.C_F) * 1
                B = (state.flag[1:config.imax+1, 1:config.jmax+1] <   0x0100) * 1
                C = A * B
                D = np.zeros((config.imax+2, config.jmax+2))
                D[1:config.imax+1, 1:config.jmax+1] = C
                state.temp[D == True] = T2[D == True]
                """
                for ii in xrange(1, config.imax+1):
                    for jj in xrange(1, config.jmax+1):
                        if ((state.flag[ii, jj] & C_F) and (state.flag[ii, jj] < 0x0100)):
                            # \* only fluid cells */
                            add = (defn.eps_E * (state.pressures[ii+1, jj]-state.pressures[ii, jj]) - 
                                   defn.eps_W * (state.pressures[ii, jj]-state.pressures[ii-1, jj])) * rdx2 +
                                  (defn.eps_N * (state.pressures[ii, jj+1]-state.pressures[ii, jj]) -
                                   defn.eps_S * (state.pressures[ii, jj]-P[ii, jj-1])) * rdy2 -  
                                   RHS[ii, jj]
                            res += add**2
                res = sqrt((res)/ifull)/p0
            if res < config.stop_toler:
                break
        elif (config.p_bound == 2):

        # \* copy values at external boundary */
        state.pressures[1:config.imax+1, 0] = state.pressures[1:config.imax+1, 1]
        state.pressures[1:config.imax+1, config.jmax+1] = state.pressures[1:config.imax+1, config.jmax]
        
        state.pressures[0, 1:config.jmax+1] = state.pressures[1, 1:config.jmax+1]
        state.pressures[config.imax+1, 1:config.jmax+1] = state.pressures[config.imax+1, 1:config.jmax+1]
        # \* and at interior boundary cells */
        for ii in xrange(1, config.imax+1):
            for jj in xrange(1, config.jmax+1):
	            if (state.flag[ii, jj] >= defn.B_N and flag[ii, jj] <= defn.B_SO):
                    if state.flag[ii, jj] == defn.B_N:
	                    state.pressures[ii, jj] = state.pressures[ii, jj+1]
	                    break
                    elif state.flag[ii, jj] == defn.B_O:
                        state.pressures[ii, jj] = state.pressures[ii+1, jj]
                        break
                    elif state.flag[ii, jj] == defn.B_S:
                        state.pressures[ii, jj] = state.pressures[ii, jj-1]
                        break
                    elif state.flag[ii, jj] == defn.B_W:
                        state.pressures[ii, jj] = state.pressures[ii-1, jj]
                        break
                    elif state.flag[ii, jj] == defn.B_NO:
                        state.pressures[ii, jj] = 0.5*(state.pressures[ii, jj+1]+state.pressures[ii+1, jj])
                        break
                    elif state.flag[ii, jj] == defn.B_SO:
                        state.pressures[ii, jj] = 0.5*(state.pressures[ii, jj-1]+state.pressures[ii+1, jj])
                        break
	                elif state.flag[ii, jj] == defn.B_SW:
                        state.pressures[ii, jj] = 0.5*(state.pressures[ii, jj-1]+state.pressures[ii-1, jj])
                        break
                    elif state.flag[ii, jj] == defn.B_NW:
                        state.pressures[ii, jj] = 0.5*(state.pressures[ii, jj+1]+state.pressures[ii-1, jj])
                        break
                    else:
                        break
                        
        # \* relaxation for fluid cells */
        for ii in xrange(1, config.imax+1):
            for jj in xrange(1, config.jmax+1):
                if ((state.flag[ii, jj] & defn.C_F) and (state.flag[ii, jj] < 0x0100)):
                    # \* only fluid cells */
                    add =  (state.pressures[ii+1, jj]-2*state.pressures[ii, jj]+state.pressures[ii-1, jj])*rdx2+
		                   (state.pressures[ii, jj+1]-2*state.pressures[ii, jj]+state.pressures[ii, jj-1])*rdy2-
		                    RHS[ii, jj]
		            res += add ** 2
		res = sqrt((res)/ifull)/p0
		# \* convergence? */
		if res < config.stop_toler:
		    break
    return state, iter

"""
\* Computation of new velocity values */
"""
def compute_new_velocity(config, state, F, G):
    for ii in xrange(1, config.imax):
        for jj in xrange(1, config.jmax+1):
        # \* only if both adjacent cells are fluid cells */
            if ( ((state.flag[ii, jj] & defn.C_F)   and (state.flag[ii, jj] < defn.C_E)) and
                ((state.flag[ii+1, jj] & defn.C_F) and (state.flag[ii+1, jj] < defn.C_E)) ):
                state.x_grid_vel[ii, jj] = F[ii, jj]-
                                           (state.pressures[ii+1, jj] - state.pressures[ii, jj])*config.delta_t/config.delta_x
    for (i=1;i<=imax;i++)
        for (j=1;j<=jmax-1;j++)
            # \* only if both adjacent cells are fluid cells */
            if( ((state.flag[ii, jj]   & defn.C_F) and (state.flag[ii, jj] < defn.C_E)) and
                ((state.flag[ii, jj+1] & defn.C_F) and (state.flag[ii, jj+1] < defn.C_E)) ):
                state.y_grid_vel[ii, jj] = G[ii, jj]-
                                           (state.pressures[ii, jj+1]-state.pressures[ii, jj])*config.delta_t/config.delta_y
    return state

"""
/*------------------------------------------------------------*/
/* Computation of adaptive time stepsize satisfying           */
/* the CFL stability criteria                                 */
/* and set the flag "write" if some data has to be written    */
/* into a file.                                               */
/*------------------------------------------------------------*/
"""
def compute_delta_t(config, state, write, delt_inj, del_streak, del_trace, del_vec):
    # \* delt satisfying CFL conditions */
    if config.safety_tau >= 10. ** (-10):
        # \* else no time stepsize control */
        umax = 10. ** (-10)
        vmax = 10. ** (-10)
        
        umax = np.max(umax, np.max(np.abs(state.x_grid_vel[1: , 1: ])))
        vmax = np.max(vmax, np.max(np.abs(state.y_grid_vel[1: , 1: ])))

        deltu = config.delta_x / umax
        deltv = config.delta_y / vmax
                 
        if(config.Pr < 1):
            deltRePr = 1 / (np.reciprocal(config.delta_x**2) + np.reciprocal(delta_y**2))*Re*Pr/2.
        else:
            deltRePr = 1/( np.reciprocal(config.delta_x**2) + np.reciprocal(delta_y**2))*Re/2.
        
        if(deltu<deltv):
          if(deltu < deltRePr):
              config.delta_t = deltu
          else:
              config.delta_t = deltRePr
        else:
            if(deltv<deltRePr):
                config.delta_t = deltv
            else:
                config.delta_t = deltRePr
        # \* multiply by safety factor */
        config.delta_t = config.safety_tau * config.delta_t
        # \* look if some data has to be written to a file in the next time step */
    write = 0
    t_neu = state.current_time + config.delta_t
    t_trace = t_inj = t_streak = t_vec = t_neu + 1.0e+10
    
    if( np.int(state.current_time / config.del_trace)!= np.int(t_neu / config.del_trace) ):
        t_trace = np.int(t_neu / config.del_trace) * config.del_trace
        write += 1
        
    if( np.int(state.current_time / config.del_inj)!= np.int(t_neu / config.del_inj) ):
        t_inj = np.int(t_neu / config.del_inj) * config.del_inj
        write += 2
        
    if( np.int( state.current_time / config.del_streak) != np.int(t_neu / config.del_streak) ):
        t_streak = np.int(t_neu / config.del_streak) * config.del_streak
        write += 4
        
    if( np.int(state.current_time / config.del_vec) != np.int(t_neu / config.del_vec) ):
        t_vec = np.int(t_neu/ config.del_vec) * config.del_vec
        write += 8
    return config, write

