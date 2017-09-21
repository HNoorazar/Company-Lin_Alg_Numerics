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
def set_boun_cond(state, config):
    # Left and right boundary
    # First Loop
	# western and eastern boundary
	if config.wW == 1:
		state.x_grid_vel[0, :] = 0.
		state.y_grid_vel[0, :] = state.y_grid_vel[1, :]
	elif config.wW == 2:
		state.x_grid_vel[0, :] = 0.
		state.y_grid_vel[0, :] = -state.y_grid_vel[1, :]
	elif config.wW == 3:
		state.x_grid_vel[0, :] = state.x_grid_vel[1, :]
		state.y_grid_vel[0, :] = state.y_grid_vel[1, :]
	elif config.wW == 4: # periodic
		state.x_grid_vel[0, :] = state.x_grid_vel[config.imax-1, :]
		state.y_grid_vel[0, :] = state.y_grid_vel[config.imax-1, :]
		state.y_grid_vel[1, :] = state.y_grid_vel[config.imax  , :]
		state.pressures[1, :] = state.pressures[config.imax, :]
	state.temp[0, :] = state.temp[1, :]

	if config.wE == 1: # free-slip
		state.x_grid_vel[config.imax, :] = 0.
		state.y_grid_vel[config.imax+1,:]= state.y_grid_vel[config.imax, :]
	elif config.wE == 2: # no-slip
		state.x_grid_vel[config.imax, :] = 0.
		state.y_grid_vel[config.imax+1,:]= -state.y_grid_vel[config.imax, :]
	elif config.wE == 3: # outflow
		state.x_grid_vel[config.imax, :] = state.x_grid_vel[config.imax-1, :]
		state.y_grid_vel[config.imax+1, :] = state.y_grid_vel[config.imax, :]
	elif config.wE == 4: # periodic
		state.x_grid_vel[config.imax, :] = state.x_grid_vel[1, :]
		state.y_grid_vel[config.imax+1, :] = state.y_grid_vel[2, :]
	state.temp[config.imax+1, :] = state.temp[config.imax, :]

	# Second Loop
	# northern and southern boundary
	if config.wN == 1:
		state.y_grid_vel[:, config.jmax] = 0.
		state.x_grid_vel[:, config.jmax+1] =  state.x_grid_vel[:, config.jmax]
	elif config.wN == 2:
		state.y_grid_vel[:, config.jmax] = 0.
		state.x_grid_vel[:, config.jmax+1] = -state.x_grid_vel[:, config.jmax]
	elif config.wN == 3:
		state.y_grid_vel[:, config.jmax]   = state.y_grid_vel[:, config.jmax-1]
		state.x_grid_vel[:, config.jmax+1] = state.x_grid_vel[:, config.jmax]
	elif config.wN == 4:
		state.y_grid_vel[:, config.jmax] = state.y_grid_vel[:, 1]
		state.x_grid_vel[:, config.jmax+1] = state.x_grid_vel[:, 2]
	state.temp[:, 0] = state.temp[:, 1]

	if config.wS == 1:
		state.y_grid_vel[:, 0] = 0.0
		state.x_grid_vel[:, 0] = state.x_grid_vel[:, 1]
	elif config.wS == 2:
		state.y_grid_vel[:, 0] = 0.0
		state.x_grid_vel[:, 0] = -state.x_grid_vel[:, 1]
	elif config.wS == 3:
		state.y_grid_vel[:, 0] = state.y_grid_vel[:, 1]
		state.x_grid_vel[:, 0] = state.x_grid_vel[:, 1]
	elif config.wS == 4:
		state.y_grid_vel[:, 0] = state.y_grid_vel[:, config.jmax-1]
		state.x_grid_vel[:, 0] = state.x_grid_vel[:, config.jmax-1]
		state.x_grid_vel[:, 1] = state.x_grid_vel[:, config.jmax]
		state.pressures[:, 1] = state.pressures[:, config.jmax]
	state.temp[:, config.jmax+1] = state.temp[:, config.jmax]
    #  /* setting the boundary values at inner obstacle cells */
    #  /*                  (only no-slip)                     */
    #  /*-----------------------------------------------------*/
    for ii in xrange(0, imax+1):
        for jj in xrange(0, jmax+1):
            if (flag[ii, jj] and 0x000f): # /* The mask 0x000f filters the */
                                          # /* obstacle cells adjacent to  */
                                          # /* fluid cells                 */
                if flag[ii, jj] == B_N:
                    state.y_grid_vel[ii, jj] = 0.
                    state.x_grid_vel[ii, jj] = -state.x_grid_vel[ii, jj+1]
                    state.x_grid_vel[ii-1, jj] = -state.x_grid_vel[ii-1, jj+1]
                    state.temp[ii, jj] = state.temp[ii, jj+1]
                    break # Do we need a break here?
                elif flag[ii, jj] == B_O:
                    state.x_grid_vel[ii, jj] = 0.
                    state.y_grid_vel[ii, jj] = -state.i_grid_vel[ii+1, jj]
                    state.y_grid_vel[ii, jj-1] = -state.i_grid_vel[ii+1, jj-1]
                    state.temp[ii, jj] = state.temp[ii+1, jj]
                    break
                elif flag[ii, jj] == B_S:
                    state.y_grid_vel[ii, jj-1] = 0.
                    state.x_grid_vel[ii, jj] = -state.x_grid_vel[ii, jj-1]
                    state.x_grid_vel[ii-1, jj] = -state.x_grid_vel[ii-1, jj-1]
                    state.temp[i][j] = state.temp[i][j-1]
                    break
                elif flag[ii, jj] == B_W:
                    state.x_grid_vel[ii-1, jj] = 0.0
                    state.y_grid_vel[ii, jj] = -state.y_grid_vel[ii-1, jj]
                    state.y_grid_vel[ii, jj-1] = -state.y_grid_vel[ii-1, jj-1]
                    state.temp[ii, jj] = state.temp[ii-1, jj]
                    break
                elif flag[ii, jj] == B_NO:
                    state.y_grid_vel[ii, jj]   = 0.;
                    state.x_grid_vel[ii, jj]   = 0.;
                    state.y_grid_vel[ii, jj-1] = -state.y_grid_vel[ii+1, jj-1];
                    state.x_grid_vel[ii-1, jj] = -state.x_grid_vel[ii-1, jj+1];
                    state.temp[ii, jj] = 0.5*(state.temp[ii, jj+1] + state.temp[ii+1, jj]);
                    break;
                elif flag[ii, jj] == B_SO:
                    state.y_grid_vel[ii, jj-1] = 0.;
                    state.x_grid_vel[ii, jj]   = 0.;
                    state.y_grid_vel[ii, jj]   = -state.y_grid_vel[ii+1, jj];
                    state.x_grid_vel[ii-1, jj] = -state.x_grid_vel[ii-1, jj-1];
                    state.temp[ii, jj] = 0.5 * (state.temp[ii, jj-1] + state.temp[ii+1, jj]);
                    break;
                elif flag[ii, jj] == B_SW:
                    state.y_grid_vel[ii, jj-1] = 0.;
                    state.x_grid_vel[ii-1, jj] = 0.;
                    state.y_grid_vel[ii, jj] = -state.y_grid_vel[ii-1, jj];
                    state.x_grid_vel[ii, jj] = -state.x_grid_vel[ii, jj-1];
                    state.temp[ii, jj] = 0.5 * (state.temp[ii, jj-1] + state.temp[ii-1, jj]);
                    break;
                elif flag[ii, jj] == B_NW:
                    state.y_grid_vel[ii, jj] = 0.;
                    state.x_grid_vel[ii-1, jj] = 0.;
                    state.y_grid_vel[ii, jj-1] = -state.y_grid_vel[ii-1, jj-1];
                    state.x_grid_vel[ii, jj] = -state.x_grid_vel[ii, jj+1];
                    state.temp[ii, jj] = 0.5 * (state.temp[ii, jj+1] + state.temp[ii-1, jj]);
                    break
    return state


"""
/*-----------------------------------------------------------------*/
/* Setting specific boundary conditions, depending on "problem"    */
/*-----------------------------------------------------------------*/
Again, this is not consistent with the book.
Problem 5 has five input, the C++ code has nine!
"""
def set_specific_conditions(state, config):
    if problem == 'drop' or problem == 'dam':
     break
     """
     /*-----------------------------------------------------------*/
     /* Driven Cavity: U = 1.0 at the upper boundary              */
     /*-----------------------------------------------------------*/
     """
    elif config.problem == 'dcavity':
        for row_count in xrange(0, config.imax+1):
            state.x_grid_vel[row_count, config.jmax+1] = 2. - 
                                                         state.x_grid_vel[row_count, jmax]
            break
    """"
    /*-----------------------------------------------------------------*/
    /* Flow past a backward facing step, with or without free boundary */
    /*                  U = 1.0 at the left boundary                   */
    /*-----------------------------------------------------------------*/
    """"
    elif (config.problem == 'backstep') or (config.problem == 'wave'):
        for col_count in xrange(1 + config.jmax/2 , config.jmax+1):
            state.x_grid_vel[0, col_count] = 1.
            break
    """
    /*--------------------------------------------------------------*/
    /* Flow past an obstacle: U = 1.0 at left boundary              */
    /*--------------------------------------------------------------*/
    """
    elif (config.problem == 'plate') or (config.problem == 'circle'):
        state.y_grid_vel[0, 0] = 2 * config.init_y_vel_scalar - state.x_grid_vel[1, 0]
        state.x_grid_vel[0, 1:config.jmax+1] = config.init_x_vel_scalar
        state.y_grid_vel[0, 1:config.jmax+1] = 2 * config.init_y_vel_scalar - 
                                                   state.y_grid_vel[0, 1:config.jmax+1]
        break
    """"
    /*---------------------------------------------------------------------*/
    /* Inflow for injection molding: U = 1.0 in the mid of left boundary   */
    /*---------------------------------------------------------------------*/
    """"
    elif (config.problem == 'molding'):
        lowe_bound = int(floor(1 + .4 * config.jmax))
        upper_bound = int(floor(.6 * config.jmax))
        state.y_grid_vel[0,lowe_bound:upper_bound+1] =  1.
        break
    """
    /*------------------------------------------------------------------*/
    /* natural convection or fluidtrap: left T = 0.5 right T = -0.5     */
    /*                          upper and lower wall adiabatic          */
    /*------------------------------------------------------------------*/
    """
    elif (config.problem == 'convection') or (config.problem == 'fluidtrap'):
        state.temp[0, 0:config.jmax+2] = 2 * 0.5 - state.temp[1, 0:config.jmax+2] # left wall heated
        # right wall heated:
        state.temp[config.imax+1, 0:config.jmax+2] = -2 * 0.5 - 
                                                       state.temp[config.imax, 0:config.jmax+2] 

        state.temp[0:config.imax+2, 0] = state.temp[0:config.imax+2, 1]
        state.temp[0:config.imax+2, 0] = state.temp[0:config.imax+2, config.jmax] # adiabatic walls
        break
    """
    /*----------------------------------------------------*/
    /* Rayleigh-Benard flow: top T = -0.5 bottom T = 0.5  */
    /*                       left and right adiabatic     */
    /*----------------------------------------------------*/
    """
    elif config.problem == 'rayleigh':
        state.temp[0, 0:config.jmax+2] = state.temp[1, 0:config.jmax+2];
        #  adiabatic walls
        state.temp[config.imax+1, 0:config.jmax+2] = state.temp[config.imax, 0:config.jmax+2]
        
        state.temp[0:config.imax+2, 0] = 2*(0.5) - state.temp[0:config.imax+2, 1] # lower wall heated
        # upper wall cooled:
        state.temp[0:config.imax+2, config.jmax+1] = 2 * (-0.5) - 
                                                     state.temp[0:config.imax+2, config.jmax]
        break
    else:
       print ('Problem {} not defined!'.format(config.problem))
       break
    return state
