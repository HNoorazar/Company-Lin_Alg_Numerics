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


    
    











