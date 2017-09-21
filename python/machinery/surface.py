import numpy as np
import definitions as defn

"""
/*------------------------------------------------------------------------*/
/* Initialize particles for free boundary problems to define fluid domain */
/*------------------------------------------------------------------------*/

Is the following a dunction or not? 
What the hell is the name of the function?
name of funcion is INIT_PARTICLES and it is of type structure named particleline?
And it is also a pointer? :'(

how in declarign a function we have pointers?! :(
struct particleline *INIT_PARTICLES (int *N,int imax,int jmax,
                                     REAL delx,REAL dely,
                                     int ppc,char *problem,REAL **U,REAL **V)
"""

def particleline(config, state, N, ppc, Particlelines):
    # Gr"o"sen f"ur Tropfen 
    height=0.
    rad=0.
    mpx=0.
    mpy=0.
    vstart=0.
          #  Initialization of some parameters
    if (config.problem == 'dam'):
        N = 1
    if (config.problem == 'dam'):
        N = 2
        height = 1./2. * config.jmax * config.delta_y  # Height of the basin
        rad    = 0.1 * config.jmax * config.delta_y    #  Radius of the drop
        mpx    = 0.5 * config.imax * config.delta_x    # Mid point of the drop
        mpy    = 2./3. * config.jmax * config.delta_y
        vstart = -2.0                                  # Initial velocity of the drop
    if What the hell is that? #line 34 of Surface.c:
        print "No Memory"
        break
    Particlelines -= 1 # Particlelines from 1 to N 
    for ii in xrange(1, int(N)+1) # Line 42 of surface.c
        Particlelines[ii].length = 0
        Particlelines[ii].Particles = PARTALLOC(-1.,-1.) # what the hell is this?


                       #  Set the particles 
    for ii in xrange(1, config.imax):
        for jj in xrange(1, config.jmax):
            for ip in xrange(1, ppc+1):
                x = (ii-1) * config.delta_x + (ip-.5)/(float(ppc)) * config.delta_x
                for jp in xrange(1, ppc+1):
                    y = (jj-1) * config.delta_y + (jp-.5)/(float(ppc)) * config.delta_y
                    if(problem == 'dam'):
                        if (x < 0.2 * config.imax * config.jmax):
                            SET_PART(Particlelines[1],x,y)
                    
                    if (problem == 'drop'):
                        if (y<height):
                            SET_PART(Particlelines[1],x,y)
                        elif ((x-mpx)*(x-mpx)+(y-mpy)*(y-mpy) <= rad**2)
                            SET_PART(Particlelines[2],x,y)
                            state.y_grid_vel[ii, jj] = vstart
    return Particlelines

def SET_PART(my_class, x, y):
    # what the hell is at the beginning of the code?
    # the pointer *part, how do you pass that to the function,
    # if it is not amongst the inputs?
    #  struct particle *part
    
    # generate a class of particle, and call it part:
    part = particle_class()
    part.x = 
    





"""
/*---------------------------------------------------------*/
/* Set boundary values at free surface                     */
/*---------------------------------------------------------*/
"""
def SET_UVP_SURFACE(condig, state):
    for jj in xrange(1, config.jmax+1):
        for ii in xrange(1, config.imax):
            if (state.flag[ii, jj]) & defn.C_E and (state.flag[ii+1, jj] & defn.C_E):
                state.x_grid_vel[ii,jj] = 0.0

    for jj in xrange(1, config.jmax)
        for ii in xrange(1, config.imax+1)
            if ((state.flag[ii, jj] & defn.C_E) && state.flag[ii, jj+1] & defn.C_E)
                state.y_grid_vel[ii,jj] = 0.0
    for jj in xrange(1, config.jmax+1):
        for ii in xrange(1, config.imax+1):
            # treat only surface cells 
            if ( (not(state.flag[ii, jj] & defn.C_E)) or (state.flag[ii, jj] < 0x0100) ):
                expression = state.flag[ii, jj] & defn.C_NSWO
                # /* mask NSWO_E=0x0f00    */
	            # /* filters surface cells */
	            dxdy = config.detla_x / config.delta_y
	            dydx = config.detla_y / config.delta_x
                if expression == defn.C_N:
                    state.y_grid_vel[ii, jj] = state.y_grid_vel[ii, jj-1] - 
                                               dydx * 
                                               (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii-1, jj])
                    if (state.flag[ii-1, jj+1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj] - 
                                                       dydx * 
                                                       (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii-1, jj])
                ##
                elif expression == defn.C_S:
                    state.y_grid_vel[ii, jj-1] = state.y_grid_vel[ii, jj] + 
                                                 dydx * 
                                                 (state.x_grid_vel[ii, jj] - state.y_grid_vel[ii-1, jj])
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj] + 
                                                       dydx * 
                                                       (state.y_grid_vel[ii, jj-1] - state.x_grid_vel[ii-1, jj-1])
                ##
                elif expression == defn.C_O:
                    state.x_grid_vel[ii, jj] = state.x_grid_vel[ii-1, jj] - 
                                               dxdy * 
                                               (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
                    if (state.flag[ii+1, jj-1] & defn.C_E):  
                        state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1] - 
                                                       dxdy * 
                                                       (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii, jj-1]);
                ##
                elif expression == defn.C_W:
                    state.x_grid_vel[ii-1, jj] = state.x_grid_vel[ii, jj] + 
                                                 dxdy * 
                                                 (state.y_grid_vel[ii, jj]-state.y_grid_vel[ii, jj-1])
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1] + 
                                                       dxdy * 
                                                       (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii-1, jj-1])
                ##
                
                elif expression == defn.C_NO:
                    state.x_grid_vel[ii, jj] = state.y_grid_vel[ii-1, jj]
                    state.y_grid_vel[ii, jj] = state.y_grid_ve[ii, jj-1]

                    if (state.flag[ii-1, jj+1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj] - 
                                                     dydx * 
                                                     (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii-1, jj])
                    if (state.flag[ii+1, jj+1] & defn.C_E):
                        state.x_grid_vel[ii, jj+1] = state.x_grid_vel[ii, jj]
                        state.y_grid_vel[ii+1, jj] = state.y_grid_vel[ii, jj]
                    if (state.flag[ii+1, jj-1] & defn.C_E):
                           state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1] - 
                                                          dxdy * 
                                                          (state.x_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
                ##
                elif expression == defn.C_NW:
                    state.x_grid_vel[ii-1, jj] = state.x_grid_vel[ii, jj];
                    state.y_grid_vel[ii, jj]   = state.y_grid_vel[ii, jj-1]
                    if (sfate.flag[ii-1, jj+1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj]
                        state.y_grid_vel[ii-1, jj]   = state.y_grid_vel[ii, jj]
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1] + 
                                                     dxdy * 
                                                     (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii-1, jj-1]);
                ##
                elif expression == defn.C_SW:
                    state.x_grid_vel[ii-1, jj] = state.x_grid_vel[ii, jj]
                    state.y_grid_vel[ii, jj-1] = state.y_grid_vel[ii, jj]
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj]
                        state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1]
                ##
                elif expression == defn.C_SO:
                    state.x_grid_vel[ii, jj]   = state.x_grid_vel[ii-1, jj]
                    state.y_grid_vel[ii, jj-1] = state.y_grid_vel[ii, jj]
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj] + 
                                                       dydx * 
                                                       (state.y_grid_vel[ii, jj-1] - state.y_grid_vel[ii-1, jj-1]);
                        if (state.flag[ii+1, jj-1] & defn.C_E):
                           state.x_grid_vel[ii, jj-1]   = state.x_grid_vel[ii, jj]
                           state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1]
                ##
                elif (expression == defn.C_WO):
                    state.x_grid_vel[ii, jj]   += config.delta_t * config.GX
                    state.x_grid_vel[ii-1, jj] += config.delta_t * config.GX
                    
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1] + 
                                                       dxdy * 
                                                       (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii-1, jj-1]);
                    if (state.flag[ii+1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1] - 
                                                       dxdy * 
                                                       (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii, jj-1])
                ##
                elif (expression == defn.C_NS):
                    state.y_grid_vel[ii, jj]   += config.delta_t * config.GY
                    state.y_grid_vel[ii, jj-1] += config.delta_t * config.GY
                    if (state.flag[ii-1, jj+1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj] - 
                                                       dydx * 
                                                       (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii-1, jj])
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.x_grid_vel[[ii-1, jj-1] = state.x_grid_vel[ii-1, jj] +
                                                        dydx *
                                                        (state.y_grid_vel[ii, jj-1] - state.y_grid_vel[ii-1, jj-1])
		        ##
                elif (expression == defn.C_NWO):
                    state.y_grid_vel[ii, jj] = state.y_grid_vel[ii, jj-1] - 
                                               dydx * 
                                               (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii-1, jj])
                    state.x_grid_vel[ii, jj]   += config.delta_t * config.GX
                    state.x_grid_vel[ii-1, jj] += config.delta_t * config.GX
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj-1]  = state.y_grid_vel[ii, jj-1] + 
                                                        dxdy * 
                                                        (state.x_grid_vel[ii-1, jj] - state.x_grid_vel[ii-1, jj-1])
                        
                    if (state.flag[ii+1, jj-1] & defn.C_E)):
                        state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1] - 
                                                       dxdy * 
                                                       (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii, jj-1])
                        
                    if (state.flag[ii-1, jj+1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj]    = state.y_grid_vel[ii, jj]
                        state.x_grid_vel[ii-1, jj+1]  = state.x_grid_vel[ii-1, jj]

                    if (state.flag[ii+1, jj+1] & defn.C_E):
                        state.y_grid_vel[ii+1, jj]  = state.y_grid_vel[ii, jj]
                        state.x_grid_vel[ii, jj+1]  = state.x_grid_vel[ii, jj]
                ##
                elif (expression == defn.C_NSW):
                    state.x_grid_vel[ii-1, jj] = state.x_grid_vel[ii, jj] + 
                                                 dxdy * 
                                                 (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
                    state.y_grid_vel[ii, jj]   += config.delta_t * config.GY;
                    state.y_grid_vel[ii, jj-1] += config.delta_t * config.GY;
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj-1]  = state.y_grid_vel[ii, jj-1];
                        state.x_grid_vel[ii-1, jj-1]  = state.x_grid_vel[ii-1, jj];
                 
                    if (state.flag[ii-1, jj+1] & defn.C_E):
                        state.y_grid_vel[ii-1, jj]    = state.y_grid_vel[ii, jj];
                        state.x_grid_vel[ii-1, jj+1]  = state.x_grid_vel[ii-1, jj];
                ##
                elif expression == defn.C_SWO:
		            state.y_grid_vel[ii, jj-1] = state.y_grid_vel[ii, jj] + 
		                                          dydx * 
		                                          (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii-1, jj])
                    state.x_grid_vel[ii, jj]   += config.delta_t * config.GX
                    state.x_grid_vel[ii-1, jj] += config.delta_t * config.GX
                    
                    if (state.flag[ii-1, jj-1] & defn.C_E):
                        state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj]
                        state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1]
                    
                    if (state.flag[ii+1, jj-1] & defn.C_E):
                        state.x_grid_vel[ii, jj-1] = state.x_grid_vel[ii, jj]
                        state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1]

                ##
                elif expression == defn.C_NSO:
		            state.x_grid_vel[ii, jj] = state.x_grid_vel[ii-1, jj] - 
		                                       dxdy *
		                                       ( state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
		            state.y_grid_vel[ii, jj]   += config.delta_t * config.GY
		            state.y_grid_vel[ii, jj-1] += config.delta_t * config.GY
		            if (state.flag[ii-1, jj+1] & defn.C_E):
		                state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj] - 
		                                               dydx *
		                                               (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii-1, jj] )
		            if (state.flag[ii-1, jj-1] & defn.C_E):
		                state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj] + 
		                                               dydx *
		                                               (state.y_grid_vel[ii, jj-1] - state.y_grid_vel[ii-1, jj-1])
		            if (state.flag[ii+1, jj-1] & defn.C_E):
		                state.x_grid_vel[ii, jj-1] = state.x_grid_vel[ii, jj]
		                state.y_grid.vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1]
		            if (state.flag[ii+1, jj+1] & defn.C_E):
		                state.x_grid_vel[ii, jj+1] = state.x_grid_vel[ii, jj]
		                state.y_grid_vel[ii+1, jj] = state.y_grid_vel[ii, jj]
                ##
                elif (expression == defn.C_NSWO):
                    state.x_grid_vel[ii, jj]   += config.delta_t * config.GX
		            state.x_grid_vel[ii-1, jj] += config.delta_t * config.GX
		            state.y_grid_vel[ii, jj]   += config.delta_t * config.GY
		            state.y_grid_vel[ii, jj-1] += config.delta_t * config.GY
		            if (state.flag[ii-1, jj+1] & defn.C_E):
		                state.x_grid_vel[ii-1, jj+1] = state.x_grid_vel[ii-1, jj]
		                state.y_grid_vel[ii-1, jj]   = state.y_grid_vel[ii, jj]
		            
		            if (state.flag[ii+1, jj+1] & defn.C_E):
		                state.x_grid_vel[ii, jj+1] = state.x_grid_vel[ii, jj]
		                state.y_grid_vel[ii+1, jj] = state.y_grid_vel[ii, jj]
		            
		            if (state.flag[ii-1, jj-1] & defn.C_E):
		                state.x_grid_vel[ii-1, jj-1] = state.x_grid_vel[ii-1, jj]
		                state.y_grid_vel[ii-1, jj-1] = state.y_grid_vel[ii, jj-1]
		            
		            if (state.flag[ii+1, jj-1] & defn.C_E):
		                state.x_grid_vel[ii, jj-1]   = state.x_grid_vel[ii, jj]
		                state.y_grid_vel[ii+1, jj-1] = state.y_grid_vel[ii, jj-1]
                ##
                else:
                    break
    #  Second loop for pressure boundary values
    for jj in xrange(1, config.jmax+1):
        for ii in xrange(1, config.imax+1):
            if (not ((state.flag[ii, jj] & defn.C_E) or (state.flag[ii, jj] < 0x0100 ))):
                expression = state.flag[ii, jj] & defn.C_NSWO
                if expression == defn.C_N:
                    state.pressures[ii, jj] = (2. / (config.Ray_no * config.delta_y)) * 
                                              (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
                elif expression == defn.C_S:
                    state.pressures[ii, jj] = (2. / (config.Ray_no * config.delta_y)) * 
                                              (state.y_grid_vel[ii, jj] - state.y_grid_vel[ii, jj-1])
                elif expression == defn.C_O:
                    state.pressures[ii, jj] = ( 2. / (config.Ray_no * config.delta_x)) *
                                              (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii-1, jj])
                elif expression == defn.C_W:
                    state.pressures[ii, jj] = ( 2. / (config.Ray_no * config.delta_x)) *
                                              (state.x_grid_vel[ii, jj] - state.x_grid_vel[ii-1, jj])
                    
                elif expression == defn.C_NO:
                    term1 = state.x_grid_vel[ii, jj]   + 
                            state.x_grid_vel[ii-1, jj] - 
                            state.x_grid_vel[ii, jj-1] - 
                            state.x_grid_vel[ii-1, jj-1]
                    term1 = term1 / config.delta_y
                    
                    term2 = state.y_grid_vel[ii, jj]   + 
                            state.y_grid_vel[ii, jj-1] -
                            state.y_grid_vel[ii-1, jj] -
                            state.y_grid_vel[ii-1, jj-1]
                    term2 = term2  / config.delta_x
                    
                    term = term1 + term2
                    state.pressures[ii, jj] = (1./ (2 * config.Ray_no)) * term
                    
                elif expression == C_NW:
                    term1 = state.x_grid_vel[ii, jj]   + 
                            state.x_grid_vel[ii-1, jj] - 
                            state.x_grid_vel[ii, jj-1] - 
                            state.x_grid_vel[ii-1, jj-1]
                    term1 = term1 / config.delta_y
                    
                    term2 = state.y_grid_vel[ii+1, jj]   + 
                            state.y_grid_vel[ii+1, jj-1] - 
                            state.y_grid_vel[ii, jj]     - 
                            state.y_grid_vel[ii, jj-1]
                    term2 = term2 / config.delta_x

                    term = term1 + term2
                    state.pressures[ii, jj] = (-1./ (2 * config.Ray_no)) * term
                    
                elif expression == C_SW:
                    term1 = state.x_grid_vel[ii, jj+1]   + 
                            state.x_grid_vel[ii-1, jj+1] - 
                            state.x_grid_vel[ii, jj]     - 
                            state.x_grid_vel[ii-1, jj]
                    term1 = term1 / config.delta_y
                    
                    term2 = state.y_grid_vel[ii+1, jj]   + 
                            state.y_grid_vel[ii+1, jj-1] - 
                            state.y_grid_vel[ii, jj]     - 
                            state.y_grid_vel[ii, jj-1]
                    term2 = term2 / config.delta_x)
                    
                    term = term1 + term2
                    state.pressures[ii, jj] = (1./ (2. * config.Ray_no)) * term
                elif expression == C_SO:
                    term1 = state.x_grid_vel[ii, jj+1]    + 
                            state.x_grid_velU[ii-1, jj+1] - 
                            state.x_grid_vel[ii, jj]      - 
                            state.x_grid_vel[ii-1, jj])
                    term1 = term1 / config.delta_y
                    
                    term2 = state.y_grid_vel[ii, jj]   + 
                            state.y_grid_vel[ii, jj-1] -
                            state.y_grid_vel[ii-1, jj] -
                            state.y_grid_vel[ii-1, jj-1]
                    term2 = term2 / config.delta_x
                    
                    term = term1 + term2
                    state.pressures[ii, jj] = (-1./ (2 * config.Ray_no)) * term

                elif expression == C_WO:
                    state.pressures[ii, jj] = 0.
                elif expression == C_NS:
                    state.pressures[ii, jj] = 0.
                elif expression == C_NWO:
                    state.pressures[ii, jj] = 0.
                elif expression == C_NSW:
                    state.pressures[ii, jj] = 0.
                elif expression == C_SWO:
                    state.pressures[ii, jj] = 0.
                elif expression == C_NSO:
                    state.pressures[ii, jj] = 0.
                elif expression == C_NSWO:
                    state.pressures[ii, jj] = 0.
                else:
                    break
                
                    
                    
                              
                              
                    























