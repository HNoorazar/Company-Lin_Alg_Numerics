class staticParameters:
    def __init__(self):
        """
        Create static parameter object with default parameters.
        """
        # Geometry data
        self.xlength = 0. # domain size in x direction
        self.ylength = 0. # domain size in y direction
        self.imax = 0 # number of interior cells in x-direction
        self.jmax = 0 # number of interior cells in y-direction
        self.delta_x = 0. # length delta_x of one cell in x-direction
        self.delta_y = 0. # length delta_x of one cell in y-direction

        # time criteria
        self.current_time = 0. # current time. [Do we need to remove this from here?]
        self.final_t = 0. # Final time
        self.deta_t = 0. # time step size
        self.safety_tau = 0. # Safety factor for time step size control \tau

        # Pressure-iteration data
        self.iteration_max = 0 # Maximal number of pressure iterations in one time step
        self.iteration_count = 0  # SOR iteration counter [Do we need to remove this?]
        self.res_norm = 0. # norm of pressure equation residual [Do we need to remove this?]
        self.stop_toler = 0. # stopping tolerance eps for pressure iteration
        self.relax_param = 0. # relaxation parameter omega for SOR iteration
        self.gamma = 0. # upwind differencing factor \gamma

        # Problem dependent data
        self.init_x_vel_scalar = 0. 
        self.init_y_vel_scalar = 0. 
        self.init_press_scalar = 0. 
        self.Ray_no = 0.
        self.wW = 0
        self.wE = 0
        self.wN = 0
        self.wS = 0
        
    def printOut(self):
        print("=================")
        print("StaticParameters:")
        print("=================")
        print("Problem = " + str(self.problem))
        print ("Geometry Data")
        print("domain size in x direction            = " + str(self.xlength))
        print("domain size in y direction            = " + str(self.ylength))
        print("No. of interior cells in x-direction  = " + str(self.imax))        
        print("No. of interior cells in y-direction  = " + str(self.jmax))
        print("length del_x of 1 cell in x-direction = " + str(self.delta_x))
        print("length del_x of 1 cell in x-direction = " + str(self.delta_y))

        print ("Time Criteria")        
        print("Final Time     = "+str(self.final_t))
        print("Time Step Size = "+str(self.deta_t))
        print("Safety Factor  = "+str(self.safety_tau))
        
        print("del_trace  = "+str(self.del_trace))
        print("del_inj    = "+str(self.del_inj))
        print("del_streak = "+str(self.del_streak))
        print("del_vec    = "+str(self.del_vec))

        print("vec_file    = "+str(self.vec_file))
        print("trace_file  = "+str(self.trace_file))
        print("streak_file = "+str(self.streak_file))

        print("in_file = "+str(self.in_file))
        print("out_file = "+str(self.out_file))
        
        print("N = "+str(self.N))
        print("pos1x = "+str(self.pos1x))
        print("pos1y = "+str(self.pos1y))
        print("pos2x = "+str(self.pos2x))
        print("pos2y = "+str(self.pos2y))

        print ("Pressure Iteration Data")
        print ("Maximal No. of pressure iterations in one time step = " + str(self.iteration_max))
        print ("stopping tolerance for pressure iteration           = " + str(self.stop_toler))  # epsilon
        print ("relaxation parameter omega for SOR iteration        = " + str(self.relax_param)) # omega
        print ("gamma   = " + str(self.gamma))
        print ("p_bound = " + str(self.p_bound))
        print ("Reynolds number = " + str(self.Ray_no))
        print ("Pr   = " + str(self.Pr))
        print ("beta = " + str(self.beta))
        print ("GX   = " + str(self.GX))
        print ("GY   = " + str(self.GY))


        print("")

    def writeToFile(self, fname):
        config = configparser.RawConfigParser()

        config.add_section('parameters')
        config.set('parameters', 'xlength', self.xlength) # 1
        config.set('parameters', 'ylength', self.ylength) # 2
        config.set('parameters', 'imax', self.imax)       # 3
        config.set('parameters', 'jmax', self.jmax)       # 4
        config.set('parameters', 'delta_x', self.delta_x) # 5
        config.set('parameters', 'delta_y', self.delta_y) # 6
        config.set('parameters', 'current_time', self.current_time) # 7
        config.set('parameters', 'final_t', self.final_t) # 8
        config.set('parameters', 'deta_t', self.deta_t)   # 9
        config.set('parameters', 'safety_tau', self.safety_tau) # 10
        config.set('parameters', 'iteration_max', self.iteration_max)     # 11
        config.set('parameters', 'iteration_count', self.iteration_count) # 12
        config.set('parameters', 'res_norm', self.res_norm)       # 13
        config.set('parameters', 'stop_toler', self.stop_toler)   # 14
        config.set('parameters', 'relax_param', self.relax_param) # 15
        config.set('parameters', 'gamma', self.gamma) # 16


        with open(fname, 'wb') as configfile:
            config.write(configfile)

    def readFromFile(self, fname):
        config = configparser.RawConfigParser()
        config.read(fname)
        self.problem = config.get('parameters', 'problem') # 1 string?
        self.xlength = config.getfloat('parameters', 'xlength') # 2
        self.ylength = config.getfloat('parameters', 'ylength') # 3
        self.imax = config.getfloat('parameters', 'imax')       # 4
        self.jmax = config.getfloat('parameters', 'jmax')       # 5
        # should not we get rid of the two followings? 
        # These can be computed by division of xlength by imax!
        self.delta_x = config.getfloat('parameters', 'xlength') / config.getfloat('parameters', 'imax') # 6 
        self.delta_y = config.getfloat('parameters', 'ylength') / config.getfloat('parameters', 'jmax') # 7
        self.final_t = config.getint('parameters', 'final_t') # 8
        self.deta_t = config.getint('parameters', 'deta_t')   # 9
        self.safety_tau = config.getint('parameters', 'safety_tau') # 10
        self.del_trace = config.getfloat('parameters', 'del_trace') # 11
        self.del_inj = config.getfloat('parameters', 'del_inj')     # 12
        self.del_streak = config.getfloat('parameters', 'del_streak') # 13
        self.del_vec = config.getfloat('parameters', 'del_vec')    # 14
        self.vec_file = config.get('parameters', 'vec_file')       # 15 string?
        self.trace_file = config.get('parameters', 'trace_file')   # 16 string?
        self.streak_file = config.get('parameters', 'streak_file') # 17 string?
        self.in_file = config.get('parameters', 'in_file')         # 18 string?
        self.out_file = config.get('parameters', 'out_file')       # 19 string?
        self.N = config.getint('parameters', 'N') # 20
        self.pos1x = config.getfloat('parameters', 'pos1x') # 21
        self.pos1y = config.getfloat('parameters', 'pos1y') # 22
        self.pos2x = config.getfloat('parameters', 'pos2x') # 23
        self.pos2y = config.getfloat('parameters', 'pos2y') # 24
        self.iteration_max = config.getint('parameters', 'iteration_max') # 25
        self.stop_toler = config.getfloat('parameters', 'stop_toler') # This is epsilon in C++
        self.relax_param = config.getfloat('parameters', 'relax_param') # This is omega in C++
        self.gamma = config.getfloat('parameters', 'gamma')     # 26 This is upwind differencing parameter.
        self.p_bound = config.getint('parameters', 'p_bound')   # 27
        # What is res_norm? It is not in the init.c reading file. 
        # (I cannot remember where I saw this first))
        self.res_norm = config.getint('parameters', 'res_norm') # 28
        self.Ray_no = config.getfloat('parameters', Ray_no)     # 29. This is called Re in C++
        self.Pr = config.getfloat('parameters', Pr)             # 30
        self.beta = config.getfloat('parameters',beta) # 31
        self.GX = config.getfloat('parameters', GX)    # 32
        self.GY = config.getfloat('parameters', GY)    # 33
        self.init_x_vel_scalar = config.getfloat('parameters', init_x_vel_scalar) # 34
        self.init_y_vel_scalar = config.getfloat('parameters', init_y_vel_scalar) # 35
        self.init_press_scalar = config.getfloat('parameters', init_press_scalar) # 36
        self.initial_temp = config.getfloat('parameters', TI) # 37
        self.wW = config.getint('parameters', wW) # 38
        self.wE = config.getint('parameters', wE) # 39
        self.wN = config.getint('parameters', wN) # 40
        self.wS = config.getint('parameters', wS) # 41


