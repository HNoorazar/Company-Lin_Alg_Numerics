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
        self.delta_x = 0.0 # length delta_x of one cell in x-direction
        self.delta_y = 0.0 # length delta_x of one cell in y-direction

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

    def printOut(self):
        print("=================")
        print("StaticParameters:")
        print("=================")
        print ("Geometry Data")
        print("domain size in x direction              = "+str(self.xlength))
        print("domain size in y direction              = "+str(self.ylength))
        print("No. of interior cells in x-direction    = " + str(self.imax))        
        print("No. of interior cells in y-direction    = "+str(self.jmax))
        print("length del_x of 1 cell in x-direction   = "+str(self.delta_x))
        print("length del_x of 1 cell in x-direction   = "+str(self.delta_y))

        print ("Time Criteria")        
        print("Current Time          = "+str(self.current_time))
        print("Final Time            = "+str(self.final_t))
        print("Time Step Size        = "+str(self.deta_t))
        print("delta_t Safety Factor = "+str(self.safety_tau))

        print ("Pressure Iteration Data")
        print ("Maximal No. of pressure iterations in one time step = " + str(self.iteration_max))
        print ("norm of pressure equation residual = " + str(self.iteration_count)) # Do we need to remove this?
        print ("SOR iteration counter = " + str(self.res_norm)) # Do we need to remove this?
        print ("stopping tolerance for pressure iteration = " + str(self.stop_toler))
        print ("relaxation parameter omega for SOR iteration = " + str(self.gamma))        
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
        config.set('parameters', 'safety_tau', self.safety_tau)           # 10
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
        self.xlength = config.getfloat('parameters', 'xlength') # 1
        self.ylength = config.getfloat('parameters', 'ylength') # 2
        self.imax = config.getfloat('parameters', 'imax')       # 3
        self.jmax = config.getfloat('parameters', 'jmax')       # 4
        self.delta_x = config.getfloat('parameters', 'xlength') / config.getfloat('parameters', 'imax') # 5
        self.delta_y = config.getfloat('parameters', 'ylength') / config.getfloat('parameters', 'jmax') # 6
        self.current_time = config.getint('parameters', 'current_time') # 7
        self.final_t = config.getint('parameters', 'final_t') # 8
        self.deta_t = config.getint('parameters', 'deta_t')   # 9
        self.safety_tau = config.getint('parameters', 'safety_tau')           # 10
        self.iteration_max = config.getint('parameters', 'iteration_max')     # 11
        self.iteration_count = config.getint('parameters', 'iteration_count') # 12
        self.res_norm = config.getint('parameters', 'res_norm')       # 13
        self.stop_toler = config.getint('parameters', 'stop_toler')   # 14
        self.relax_param = config.getint('parameters', 'relax_param') # 15
        self.gamma = config.getint('parameters', 'gamma')             # 16

        








