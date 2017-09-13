import numpy as np
import machinery.core as my_core
import machinery.config as my_cfg
import machinery.surface as my_surf
import machinery.boundary as my_bnd
import machinery.definitions as my_def


config = my_cfg.staticParameters()
config.readFromFile('staticParameters.cfg')
config.printOut()

