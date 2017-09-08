import coupling as coup
import numpy as np

class GridState:
    def __init__(self, x_velocities, y_grid_vel, pressures, ibound):
        self.x_grid_vel = x_velocities
        self.y_grid_vel = y_grid_vel
        self.pressures = pressures
        self.temp = temp
        self.flag = flag
        self.ibound = ibound
        self.F = F
        self.G = G
        self.current_time = current_time