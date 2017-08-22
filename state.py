import coupling as coup
import numpy as np

class GridState:
    def __init__(self, x_velocities, y_grid_vel, pressures):
        self.x_grid_vel = x_velocities
        self.y_grid_vel = y_grid_vel
        self.pressures = pressures