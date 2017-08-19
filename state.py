import coupling as coup
import numpy as np

class GridState:
    def __init__(self, x_velocities, y_velocities, pressures):
        self.x_velocities = x_velocities
        self.y_velocities = y_velocities
        self.pressures = pressures