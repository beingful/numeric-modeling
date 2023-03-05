from parameters import *

class Conditions:
    __parameters = Parameters()

    def __init__(self, function_w, grid_nodes_values):
        self.__w = function_w
        self.__grid_nodes_values = grid_nodes_values

    def set(self):
        self.__set_initial_conditions()
        self.__set_boundary_conditions()

    def __set_initial_conditions(self):
        for i in range(0, self.__parameters.N + 1, 1):
            x = self.__parameters.x0 + i * self.__parameters.hx
            for j in range(0, self.__parameters.N + 1, 1):
                y = self.__parameters.y0 + j * self.__parameters.hy
                self.__grid_nodes_values[0, i, j] = self.__w.calculate(x, y, self.__parameters.t0)

    def __set_boundary_conditions(self):
        for i in range(1, self.__parameters.M + 1, 1):
            t = self.__parameters.t0 + i * self.__parameters.tau
            for x in (self.__parameters.x0, self.__parameters.l1):
                x_index = int(x / self.__parameters.hx)
                for k in range(0, self.__parameters.N + 1, 1):
                    y = self.__parameters.y0 + k * self.__parameters.hy
                    self.__grid_nodes_values[i, x_index, k] = self.__w.calculate(x, y, t)
            for j in range(0, self.__parameters.N + 1, 1):
                x = self.__parameters.x0 + j * self.__parameters.hx
                for y in (self.__parameters.y0, self.__parameters.l2):
                    y_index = int(y / self.__parameters.hy)
                    self.__grid_nodes_values[i, j, y_index] = self.__w.calculate(x, y, t)