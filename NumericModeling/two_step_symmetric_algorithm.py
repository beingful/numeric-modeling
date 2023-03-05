from conditions import *
from numpy import zeros
from parameters import *

class TwoStepSymmetricAlgorithm:
    __parameters = Parameters()
    __grid_nodes_values = zeros((__parameters.M + 1, __parameters.N + 1, __parameters.N + 1))
    __METHOD_NAME = 'Two step symmetric algorithm'

    def __init__(self, function_w, function_f):
        self.__w = function_w
        self.__f =  function_f

    def use_algorithm(self):
        conditions = Conditions(self.__w, self.__grid_nodes_values)
        conditions.set()

        for i in range(0, self.__parameters.M, 1):
            for j in range(1, self.__parameters.N, 1):
                for k in range(1, self.__parameters.N, 1):
                    if (i + j + k) % 2 == 0:
                        self.__grid_nodes_values[i + 1, j, k] = \
                            self.__calculate_explicit(i, j, k)
            for j in range(1, self.__parameters.N, 1):
                for k in range(1, self.__parameters.N, 1):
                    if (i + j + k) % 2 == 1:
                        self.__grid_nodes_values[i + 1, j, k] = \
                            self.__calculate_implicit(i, j, k)

        return {'method_name': self.__METHOD_NAME,
                'method_results': self.__grid_nodes_values}

    def __calculate_explicit(self, t_step, x_step, y_step):
        w_x_y_t = self.__grid_nodes_values[t_step, x_step, y_step]
        w_x_next_y_t = self.__grid_nodes_values[t_step, x_step + 1, y_step]
        w_x_prev_y_t = self.__grid_nodes_values[t_step, x_step - 1, y_step]
        w_x_y_next_t = self.__grid_nodes_values[t_step, x_step, y_step + 1]
        w_x_y_prev_t = self.__grid_nodes_values[t_step, x_step, y_step - 1]

        f_x_y_t = self.__f.calculate(self.__parameters.x0 + x_step * self.__parameters.hx,
                                     self.__parameters.y0 + y_step * self.__parameters.hy,
                                     self.__parameters.t0 + t_step * self.__parameters.tau)

        return (w_x_y_t + self.__parameters.tau *
                (f_x_y_t + self.__w.a *
                 ((w_x_next_y_t - 2 * w_x_y_t + w_x_prev_y_t) / self.__parameters.hx ** 2
                 + (w_x_y_next_t - 2 * w_x_y_t + w_x_y_prev_t) / self.__parameters.hy ** 2)))

    def __calculate_implicit(self, t_step, x_step, y_step):
        w_x_y_t = self.__grid_nodes_values[t_step, x_step, y_step]
        w_x_next_y_t_next = self.__grid_nodes_values[t_step + 1, x_step + 1, y_step]
        w_x_prev_y_t_next = self.__grid_nodes_values[t_step + 1, x_step - 1, y_step]
        w_x_y_next_t_next = self.__grid_nodes_values[t_step + 1, x_step, y_step + 1]
        w_x_y_prev_t_next = self.__grid_nodes_values[t_step + 1, x_step, y_step - 1]

        f_x_y_t_next = self.__f.calculate(self.__parameters.x0 + x_step * self.__parameters.hx,
                                          self.__parameters.y0 + y_step * self.__parameters.hy,
                                          self.__parameters.t0 + (t_step + 1) * self.__parameters.tau)

        return ((w_x_y_t + self.__parameters.tau *
                 (f_x_y_t_next + self.__w.a *
                  ((w_x_next_y_t_next + w_x_prev_y_t_next) / self.__parameters.hx ** 2
                   + (w_x_y_next_t_next + w_x_y_prev_t_next) / self.__parameters.hy ** 2)))
                / (1 + 2 * self.__parameters.tau * (1 / self.__parameters.hx ** 2 + 1 / self.__parameters.hy ** 2)))
