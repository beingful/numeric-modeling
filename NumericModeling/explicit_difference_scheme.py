from conditions import *
from function import *
from numpy import zeros

class ExplicitDifferenceScheme:
    __N, __M = 5, 20
    __x0 = __y0 = __t0 = 0
    __l1, __l2, __T = 5, 5, 1
    __hx, __hy, __tau = __l1/__N, __l2/__N, __T/__M
    __METHOD_NAME = 'Explicit difference scheme'

    def __init__(self, function_w, function_f):
        self.__w = function_w
        self.__f =  function_f
        self.__grid_nodes_values = zeros((self.__M + 1, self.__N + 1, self.__N + 1))

    def use_algorithm(self):
        conditions = Conditions(self.__w, self.__grid_nodes_values, self.__x0, self.__y0, self.__t0,
                                self.__l1, self.__l2, self.__hx, self.__hy, self.__tau)
        conditions.set()

        for i in range(0, self.__M, 1):
            for j in range(1, self.__N, 1):
                for k in range(1, self.__N, 1):
                    result = self.__calculate(i, j, k)
                    self.__grid_nodes_values[i + 1, j, k] = result

        return {'method_name': self.__METHOD_NAME,
                'method_results': self.__grid_nodes_values}

    def __calculate(self, t_step, x_step, y_step):
        w_x_y_t = self.__grid_nodes_values[t_step, x_step, y_step]
        w_x_next_y_t = self.__grid_nodes_values[t_step, x_step + 1, y_step]
        w_x_prev_y_t = self.__grid_nodes_values[t_step, x_step - 1, y_step]
        w_x_y_next_t = self.__grid_nodes_values[t_step, x_step, y_step + 1]
        w_x_y_prev_t = self.__grid_nodes_values[t_step, x_step, y_step - 1]

        f_x_y_t = self.__f.calculate(self.__x0 + x_step * self.__hx,
                                     self.__y0 + y_step * self.__hy,
                                     self.__t0 + t_step * self.__tau)

        return (w_x_y_t + self.__tau *
                (f_x_y_t + self.__w.a * 
                 ((w_x_next_y_t - 2 * w_x_y_t + w_x_prev_y_t) / self.__hx ** 2
                  + (w_x_y_next_t - 2 * w_x_y_t + w_x_y_prev_t) / self.__hy ** 2)))
