from matplotlib.pyplot import figure, show
from numpy import zeros
from parameters import *

class HeatConductionPlot:
    __parameters = Parameters()
    __x_vector = zeros((__parameters.M + 1) * (__parameters.N + 1) ** 2)
    __y_vector = zeros((__parameters.M + 1) * (__parameters.N + 1) ** 2)
    __t_vector = zeros((__parameters.M + 1) * (__parameters.N + 1) ** 2)
    __w_vector = zeros((__parameters.M + 1) * (__parameters.N + 1) ** 2)

    def __init__(self, method_results, method_name):
        self.__method_results = method_results
        self.__method_name = method_name

    def build_plot(self):
        self.__set_vectors()
        self.__set_plot()

    def __set_vectors(self):
        counter = 0

        for i in range(0, self.__parameters.M + 1, 1):
            self.__t_vector[counter:counter + (self.__parameters.N + 1) ** 2] = \
                self.__parameters.t0 + i * self.__parameters.tau
            for j in range(0, self.__parameters.N + 1, 1):
                self.__x_vector[counter:counter + self.__parameters.N + 1] = \
                    self.__parameters.x0 + j * self.__parameters.hx
                for k in range(0, self.__parameters.N + 1, 1):
                    self.__y_vector[counter] = \
                        self.__parameters.y0 + k * self.__parameters.hy
                    self.__w_vector[counter] = self.__method_results[i][j][k]

                    counter += 1

    def __set_plot(self):
        canvas = figure()

        plot = canvas.add_subplot(111, projection='3d')
        plot.scatter(self.__x_vector, self.__y_vector, self.__t_vector, c = self.__w_vector, cmap='coolwarm')
        plot.set_xlabel('x')
        plot.set_ylabel('y')
        plot.set_zlabel('t')
        plot.set_title(self.__method_name)

        show()
