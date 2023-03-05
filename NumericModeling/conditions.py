class Conditions:
    def __init__(self, function_w, grid_nodes_values, x0, y0, t0, x_last, y_last, x_step, y_step, t_step):
        self.__w = function_w
        self.__grid_nodes_values = grid_nodes_values
        self.__x0, self.__y0, self.__t0 = x0, y0, t0
        self.__x_last, self.__y_last = x_last, y_last
        self.__x_step, self.__y_step, self.__t_step = x_step, y_step, t_step

    def set(self):
        self.__set_initial_conditions()
        self.__set_boundary_conditions()

    def __set_initial_conditions(self):
        for i in range(0, len(self.__grid_nodes_values[0]), 1):
            x = i * self.__x_step
            for j in range(0, len(self.__grid_nodes_values[i]), 1):
                y = j * self.__y_step
                self.__grid_nodes_values[0, i, j] = self.__w.calculate(x, y, self.__t0)

    def __set_boundary_conditions(self):
        for i in range(1, len(self.__grid_nodes_values), 1):
            t = i * self.__t_step
            for x in (self.__x0,  self.__x_last):
                x_index = int(x / self.__x_step)
                for k in range(0, len(self.__grid_nodes_values[i]), 1):
                    y = k * self.__y_step
                    self.__grid_nodes_values[i, x_index, k] = self.__w.calculate(x, y, t)
            for j in range(0, len(self.__grid_nodes_values[i]), 1):
                x = j * self.__x_step
                for y in (self.__y0, self.__y_last):
                    y_index = int(y / self.__y_step)
                    self.__grid_nodes_values[i, j, y_index] = self.__w.calculate(x, y, t)