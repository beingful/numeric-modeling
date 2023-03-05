from function import *
from conditions import *
from numpy import zeros
from tridiagonal_matrix_algorithm import *
from parameters import *

class AlternatingDirectionImplicitMethod:
    __parameters = Parameters()
    __grid_nodes_values = zeros((__parameters.M + 1, __parameters.N + 1, __parameters.N + 1))
    __METHOD_NAME = 'Alternating direction implicit method'

    def __init__(self, function_w, function_f):
        self.__w = function_w
        self.__f =  function_f
        self.__alpha_x = self.__gamma_x = self.__w.a * self.__parameters.tau / (2 * self.__parameters.hx ** 2)
        self.__beta_x = - (1 + self.__w.a * self.__parameters.tau / self.__parameters.hx ** 2)
        self.__alpha_y = self.__gamma_y = self.__w.a * self.__parameters.tau / (2 * self.__parameters.hy ** 2)
        self.__beta_y = - (1 + self.__w.a * self.__parameters.tau / self.__parameters.hy ** 2)

    def use_algorithm(self):
        conditions = Conditions(self.__w, self.__grid_nodes_values)
        conditions.set()

        tridiagonal_matrix = zeros((self.__parameters.N - 1, self.__parameters.N - 1))
        solution_vector = zeros(self.__parameters.N - 1)
        tridiagonal_matrix_results = zeros((self.__parameters.N - 1, self.__parameters.N - 1))

        for i in range(0, self.__parameters.M, 1):
            for j in range(1, self.__parameters.N, 1):
                for k in range(1, self.__parameters.N, 1):
                    self.__calculate_for_half_step(i, k, j, tridiagonal_matrix, solution_vector)
                tridiagonal_matrix_results[j - 1] = \
                    self.__use_tridiagonal_matrix_algorithm(tridiagonal_matrix, solution_vector)
            for j in range(1, self.__parameters.N, 1):
                for k in range(1, self.__parameters.N, 1):
                    self.__calculate_for_full_step(i, j, k, tridiagonal_matrix, solution_vector, tridiagonal_matrix_results[k - 1])
                self.__grid_nodes_values[i + 1][j][1:self.__parameters.N] = \
                    self.__use_tridiagonal_matrix_algorithm(tridiagonal_matrix, solution_vector)

        return {'method_name': self.__METHOD_NAME,
                'method_results': self.__grid_nodes_values}

    def __calculate_for_half_step(self, t_step, x_step, y_step, tridiagonal_matrix, solutions_vector):
        parameters_vector = zeros(self.__parameters.N - 1)

        if x_step != 1:
            parameters_vector[x_step - 2] = self.__alpha_x
        if x_step != self.__parameters.N - 1:
            parameters_vector[x_step] = self.__gamma_x

        parameters_vector[x_step - 1] = self.__beta_x

        w_x_y_t = self.__grid_nodes_values[t_step, x_step, y_step]
        w_x_y_next_t = self.__grid_nodes_values[t_step, x_step, y_step + 1]
        w_x_y_prev_t = self.__grid_nodes_values[t_step, x_step, y_step - 1]

        x = self.__parameters.x0 + x_step * self.__parameters.hx
        y = self.__parameters.y0 + y_step * self.__parameters.hy
        t_half = self.__parameters.t0 + (t_step + 1/2) * self.__parameters.tau

        f_x_y_t_half = self.__f.calculate(x, y, t_half)

        solution = (self.__parameters.tau / 2 *
                    (- f_x_y_t_half - self.__w.a *
                     ((w_x_y_next_t - 2 * w_x_y_t + w_x_y_prev_t) / self.__parameters.hy ** 2))
                    - w_x_y_t)

        if x_step == 1:
            w_in_left_bound = self.__w.calculate(self.__parameters.x0, y, t_half)
            solution -= self.__alpha_x * w_in_left_bound
        elif x_step == self.__parameters.N - 1:
            w_in_right_bound = self.__w.calculate(self.__parameters.l1, y, t_half)
            solution -= self.__gamma_x * w_in_right_bound
        
        tridiagonal_matrix[x_step - 1] = parameters_vector
        solutions_vector[x_step - 1] = solution

    def __calculate_for_full_step(self, t_step, x_step, y_step, tridiagonal_matrix,
                                  solutions_vector, half_step_tridiagonal_matrix_results):
        vector_with_parameters = zeros(self.__parameters.N - 1)

        if y_step != 1:
            vector_with_parameters[y_step - 2] = self.__alpha_y
        if y_step != self.__parameters.N - 1:
            vector_with_parameters[y_step] = self.__gamma_y

        vector_with_parameters[y_step - 1] = self.__beta_y

        x = self.__parameters.x0 + x_step * self.__parameters.hx
        y = self.__parameters.y0 + y_step * self.__parameters.hy
        t_half = self.__parameters.t0 + (t_step + 1/2) * self.__parameters.tau

        if x_step == 1:
            w_x_prev_y_t_half = self.__w.calculate(self.__parameters.x0, y, t_half)
            w_x_next_y_t_half = half_step_tridiagonal_matrix_results[x_step]
        elif x_step == self.__parameters.N - 1:
            w_x_prev_y_t_half = half_step_tridiagonal_matrix_results[x_step - 2]
            w_x_next_y_t_half = self.__w.calculate(self.__parameters.l1, y, t_half)
        else:
            w_x_prev_y_t_half = half_step_tridiagonal_matrix_results[x_step - 2]
            w_x_next_y_t_half = half_step_tridiagonal_matrix_results[x_step]

        w_x_y_t_half = half_step_tridiagonal_matrix_results[x_step - 1]
        f_x_y_t_half = self.__f.calculate(x, y, t_half)

        solution = (self.__parameters.tau / 2 *
                    (- f_x_y_t_half - self.__w.a *
                     ((w_x_next_y_t_half - 2 * w_x_y_t_half + w_x_prev_y_t_half)
                      / self.__parameters.hx ** 2)) - w_x_y_t_half)

        if y_step == 1:
            w_in_left_bound = self.__grid_nodes_values[t_step, x_step, y_step - 1]
            solution -= self.__alpha_y * w_in_left_bound
        elif y_step == self.__parameters.N - 1:
            w_in_right_bound = self.__grid_nodes_values[t_step, x_step, y_step + 1]
            solution -= self.__gamma_y * w_in_right_bound

        tridiagonal_matrix[y_step - 1] = vector_with_parameters
        solutions_vector[y_step - 1] = solution

    def __use_tridiagonal_matrix_algorithm(self, tridiagonal_matrix, solutions):
        main_diag = zeros(self.__parameters.N - 1)
        lower_diag = upper_diag = zeros(self.__parameters.N - 2)

        for i in range(0, self.__parameters.N - 1, 1):
            for j in range(0, self.__parameters.N - 1, 1):
                if i == j:
                    main_diag[i] = tridiagonal_matrix[i][j]
                elif i == j - 1:
                    upper_diag[i] = tridiagonal_matrix[i][j]
                elif i == j + 1:
                    lower_diag[j] = tridiagonal_matrix[i][j]
        
        algorithm = TridiagonalMatrixAlgorithm(main_diag, lower_diag, upper_diag, solutions)

        return algorithm.use_algorithm()