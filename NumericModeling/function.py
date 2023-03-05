from numpy import zeros
from partial_derivative import PartialDerivatives
from sympy import symbols, cos, exp, N

class FunctionW:
    a = 1
    __N, __M = 5, 20
    __x0 = __y0 = __t0 = 0
    __l1, __l2, __T = 5, 5, 1
    __hx, __hy, __tau = __l1/__N, __l2/__N, __T/__M
    __variables = symbols('x y t')

    def __init__(self):
        x, y, t, A, a, k1, k2, C1, C2 = symbols('x y t A a k1 k2 C1 C2')
        w = A * cos(k1 * x + C1) * cos(k2 * y + C2) * exp(-(k1 ** 2 + k2 ** 2) * a * t)
        constants = {A: 1, a: self.a, k1: 1, k2: 1, C1: 0, C2: 0}
        self.__w = w.subs(constants)

    def get_function(self):
        return self.__w

    def get_variables(self):
        return self.__variables

    def calculate(self, x_value, y_value, t_value):
        x, y, t = self.__variables
        return N(self.__w.subs({x: x_value, y: y_value, t: t_value}))

    def calculate_all(self):
        results = zeros((self.__M + 1, self.__N + 1, self.__N + 1))

        for i in range(0, self.__M + 1, 1):
            for j in range(0, self.__N + 1, 1):
                for k in range(0, self.__N + 1, 1):
                    results[i][j][k] = self.calculate(self.__x0 + j * self.__hx, \
                        self.__y0 + k * self.__hy, self.__t0 + i * self.__tau)

        return results

class FunctionF:
    def __init__(self, w_function):
        self.__variables = w_function.get_variables()
        self.__a = w_function.a
        self.__f = self.__restore_function_f(w_function.get_function(), self.__variables)

    def calculate(self, x_value, y_value, t_value):
        x, y, t = self.__variables
        return self.__f.subs({x: x_value, y: y_value, t: t_value})

    def __restore_function_f(self, w, variables):
        partialDerivative = PartialDerivatives(w, variables)

        dw_dt = partialDerivative.calculate_partial_derivative('t')
        d2w_dx2 = partialDerivative.calculate_second_partial_derivative('x')
        d2w_dy2 = partialDerivative.calculate_second_partial_derivative('y')

        return dw_dt - self.__a * (d2w_dx2 + d2w_dy2)