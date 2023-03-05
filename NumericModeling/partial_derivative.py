from sympy import symbols, diff

class PartialDerivatives:
    def __init__(self, w_function, variables):
        self.__w = w_function
        self.__variables = variables

    def calculate_partial_derivative(self, variable):
        if symbols(variable) in self.__variables:
            return diff(self.__w, variable)

    def calculate_second_partial_derivative(self, variable):
        if symbols(variable) in self.__variables:
            return diff(self.__w, variable, variable)
