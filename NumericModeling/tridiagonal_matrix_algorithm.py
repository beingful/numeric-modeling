from numpy import zeros

class TridiagonalMatrixAlgorithm:
    def __init__(self, main_diag, lower_diag, upper_diag, solutions):
        self.__main_diag = main_diag
        self.__lower_diag = lower_diag
        self.__upper_diag = upper_diag
        self.__solutions = solutions

    def use_algorithm(self):
        n = len(self.__solutions)

        alpha = zeros(n - 1, float)
        betha = zeros(n, float)
        result = zeros(n, float)

        alpha[0] = self.__upper_diag[0] / self.__main_diag[0]
        betha[0] = self.__solutions[0] / self.__main_diag[0]

        for i in range(1, n - 1):
            alpha[i] = (self.__upper_diag[i]
                        / (self.__main_diag[i] - self.__lower_diag[i - 1] * alpha[i - 1]))

        for i in range(1, n):
            betha[i] = ((self.__solutions[i] - self.__lower_diag[i - 1] * betha[i - 1]) 
                        / (self.__main_diag[i] - self.__lower_diag[i - 1] * alpha[i - 1]))

        result[n - 1] = betha[n - 1]

        for i in range(n - 1, 0, -1):
            result[i - 1] = betha[i - 1] - alpha[i - 1] * result[i]

        return result
