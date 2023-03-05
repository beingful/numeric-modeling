from matplotlib.pyplot import xlabel, ylabel, title, plot, show
from numpy import zeros

class AbsoluteError:
    __N, __M = 5, 20

    def __init__(self, true_function_results, method_results, method_name):
        self.__method_results = method_results
        self.__true_function_results = true_function_results
        self.__method_name = method_name
        self.__steps = zeros((self.__M + 1) * (self.__N + 1) ** 2)
        self.__method_abs_error = zeros((self.__M + 1) * (self.__N + 1) ** 2)

    def build_plot(self):
        self.__calculate_abs_errors()
        self.__setup_plot()

    def __calculate_abs_errors(self):
        counter = 0

        for i in range(0, self.__M + 1, 1):
            for j in range(0, self.__N + 1, 1):
                for k in range(0, self.__N + 1, 1):
                    self.__method_abs_error[counter] = \
                        abs(self.__true_function_results[i][j][k] - self.__method_results[i][j][k])
                    self.__steps[counter] = counter
                    counter += 1

    def __setup_plot(self):
        xlabel('step')
        ylabel('results')
        title(f'Absolute error for {self.__method_name} on each step')
        plot(self.__steps, self.__method_abs_error, 'r.', markersize=4)
        show()
