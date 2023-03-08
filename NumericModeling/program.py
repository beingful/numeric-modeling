from absolute_error import *
from alternating_direction_implicit_method import *
from concurrent.futures import ThreadPoolExecutor
from explicit_difference_scheme import *
from function import *
from heat_conduction_plot import *
from two_step_symmetric_algorithm import *

function_w = FunctionW()

with ThreadPoolExecutor() as true_function_executor:
    true_function_future = true_function_executor.submit(lambda: function_w.calculate_all())

function_f = FunctionF(function_w)

methods = [ExplicitDifferenceScheme(function_w, function_f),
           AlternatingDirectionImplicitMethod(function_w, function_f),
           TwoStepSymmetricAlgorithm(function_w, function_f)]

with ThreadPoolExecutor() as methods_executor:
    methods_results = methods_executor.map(lambda method: method.use_algorithm(), methods)

true_function_results = true_function_future.result()

heat_conduction_plot = HeatConductionPlot(true_function_results, 'True function results')
heat_conduction_plot.build_plot()

for result in methods_results:
    method_name = result['method_name']
    method_results = result['method_results']

    abs_error = AbsoluteError(true_function_results, method_results, method_name)
    abs_error.build_plot()

    heat_conduction_plot = HeatConductionPlot(method_results, method_name)
    heat_conduction_plot.build_plot()