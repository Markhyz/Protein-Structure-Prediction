import sys
import numpy as np
from pymoo.factory import get_problem, get_performance_indicator

fitness_data_file = sys.argv[1]

def read_fitness_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            values.append(list(map(float, line.split())))
            values[-1][1] *= -1
            values[-1][2] *= -1
    return values

file_path = f'{fitness_data_file}.norm_fitness'
fitness_data = np.array(read_fitness_file(file_path))

hypervolume = get_performance_indicator("hv", ref_point=np.array([1.0, 0.0, 0.0])).calc(fitness_data)
print(hypervolume)
