import sys
import numpy as np
from pymoo.factory import get_problem, get_performance_indicator

fitness_data_file = sys.argv[1]
gen_num = int(sys.argv[2])

def read_fitness_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            values.append(list(map(float, line.split())))
            values[-1][1] *= -1
            values[-1][2] *= -1
    return values

for gen in range(gen_num):
    file_path = f'{fitness_data_file}{gen + 1}.norm_fitness'
    fitness_data = np.array(read_fitness_file(file_path))
    hypervolume = get_performance_indicator("hv", ref_point=np.array([1.0, 0.0, 0.0])).calc(fitness_data)
    print(gen + 1, hypervolume)

# fitness_data = [[0.5, 0.3, 0.7], [0.2, 0.8, 0.4], [0.1, 0.8, 0.1], [0.3, 0.5, 0.6], [0.8, 1.0, 0.9]]
# fitness_2 = np.array([[x[0], -x[1], -x[2]] for x in fitness_data])
# print(get_performance_indicator("hv", ref_point=np.array([1.0, 0.0, 0.0])).calc(fitness_2))
