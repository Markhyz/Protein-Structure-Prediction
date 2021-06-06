import sys
import numpy as np
from pymoo.visualization.radviz import Radviz

fitness_data_file = sys.argv[1]
output_file = sys.argv[2]
graph_title = sys.argv[3]

def read_fitness_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            values.append(list(map(float, line.split())))
    return values

file_path = f'{fitness_data_file}.norm_fitness'
fitness_data = np.array(read_fitness_file(file_path))

plt = Radviz(title=graph_title, labels=["Energy Function", "Secondary Structure", "Contact Map"]).add(fitness_data, s = 5)
plt.save(output_file, dpi=1200)