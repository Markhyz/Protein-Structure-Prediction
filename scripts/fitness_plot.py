import sys
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('GTK3Agg')

fitness_data_dir = sys.argv[1]
fitness_data_files = sys.argv[2:]

def read_fitness_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            values.append(list(map(float, line.split())))
    return values

fitness_data = []
for file in fitness_data_files:
    file_path = f'{fitness_data_dir}/{file}.fitness'
    fitness_by_type = list(zip(*read_fitness_file(file_path)))
    fitness_data.append(fitness_by_type)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

for fitness_graph in fitness_data:
    ax.scatter(fitness_graph[1], fitness_graph[2], fitness_graph[0])

plt.show()