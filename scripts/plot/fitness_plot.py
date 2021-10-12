import sys
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('GTK3Agg')

output_file = sys.argv[1]
graph_title = sys.argv[2]
is_full = int(sys.argv[3])
fitness_data_dir = sys.argv[4]
fitness_data_files = sys.argv[5:]

def read_fitness_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            values.append(list(map(float, line.split())))
    return values

fitness_data = []
for file in fitness_data_files:
    file_path = f'{fitness_data_dir}/{file}.norm_fitness'
    fitness_by_type = list(zip(*read_fitness_file(file_path)))
    fitness_data.append(fitness_by_type)

markers = ['o', '^', 's', 'x', 'd']

plt.rc('axes', axisbelow=True)

plt.title(graph_title)
plt.xlabel("Secondary Structure", labelpad=5)
plt.ylabel("Contact Map", labelpad=10)

plt.grid(True)

for (index, fitness) in enumerate(fitness_data):
    plt.scatter(fitness[1], fitness[2], 30, marker=markers[index], c=fitness[0], cmap=plt.cm.jet, label=fitness_data_files[index])


if is_full:
    plt.xlim([0, 1])
    plt.ylim([0, 1])

plt.clim(0, 1)

color_bar = plt.colorbar()
color_bar.set_label("Energy Function", labelpad=10)

legend = plt.legend()

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# for fitness in fitness_data:
#     ax.scatter(fitness[1], fitness[2], fitness[0])

plt.savefig(output_file, dpi=1200)