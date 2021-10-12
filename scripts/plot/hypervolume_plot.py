import sys
import matplotlib
import matplotlib.pyplot as plt

output_file = sys.argv[1]
graph_title = sys.argv[2]
data_files = sys.argv[3:]


def read_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            line_values = line.split()
            values.append([int(line_values[0]), float(line_values[1])])
    return values


hypervolumes = []
for file in data_files:
    hypervolume_data = read_file(file)
    hypervolumes.append(list(zip(*hypervolume_data)))

hypervolumes_x = []
hypervolumes_y = []
cnt = 1
for hv in hypervolumes:
    for value in hv[1]:
        hypervolumes_x.append(cnt)
        hypervolumes_y.append(value)
        cnt += 1

plt.rc("axes", axisbelow=True)

plt.title(graph_title)
plt.xlabel("Generations", labelpad=5)
plt.ylabel("Hypervolume", labelpad=10)

plt.grid(True)

plt.plot(hypervolumes_x, hypervolumes_y)

plt.ylim([0.0, 1.0])

plt.savefig(output_file, dpi=1200)
