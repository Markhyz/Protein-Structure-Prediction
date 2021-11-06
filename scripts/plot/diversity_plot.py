import sys
import matplotlib
import matplotlib.pyplot as plt

output_file = sys.argv[1]
data_file = sys.argv[2]
graph_title = sys.argv[3]


def read_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            line_values = line.split()
            values.append([int(line_values[0]), float(line_values[1])])
    return values


diversity_data = read_file(data_file)

second_phase_first_index = len(diversity_data) // 2

diversity_data[second_phase_first_index][1] = diversity_data[
    second_phase_first_index + 1
][1]

diversity_data = list(zip(*diversity_data))

plt.rc("axes", axisbelow=True)

plt.title(graph_title)
plt.xlabel("Generations", labelpad=5)
plt.ylabel("Diversity", labelpad=10)

plt.grid(True)

plt.plot(diversity_data[0], diversity_data[1])

plt.ylim([0, 0.5])

plt.savefig(output_file, dpi=1200)
