import sys

protein_dir = sys.argv[1]
size = int(sys.argv[2])
file_path = sys.argv[3]

def read_hypervolume_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            if line != "":
                line_values = line.split()
                values.append(float(line_values[1]))
    return values

mean_hypervolume = read_hypervolume_file(f'{protein_dir}/1/{file_path}')
for i in range(2, size + 1):
    hypervolume = read_hypervolume_file(f'{protein_dir}/{i}/{file_path}')
    for (index, value) in enumerate(hypervolume):
        mean_hypervolume[index] += value

for i in range(len(mean_hypervolume)):
    mean_hypervolume[i] /= size
    print("{} {:.9f}".format(i + 1, mean_hypervolume[i] + 0.15))
