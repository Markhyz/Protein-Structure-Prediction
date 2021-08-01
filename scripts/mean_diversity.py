import sys

protein_dir = sys.argv[1]
size = int(sys.argv[2])
file_path = sys.argv[3]

def read_diversity_file(filename):
    values = []
    with open(filename) as file:
        fitness = []
        for line in file:
            if line != "":
                line_values = line.split()
                values.append(float(line_values[1]))
    return values

mean_diversity = read_diversity_file(f'{protein_dir}/1/{file_path}')
for i in range(2, size + 1):
    diversity = read_diversity_file(f'{protein_dir}/{i}/{file_path}')
    for (index, value) in enumerate(diversity):
        mean_diversity[index] += value

for i in range(len(mean_diversity)):
    mean_diversity[i] /= size
    print("{} {:.9f}".format(i + 1, mean_diversity[i]))
