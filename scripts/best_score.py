import sys

protein = sys.argv[1]

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values

mufold_rmsd = read_file(f'statistics/predictor/{protein}/mufold_rmsd')

print(min([(y, x) for (x, y) in enumerate(mufold_rmsd)])[1] + 1)
