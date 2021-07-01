from statistics import mean, stdev
import sys

def stats(values, type):
    return [max(values) if type == 1 else min(values), mean(values) * 100 if type == 1 else mean(values), stdev(values) * 100 if type == 1 else stdev(values)]

protein = sys.argv[1]

def read_external_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            line_values = line.split(' ')
            if line_values[0] == protein:
                values.append(line_values[1:])
    return values

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values

scores = []

scores.append(read_external_file('statistics/chen2020'))
scores.append(read_external_file('statistics/narloch2020'))
scores.append(read_external_file('statistics/song2018'))
scores.append(read_external_file('statistics/zhang2018'))

rosetta_rmsd = read_file(f"statistics/predictor/{protein}/rosetta_rmsd")
rosetta_gdt  = read_file(f"statistics/predictor/{protein}/rosetta_gdt")
scores.append(["rosetta"] + stats(rosetta_rmsd, 0) + stats(rosetta_gdt, 1))

best_rmsd = read_file(f"statistics/predictor/{protein}/best_rmsd")
best_gdt  = read_file(f"statistics/predictor/{protein}/best_gdt")
scores.append(["MO-BRKGA (best)"] + stats(best_rmsd, 0) + stats(best_gdt, 1))

mean_rmsd = read_file(f"statistics/predictor/{protein}/mean_rmsd")
mean_gdt  = read_file(f"statistics/predictor/{protein}/mean_gdt")
scores.append(["MO-BRKGA (mean)"] + stats(mean_rmsd, 0) + stats(mean_gdt, 1))

mufold_rmsd = read_file(f"statistics/predictor/{protein}/mufold_rmsd")
mufold_gdt  = read_file(f"statistics/predictor/{protein}/mufold_gdt")
scores.append(["MO-BRKGA (MUFOLD-CL)"] + stats(mufold_rmsd, 0) + stats(mufold_gdt, 1))

scores = list(filter(lambda x: x != [], scores))

print("\\begin{table}[ht]")
print("\\center")
print("\\caption{", protein.upper(), "prediction comparison}")
print("\\begin{tabular}{|c|c|c|c|c|}")
print("\\hline")
print("\\textbf{Method} & \\textbf{Best RMSD} & \\textbf{Average RMSD} & \\textbf{Best GDT\\_TS} & \\textbf{Average GDT\\_TS} \\\\ \\hline")

for score in scores:
    print(score)
    print(score[0], end="")
    
    if score[1] != "-":
        print("& {:.1f}".format(float(score[1])), end="")
    else:
        print("& -", end="")

    if score[2] != "-":
        if score[3] != "-":
            print("& {:.1f} $\\pm$ {:.1f}".format(score[2], score[3]), end="")
        else:
            print("& {:.1f}".format(score[2]), end="")
    else:
        print("& -", end="")

    if score[4] != "-":
        print("& {:.1f}".format(float(score[4])), end="")
    else:
        print("& -", end="")

    if score[5] != "-":
        if score[6] != "-":
            print("& {:.1f} $\\pm$ {:.1f}".format(score[5], score[6]), end="")
        else:
            print("& {:.1f}".format(score[5]), end="")
    else:
        print("& -", end="")
    
    print("\\\\ \\hline")

print("\\end{tabular}")
print("\\label{tab:", protein, "_pred}")
print("\\end{table}")


