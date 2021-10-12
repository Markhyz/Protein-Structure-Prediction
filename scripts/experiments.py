from statistics import mean, stdev
import sys

def stats(values, type):
    return [max(values) * 100 if type == 1 else min(values), mean(values) * 100 if type == 1 else mean(values), stdev(values) * 100 if type == 1 else stdev(values)]

metric = sys.argv[1]
# offset = int(sys.argv[2])

def read_external_file(values, filename):
    with open(filename) as file:
        for line in file:
            line_values = list(filter(lambda x: x != '', line.split(' ')))
            if line_values[1] not in values:
                values[line_values[1]] = { "bb": 0, "wb": 0, "ba": 0, "wa": 0}
            values[line_values[1]][line_values[0].upper()] = line_values[2:]
    return values

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values

print("\\begin{table}[H]")
print("\\center")
print("\\begin{tabular}{|c|ccc|}")
print("\\hline")
print("\\textbf{Protein} & \\multicolumn{1}{c|}{\\textbf{Method}} & \\multicolumn{1}{c|}{\\textbf{Best}} & \\textbf{Average} \\\\ \\hline")

proteins = ["1ab1", "1acw", "1ail", "1aly", "1bdd", "1crn", "1dfn", "1enh", "1gb1", "1hhp", "1i6c", "1rop", "1zdd", "2kdl", "2mr9", "2p81", "T0868", "T0900", "T0968s1", "T1010"]
# proteins = ["1ab1", "1acw"]

scores = {}

read_external_file(scores, 'statistics/chen2020')
read_external_file(scores, 'statistics/narloch2020')
read_external_file(scores, 'statistics/song2018')
read_external_file(scores, 'statistics/zhang2018')

scores["Rosetta"] = { "bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (best)"] = { "bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (mean)"] = { "bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (MUFOLD-CL)"] = { "bb": 0, "wb": 0, "ba": 0, "wa": 0}

for protein in proteins:
    rosetta_rmsd = read_file(f"statistics/predictor/{protein}/rosetta_rmsd")
    rosetta_gdt  = read_file(f"statistics/predictor/{protein}/rosetta_gdt")
    scores["Rosetta"][protein.upper()] = stats(rosetta_rmsd, 0) + stats(rosetta_gdt, 1)

    best_rmsd = read_file(f"statistics/predictor/{protein}/best_rmsd")
    best_gdt  = read_file(f"statistics/predictor/{protein}/best_gdt")
    scores["MO-BRKGA (best)"][protein.upper()] = stats(best_rmsd, 0) + stats(best_gdt, 1)

    mean_rmsd = read_file(f"statistics/predictor/{protein}/mean_rmsd")
    mean_gdt  = read_file(f"statistics/predictor/{protein}/mean_gdt")
    scores["MO-BRKGA (mean)"][protein.upper()] = stats(mean_rmsd, 0) + stats(mean_gdt, 1)

    mufold_rmsd = read_file(f"statistics/predictor/{protein}/mufold_rmsd")
    mufold_gdt  = read_file(f"statistics/predictor/{protein}/mufold_gdt")
    scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()] = stats(mufold_rmsd, 0) + stats(mufold_gdt, 1)

    methods = ["NSGA2", "GDE3", "DEMO", "MOPSO", "MODE-K", "SCDE", "Rosetta", "MO-BRKGA (best)", "MO-BRKGA (mean)", "MO-BRKGA (MUFOLD-CL)"]

    best_best = -1
    best_idx = -1
    best_avg = -1
    avg_idx = -1
    

    for (index, method) in enumerate(methods):
        best_mu = float(scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][0] if metric == "RMSD" else scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][3])
        avg_mu = float(scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][1] if metric == "RMSD" else scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][4])
        if protein.upper() in scores[method]:
            best = scores[method][protein.upper()][0] if metric == "RMSD" else scores[method][protein.upper()][3]
            avg = scores[method][protein.upper()][1] if metric == "RMSD" else scores[method][protein.upper()][4]
            if metric == "RMSD":
                if best != "-":
                    best = float(best)
                    if best < best_best or best_best == -1:
                        best_best = best
                        best_idx = index
                    if best < best_mu:
                        scores[method]["bb"] += 1
                    else:
                        scores[method]["wb"] += 1 
                if avg != "-":
                    avg = float(avg)
                    if avg < best_avg or best_avg == -1:
                        best_avg = avg
                        avg_idx = index
                    if avg < avg_mu:
                        scores[method]["ba"] += 1
                    else:
                        scores[method]["wa"] += 1
            else:
                if best != "-":
                    best = float(best)
                    if best > best_best or best_best == -1:
                        best_best = best
                        best_idx = index
                    if best > best_mu:
                        scores[method]["bb"] += 1
                    else:
                        scores[method]["wb"] += 1 
                if avg != "-":
                    avg = float(avg)
                    if avg > best_avg or best_avg == -1:
                        best_avg = avg
                        avg_idx = index
                    if avg > avg_mu:
                        scores[method]["ba"] += 1
                    else:
                        scores[method]["wa"] += 1 

    for (index, method) in enumerate(methods):
        color = "{FFFFFF}"
        if index % 2 != 0:
            color = "{EFEFEF}"
        if index == len(methods) - 1:
            print("\\multirow{-10}{*}{", protein.upper(), "} ", sep="")
        print(f'& \cellcolor[HTML]{color} {method}', end="") 
        
        best = "-"
        avg = "-"
        std = "-"

        if protein.upper() in scores[method]:
            best = scores[method][protein.upper()][0] if metric == "RMSD" else scores[method][protein.upper()][3]
            avg = scores[method][protein.upper()][1] if metric == "RMSD" else scores[method][protein.upper()][4]
            std = scores[method][protein.upper()][2] if metric == "RMSD" else scores[method][protein.upper()][5]

        if best != "-":
            if best_idx == index:
                print(" & \\cellcolor[HTML]{} \\textbf".format(color), "{", "{:.2f}".format(float(best)), "}", sep="", end="")
            else:
                print(" & \\cellcolor[HTML]{} {:.2f}".format(color, float(best)), end="")
        else:
            print(" & \\cellcolor[HTML]{} -".format(color), end="")

        print("&\\multicolumn{1}{c|}{", end="")
        if avg != "-":
            if std != "-":
                if avg_idx == index:
                    print("\\cellcolor[HTML]{} \\textbf".format(color), "{", "{:.2f} $\\pm$ {:.2f}".format(float(avg), float(std)), "}", end="")
                else:
                    print("\\cellcolor[HTML]{} {:.2f} $\\pm$ {:.2f}".format(color, float(avg), float(std)), end="")
            else:
                if avg_idx == index:
                    print("\\cellcolor[HTML]{} \\textbf".format(color), "{", "{:.2f}".format(float(avg)), "}", end="")
                else:
                    print("\\cellcolor[HTML]{} {:.2f}".format(color, float(avg)), end="")
        else:
            print("\\cellcolor[HTML]{} -".format(color), end="")
        
        print("}\\\\")

    print("\\hline", sep="")

for (index, method) in enumerate(methods):
    color = "{FFFFFF}"
    if index % 2 != 0:
        color = "{EFEFEF}"
    if index == len(methods) - 1:
        print("\\multirow{-10}{*}{ \\textbf{B/S/W} } ", sep="")
    print(f'& \cellcolor[HTML]{color} {method}', end="") 

    print(" & \\cellcolor[HTML]", color, "{}/0/{}".format(scores[method]["bb"], scores[method]["wb"]), end = "")
    print("&\\multicolumn{1}{c|}{\\cellcolor[HTML]", color, "{}/0/{}".format(scores[method]["ba"], scores[method]["wa"]), "}", end="")
    print("\\\\")

print("\\hline", sep="")


print("\\end{tabular}")
print("\\\\ Source: Author")
print("\\end{table}")


