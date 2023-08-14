from statistics import mean, stdev
import sys


def stats(values, type):
    return [
        max(values) * 100 if type == 1 else min(values),
        mean(values) * 100 if type == 1 else mean(values),
        stdev(values) * 100 if type == 1 else stdev(values),
    ]


metric = sys.argv[1]
# offset = int(sys.argv[2])


def read_external_file(values, filename):
    with open(filename) as file:
        for line in file:
            line_values = list(filter(lambda x: x != "", line.split(" ")))
            if line_values[1] not in values:
                values[line_values[1]] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
            values[line_values[1]][line_values[0].upper()] = line_values[2:]
    return values


def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values


print("\\begin{landscape}")
print("\\centering")
print("\\small")
print("\\begin{longtable}{|cc|c|c|c|c|c|c|c|c|c|c|}")
print("\\hline")
print(
    "\\multicolumn{2}{|c|}{\\textbf{Protein}} & \\textbf{NSGA2} & \\textbf{GDE3} & \\textbf{DEMO} & \\textbf{MOPSO} & \\textbf{MODE-K} & \\textbf{SCDE} & \\textbf{Rosetta} & \\textbf{\\begin{tabular}[c]{@{}c@{}}MO-BRKGA\\\\ Best\\end{tabular}} & \\textbf{\\begin{tabular}[c]{@{}c@{}}MO-BRKGA\\\\ Mean\\end{tabular}} & \\textbf{\\begin{tabular}[c]{@{}c@{}}MO-BRKGA\\\\ MUFOLD-CL\\end{tabular}} \\\\ \\hline"
)


proteins = [
    "1ab1",
    "1acw",
    "1ail",
    "1aly",
    "1bdd",
    "1crn",
    "1dfn",
    "1enh",
    "1gb1",
    "1hhp",
    "1i6c",
    "1rop",
    "1zdd",
    "2kdl",
    "2mr9",
    "2p81",
    "T0868",
    "T0900",
    "T0968s1",
    "T1010",
]
# proteins = ["1ab1", "1acw"]

scores = {}

read_external_file(scores, "statistics/chen2020")
read_external_file(scores, "statistics/narloch2020")
read_external_file(scores, "statistics/song2018")
read_external_file(scores, "statistics/zhang2018")

scores["Rosetta"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (best)"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (mean)"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MO-BRKGA (MUFOLD-CL)"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}

for protein in proteins:
    rosetta_rmsd = read_file(f"statistics/predictor/{protein}/rosetta_rmsd")
    rosetta_gdt = read_file(f"statistics/predictor/{protein}/rosetta_gdt")
    scores["Rosetta"][protein.upper()] = stats(rosetta_rmsd, 0) + stats(rosetta_gdt, 1)

    best_rmsd = read_file(f"statistics/final_2/{protein}/best.rmsd")
    best_gdt = read_file(f"statistics/final_2/{protein}/best.gdt")
    scores["MO-BRKGA (best)"][protein.upper()] = stats(best_rmsd, 0) + stats(
        best_gdt, 1
    )

    mean_rmsd = read_file(f"statistics/final_2/{protein}/mean.rmsd")
    mean_gdt = read_file(f"statistics/final_2/{protein}/mean.gdt")
    scores["MO-BRKGA (mean)"][protein.upper()] = stats(mean_rmsd, 0) + stats(
        mean_gdt, 1
    )

    mufold_rmsd = read_file(f"statistics/final_2/{protein}/mufold.rmsd")
    mufold_gdt = read_file(f"statistics/final_2/{protein}/mufold.gdt")
    scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()] = stats(mufold_rmsd, 0) + stats(
        mufold_gdt, 1
    )

    methods = [
        "NSGA2",
        "GDE3",
        "DEMO",
        "MOPSO",
        "MODE-K",
        "SCDE",
        "Rosetta",
        "MO-BRKGA (best)",
        "MO-BRKGA (mean)",
        "MO-BRKGA (MUFOLD-CL)",
    ]

    best_best = -1
    best_idx = -1
    best_avg = -1
    avg_idx = -1

    for (index, method) in enumerate(methods):
        best_mu = float(
            scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][0]
            if metric == "RMSD"
            else scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][3]
        )
        avg_mu = float(
            scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][1]
            if metric == "RMSD"
            else scores["MO-BRKGA (MUFOLD-CL)"][protein.upper()][4]
        )
        if protein.upper() in scores[method]:
            best = (
                scores[method][protein.upper()][0]
                if metric == "RMSD"
                else scores[method][protein.upper()][3]
            )
            avg = (
                scores[method][protein.upper()][1]
                if metric == "RMSD"
                else scores[method][protein.upper()][4]
            )
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

    print(
        f"{protein.upper()}",
        "& \\begin{tabular}[c]{@{}c@{}}$f^*$ \\\\ $\\overline{x}$\\\\ $s$\\end{tabular}",
        end="",
    )

    for (index, method) in enumerate(methods):
        print("& \\begin{tabular}[c]{@{}c@{}}", end="")

        best = "-"
        avg = "-"
        std = "-"

        if protein.upper() in scores[method]:
            best = (
                scores[method][protein.upper()][0]
                if metric == "RMSD"
                else scores[method][protein.upper()][3]
            )
            avg = (
                scores[method][protein.upper()][1]
                if metric == "RMSD"
                else scores[method][protein.upper()][4]
            )
            std = (
                scores[method][protein.upper()][2]
                if metric == "RMSD"
                else scores[method][protein.upper()][5]
            )

        if best != "-":
            if best_idx == index:
                print(
                    "\\textbf{", "{:.2f}".format(float(best)), "} \\\\", sep="", end=""
                )
            else:
                print("{:.2f} \\\\".format(float(best)), end="")
        else:
            print("- \\\\", end="")

        if avg != "-":
            print("{:.2f} \\\\ {:.2f}".format(float(avg), float(std)), end="")
        else:
            print("- \\\\ -", end="")

        print("\\end{tabular}", end="")

    print("\\\\ \\hline")

print("\\multicolumn{2}{|c|}{\\textbf{B/S/W}}", end="")
for (index, method) in enumerate(methods):

    print("& {}/0/{}".format(scores[method]["ba"], scores[method]["wa"]), end="")

print("\\\\ \\hline")


print("\\end{longtable}")
print("\\\\ Source: Author")
print("\\end{landscape}")
