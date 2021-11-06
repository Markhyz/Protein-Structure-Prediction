from statistics import mean, stdev
import sys


def stats(values, type):
    return [
        max(values) * 100 if type == 1 else min(values),
        mean(values) * 100 if type == 1 else mean(values),
        stdev(values) * 100 if type == 1 else stdev(values),
    ]


metric = sys.argv[1]


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

scores = {}

read_external_file(scores, "statistics/chen2020")
read_external_file(scores, "statistics/narloch2020")
read_external_file(scores, "statistics/song2018")

if metric == "RMSD":
    read_external_file(scores, "statistics/zhang2018")

scores["Rosetta"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["best"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["mean"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MUFOLD"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}

methods = [
    "NSGA2",
    "GDE3",
    "DEMO",
    "SCDE",
    "Rosetta",
    "best",
    "mean",
    "MUFOLD",
]

if metric == "GDT":
    methods = list(filter(lambda x: x != "SCDE", methods))

print("protein", " ".join(methods))

for protein in proteins:
    rosetta_rmsd = read_file(f"statistics/predictor/{protein}/rosetta_rmsd")
    rosetta_gdt = read_file(f"statistics/predictor/{protein}/rosetta_gdt")
    scores["Rosetta"][protein.upper()] = stats(rosetta_rmsd, 0) + stats(rosetta_gdt, 1)

    best_rmsd = read_file(f"statistics/final_2/{protein}/best.rmsd")
    best_gdt = read_file(f"statistics/final_2/{protein}/best.gdt")
    scores["best"][protein.upper()] = stats(best_rmsd, 0) + stats(best_gdt, 1)

    mean_rmsd = read_file(f"statistics/final_2/{protein}/mean.rmsd")
    mean_gdt = read_file(f"statistics/final_2/{protein}/mean.gdt")
    scores["mean"][protein.upper()] = stats(mean_rmsd, 0) + stats(mean_gdt, 1)

    mufold_rmsd = read_file(f"statistics/final_2/{protein}/mufold.rmsd")
    mufold_gdt = read_file(f"statistics/final_2/{protein}/mufold.gdt")
    scores["MUFOLD"][protein.upper()] = stats(mufold_rmsd, 0) + stats(mufold_gdt, 1)

    print(protein, "20 " * len(methods))

    print(protein, end=" ")
    for (index, method) in enumerate(methods):
        avg = 0.0

        if protein.upper() in scores[method]:
            avg = (
                scores[method][protein.upper()][1]
                if metric == "RMSD"
                else scores[method][protein.upper()][4]
            )

        print(avg, end=" ")
    print("")
    print(protein, end=" ")
    for (index, method) in enumerate(methods):
        std = 0.0

        if protein.upper() in scores[method]:
            std = (
                scores[method][protein.upper()][2]
                if metric == "RMSD"
                else scores[method][protein.upper()][5]
            )

        print(float(std), end=" ")
    print("")
