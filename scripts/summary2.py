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

scores["best"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["MUFOLD"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}

scores["old_best"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}
scores["old_MUFOLD"] = {"bb": 0, "wb": 0, "ba": 0, "wa": 0}

methods = [
    "best",
    "MUFOLD",
    "old_best",
    "old_MUFOLD",
]

print("protein", " ".join(methods))

for protein in proteins:
    old_best_rmsd = read_file(f"statistics/predictor/{protein}/best_rmsd")
    old_best_gdt = read_file(f"statistics/predictor/{protein}/best_gdt")
    scores["old_best"][protein.upper()] = stats(old_best_rmsd, 0) + stats(
        old_best_gdt, 1
    )

    old_mufold_rmsd = read_file(f"statistics/predictor/{protein}/mufold_rmsd")
    old_mufold_gdt = read_file(f"statistics/predictor/{protein}/mufold_gdt")
    scores["old_MUFOLD"][protein.upper()] = stats(old_mufold_rmsd, 0) + stats(
        old_mufold_gdt, 1
    )

    best_rmsd = read_file(f"statistics/final_2/{protein}/best.rmsd")
    best_gdt = read_file(f"statistics/final_2/{protein}/best.gdt")
    scores["best"][protein.upper()] = stats(best_rmsd, 0) + stats(best_gdt, 1)

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
