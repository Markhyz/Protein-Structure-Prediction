from statistics import mean, stdev
import sys


def stats(values, type):
    return [
        max(values) * 100 if type == 1 else min(values),
        mean(values) * 100 if type == 1 else mean(values),
        stdev(values) * 100 if type == 1 else stdev(values),
    ]


metric = sys.argv[1]
s_method = sys.argv[2]


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

scores["old_best"] = {}
scores["old_MUFOLD"] = {}
scores["best"] = {}
scores["MUFOLD"] = {}

methods = [
    "best",
    "MUFOLD",
    "old_best",
    "old_MUFOLD",
]

print("protein method value")

for protein in proteins:
    scores["old_best"]["RMSD"] = read_file(f"statistics/predictor/{protein}/best_rmsd")
    scores["old_best"]["GDT"] = read_file(f"statistics/predictor/{protein}/best_gdt")
    scores["old_MUFOLD"]["RMSD"] = read_file(
        f"statistics/predictor/{protein}/mufold_rmsd"
    )
    scores["old_MUFOLD"]["GDT"] = read_file(
        f"statistics/predictor/{protein}/mufold_gdt"
    )

    scores["best"]["RMSD"] = read_file(f"statistics/final_2/{protein}/best.rmsd")
    scores["best"]["GDT"] = read_file(f"statistics/final_2/{protein}/best.gdt")

    scores["MUFOLD"]["RMSD"] = read_file(f"statistics/final_2/{protein}/mufold.rmsd")
    scores["MUFOLD"]["GDT"] = read_file(f"statistics/final_2/{protein}/mufold.gdt")

    for (index, method) in enumerate(methods):
        if method != "best" and method != s_method:
            continue
        for value in scores[method][metric]:
            print(protein, method, value)
