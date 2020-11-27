from statistics import mean, stdev
import sys

protein_list = ["1crn", "1gb1", "1hhp", "1i6c", "1vii", "T0868", "T0900"]

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values
        
def stats(values):
    return (mean(values), stdev(values))

for protein in protein_list:
    ca_rmsd      = read_file(f"results/{protein}/best_ca_rmsd2")
    ca_gdt       = read_file(f"results/{protein}/best_ca_gdt2")

    mean_rmsd    = read_file(f"results/{protein}/mean_rmsd")
    mean_gdt     = read_file(f"results/{protein}/mean_gdt")

    mufold_rmsd  = read_file(f"results/{protein}/mufold_rmsd")
    mufold_gdt   = read_file(f"results/{protein}/mufold_gdt")

    rosetta_rmsd = read_file(f"results/{protein}/rosetta_rmsd")
    rosetta_gdt  = read_file(f"results/{protein}/rosetta_gdt")

    best_rmsd = min(ca_rmsd)
    best_gdt = max(ca_gdt)

    stat_rmsd = stats(ca_rmsd)
    stat_gdt = stats(ca_gdt)

    stat_mean_rmsd = stats(mean_rmsd)
    stat_mean_gdt = stats(mean_gdt)

    stat_mufold_rmsd = stats(mufold_rmsd)
    stat_mufold_gdt = stats(mufold_gdt)

    stat_rosetta_rmsd = stats(rosetta_rmsd)
    stat_rosetta_gdt = stats(rosetta_gdt)

    print(protein.upper(), ": {.3f} {.3f} {.3f}/{.3f} {.3f}/{.3f} {.3f}/{.3f} {.3f}/{.3f}", 
          best_rmsd, best_gdt, stat_rmsd[0], stat_rmsd[1], stat_gdt[0], stat_gdt[1], stat_mean_rmsd[0], stat_mean_rmsd[1],
          stat_mufold_rmsd[0], stat_mufold_rmsd[1]. stat_rosetta_rmsd[0], stat_rosetta_rmsd[1])