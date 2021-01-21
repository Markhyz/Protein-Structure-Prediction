protein_list = ["1vii", "1i6c", "1crn", "1gb1", "1hhp", "T0868", "T0900"]

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values
        
def stats(values):
    return (mean(values), stdev(values))

rmsds = []
gdts = []

for protein in protein_list:
    ca_rmsd      = read_file(f"results/{protein}/best_ca_rmsd2")
    ca_gdt       = read_file(f"results/{protein}/best_ca_gdt2")

    mufold_rmsds  = read_file(f"results/{protein}/mufold_rmsd")
    mufold_gdts   = read_file(f"results/{protein}/mufold_gdt")

    rosetta_rmsds = read_file(f"results/{protein}/rosetta_rmsd")
    rosetta_gdts  = read_file(f"results/{protein}/rosetta_gdt")
    
    for index in range(len(ca_rmsd)):
        rmsds.append([protein, ca_rmsd[index], "Melhor"])
        rmsds.append([protein, mufold_rmsds[index], "MUFOLD-CL"])
        rmsds.append([protein, rosetta_rmsds[index], "Rosetta"])
        gdts.append([protein, ca_gdt[index], "Melhor"])
        gdts.append([protein, mufold_gdts[index], "MUFOLD-CL"])
        gdts.append([protein, rosetta_gdts[index], "Rosetta"])

print("protein value method")
for values in rmsds:
    print(values[0], values[1], values[2])

print("\nprotein value method")
for values in gdts:
    print(values[0], values[1], values[2])