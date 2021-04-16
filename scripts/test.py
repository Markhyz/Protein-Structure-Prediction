from statistics import mean, stdev
import sys

#protein_list = ["1crn", "1gb1", "1hhp", "1i6c", "1vii", "T0868", "T0900"]
protein_list = ["1acw", "1ail", "1crn", "1enh", "1rop", "1zdd", "2mr9", "2p81"]

def read_file(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append(float(line))
    return values

def read_file2(filename):
    values = []
    with open(filename) as file:
        for line in file:
            values.append([float(x) for x in line.split(" ")])
    return values
        
def stats(values):
    return (mean(values), stdev(values))

header = """\\begin{table}[ht]
\\center
\\begin{tabular}{|l|l|l|l|l|l|}
\\hline
\\multicolumn{1}{|c|}{\\multirow{2}{*}{Protein}} & \\multicolumn{1}{c|}{\\multirow{2}{*}{Algorithm}} & \\multicolumn{2}{c|}{RMSD} & \\multicolumn{2}{c|}{GDT} \\\\ \\cline{3-6} 
\\multicolumn{1}{|c|}{} & \\multicolumn{1}{c|}{} & \\textbf{Best} & \\textbf{Mean} & \\textbf{Best} & \\textbf{Mean} \\\\ \\hline"""

footer = """\\end{tabular}
   \\end{table}"""


print(header)

narloch = read_file2("narloch2020")
n_index = 0

for protein in protein_list:
    mufold_rmsds  = read_file(f"results/{protein}/mufold_rmsd")
    mufold_gdts   = read_file(f"results/{protein}/mufold_gdt")

    mufold_rmsd = min(mufold_rmsds)
    mufold_gdt = max(mufold_gdts)

    stat_mufold_rmsd = stats(mufold_rmsds)
    stat_mufold_gdt = stats(mufold_gdts)
    
    print("\\multirow{4}{*}{", protein.upper(), "} & NSGA-II &", " {:.2f} & {:.2f} $\\pm$ {:.2f}  & {:.2f}  & {:.2f} $\\pm$ {:.2f}".format(
            narloch[n_index][0],  narloch[n_index][1], narloch[n_index][2], narloch[n_index][3], narloch[n_index][4], narloch[n_index][5]), " \\\\ \\cline{2-6}")
    print("& GDE3 & {:.2f}  & {:.2f} $\\pm$ {:.2f} & {:.2f} & {:.2f} $\\pm$ {:.2f} \\\\ \\cline\{2-6\}".format(
            narloch[n_index + 1][0],  narloch[n_index + 1][1], narloch[n_index + 1][2], narloch[n_index + 1][3], narloch[n_index + 1][4], narloch[n_index + 1][5]))
    print("& DEMO & {:.2f}  & {:.2f} $\\pm$ {:.2f} & {:.2f} & {:.2f} $\\pm$ {:.2f} \\\\ \\cline\{2-6\}".format(
            narloch[n_index + 2][0],  narloch[n_index + 2][1], narloch[n_index + 2][2], narloch[n_index + 2][3], narloch[n_index + 2][4], narloch[n_index + 2][5]))

    n_index += 3

    print("   & \\textbf{MO-BRKGA} & ", "{:.2f} & {:.2f}$\\pm${:.2f} & {:.2f} & {:.2f}$\\pm${:.2f}  \\\\ \\hline".format(
          mufold_rmsd,
          stat_mufold_rmsd[0],
          stat_mufold_rmsd[1],
          mufold_gdt,
          stat_mufold_gdt[0],
          stat_mufold_gdt[1]))

print(footer)
