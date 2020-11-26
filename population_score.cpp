#include "rosetta.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    ofstream out("/dev/null");
    auto* coutbuf = cout.rdbuf();
    cout.rdbuf(out.rdbuf());

    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string structure_dir = argv[1];
    string native_structure_path = argv[2];
    int num_decoys = stoi(argv[3]);
    string type = argv[4];
    size_t pose_start = argc > 5 ? stoi(argv[5]) : 1;

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, native_structure_path, false, core::import_pose::PDB_file);

    double best_rmsd = 1e9, mean_rmsd = 0, mean_gdt = 0, best_gdt = 0;

    for (int i = 1; i <= num_decoys; ++i) {
        core::pose::Pose structure;
        core::import_pose::pose_from_file(structure, structure_dir + string("/decoy_") + to_string(i) + string(".pdb"),
                                          false, core::import_pose::PDB_file);

        if (pose_start - 1 > 0) {
            structure.delete_residue_range_slow(1, pose_start - 1);
        }

        double rmsd = core::scoring::CA_rmsd(structure, native_structure);
        double gdt = core::scoring::CA_gdtmm(structure, native_structure);

        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
        }

        if (gdt > best_gdt) {
            best_gdt = gdt;
        }

        mean_rmsd += rmsd;
        mean_gdt += gdt;
    }

    mean_rmsd /= num_decoys;
    mean_gdt /= num_decoys;

    if (type == "rmsd") {
        printf("%.9f\n", best_rmsd);
    }

    if (type == "gdt") {
        printf("%.9f\n", best_gdt);
    }

    if (type == "mean_rmsd") {
        printf("%.9f\n", mean_rmsd);
    }

    if (type == "mean_gdt") {
        printf("%.9f\n", mean_gdt);
    }

    cout.rdbuf(coutbuf);

    return 0;
}