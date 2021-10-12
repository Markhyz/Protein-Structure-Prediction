#include "../rosetta.hpp"

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
    string structure_file_1 = argv[1];
    string structure_file_2 = argv[2];
    string type = argv[3];
    size_t pose_start = argc > 4 ? stoi(argv[4]) : 1;

    core::pose::Pose structure_1;
    core::import_pose::pose_from_file(structure_1, structure_file_1, false, core::import_pose::PDB_file);

    if (pose_start - 1 > 0) {
        structure_1.delete_residue_range_slow(1, pose_start - 1);
    }

    core::pose::Pose structure_2;
    core::import_pose::pose_from_file(structure_2, structure_file_2, false, core::import_pose::PDB_file);

    cout.rdbuf(coutbuf);

    cout.precision(9);
    cout << fixed;

    if (type == "rmsd") {
        cout << core::scoring::CA_rmsd(structure_1, structure_2) << endl;
    }

    if (type == "gdt") {
        cout << core::scoring::CA_gdtmm(structure_1, structure_2) << endl;
    }

    return 0;
}