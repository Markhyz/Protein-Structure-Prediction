#include "rosetta.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string structure_file_1 = argv[1];
    string structure_file_2 = argv[2];

    core::pose::Pose structure_1;
    core::import_pose::pose_from_file(structure_1, structure_file_1, false, core::import_pose::PDB_file);

    core::pose::Pose structure_2;
    core::import_pose::pose_from_file(structure_2, structure_file_2, false, core::import_pose::PDB_file);

    cout.precision(9);
    cout << fixed << endl;

    cout << "CA RMSD: " << core::scoring::CA_rmsd(structure_1, structure_2) << endl;
    cout << "AA RMSD: " << core::scoring::all_atom_rmsd(structure_1, structure_2) << endl;
    cout << "GDT-TM: " << core::scoring::CA_gdtmm(structure_1, structure_2) << endl;

    return 0;
}