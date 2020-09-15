#include "rosetta.hpp"

#include <evo_alg/algorithms/ga.hpp>
#include <evo_alg/core.hpp>
#include <evo_alg/operators.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

using namespace std;

namespace angle_type {
    constexpr uint8_t phi = 1;
    constexpr uint8_t psi = 2;
    constexpr uint8_t chi = 3;
}

vector<uint8_t> angle_types;
vector<size_t> residue_indexes;

string ss;
vector<vector<double>> ss_prob;

vector<vector<vector<vector<double>>>> frag3, frag9;

double max_energy = 0;

class RosettaEnergyFunction : public evo_alg::FitnessFunction<double> {
  public:
    static RosettaEnergyFunction* create(string sequence) {
        core::pose::Pose protein_structure;
        core::pose::make_pose_from_sequence(
            protein_structure, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));
        size_t num_angles = 0;

        for (size_t residue_idx = 1; residue_idx <= protein_structure.total_residue(); ++residue_idx) {
            protein_structure.set_omega(residue_idx, 180);
            num_angles += 2; // + protein_structure.residue(residue_idx).nchi();
        }

        return new RosettaEnergyFunction(num_angles, protein_structure);
    }

    size_t getDimension() const override {
        return dimension_;
    }

    RosettaEnergyFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    double scoreProtein(core::pose::Pose& pose) {
        return energy_score_->score(pose);
    }

    fitness_t operator()(evo_alg::Genotype<double> const& genotype) override {
        vector<double> chromosome = genotype.getChromosome();
        size_t angle_idx = 0;
        size_t residue_num = protein_structure_.total_residue();

        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            protein_structure_.set_phi(residue_idx, chromosome[angle_idx++]);
            protein_structure_.set_psi(residue_idx, chromosome[angle_idx++]);

            // size_t chi_num = protein_structure_.residue(residue_idx).nchi();
            // for (size_t chi_idx = 1; chi_idx <= chi_num; ++chi_idx)
            //     protein_structure_.set_chi(chi_idx, residue_idx, chromosome[angle_idx++]);
        }

        core::scoring::dssp::Dssp dssp(protein_structure_);

        string protein_ss = dssp.get_dssp_secstruct();

        double ss_score = 0, ss_total = 0;
        for (size_t ss_index = 0; ss_index < protein_ss.size(); ++ss_index) {
            // if (protein_ss[ss_index] != ss[ss_index])
            //     ++ss_score;
            if (protein_ss[ss_index] == 'L')
                ss_score += ss_prob[ss_index][0];
            if (protein_ss[ss_index] == 'H')
                ss_score += ss_prob[ss_index][1];
            if (protein_ss[ss_index] == 'E')
                ss_score += ss_prob[ss_index][2];

            if (ss[ss_index] == 'L')
                ss_total += ss_prob[ss_index][0];
            if (ss[ss_index] == 'H')
                ss_total += ss_prob[ss_index][1];
            if (ss[ss_index] == 'E')
                ss_total += ss_prob[ss_index][2];
        }
        ss_score = ss_score / ss_total;

        double rosetta_score = energy_score_->score(protein_structure_);
        max_energy = max(max_energy, rosetta_score);
        rosetta_score = 1 - rosetta_score / max_energy;

        double result = (rosetta_score + pow(ss_score, 1)) / 2;

        return {result};
    }

    const core::pose::Pose& getProteinStructure() const {
        return protein_structure_;
    }

  private:
    RosettaEnergyFunction(size_t num_angles, core::pose::Pose& protein_structure)
        : FitnessFunction<double>({num_angles, {-180, 180}}), protein_structure_{protein_structure},
          energy_score_{core::scoring::ScoreFunctionFactory::create_score_function("score3")} {};

    size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

evo_alg::real_individual_t polynomialMutation(evo_alg::real_individual_t const& individual, double const pr) {
    return evo_alg::mutator::polynomial(individual, pr, 60);
}

pair<evo_alg::real_individual_t, evo_alg::real_individual_t> sbxCrossover(evo_alg::real_individual_t const& parent_1,
                                                                          evo_alg::real_individual_t const& parent_2) {
    return evo_alg::recombinator::sbx(parent_1, parent_2, 2, 0.5);
}

size_t tournamentSelection(std::vector<double> const& individuals_fit) {
    return evo_alg::selector::tournament(individuals_fit, 2);
}

evo_alg::real_individual_t fragmentMutation(evo_alg::real_individual_t const& individual, double const pr) {
    evo_alg::real_individual_t mutated_individual(individual);
    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = mutated_individual.getChromosome();

    auto mutateFrag3 = [&mutated_chromosome, pr]() {
        for (size_t index = 0; index < residue_indexes.size() - 2; ++index) {
            double const cur_pr = evo_alg::utils::uniformProbGen();
            if (cur_pr < pr) {
                // cout << "Doing mutation: " << index << endl;
                // cout << "\t" << frag3[index].size() << " fragments" << endl;
                uniform_int_distribution<size_t> index_gen(0, frag3[index].size() - 1);
                size_t selected_frag = index_gen(evo_alg::utils::rng);
                // cout << "\t"
                //      << "Selected " << selected_frag << endl;
                vector<vector<double>> fragment = frag3[index][selected_frag];
                // cout << "\t"
                //      << "Fragment: ";
                // for (auto v : fragment) {
                //     cout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ") ";
                // }
                // cout << endl;
                for (size_t frag_index = 0; frag_index < 3; ++frag_index) {
                    size_t res_index = residue_indexes[index + frag_index];
                    // cout << "\t\tChanging residue " << res_index << endl;
                    mutated_chromosome[res_index] = fragment[frag_index][0];
                    mutated_chromosome[res_index + 1] = fragment[frag_index][1];
                }
                // cout << endl;
            }
        }
    };

    auto mutateFrag9 = [&mutated_chromosome, pr]() {
        for (size_t index = 0; index < residue_indexes.size() - 8; ++index) {
            double const cur_pr = evo_alg::utils::uniformProbGen();
            if (cur_pr < pr) {
                uniform_int_distribution<size_t> index_gen(0, frag9[index].size() - 1);
                vector<vector<double>> fragment = frag9[index][index_gen(evo_alg::utils::rng)];
                for (size_t frag_index = 0; frag_index < 9; ++frag_index) {
                    size_t res_index = residue_indexes[index + frag_index];
                    mutated_chromosome[res_index] = fragment[frag_index][0];
                    mutated_chromosome[res_index + 1] = fragment[frag_index][1];
                }
            }
        }
    };

    mutateFrag3();

    mutated_individual.setChromosome(mutated_chromosome);

    return mutated_individual;
}

evo_alg::Population<evo_alg::Individual<double>>
proteinInit(typename evo_alg::FitnessFunction<double>::const_shared_ptr fitness, size_t const size) {
    evo_alg::Population<evo_alg::Individual<double>> population(size);
    vector<pair<double, double>> bounds = fitness->getBounds();
    vector<pair<double, double>> bounds2(bounds);
    size_t const chromosome_size = bounds.size();

    size_t residue_idx = -1;
    for (size_t index = 0; index < chromosome_size; ++index) {
        if (angle_types[index] == angle_type::phi)
            ++residue_idx;

        if (ss[residue_idx] == 'H') {
            if (angle_types[index] == angle_type::phi) {
                bounds2[index] = {-67, -47};
            } else if (angle_types[index] == angle_type::psi) {
                bounds2[index] = {-57, -37};
            }
        } else if (ss[residue_idx] == 'E') {
            if (angle_types[index] == angle_type::phi) {
                bounds2[index] = {-130, -110};
            } else if (angle_types[index] == angle_type::psi) {
                bounds2[index] = {110, 130};
            }
        }
    }

    // for (size_t index = 0; index < size; ++index) {
    //     vector<double> new_chromosome(chromosome_size);

    //     for (size_t gene_index = 0; gene_index < chromosome_size; ++gene_index) {
    //         std::uniform_real_distribution<double> geneGenerator(bounds[gene_index].first,
    //         bounds[gene_index].second); new_chromosome[gene_index] = geneGenerator(evo_alg::utils::rng);
    //     }

    //     population[index] = {fitness, new_chromosome};
    // }

    for (size_t index = 0; index < size; ++index) {
        vector<double> new_chromosome(chromosome_size);

        for (size_t gene_index = 0; gene_index < chromosome_size; ++gene_index) {
            std::uniform_real_distribution<double> geneGenerator(bounds2[gene_index].first, bounds2[gene_index].second);
            new_chromosome[gene_index] = geneGenerator(evo_alg::utils::rng);
        }

        population[index] = {fitness, new_chromosome};
    }

    return population;
}

void getSS(string protein_name) {
    ifstream ss_in("proteins/" + protein_name + ".ss2");
    string line;

    getline(ss_in, line);
    getline(ss_in, line);

    int index;
    string residue, res_ss;
    double c_prob, h_prob, e_prob;
    while (ss_in >> index >> residue >> res_ss >> c_prob >> h_prob >> e_prob) {
        if (res_ss == "C")
            res_ss = "L";
        ss += res_ss;
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
}

void getFragments3(string protein_name) {
    ifstream frag3_in("proteins/" + protein_name + ".frag3");

    string _, line;
    size_t pos, n;
    while (getline(frag3_in, line)) {
        istringstream iss(line);
        iss >> _ >> pos >> _ >> n;
        getline(frag3_in, line);
        vector<vector<vector<double>>> frags;
        for (size_t index = 0; index < n; ++index) {
            vector<vector<double>> frag;
            for (size_t index2 = 0; index2 < 3; ++index2) {
                getline(frag3_in, line);
                istringstream iss2(line);
                double phi, psi, omega;

                iss2 >> _ >> _ >> _ >> _ >> _;
                iss2 >> phi >> psi >> omega;

                frag.push_back({phi, psi, omega});
            }
            getline(frag3_in, line);
            frags.push_back(frag);
        }
        frag3.push_back(frags);
    }
}

void getFragments9(string protein_name) {
    ifstream frag9_in("proteins/" + protein_name + ".frag9");

    string _, line;
    size_t pos, n;
    while (getline(frag9_in, line)) {
        istringstream iss(line);
        iss >> _ >> pos >> _ >> n;
        getline(frag9_in, line);
        vector<vector<vector<double>>> frags;
        for (size_t index = 0; index < n; ++index) {
            vector<vector<double>> frag;
            for (size_t index2 = 0; index2 < 9; ++index2) {
                getline(frag9_in, line);
                istringstream iss2(line);
                double phi, psi, omega;

                iss2 >> _ >> _ >> _ >> _ >> _;
                iss2 >> phi >> psi >> omega;

                frag.push_back({phi, psi, omega});
            }
            getline(frag9_in, line);
            frags.push_back(frag);
        }
        frag9.push_back(frags);
    }
}

int main(int argc, char** argv) {
    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string protein_name = argv[1];

    ifstream fasta_in("proteins/" + protein_name + ".fasta");

    string protein_desc, sequence;
    getline(fasta_in, protein_desc);
    getline(fasta_in, sequence);

    getSS(protein_name);
    getFragments3(protein_name);
    getFragments9(protein_name);

    cout << endl << "Protein: " << endl;
    cout << protein_desc << endl;
    cout << sequence << endl << endl;

    cout << "Sec. Struct: " << ss << endl << endl;

    core::pose::Pose pose;

    core::pose::make_pose_from_sequence(
        pose, sequence, *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));

    size_t residue_num = pose.total_residue();
    size_t res_index = 0;
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        residue_indexes.push_back(res_index);
        angle_types.push_back(angle_type::phi);
        angle_types.push_back(angle_type::psi);

        pose.set_phi(residue_idx, 0);
        pose.set_psi(residue_idx, 0);
        pose.set_omega(residue_idx, 180);

        // size_t chi_num = pose.residue(residue_idx).nchi();
        // for (size_t chi_idx = 1; chi_idx <= chi_num; ++chi_idx)
        //     angle_types.push_back(angle_type::chi);

        res_index += 2; // + chi_num;
    }

    evo_alg::FitnessFunction<double>::shared_ptr fit(RosettaEnergyFunction::create(sequence));
    max_energy = dynamic_pointer_cast<RosettaEnergyFunction>(fit)->scoreProtein(pose);

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, "proteins/" + protein_name + ".pdb", false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, "proteins/" + protein_name + ".pdb", false,
                                      core::import_pose::PDB_file);

    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::ScoreFunctionOP energy_score = core::scoring::ScoreFunctionFactory::create_score_function("score3");

    evo_alg::Population<evo_alg::Individual<double>> pop;
    evo_alg::Individual<double> best_ind;
    vector<double> best_fit, mean_fit, diversity;
    tie(best_ind, pop, best_fit, mean_fit, diversity) =
        evo_alg::ga<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(
            1000, 100, 0, fit, proteinInit, evo_alg::selector::roulette, evo_alg::recombinator::onePoint<double>, 0.9,
            polynomialMutation, 1.0 / residue_num, 1);

    shared_ptr<const RosettaEnergyFunction> fit2 =
        dynamic_pointer_cast<const RosettaEnergyFunction>(best_ind.getFitnessFunction());

    vector<double> chromosome = best_ind.getChromosome();
    size_t angle_idx = 0;
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        pose.set_phi(residue_idx, chromosome[angle_idx++]);
        pose.set_psi(residue_idx, chromosome[angle_idx++]);

        // size_t chi_num = pose.residue(residue_idx).nchi();
        // for (size_t chi_idx = 1; chi_idx <= chi_num; ++chi_idx)
        //     pose.set_chi(chi_idx, residue_idx, chromosome[angle_idx++]);
    }

    cout << fixed << endl;
    cout.precision(9);
    cout << energy_score->score(native_structure2) << " " << energy_score->score(pose) << endl << endl;

    core::util::switch_to_residue_type_set(pose, core::chemical::FA_STANDARD);

    core::scoring::dssp::Dssp dssp(pose);
    core::scoring::dssp::Dssp dssp_native(native_structure);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp.get_dssp_secstruct() << endl << endl;

    cout << core::scoring::native_CA_rmsd(native_structure, pose) << endl << endl;

    cout.precision(2);

    for (size_t ind_index = 0; ind_index < pop.getSize(); ++ind_index) {
        vector<double> chromosome = pop[ind_index].getChromosome();
        size_t angle_idx = 0;
        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            pose.set_phi(residue_idx, chromosome[angle_idx++]);
            pose.set_psi(residue_idx, chromosome[angle_idx++]);

            // size_t chi_num = pose.residue(residue_idx).nchi();
            // for (size_t chi_idx = 1; chi_idx <= chi_num; ++chi_idx)
            //     pose.set_chi(chi_idx, residue_idx, chromosome[angle_idx++]);
        }

        cout << core::scoring::native_CA_rmsd(native_structure, pose) << " ";
    }
    cout << endl;

    pose.dump_pdb("res.pdb");

    ofstream fit_out("res.fit"), diver_out("res.diver");
    diver_out << fixed;

    fit_out << fixed;
    fit_out.precision(9);
    for (size_t index = 0; index < best_fit.size(); ++index) {
        fit_out << index << " " << best_fit[index] << " " << mean_fit[index] << endl;
    }

    diver_out << fixed;
    diver_out.precision(9);
    for (size_t index = 0; index < diversity.size(); ++index) {
        diver_out << index << " " << diversity[index] << endl;
    }

    return 0;
}