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
    constexpr uint8_t omega = 3;
    constexpr uint8_t chi = 4;
}

vector<uint8_t> angle_types;
vector<size_t> residue_indexes;
vector<size_t> residue_chi_num;

string ss;
vector<vector<double>> ss_prob;

vector<vector<vector<vector<double>>>> frag3, frag9;

double max_energy = 0;

void getSS(string protein_name) {
    ifstream ss_in("proteins/" + protein_name + "/ss2");
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

double fix_angle(double angle) {
    return angle + (angle < -180 ? 360 : (angle > 180 ? -360 : 0));
}

void getFragments3(string protein_name) {
    ifstream frag3_in("proteins/" + protein_name + "/frag3");

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

                frag.push_back({fix_angle(phi), fix_angle(psi), fix_angle(omega)});
            }
            getline(frag3_in, line);
            frags.push_back(frag);
        }
        frag3.push_back(frags);
    }
}

void getFragments9(string protein_name) {
    ifstream frag9_in("proteins/" + protein_name + "/frag9");

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

                frag.push_back({fix_angle(phi), fix_angle(psi), fix_angle(omega)});
            }
            getline(frag9_in, line);
            frags.push_back(frag);
        }
        frag9.push_back(frags);
    }
}

class RosettaFaEnergyFunction : public evo_alg::FitnessFunction<double> {
  public:
    static RosettaFaEnergyFunction* create(string sequence) {
        core::pose::Pose protein_structure;
        core::pose::make_pose_from_sequence(
            protein_structure, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        size_t num_angles = 0;

        for (size_t residue_idx = 1; residue_idx <= protein_structure.total_residue(); ++residue_idx) {
            num_angles += 3 + residue_chi_num[residue_idx - 1];
        }

        return new RosettaFaEnergyFunction(num_angles, protein_structure);
    }

    size_t getDimension() const override {
        return dimension_;
    }

    RosettaFaEnergyFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    fitness_t operator()(evo_alg::Genotype<double> const& genotype) override {
        vector<double> chromosome = genotype.getChromosome();
        size_t angle_idx = 0;
        size_t residue_num = protein_structure_.total_residue();

        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            protein_structure_.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            protein_structure_.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            protein_structure_.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                protein_structure_.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
        }

        double rosetta_score = energy_score_->score(protein_structure_);

        return {-rosetta_score};
    }

    const core::pose::Pose& getProteinStructure() const {
        return protein_structure_;
    }

  private:
    RosettaFaEnergyFunction(size_t num_angles, core::pose::Pose& protein_structure)
        : FitnessFunction<double>({num_angles, {-180, 180}}), protein_structure_{protein_structure},
          energy_score_{core::scoring::ScoreFunctionFactory::create_score_function("ref2015")} {};

    size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

class RosettaCentroidEnergyFunction : public evo_alg::FitnessFunction<double> {
  public:
    static RosettaCentroidEnergyFunction* create(string sequence) {
        core::pose::Pose protein_structure;
        core::pose::make_pose_from_sequence(
            protein_structure, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));
        size_t num_angles = 0;

        for (size_t residue_idx = 1; residue_idx <= protein_structure.total_residue(); ++residue_idx) {
            num_angles += 3 + residue_chi_num[residue_idx - 1];
        }

        return new RosettaCentroidEnergyFunction(num_angles, protein_structure);
    }

    size_t getDimension() const override {
        return dimension_;
    }

    RosettaCentroidEnergyFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    fitness_t operator()(evo_alg::Genotype<double> const& genotype) override {
        vector<double> chromosome = genotype.getChromosome();
        size_t angle_idx = 0;
        size_t residue_num = protein_structure_.total_residue();

        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t res_chromosome_index = residue_indexes[residue_idx - 1];
            protein_structure_.set_phi(residue_idx, chromosome[res_chromosome_index]);
            protein_structure_.set_psi(residue_idx, chromosome[res_chromosome_index + 1]);
            protein_structure_.set_omega(residue_idx, chromosome[res_chromosome_index + 2]);
        }

        core::scoring::dssp::Dssp dssp(protein_structure_);

        string protein_ss = dssp.get_dssp_secstruct();

        double ss_score = 0, ss_total = 0;
        for (size_t ss_index = 0; ss_index < protein_ss.size(); ++ss_index) {
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

        double result = (rosetta_score + pow(ss_score, 3)) / 2;

        return {result};
    }

    const core::pose::Pose& getProteinStructure() const {
        return protein_structure_;
    }

  private:
    RosettaCentroidEnergyFunction(size_t num_angles, core::pose::Pose& protein_structure)
        : FitnessFunction<double>({num_angles, {-180, 180}}), protein_structure_{protein_structure},
          energy_score_{core::scoring::ScoreFunctionFactory::create_score_function("score3")} {};

    size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

size_t tournamentSelection(std::vector<double> const& individuals_fit) {
    return evo_alg::selector::tournament(individuals_fit, 2);
}

evo_alg::real_individual_t fragPolynomialMutation(evo_alg::real_individual_t const& individual, double const pr) {
    evo_alg::real_individual_t mutated_individual(individual);

    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = mutated_individual.getChromosome();
    uint32_t n = 20;

    std::vector<std::pair<double, double>> bounds = individual.getBounds();
    for (size_t index = 0; index < mutated_chromosome.size(); ++index) {
        if (angle_types[index] != angle_type::chi)
            continue;

        double const cur_pr = evo_alg::utils::uniformProbGen();
        if (cur_pr < pr) {
            double const u = evo_alg::utils::uniformProbGen();
            double const delta = u < 0.5 ? pow(2 * u, 1.0 / (n + 1)) - 1 : 1 - pow((2 * (1 - u)), 1 / (n + 1));
            mutated_chromosome[index] += delta * (u < 0.5 ? individual_chromosome[index] - bounds[index].first
                                                          : bounds[index].second - individual_chromosome[index]);
        }

        assert(evo_alg::utils::numericGreaterEqual(mutated_chromosome[index], bounds[index].first));
        assert(evo_alg::utils::numericLowerEqual(mutated_chromosome[index], bounds[index].second));
    }
    mutated_individual.setChromosome(mutated_chromosome);

    return mutated_individual;
}

void fragmentInsertion(size_t index, vector<double>& residue_angles, vector<vector<vector<double>>>& fragment_list) {
    uniform_int_distribution<size_t> index_gen(0, fragment_list.size() - 1);
    size_t selected_frag = index_gen(evo_alg::utils::rng);
    vector<vector<double>> fragment = fragment_list[selected_frag];
    for (size_t frag_index = 0; frag_index < fragment.size(); ++frag_index) {
        size_t res_index = residue_indexes[index + frag_index];
        residue_angles[res_index] = fragment[frag_index][0];
        residue_angles[res_index + 1] = fragment[frag_index][1];
        residue_angles[res_index + 2] = fragment[frag_index][2];
    }
}

evo_alg::real_individual_t frag3Mutation(evo_alg::real_individual_t const& individual, double const pr) {
    evo_alg::real_individual_t mutated_individual(individual);
    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = mutated_individual.getChromosome();

    for (size_t index = 0; index < residue_indexes.size() - 2; ++index) {
        double const cur_pr = evo_alg::utils::uniformProbGen();
        if (cur_pr < pr) {
            fragmentInsertion(index, mutated_chromosome, frag3[index]);
        }
    }

    mutated_individual.setChromosome(mutated_chromosome);

    return mutated_individual;
}

evo_alg::real_individual_t frag9Mutation(evo_alg::real_individual_t const& individual, double const pr) {
    evo_alg::real_individual_t mutated_individual(individual);
    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = mutated_individual.getChromosome();

    for (size_t index = 0; index < residue_indexes.size() - 8; ++index) {
        double const cur_pr = evo_alg::utils::uniformProbGen();
        if (cur_pr < pr) {
            fragmentInsertion(index, mutated_chromosome, frag9[index]);
        }
    }

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

        string ss_types = "LHE";
        char ss_type = 'E';
        double pr = evo_alg::utils::uniformProbGen();
        double cur_pr = 0;
        for (size_t ss_index = 0; ss_index < 3; ++ss_index) {
            cur_pr += ss_prob[residue_idx][ss_index];
            if (cur_pr > pr) {
                ss_type = ss_types[ss_index];
                break;
            }
        }

        if (ss_type == 'H') {
            if (angle_types[index] == angle_type::phi) {
                bounds2[index] = {-67, -47};
            } else if (angle_types[index] == angle_type::psi) {
                bounds2[index] = {-57, -37};
            }
        } else if (ss_type == 'E') {
            if (angle_types[index] == angle_type::phi) {
                bounds2[index] = {-130, -110};
            } else if (angle_types[index] == angle_type::psi) {
                bounds2[index] = {110, 130};
            }
        }
    }

    auto generateAngles = [&population, &fitness, &chromosome_size](vector<pair<double, double>> const& bounds,
                                                                    size_t start, size_t end) {
        for (size_t index = start; index < end; ++index) {
            vector<double> new_chromosome(chromosome_size);

            for (size_t gene_index = 0; gene_index < chromosome_size; ++gene_index) {
                if (angle_types[gene_index] == angle_type::omega) {
                    new_chromosome[gene_index] = 180;
                } else {
                    std::uniform_real_distribution<double> geneGenerator(bounds[gene_index].first,
                                                                         bounds[gene_index].second);
                    new_chromosome[gene_index] = geneGenerator(evo_alg::utils::rng);
                }
            }

            population[index] = {fitness, new_chromosome};
        }
    };

    generateAngles(bounds, 0, (2 * size) / 3);
    generateAngles(bounds2, (2 * size) / 3, size);

    return population;
}

void repackProtein(core::pose::Pose& pose, core::scoring::ScoreFunctionOP& scorefnx) {
    core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(pose);
    task->restrict_to_repacking();
    core::pack::pack_rotamers(pose, *scorefnx, task);
}

int main(int argc, char** argv) {
    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string protein_name = argv[1];

    ifstream fasta_in("proteins/" + protein_name + "/fasta");

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

    core::pose::Pose pose, faPose;

    core::pose::make_pose_from_sequence(
        pose, sequence, *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));

    core::pose::make_pose_from_sequence(
        faPose, sequence,
        *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));

    size_t residue_num = pose.total_residue();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx)
        residue_chi_num.push_back(0); // faPose.residue(residue_idx).nchi());

    size_t res_index = 0;
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        residue_indexes.push_back(res_index);
        angle_types.push_back(angle_type::phi);
        angle_types.push_back(angle_type::psi);
        angle_types.push_back(angle_type::omega);

        pose.set_phi(residue_idx, 0);
        pose.set_psi(residue_idx, 0);
        pose.set_omega(residue_idx, 180);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            angle_types.push_back(angle_type::chi);

        res_index += 3 + residue_chi_num[residue_idx - 1];
    }

    evo_alg::FitnessFunction<double>::shared_ptr fit_centroid(RosettaCentroidEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_centroid =
        core::scoring::ScoreFunctionFactory::create_score_function("score3");
    max_energy = energy_score_centroid->score(pose);

    evo_alg::FitnessFunction<double>::shared_ptr fit_fa(RosettaFaEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_fa =
        core::scoring::ScoreFunctionFactory::create_score_function("ref2015");

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    size_t iteration_num = 1000;
    size_t it_explore = iteration_num * 0.7, it_exploit = (iteration_num - it_explore) * 1,
           it_refine = iteration_num - it_explore - it_exploit;
    size_t pop_size = 200;
    double cross_prob = 0.95;
    double mut_prob = 1.0 / residue_num;
    double gen_gap = 0.5;

    evo_alg::Population<evo_alg::Individual<double>> pop;
    evo_alg::Individual<double> best_ind;
    vector<double> best_fit, mean_fit, diversity;
    tie(best_ind, pop, best_fit, mean_fit, diversity) =
        evo_alg::ga<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(
            it_explore, pop_size, 0, fit_centroid, proteinInit, evo_alg::selector::roulette,
            evo_alg::recombinator::onePoint<double>, cross_prob, frag9Mutation, mut_prob, gen_gap, pop, 1, NAN, 0);

    const std::vector<size_t> sorted_individuals = pop.getSortedIndividuals();
    pop[sorted_individuals[pop_size - 1]] = best_ind;

    vector<double> best_fit_2, mean_fit_2, diversity_2;
    tie(best_ind, pop, best_fit_2, mean_fit_2, diversity_2) =
        evo_alg::ga<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(
            it_exploit, pop_size, 1, fit_centroid, proteinInit, tournamentSelection,
            evo_alg::recombinator::onePoint<double>, cross_prob, frag3Mutation, mut_prob / 3, 1, pop, 1, NAN, 0);

    for (size_t index = 0; index < best_fit_2.size(); ++index) {
        best_fit.push_back(best_fit_2[index]);
        mean_fit.push_back(mean_fit_2[index]);
        diversity.push_back(diversity_2[index]);
    }

    vector<double> chromosome = best_ind.getChromosome();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
        pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
        pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

        faPose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
        faPose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
        faPose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            faPose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
    }

    cout.precision(9);
    cout << fixed << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << " "
         << energy_score_centroid->score(pose) << endl
         << endl;

    size_t pose_start = argc > 2 ? stoi(argv[2]) : 1;
    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
        repackProtein(faPose, energy_score_fa);
        cout << endl << core::scoring::all_atom_rmsd(faPose, native_structure) << endl << endl;
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
    } else {
        repackProtein(faPose, energy_score_fa);
        cout << endl << core::scoring::all_atom_rmsd(faPose, native_structure) << endl << endl;
    }

    double best_rmsd = 1e9;
    size_t best_structure = -1;
    for (size_t ind_index = 0; ind_index < pop.getSize(); ++ind_index) {
        vector<double> chromosome = pop[ind_index].getChromosome();
        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            faPose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            faPose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            faPose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                faPose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
        }
        double rmsd;

        if (pose_start - 1 > 0) {
            repackProtein(faPose, energy_score_fa);
            faPose.delete_residue_range_slow(1, pose_start - 1);
            rmsd = core::scoring::all_atom_rmsd(faPose, native_structure);
            core::pose::make_pose_from_sequence(
                faPose, sequence,
                *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        } else {
            repackProtein(faPose, energy_score_fa);
            rmsd = core::scoring::all_atom_rmsd(faPose, native_structure);
        }

        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
            best_structure = ind_index;
        }
    }

    chromosome = pop[best_structure].getChromosome();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        faPose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
        faPose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
        faPose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            faPose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
    }

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    core::scoring::dssp::Dssp dssp_native(native_structure);
    core::scoring::dssp::Dssp dssp(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp.get_dssp_secstruct() << endl << endl;

    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl;

    cout << "Full Atom energy: " << energy_score_fa->score(native_structure) << " " << energy_score_fa->score(faPose)
         << endl
         << endl;

    repackProtein(faPose, energy_score_fa);

    core::scoring::dssp::Dssp dssp2(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp2.get_dssp_secstruct() << endl << endl;

    cout << "Full Atom energy: " << energy_score_fa->score(native_structure) << " " << energy_score_fa->score(faPose)
         << endl
         << endl;

    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl;

    faPose.dump_pdb("res.pdb");

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