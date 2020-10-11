#include "protein_commons.hpp"

#include <evo_alg/algorithms/ga.hpp>

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

        double ss_score, ss_total;
        tie(ss_score, ss_total) = ssScore(protein_structure_);
        ss_score = ss_score / ss_total;

        double cm_score, cm_total;
        tie(cm_score, cm_total) = cmScore(protein_structure_);
        cm_score = cm_score / cm_total;

        double rosetta_score = energy_score_->score(protein_structure_);
        rosetta_score = 1 - rosetta_score / max_energy;

        double result = (rosetta_score + pow(ss_score, 10)) / 2;

        return {result};
    }

    vector<double> getScores(core::pose::Pose& pose) {
        double ss_score, ss_total;
        tie(ss_score, ss_total) = ssScore(protein_structure_);

        double cm_score, cm_total;
        tie(cm_score, cm_total) = cmScore(protein_structure_);

        double rosetta_score = energy_score_->score(pose);

        return {rosetta_score, ss_score, ss_total, cm_score, cm_total};
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

void residueOnePointCrossover(evo_alg::real_individual_t const& parent_1, evo_alg::real_individual_t const& parent_2,
                              evo_alg::real_individual_t& child_1, evo_alg::real_individual_t& child_2) {
    evo_alg::real_chromosome_t parent_1_chromosome = parent_1.getChromosome();
    evo_alg::real_chromosome_t parent_2_chromosome = parent_2.getChromosome();
    evo_alg::real_chromosome_t child_1_chromosome = parent_1_chromosome;
    evo_alg::real_chromosome_t child_2_chromosome = parent_2_chromosome;

    size_t const chromosome_size = parent_1_chromosome.size();
    size_t const residue_num = residue_indexes.size();
    std::uniform_int_distribution<size_t> point_dist(1, residue_num - 1);
    size_t const cut_point = residue_indexes[point_dist(evo_alg::utils::rng)];
    for (size_t index = cut_point; index < chromosome_size; ++index) {
        child_1_chromosome[index] = parent_2_chromosome[index];
        child_2_chromosome[index] = parent_1_chromosome[index];
    }
    child_1.setChromosome(child_1_chromosome);
    child_2.setChromosome(child_2_chromosome);
}

void fragmentInsertion(size_t index, vector<double>& residue_angles, vector<vector<vector<double>>>& fragment_list,
                       vector<double>& fragment_prob_list) {
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

void frag3Mutation(evo_alg::real_individual_t& individual, double const pr) {
    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = individual_chromosome;

    for (size_t index = 0; index < residue_indexes.size() - 2; ++index) {
        double const cur_pr = evo_alg::utils::uniformProbGen();
        if (cur_pr < pr) {
            fragmentInsertion(index, mutated_chromosome, frag3[index], frag3_prob[index]);
        }
    }

    individual.setChromosome(mutated_chromosome);
}

void frag9Mutation(evo_alg::real_individual_t& individual, double const pr) {
    evo_alg::real_chromosome_t individual_chromosome = individual.getChromosome();
    evo_alg::real_chromosome_t mutated_chromosome = individual_chromosome;

    for (size_t index = 0; index < residue_indexes.size() - 8; ++index) {
        double const cur_pr = evo_alg::utils::uniformProbGen();
        if (cur_pr < pr) {
            fragmentInsertion(index, mutated_chromosome, frag9[index], frag9_prob[index]);
        }
    }

    individual.setChromosome(mutated_chromosome);
}

void proteinInit(evo_alg::Population<evo_alg::Individual<double>>& population,
                 typename evo_alg::FitnessFunction<double>::const_shared_ptr fitness, size_t const size) {
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

    auto generateAngles = [&population, &chromosome_size](vector<pair<double, double>> const& bounds, size_t start,
                                                          size_t end) {
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
            population[index].setChromosome(new_chromosome);
        }
    };

    generateAngles(bounds2, 0, size * 0.5);
    generateAngles(bounds, size * 0.5, size);
}

void repackProtein(core::pose::Pose& pose, core::scoring::ScoreFunctionOP& scorefnx) {
    core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(pose);
    task->restrict_to_repacking();
    core::pack::pack_rotamers(pose, *scorefnx, task);
}

// argv -> 1: protein name | 2: start offset of pdb structure (default 1) | 3: rmsd file name
int main(int argc, char** argv) {
    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string protein_name = argv[1];
    size_t pose_start = argc > 2 ? stoi(argv[2]) : 1;

    ifstream fasta_in("proteins/" + protein_name + "/fasta");

    string protein_desc, sequence;
    getline(fasta_in, protein_desc);
    getline(fasta_in, sequence);

    core::pose::Pose pose, faPose;

    core::pose::make_pose_from_sequence(
        pose, sequence, *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));

    core::pose::make_pose_from_sequence(
        faPose, sequence,
        *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));

    size_t residue_num = pose.total_residue();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx)
        residue_chi_num.push_back(faPose.residue(residue_idx).nchi());

    size_t res_index = 0;
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        residue_indexes.push_back(res_index);
        angle_types.push_back(angle_type::phi);
        angle_types.push_back(angle_type::psi);
        angle_types.push_back(angle_type::omega);

        pose.set_phi(residue_idx, 0);
        pose.set_psi(residue_idx, 0);
        pose.set_omega(residue_idx, 180);

        faPose.set_phi(residue_idx, 0);
        faPose.set_psi(residue_idx, 0);
        faPose.set_omega(residue_idx, 180);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            angle_types.push_back(angle_type::chi);

        res_index += 3 + residue_chi_num[residue_idx - 1];
    }

    getSS(protein_name);
    getFragments(protein_name, 3, frag3, frag3_prob);
    getFragments(protein_name, 9, frag9, frag9_prob);
    getContactMap(protein_name);

    evo_alg::FitnessFunction<double>::shared_ptr fit_centroid(RosettaCentroidEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_centroid =
        core::scoring::ScoreFunctionFactory::create_score_function("score3");
    max_energy = energy_score_centroid->score(pose);

    evo_alg::FitnessFunction<double>::shared_ptr fit_fa(RosettaFaEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_fa =
        core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    max_energy_fa = energy_score_fa->score(faPose);

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::dssp::Dssp dssp_native(native_structure);

    cout << endl << "Protein: " << endl;
    cout << protein_desc << endl;
    cout << sequence << endl << endl;

    cout << "Sec. Struct: " << ss << endl << endl;

    size_t iteration_num = 1000;
    size_t it_explore = iteration_num * 0.5, it_exploit = (iteration_num - it_explore) * 1.0,
           it_refine = iteration_num - it_explore - it_exploit;
    size_t pop_size = 500;
    double cross_prob = 0.95;
    double mut_prob = 1.0 / residue_num;
    double gen_gap = 0.5, gap_inc = 0.5 / it_explore;
    double scale_c = 1.2, scale_inc = 0.8 / it_explore;

    evo_alg::Population<evo_alg::Individual<double>> pop;
    evo_alg::Individual<double> best_ind;
    vector<double> best_fit, mean_fit, diversity;
    tie(best_ind, pop, best_fit, mean_fit, diversity) =
        evo_alg::ga<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(
            it_explore, pop_size, 0, fit_centroid, proteinInit, evo_alg::selector::roulette, residueOnePointCrossover,
            cross_prob, frag9Mutation, mut_prob, pop, gen_gap, 0, scale_c, scale_inc, 1, NAN, 0);

    const std::vector<size_t> sorted_individuals = pop.getSortedIndividuals();
    pop[sorted_individuals[pop_size - 1]] = best_ind;

    vector<double> best_fit_2, mean_fit_2, diversity_2;
    tie(best_ind, pop, best_fit_2, mean_fit_2, diversity_2) =
        evo_alg::ga<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(
            it_exploit, pop_size, 1, fit_centroid, proteinInit, tournamentSelection, residueOnePointCrossover,
            cross_prob, frag3Mutation, mut_prob / 3, pop, 1, 0, 0, 0, 1, NAN, 0);

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

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    core::scoring::dssp::Dssp dssp0(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp0.get_dssp_secstruct() << endl << endl;

    vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);

    cout << "Energy: " << scores[0] << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    repackProtein(faPose, energy_score_fa);
    cout << endl
         << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
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

    ofstream ca_rmsd_out(argc > 3 ? string(argv[3]) + string(".ca_rmsd") : "res.ca_rmsd", ios::app);
    ca_rmsd_out << fixed;
    ca_rmsd_out.precision(9);
    ca_rmsd_out << core::scoring::CA_rmsd(faPose, native_structure) << endl;

    ofstream aa_rmsd_out(argc > 3 ? string(argv[3]) + string(".aa_rmsd") : "res.aa_rmsd", ios::app);
    aa_rmsd_out << fixed;
    aa_rmsd_out.precision(9);
    aa_rmsd_out << core::scoring::all_atom_rmsd(faPose, native_structure) << endl;

    return 0;
}