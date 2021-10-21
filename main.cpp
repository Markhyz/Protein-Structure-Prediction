#include "protein_brkga.hpp"

#include <sstream>

evo_alg::utils::Timer timer;

string leftPadding(size_t size) {
    return string(size, '.');
}

void printProteinInfo(vector<double> chromosome, string name, string sequence, size_t residue_num, size_t pose_start,
                      core::pose::Pose& fa_pose, core::pose::Pose& pose, core::pose::Pose& native_structure,
                      core::pose::Pose& native_structure2, evo_alg::FitnessFunction<double>::shared_ptr& fit_fa,
                      core::scoring::ScoreFunctionOP& energy_score_fa,
                      core::scoring::ScoreFunctionOP& energy_score_centroid) {
    cout << leftPadding(4) << name << endl << endl;

    if (fa_pose.total_residue() < residue_num) {
        core::pose::make_pose_from_sequence(
            fa_pose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
    }

    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
        pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
        pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

        fa_pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
        fa_pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
        fa_pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            fa_pose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
    }

    auto print_phase = [&](string phase) {
        cout << leftPadding(8) << phase << endl << endl;

        core::scoring::sasa::SasaCalc sasa_calc;
        vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(fa_pose);

        cout << leftPadding(12) << "Predicted" << endl;
        cout << leftPadding(16) << "FA Energy: " << scores[0] << endl;
        cout << leftPadding(16) << "Centroid energy: " << energy_score_centroid->score(pose) << endl;
        cout << leftPadding(16) << "SS: " << scores[1] << " / " << scores[2] << endl;
        cout << leftPadding(16) << "CM: " << scores[3] << " / " << scores[4] << endl;
        cout << leftPadding(16) << "SASA: " << sasa_calc.calculate(fa_pose) << endl << endl;

        if (pose_start - 1 > 0) {
            fa_pose.delete_residue_range_slow(1, pose_start - 1);
        }

        if (fa_pose.total_residue() > native_structure.total_residue()) {
            fa_pose.delete_residue_range_slow(native_structure.total_residue() + 1, fa_pose.total_residue());
        }

        cout << leftPadding(16) << "AA_RMSD: " << core::scoring::all_atom_rmsd(fa_pose, native_structure) << endl;
        cout << leftPadding(16) << "CA_RMSD: " << core::scoring::CA_rmsd(fa_pose, native_structure) << endl;
        cout << leftPadding(16) << "GDT-TS: " << core::scoring::CA_gdtmm(fa_pose, native_structure) << endl << endl;
    };

    print_phase("Before repack");

    if (fa_pose.total_residue() < residue_num) {
        core::pose::make_pose_from_sequence(
            fa_pose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));

        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            fa_pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            fa_pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            fa_pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                fa_pose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
        }
    }

    repackProtein(fa_pose, energy_score_fa);

    print_phase("After repack");
}

void outputProteinResults(config_t config, evo_alg::Population<evo_alg::Individual<double>>& population,
                          string method_name, string sequence, size_t residue_num, size_t pose_start,
                          core::pose::Pose& fa_pose, core::pose::Pose& pose, core::pose::Pose& native_structure,
                          core::pose::Pose& native_structure2, evo_alg::FitnessFunction<double>::shared_ptr& fit_fa,
                          core::scoring::ScoreFunctionOP& energy_score_fa,
                          core::scoring::ScoreFunctionOP& energy_score_centroid,
                          core::scoring::dssp::Dssp& dssp_native) {
    vector<size_t> best_frontier_dirty = population.getBestIndividuals();
    vector<size_t> best_frontier;
    for (size_t index_1 : best_frontier_dirty) {
        bool unique = true;
        for (size_t index_2 : best_frontier) {
            evo_alg::fitness::FitnessValue fitness_1 = population[index_1].getFitnessValue();
            evo_alg::fitness::FitnessValue fitness_2 = population[index_2].getFitnessValue();

            for (size_t fit_index = 0; fit_index < fitness_1.getDimension(); ++fit_index) {
                if (fabs(fitness_1[fit_index] - fitness_2[fit_index]) < 1e-9) {
                    unique = false;
                    break;
                }
            }
            if (!unique)
                break;
        }
        if (unique) {
            best_frontier.push_back(index_1);
        }
    }

    auto set_pose = [&](vector<double> chromosome) {
        if (fa_pose.total_residue() < residue_num) {
            core::pose::make_pose_from_sequence(
                fa_pose, sequence,
                *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        }

        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            fa_pose.set_phi(residue_idx, chromosome[chromosome_res_idx]);
            fa_pose.set_psi(residue_idx, chromosome[chromosome_res_idx + 1]);
            fa_pose.set_omega(residue_idx, chromosome[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                fa_pose.set_chi(chi_idx + 1, residue_idx, chromosome[chromosome_res_idx + 3 + chi_idx]);
        }

        if (pose_start - 1 > 0) {
            fa_pose.delete_residue_range_slow(1, pose_start - 1);
        }

        if (fa_pose.total_residue() > native_structure.total_residue()) {
            fa_pose.delete_residue_range_slow(native_structure.total_residue() + 1, fa_pose.total_residue());
        }
    };

    cout << "\nMO-BRKGA/" << method_name << endl << endl;
    cout << leftPadding(4) << "Total individuals: " << best_frontier_dirty.size() << endl;
    cout << leftPadding(4) << "Unique individuals: " << best_frontier.size() << endl << endl;

    double mean_rmsd = 0, mean_gdt = 0;
    double best_gdt = 0, best_rmsd = 1e9;
    size_t best_structure_gdt = -1, best_structure_rmsd = -1;
    for (size_t ind_index = 0; ind_index < population.getSize(); ++ind_index) {
        set_pose(population[ind_index].getChromosome());

        double gdt, rmsd;

        gdt = core::scoring::CA_gdtmm(fa_pose, native_structure);
        rmsd = core::scoring::CA_rmsd(fa_pose, native_structure);

        fa_pose.dump_pdb(config.decoy_output_dir + method_name + string("_decoy_") + to_string(ind_index + 1) +
                         string(".pdb"));

        mean_gdt += gdt;
        mean_rmsd += rmsd;

        if (gdt > best_gdt) {
            best_gdt = gdt;
            best_structure_gdt = ind_index;
        }
        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
            best_structure_rmsd = ind_index;
        }
    }

    mean_gdt /= population.getSize();
    mean_rmsd /= population.getSize();

    core::scoring::sasa::SasaCalc sasa_calc;
    vector<double> scores_native;
    if (native_structure.total_residue() == residue_num) {
        scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);
        cout << leftPadding(4) << "Native" << endl;
        cout << leftPadding(8) << "FA Energy: " << scores_native[0] << endl;
        cout << leftPadding(8) << "Centroid energy: " << energy_score_centroid->score(native_structure2) << endl;
        cout << leftPadding(8) << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
        cout << leftPadding(8) << "CM: " << scores_native[3] << " / " << scores_native[4] << endl;
        cout << leftPadding(8) << "SASA: " << sasa_calc.calculate(native_structure) << endl << endl;
    }

    printProteinInfo(population[best_structure_gdt].getChromosome(), "Best GDT_TS", sequence, residue_num, pose_start,
                     fa_pose, pose, native_structure, native_structure2, fit_fa, energy_score_fa,
                     energy_score_centroid);
    printProteinInfo(population[best_structure_rmsd].getChromosome(), "Best RMSD", sequence, residue_num, pose_start,
                     fa_pose, pose, native_structure, native_structure2, fit_fa, energy_score_fa,
                     energy_score_centroid);

    cout << leftPadding(4) << "Mean RMSD: " << mean_rmsd << endl;
    cout << leftPadding(4) << "Mean GDT_TS: " << mean_gdt << endl << endl;

    core::scoring::dssp::Dssp dssp(fa_pose);

    cout << leftPadding(4) << "Secondary structure" << endl;
    cout << leftPadding(8) << "Native:    " << dssp_native.get_dssp_secstruct() << endl;
    cout << leftPadding(8) << "Input:     " << ss.substr(pose_start - 1) << endl;
    cout << leftPadding(8) << "Predicted: " << dssp.get_dssp_secstruct() << endl << endl;

    ofstream ca_rmsd_out(config.score_output_dir + method_name + string(".ca_rmsd"));
    ofstream aa_rmsd_out(config.score_output_dir + method_name + string(".aa_rmsd"));
    ofstream gdtts_out(config.score_output_dir + method_name + string(".gdtts"));
    ofstream mean_rmsd_out(config.score_output_dir + method_name + string(".mean_rmsd"));
    ofstream mean_gdt_out(config.score_output_dir + method_name + string(".mean_gdt"));

    set_pose(population[best_structure_gdt].getChromosome());

    gdtts_out << fixed;
    gdtts_out.precision(9);
    gdtts_out << core::scoring::CA_gdtmm(fa_pose, native_structure) << endl;

    set_pose(population[best_structure_rmsd].getChromosome());

    ca_rmsd_out << fixed;
    ca_rmsd_out.precision(9);
    ca_rmsd_out << core::scoring::CA_rmsd(fa_pose, native_structure) << endl;

    aa_rmsd_out << fixed;
    aa_rmsd_out.precision(9);
    aa_rmsd_out << core::scoring::all_atom_rmsd(fa_pose, native_structure) << endl;

    mean_rmsd_out << fixed;
    mean_rmsd_out.precision(9);
    mean_rmsd_out << mean_rmsd << endl;

    mean_gdt_out << fixed;
    mean_gdt_out.precision(9);
    mean_gdt_out << mean_gdt << endl;
}

void outputFitness(config_t config, string output_file_name, evo_alg::fitness::frontier_t frontier) {
    ofstream fitness_out(config.algorithm_output_dir + output_file_name + string(".fitness"));
    ofstream normalized_fitness_out(config.algorithm_output_dir + output_file_name + string(".norm_fitness"));

    fitness_out << fixed;
    normalized_fitness_out << fixed;
    fitness_out.precision(9);
    normalized_fitness_out.precision(9);
    for (auto fitness : frontier) {
        vector<double> values = fitness.getValues();
        double ss_value = values[0];
        double cm_value = values[1];
        double energy_value = -values[2];
        fitness_out << energy_value << " " << ss_value << " " << cm_value << endl;

        double norm_energy_value = max(0.0, min(1.0, (energy_value - energy_min) / (energy_max - energy_min)));
        double norm_ss_value = ss_value / ss_total;
        double norm_cm_value = cm_value / cm_total;
        normalized_fitness_out << norm_energy_value << " " << norm_ss_value << " " << norm_cm_value << endl;
    }
}

void outputAlgorithmResults(config_t config,
                            vector<tuple<vector<vector<long double>>, evo_alg::fitness::frontier_t>>& best_frontiers,
                            vector<double>& diversity, string method_name) {
    for (size_t index = 0; index < best_frontiers.size(); ++index) {
        evo_alg::fitness::frontier_t best_frontier = get<1>(best_frontiers[index]);
        string fitness_file_name = method_name + "_gen_" + to_string(index + 1);
        outputFitness(config, fitness_file_name, best_frontier);
    }

    ofstream diversity_out(config.algorithm_output_dir + method_name + ".diversity");
    diversity_out << fixed;
    diversity_out.precision(9);
    for (size_t index = 0; index < diversity.size(); ++index) {
        diversity_out << index << " " << diversity[index] << endl;
    }

    ofstream time_out(config.algorithm_output_dir + method_name + ".time");
    time_out << (int)timer.getTime(method_name) << endl;
}

int main(int argc, char** argv) {
    char* rosetta_argv[] = {argv[0]};
    if (argc < 2) {
        cout << "./mobrkga_pred CONFIG_FILE" << endl;

        return EXIT_FAILURE;
    }

    string config_file_name = argv[1];
    config_t config;

    setConfig(config_file_name, config);

    core::init::init(1, rosetta_argv);

    ifstream fasta_in(config.protein_fasta_path);

    string sequence;
    string protein_desc;
    getline(fasta_in, protein_desc);
    getline(fasta_in, sequence);

    core::pose::Pose pose, fa_pose;

    core::pose::make_pose_from_sequence(
        pose, sequence, *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));

    core::pose::make_pose_from_sequence(
        fa_pose, sequence,
        *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));

    size_t pose_start = config.pose_start;
    residue_num = pose.total_residue();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx)
        residue_chi_num.push_back(fa_pose.residue(residue_idx).nchi());

    size_t res_index = 0;
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        residue_indexes.push_back(res_index);
        angle_types.push_back(angle_type::phi);
        angle_types.push_back(angle_type::psi);
        angle_types.push_back(angle_type::omega);

        pose.set_phi(residue_idx, 0);
        pose.set_psi(residue_idx, 0);
        pose.set_omega(residue_idx, 180);

        fa_pose.set_phi(residue_idx, 0);
        fa_pose.set_psi(residue_idx, 0);
        fa_pose.set_omega(residue_idx, 180);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            angle_types.push_back(angle_type::chi);

        res_index += 3 + residue_chi_num[residue_idx - 1];
    }

    core::pose::Pose native_structure, native_structure2;

    core::import_pose::pose_from_file(native_structure, config.protein_pdb_path, false, core::import_pose::PDB_file);

    core::import_pose::pose_from_file(native_structure2, config.protein_pdb_path, false, core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::dssp::Dssp dssp_native(native_structure);

    if (config.ss_w_coil) {
        ss_w_coil = config.ss_w_coil;
    }
    if (config.ss_w_helix) {
        ss_w_helix = config.ss_w_helix;
    }
    if (config.ss_w_sheet) {
        ss_w_sheet = config.ss_w_sheet;
    }

    getFragments(config.protein_frag3_path, 3, frag3);
    getFragments(config.protein_frag9_path, 9, frag9);
    getSS(config.protein_ss_path);
    getContactMap(config.protein_cm_path);

    evo_alg::FitnessFunction<double>::shared_ptr fit_centroid(RosettaCentroidEnergyMOFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_centroid =
        core::scoring::ScoreFunctionFactory::create_score_function("score3", "score4L");
    max_energy = energy_score_centroid->score(pose);

    evo_alg::FitnessFunction<double>::shared_ptr fit_fa(RosettaFaEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_fa =
        core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    max_energy_fa = energy_score_fa->score(fa_pose);

    cout << endl << "Protein: " << endl;
    cout << protein_desc << endl;
    cout << sequence << endl << endl;

    cout << "Sec. Struct: " << ss << endl << endl;

    energy_max = 50 * residue_num;

    cout.precision(9);
    cout << fixed;

    evo_alg::Population<evo_alg::Individual<double>> frag_pop, frag_archive, res_pop, res_archive;
    vector<tuple<vector<vector<long double>>, evo_alg::fitness::frontier_t>> frag_best_frontiers, res_best_frontiers;
    vector<double> frag_diversity, res_diversity;
    vector<size_t> best_individuals;
    evo_alg::fitness::frontier_t best_frontier;

    size_t frag3_chromosome_size = ceil(residue_indexes.size() / 3.0);
    size_t frag9_chromosome_size = ceil(residue_indexes.size() / 9.0);
    size_t residue_chromosome_size = residue_indexes.size();

    size_t iteration_num = config.iteration_num;
    size_t pop_size = config.population_size;

    evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>> brkga_config;

    if (config.frag_type == "frag3") {
        brkga_config = evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>>(
            iteration_num, pop_size, frag3_chromosome_size, fit_centroid, frag3Decoder);
    } else if (config.frag_type == "frag9") {
        brkga_config = evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>>(
            iteration_num, pop_size, frag9_chromosome_size, fit_centroid, frag9Decoder);
    }

    brkga_config.log_step = config.output_level == "none" ? 0 : 1;

    timer.startTimer("frag");
    tie(frag_pop, frag_archive, frag_best_frontiers, frag_diversity) =
        evo_alg::brkga::runMultiObjective<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);
    timer.stopTimer("frag");

    vector<vector<long double>> initial_pop;

    evo_alg::Population<evo_alg::Individual<double>> transcoded_pop;
    evo_alg::Population<evo_alg::Individual<double>> rough_pop;
    for (size_t ind_index = 0; ind_index < min(frag_archive.getSize(), frag_pop.getSize()); ++ind_index) {
        vector<double> chromosome = frag_archive[ind_index].getChromosome();
        vector<long double> coded_chromosome;
        for (size_t res_index = 0; res_index < residue_num; ++res_index) {
            size_t chromosome_res_idx = residue_indexes[res_index];
            uint64_t k = 1e9;
            uint64_t coded_phi = trunc((chromosome[chromosome_res_idx] + 180) / 360 * k);
            uint64_t coded_psi = trunc((chromosome[chromosome_res_idx + 1] + 180) / 360 * k);

            long double coded_residue = (coded_phi * k + coded_psi) / 1e18;

            coded_chromosome.push_back(coded_residue);
        }

        evo_alg::Individual<double> rough_ind(fit_centroid, chromosome);
        rough_pop.appendIndividual(rough_ind);

        vector<double> chromosome_2 = residueDecoder(coded_chromosome);

        for (size_t index = 0; index < chromosome.size(); ++index) {
            if (angle_types[index] != angle_type::omega &&
                !evo_alg::utils::numericEqual(chromosome[index], chromosome_2[index], 1e-4)) {
                double gene_1 = chromosome[index];
                double gene_2 = chromosome_2[index];
                if (gene_1 < 0 && gene_2 > 0) {
                    gene_1 += 360;
                }
                if (gene_1 > 0 && gene_2 < 0) {
                    gene_2 += 360;
                }
                if (!evo_alg::utils::numericEqual(gene_1, gene_2, 1e-4)) {
                    cout << "mapping error" << endl;
                    cout << index << " " << (int)angle_types[index] << endl;
                    cout << gene_1 << " " << gene_2 << endl;
                    exit(-1);
                }
            }
        }

        evo_alg::Individual<double> transcoded_ind(fit_centroid, chromosome_2);
        transcoded_pop.appendIndividual(transcoded_ind);

        initial_pop.push_back(coded_chromosome);
    }

    brkga_config.chromosome_size = residue_chromosome_size;
    brkga_config.exploration_diversity = 0.25;
    brkga_config.decoder = residueDecoder;
    brkga_config.initial_pop = initial_pop;

    timer.startTimer("res");
    tie(res_pop, res_archive, res_best_frontiers, res_diversity) =
        evo_alg::brkga::runMultiObjective<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);
    timer.stopTimer("res");

    outputProteinResults(config, frag_archive, "frag", sequence, residue_num, pose_start, fa_pose, pose,
                         native_structure, native_structure2, fit_fa, energy_score_fa, energy_score_centroid,
                         dssp_native);
    outputProteinResults(config, res_archive, "res", sequence, residue_num, pose_start, fa_pose, pose, native_structure,
                         native_structure2, fit_fa, energy_score_fa, energy_score_centroid, dssp_native);

    outputAlgorithmResults(config, frag_best_frontiers, frag_diversity, "frag");
    outputAlgorithmResults(config, res_best_frontiers, res_diversity, "res");

    rough_pop.evaluateFitness();
    transcoded_pop.evaluateFitness();

    best_individuals = rough_pop.getBestIndividuals();
    best_frontier.clear();
    for (size_t ind : best_individuals) {
        best_frontier.push_back(rough_pop[ind].getFitnessValue());
    }
    outputFitness(config, "rough", best_frontier);
    ofstream rough_chromosome_out(config.algorithm_output_dir + "rough.chromosome");
    rough_chromosome_out << fixed;
    rough_chromosome_out.precision(5);
    for (size_t index = 0; index < rough_pop.getSize(); ++index) {
        vector<double> chromosome = rough_pop[index].getChromosome();
        rough_chromosome_out << chromosome[0];
        for (size_t i = 1; i < chromosome.size(); ++i)
            rough_chromosome_out << " " << chromosome[i];
        rough_chromosome_out << endl;
    }

    best_individuals = transcoded_pop.getBestIndividuals();
    best_frontier.clear();
    for (size_t ind : best_individuals) {
        best_frontier.push_back(transcoded_pop[ind].getFitnessValue());
    }
    outputFitness(config, "transcoded", best_frontier);
    ofstream transcoded_chromosome_out(config.algorithm_output_dir + "transcoded.chromosome");
    transcoded_chromosome_out << fixed;
    transcoded_chromosome_out.precision(5);
    for (size_t index = 0; index < transcoded_pop.getSize(); ++index) {
        vector<double> chromosome = transcoded_pop[index].getChromosome();
        transcoded_chromosome_out << chromosome[0];
        for (size_t i = 1; i < chromosome.size(); ++i)
            transcoded_chromosome_out << " " << chromosome[i];
        transcoded_chromosome_out << endl;
    }

    best_individuals = frag_archive.getBestIndividuals();
    best_frontier.clear();
    for (size_t ind : best_individuals) {
        best_frontier.push_back(frag_archive[ind].getFitnessValue());
    }
    outputFitness(config, "first_phase", best_frontier);

    best_individuals = res_archive.getBestIndividuals();
    best_frontier.clear();
    for (size_t ind : best_individuals) {
        best_frontier.push_back(res_archive[ind].getFitnessValue());
    }
    outputFitness(config, "second_phase", best_frontier);

    vector<double> total_diversity;
    total_diversity.insert(total_diversity.end(), frag_diversity.begin(), frag_diversity.end());
    total_diversity.insert(total_diversity.end(), res_diversity.begin(), res_diversity.end());

    ofstream diversity_out(config.algorithm_output_dir + "total.diversity");
    diversity_out << fixed;
    diversity_out.precision(9);
    for (size_t index = 0; index < total_diversity.size(); ++index) {
        diversity_out << index << " " << total_diversity[index] << endl;
    }

    ofstream time_out(config.algorithm_output_dir + "total.time");
    time_out << (int)(timer.getTime("frag") + timer.getTime("res")) << endl;

    return EXIT_SUCCESS;
}
