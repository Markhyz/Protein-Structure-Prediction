#include "protein_brkga.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "./mobrkga_pred CONFIG_FILE" << endl;

        return EXIT_FAILURE;
    }

    string config_file_name = argv[1];
    config_t config;

    setConfig(config_file_name, config);

    ofstream ca_rmsd_out(config.score_output_dir + config.name + string(".ca_rmsd"));
    ofstream aa_rmsd_out(config.score_output_dir + config.name + string(".aa_rmsd"));
    ofstream gdttm_out(config.score_output_dir + config.name + string(".gdttm"));
    ofstream mean_rmsd_out(config.score_output_dir + config.name + string(".mean_rmsd"));
    ofstream mean_gdt_out(config.score_output_dir + config.name + string(".mean_gdt"));

    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    size_t pose_start = config.pose_start;

    ifstream fasta_in(config.protein_fasta_path);

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

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, config.protein_pdb_path, false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, config.protein_pdb_path, false,
                                      core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::sasa::SasaCalc sasa_calc;

    core::scoring::dssp::Dssp dssp_native(native_structure);
    getSS(config.protein_ss_path);
    getFragments(config.protein_frag3_path, 3, frag3, frag3_prob);
    getFragments(config.protein_frag9_path, 9, frag9, frag9_prob);
    getContactMap(config.protein_cm_path);

    evo_alg::FitnessFunction<double>::shared_ptr fit_centroid(RosettaCentroidEnergyMOFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_centroid =
        core::scoring::ScoreFunctionFactory::create_score_function("score3", "score4L");
    max_energy = energy_score_centroid->score(pose);

    evo_alg::FitnessFunction<double>::shared_ptr fit_fa(RosettaFaEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_fa =
        core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    max_energy_fa = energy_score_fa->score(faPose);

    cout << endl << "Protein: " << endl;
    cout << protein_desc << endl;
    cout << sequence << endl << endl;

    cout << "Sec. Struct: " << ss << endl << endl;

    evo_alg::Population<evo_alg::Individual<double>> pop, archive;
    vector<evo_alg::fitness::frontier_t> best_frontiers;
    vector<double> diversity;

    size_t frag3_chromosome_size = ceil(residue_indexes.size() / 3.0);
    size_t frag9_chromosome_size = ceil(residue_indexes.size() / 9.0);
    size_t residue_chromosome_size = residue_indexes.size();

    size_t iteration_num = config.iteration_num;
    size_t pop_size = config.population_size;
    double elite_x = config.elite_fraction;
    double mut_x = config.mutant_fraction;
    double cr_x = config.crossover_prob;
    double diversity_threshold = config.diversity_threshold;
    double diversity_enforcement = config.diversity_enforcement;

    auto update_fn = [&](evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>>& config, size_t it) {
        if (it > iteration_num / 2) {
            double progress = 2 * (2.0 * it / (double)iteration_num - 1);

            config.diversity_enforcement = max(0.0, diversity_enforcement * (1 - progress));
        }
    };

    evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>> brkga_config(
        iteration_num, pop_size, frag9_chromosome_size, elite_x, mut_x, cr_x, fit_centroid, frag9Decoder);

    brkga_config.log_step = 1;
    brkga_config.diversity_threshold = diversity_threshold;
    brkga_config.diversity_enforcement = diversity_enforcement;
    brkga_config.update_fn = update_fn;

    tie(pop, archive, best_frontiers, diversity) =
        evo_alg::brkga::runMultiObjective<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);

    vector<size_t> best_frontier_dirty = archive.getBestIndividuals();
    vector<size_t> best_frontier;
    for (size_t index_1 : best_frontier_dirty) {
        bool unique = true;
        for (size_t index_2 : best_frontier) {
            evo_alg::fitness::FitnessValue fitness_1 = archive[index_1].getFitnessValue();
            evo_alg::fitness::FitnessValue fitness_2 = archive[index_2].getFitnessValue();

            for (size_t fit_index = 0; fit_index < fitness_1.getDimension(); ++fit_index) {
                if (fabs(fitness_1[fit_index] - fitness_2[fit_index]) < 1e-9) {
                    unique = false;
                    break;
                }
            }
            if (!unique) {
                break;
            }
        }
        if (unique) {
            best_frontier.push_back(index_1);
        }
    }

    cout.precision(9);
    cout << fixed << endl;
    cout << "Pre-refine" << endl << endl;
    cout << "Total: " << best_frontier_dirty.size() << endl;
    cout << "True: " << best_frontier.size() << endl << endl;

    double best_gdt = 0;
    size_t best_structure = -1;
    for (size_t ind_index = 0; ind_index < archive.getSize(); ++ind_index) {
        vector<double> chromosome = archive[ind_index].getChromosome();
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
        double gdt;

        if (pose_start - 1 > 0) {
            faPose.delete_residue_range_slow(1, pose_start - 1);
        }

        if (faPose.total_residue() > native_structure.total_residue()) {
            faPose.delete_residue_range_slow(native_structure.total_residue() + 1, faPose.total_residue());
        }

        gdt = core::scoring::CA_gdtmm(faPose, native_structure);

        if (gdt > best_gdt) {
            best_gdt = gdt;
            best_structure = ind_index;
        }

        if (faPose.total_residue() < residue_num) {
            core::pose::make_pose_from_sequence(
                faPose, sequence,
                *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        }
    }

    vector<double> chromosome_pre = archive[best_structure].getChromosome();

    vector<vector<double>> initial_pop;

    for (size_t ind_index = 0; ind_index < min(archive.getSize(), pop.getSize()); ++ind_index) {
        vector<double> chromosome = archive[ind_index].getChromosome();
        vector<double> coded_chromosome;
        for (size_t res_index = 0; res_index < residue_num; ++res_index) {
            size_t chromosome_res_idx = residue_indexes[res_index];
            double coded_phi = trunc((chromosome[chromosome_res_idx] + 180) / 360 * 1e7);
            double coded_psi = trunc((chromosome[chromosome_res_idx + 1] + 180) / 360 * 1e7);

            double coded_residue = (coded_phi * 1e7 + coded_psi) / 1e14;

            coded_chromosome.push_back(coded_residue);
        }

        vector<double> chromosome_2 = residueDecoder(coded_chromosome);
        for (size_t index = 0; index < chromosome.size(); ++index) {
            if (angle_types[index] != angle_type::omega &&
                !evo_alg::utils::numericEqual(chromosome[index], chromosome_2[index], 1e-2)) {
                double gene_1 = chromosome[index];
                double gene_2 = chromosome_2[index];
                if (gene_1 < 0 && gene_2 > 0) {
                    gene_1 += 360;
                }
                if (gene_1 > 0 && gene_2 < 0) {
                    gene_2 += 360;
                }
                if (!evo_alg::utils::numericEqual(gene_1, gene_2, 1e-2)) {
                    cout << "mapping error" << endl;
                    cout << index << " " << (int)angle_types[index] << endl;
                    cout << gene_1 << " " << gene_2 << endl;
                    exit(-1);
                }
            }
        }

        initial_pop.push_back(coded_chromosome);
    }

    brkga_config.chromosome_size = residue_chromosome_size;
    brkga_config.decoder = residueDecoder;
    brkga_config.initial_pop = initial_pop;

    vector<double> diversity_2;
    tie(pop, archive, best_frontiers, diversity_2) =
        evo_alg::brkga::runMultiObjective<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);

    for (size_t index = 0; index < diversity_2.size(); ++index) {
        diversity.push_back(diversity_2[index]);
    }

    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        pose.set_phi(residue_idx, chromosome_pre[chromosome_res_idx]);
        pose.set_psi(residue_idx, chromosome_pre[chromosome_res_idx + 1]);
        pose.set_omega(residue_idx, chromosome_pre[chromosome_res_idx + 2]);

        faPose.set_phi(residue_idx, chromosome_pre[chromosome_res_idx]);
        faPose.set_psi(residue_idx, chromosome_pre[chromosome_res_idx + 1]);
        faPose.set_omega(residue_idx, chromosome_pre[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            faPose.set_chi(chi_idx + 1, residue_idx, chromosome_pre[chromosome_res_idx + 3 + chi_idx]);
    }

    vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);

    cout << "FA Energy: " << scores[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(pose) << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    vector<double> scores_native;
    if (native_structure.total_residue() == residue_num) {
        scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);
        cout << "FA Energy: " << scores_native[0] << endl;
        cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << endl;
        cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
        cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;
    }

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    if (faPose.total_residue() > native_structure.total_residue()) {
        faPose.delete_residue_range_slow(native_structure.total_residue() + 1, faPose.total_residue());
    }

    cout << "RMSD:" << endl;
    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TM:" << endl;
    cout << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    core::scoring::dssp::Dssp dssp1(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "SS:     " << ss.substr(pose_start - 1) << endl;
    cout << "Found:  " << dssp1.get_dssp_secstruct() << endl << endl;

    if (faPose.total_residue() < residue_num) {
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
    }

    best_frontier.clear();
    best_frontier_dirty = archive.getBestIndividuals();
    for (size_t index_1 : best_frontier_dirty) {
        bool unique = true;
        for (size_t index_2 : best_frontier) {
            evo_alg::fitness::FitnessValue fitness_1 = archive[index_1].getFitnessValue();
            evo_alg::fitness::FitnessValue fitness_2 = archive[index_2].getFitnessValue();

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

    cout << "Post-refine" << endl << endl;
    cout << "Total: " << best_frontier_dirty.size() << endl;
    cout << "True: " << best_frontier.size() << endl << endl;

    ofstream frontier_info("frontier.info");
    frontier_info << fixed;
    frontier_info.precision(9);

    best_gdt = 0;
    best_structure = -1;
    double best_rmsd = 1e9;
    int best_structure_2 = -1;
    double mean_rmsd = 0, mean_gdt = 0;
    for (size_t ind_index = 0; ind_index < archive.getSize(); ++ind_index) {
        vector<double> chromosome = archive[ind_index].getChromosome();
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

        if (pose_start - 1 > 0) {
            faPose.delete_residue_range_slow(1, pose_start - 1);
        }

        if (faPose.total_residue() > native_structure.total_residue()) {
            faPose.delete_residue_range_slow(native_structure.total_residue() + 1, faPose.total_residue());
        }

        double gdt, rmsd;

        frontier_info << "Index: " << ind_index << endl << endl;

        gdt = core::scoring::CA_gdtmm(faPose, native_structure);
        rmsd = core::scoring::CA_rmsd(faPose, native_structure);

        frontier_info << "GDT-TS: " << gdt << endl;
        frontier_info << "CA RMSD: " << rmsd << endl;
        frontier_info << "Full RMSD: " << core::scoring::all_atom_rmsd(faPose, native_structure) << endl;

        frontier_info << "Centroid Energy: " << energy_score_centroid->score(pose) << endl;

        if (native_structure.total_residue() == residue_num) {
            vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);
            frontier_info << "Full Atoms Energy: " << scores[0] << endl;
            frontier_info << "SS: " << scores[1] << " / " << scores[2] << endl;
            frontier_info << "CM: " << scores[3] << " / " << scores[4] << endl << endl;
        }

        faPose.dump_pdb(config.decoy_output_dir + string("decoy_") + to_string(ind_index + 1) + string(".pdb"));

        mean_gdt += gdt;
        mean_rmsd += rmsd;

        if (gdt > best_gdt) {
            best_gdt = gdt;
            best_structure = ind_index;
        }
        if (rmsd < best_rmsd) {
            best_rmsd = rmsd;
            best_structure_2 = ind_index;
        }

        if (faPose.total_residue() < residue_num) {
            core::pose::make_pose_from_sequence(
                faPose, sequence,
                *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        }
    }

    mean_gdt /= archive.getSize();
    mean_rmsd /= archive.getSize();

    vector<double> chromosome_post = archive[best_structure].getChromosome();

    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        pose.set_phi(residue_idx, chromosome_post[chromosome_res_idx]);
        pose.set_psi(residue_idx, chromosome_post[chromosome_res_idx + 1]);
        pose.set_omega(residue_idx, chromosome_post[chromosome_res_idx + 2]);

        faPose.set_phi(residue_idx, chromosome_post[chromosome_res_idx]);
        faPose.set_psi(residue_idx, chromosome_post[chromosome_res_idx + 1]);
        faPose.set_omega(residue_idx, chromosome_post[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            faPose.set_chi(chi_idx + 1, residue_idx, chromosome_post[chromosome_res_idx + 3 + chi_idx]);
    }

    scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);

    cout << "FA Energy: " << scores[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(pose) << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    if (native_structure.total_residue() == residue_num) {
        scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);
        cout << "FA Energy: " << scores_native[0] << endl;
        cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << endl;
        cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
        cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;
    }

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    if (faPose.total_residue() > native_structure.total_residue()) {
        faPose.delete_residue_range_slow(native_structure.total_residue() + 1, faPose.total_residue());
    }

    cout << "RMSD:" << endl;
    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TS:" << endl;
    cout << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    core::scoring::dssp::Dssp dssp2(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "SS:     " << ss.substr(pose_start - 1) << endl;
    cout << "Found:  " << dssp2.get_dssp_secstruct() << endl << endl;

    if (pose_start - 1 > 0) {
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            pose.set_phi(residue_idx, chromosome_post[chromosome_res_idx]);
            pose.set_psi(residue_idx, chromosome_post[chromosome_res_idx + 1]);
            pose.set_omega(residue_idx, chromosome_post[chromosome_res_idx + 2]);

            faPose.set_phi(residue_idx, chromosome_post[chromosome_res_idx]);
            faPose.set_psi(residue_idx, chromosome_post[chromosome_res_idx + 1]);
            faPose.set_omega(residue_idx, chromosome_post[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                faPose.set_chi(chi_idx + 1, residue_idx, chromosome_post[chromosome_res_idx + 3 + chi_idx]);
        }
    }

    repackProtein(faPose, energy_score_fa);

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    cout << "Full Atom energy: " << energy_score_fa->score(native_structure) << " " << energy_score_fa->score(faPose)
         << endl
         << endl;

    cout << "RMSD: " << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TM: " << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    cout << "Mean RMSD: " << mean_rmsd << endl;
    cout << "Mean GDT-TM: " << mean_gdt << endl << endl;

    cout << sasa_calc.calculate(native_structure) << endl;
    cout << sasa_calc.calculate(faPose) << endl;

    gdttm_out << fixed;
    gdttm_out.precision(9);
    gdttm_out << core::scoring::CA_gdtmm(faPose, native_structure) << endl;

    if (faPose.total_residue() < residue_num) {
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
    }

    vector<double> chromosome_rmsd = archive[best_structure_2].getChromosome();
    for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
        size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
        pose.set_phi(residue_idx, chromosome_rmsd[chromosome_res_idx]);
        pose.set_psi(residue_idx, chromosome_rmsd[chromosome_res_idx + 1]);
        pose.set_omega(residue_idx, chromosome_rmsd[chromosome_res_idx + 2]);

        faPose.set_phi(residue_idx, chromosome_rmsd[chromosome_res_idx]);
        faPose.set_psi(residue_idx, chromosome_rmsd[chromosome_res_idx + 1]);
        faPose.set_omega(residue_idx, chromosome_rmsd[chromosome_res_idx + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
            faPose.set_chi(chi_idx + 1, residue_idx, chromosome_rmsd[chromosome_res_idx + 3 + chi_idx]);
    }

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    if (faPose.total_residue() > native_structure.total_residue()) {
        faPose.delete_residue_range_slow(native_structure.total_residue() + 1, faPose.total_residue());
    }

    ca_rmsd_out << fixed;
    ca_rmsd_out.precision(9);
    ca_rmsd_out << core::scoring::CA_rmsd(faPose, native_structure) << endl;

    aa_rmsd_out << fixed;
    aa_rmsd_out.precision(9);
    aa_rmsd_out << core::scoring::all_atom_rmsd(faPose, native_structure) << endl;

    mean_rmsd_out << fixed;
    mean_rmsd_out.precision(9);
    mean_rmsd_out << mean_rmsd << endl;

    mean_gdt_out << fixed;
    mean_gdt_out.precision(9);
    mean_gdt_out << mean_gdt << endl;

    // for (size_t index = 0; index < 5; ++index) {
    //     evo_alg::fitness::frontier_t cur_frontier =
    //         best_frontiers[index ? (size_t)((index / 4.0) * best_frontiers.size()) - 1 : index];
    //     ofstream fit_out(string("res_") + to_string(index) + string(".fit"));
    //     fit_out << fixed;
    //     fit_out.precision(9);
    //     for (auto fit : cur_frontier) {
    //         for (double fit_value : fit.getValues())
    //             fit_out << fit_value << " ";
    //         fit_out << endl;
    //     }
    // }

    // ofstream diver_out("res.diver");
    // diver_out << fixed;
    // diver_out.precision(9);
    // for (size_t index = 0; index < diversity.size(); ++index) {
    //     diver_out << index << " " << diversity[index] << endl;
    // }

    return EXIT_SUCCESS;
}
