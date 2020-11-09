#include "brkga_commons.hpp"

class RosettaCentroidEnergyMOFunction : public evo_alg::FitnessFunction<double> {
  public:
    static RosettaCentroidEnergyMOFunction* create(string sequence) {
        core::pose::Pose protein_structure;
        core::pose::make_pose_from_sequence(
            protein_structure, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID)));
        size_t num_angles = 0;

        for (size_t residue_idx = 1; residue_idx <= protein_structure.total_residue(); ++residue_idx) {
            num_angles += 3 + residue_chi_num[residue_idx - 1];
        }

        return new RosettaCentroidEnergyMOFunction(num_angles, protein_structure);
    }

    size_t getDimension() const override {
        return dimension_;
    }

    RosettaCentroidEnergyMOFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    evo_alg::fitness::FitnessValue operator()(evo_alg::Genotype<double> const& genotype) override {
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

        double cm_score, cm_total;
        tie(cm_score, cm_total) = cmScore(protein_structure_);

        double rosetta_score = energy_score_->score(protein_structure_);

        return {-rosetta_score, ss_score, cm_score};
    }

    vector<double> getScores(core::pose::Pose& pose) {
        double ss_score, ss_total;
        tie(ss_score, ss_total) = ssScore(pose);

        double cm_score, cm_total;
        tie(cm_score, cm_total) = cmScore(pose);

        double rosetta_score = energy_score_->score(pose);

        return {rosetta_score, ss_score, ss_total, cm_score, cm_total};
    }

  private:
    RosettaCentroidEnergyMOFunction(size_t num_angles, core::pose::Pose& protein_structure)
        : FitnessFunction<double>({num_angles, {-180, 180}}), protein_structure_{protein_structure},
          energy_score_{core::scoring::ScoreFunctionFactory::create_score_function("score3", "score4L")} {};

    static constexpr size_t dimension_ = 3;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

// argv -> 1: protein name | 2: result score directory | 3: result decoy directory
int main(int argc, char** argv) {
    ofstream ca_rmsd_out(argc > 2 ? string(argv[2]) + string("ca_rmsd") : "ca_rmsd", ios::app);
    ofstream aa_rmsd_out(argc > 2 ? string(argv[2]) + string("aa_rmsd") : "aa_rmsd", ios::app);
    ofstream gdttm_out(argc > 2 ? string(argv[2]) + string("gdttm") : "gdttm", ios::app);

    char* argv2[] = {argv[0]};
    core::init::init(1, argv2);
    string protein_name = argv[1];

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

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::sasa::SasaCalc sasa_calc;

    core::scoring::dssp::Dssp dssp_native(native_structure);

    getSS(protein_name);
    getFragments(protein_name, 3, frag3, frag3_prob);
    getFragments(protein_name, 9, frag9, frag9_prob);
    getContactMap(protein_name);

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

    size_t iteration_num = 500;
    size_t pop_size = 100;
    double elite_x = 0.5;
    double mut_x = 0.2;
    double cr_x = 0.55;
    double diversity_threshold = 0.2;
    double diversity_enforcement = 0.5;

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

        gdt = core::scoring::CA_gdtmm(faPose, native_structure);

        if (gdt > best_gdt) {
            best_gdt = gdt;
            best_structure = ind_index;
        }
    }

    vector<double> chromosome_pre = archive[best_structure].getChromosome();

    vector<vector<double>> initial_pop;

    for (size_t ind_index = 0; ind_index < min(archive.getSize(), pop.getSize()); ++ind_index) {
        vector<double> chromosome = archive[ind_index].getChromosome();
        vector<double> coded_chromosome;
        for (size_t res_index = 0; res_index < residue_num; ++res_index) {
            size_t chromosome_res_idx = residue_indexes[res_index];
            double coded_phi = trunc((chromosome[chromosome_res_idx] + 180) / 360 * 1e5);
            double coded_psi = trunc((chromosome[chromosome_res_idx + 1] + 180) / 360 * 1e5);
            double coded_omega = trunc((chromosome[chromosome_res_idx + 2] + 180) / 360 * 1e5);

            double wtf = (coded_phi * 1e10 + coded_psi * 1e5 + coded_omega) / 1e15;

            coded_chromosome.push_back(wtf);
        }

        vector<double> chromosome_2 = residueDecoder(coded_chromosome);
        for (size_t index = 0; index < chromosome.size(); ++index) {
            if (!evo_alg::utils::numericEqual(chromosome[index], chromosome_2[index], 1e-2)) {
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
    vector<double> scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);

    cout << "FA Energy: " << scores[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(pose) << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    cout << "FA Energy: " << scores_native[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << endl;
    cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
    cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;

    cout << "RMSD:" << endl;
    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TM:" << endl;
    cout << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    core::scoring::dssp::Dssp dssp1(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "SS:     " << ss << endl;
    cout << "Found:  " << dssp1.get_dssp_secstruct() << endl << endl;

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

        frontier_info << "Index: " << ind_index << endl << endl;

        gdt = core::scoring::CA_gdtmm(faPose, native_structure);

        frontier_info << "GDT-TM: " << gdt << endl;
        frontier_info << "CA RMSD: " << core::scoring::CA_rmsd(faPose, native_structure) << endl;
        frontier_info << "Full RMSD: " << core::scoring::all_atom_rmsd(faPose, native_structure) << endl;

        frontier_info << "Centroid Energy: " << energy_score_centroid->score(pose) << endl;
        vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);
        frontier_info << "Full Atoms Energy: " << scores[0] << endl;
        frontier_info << "SS: " << scores[1] << " / " << scores[2] << endl;
        frontier_info << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

        faPose.dump_pdb(string(argc > 3 ? argv[3] : "") + string("decoy_") + to_string(ind_index + 1) + string(".pdb"));

        if (gdt > best_gdt) {
            best_gdt = gdt;
            best_structure = ind_index;
        }
    }

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
    scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);

    cout << "FA Energy: " << scores[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(pose) << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    cout << "FA Energy: " << scores_native[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << endl;
    cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
    cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;

    cout << "RMSD:" << endl;
    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TM:" << endl;
    cout << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    core::scoring::dssp::Dssp dssp2(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "SS:     " << ss << endl;
    cout << "Found:  " << dssp2.get_dssp_secstruct() << endl << endl;

    repackProtein(faPose, energy_score_fa);

    scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);

    cout << "Full Atom energy: " << energy_score_fa->score(native_structure) << " " << energy_score_fa->score(faPose)
         << endl
         << endl;

    cout << "RMSD:" << endl;
    cout << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << "GDT-TM:" << endl;
    cout << core::scoring::CA_gdtmm(faPose, native_structure) << endl << endl;

    cout << sasa_calc.calculate(native_structure) << endl;
    cout << sasa_calc.calculate(faPose) << endl;

    faPose.dump_pdb("res.pdb");

    for (size_t index = 0; index < 5; ++index) {
        evo_alg::fitness::frontier_t cur_frontier =
            best_frontiers[index ? (size_t)((index / 4.0) * best_frontiers.size()) - 1 : index];
        ofstream fit_out(string("res_") + to_string(index) + string(".fit"));
        fit_out << fixed;
        fit_out.precision(9);
        for (auto fit : cur_frontier) {
            for (double fit_value : fit.getValues())
                fit_out << fit_value << " ";
            fit_out << endl;
        }
    }

    ofstream diver_out("res.diver");
    diver_out << fixed;
    diver_out.precision(9);
    for (size_t index = 0; index < diversity.size(); ++index) {
        diver_out << index << " " << diversity[index] << endl;
    }

    gdttm_out << fixed;
    gdttm_out.precision(9);
    gdttm_out << core::scoring::CA_gdtmm(faPose, native_structure) << endl;

    ca_rmsd_out << fixed;
    ca_rmsd_out.precision(9);
    ca_rmsd_out << core::scoring::CA_rmsd(faPose, native_structure) << endl;

    aa_rmsd_out << fixed;
    aa_rmsd_out.precision(9);
    aa_rmsd_out << core::scoring::all_atom_rmsd(faPose, native_structure) << endl;

    return 0;
}