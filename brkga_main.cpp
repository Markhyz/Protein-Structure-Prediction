#include "protein_commons.hpp"

#include <evo_alg/algorithms/brkga.hpp>

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
        if (rosetta_score < 0)
            rosetta_score = 0;

        double result = (rosetta_score + ss_score) / 2;

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

    static constexpr size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

vector<double> frag3Decoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t gene_index = 0; gene_index < encoded_chromosome.size(); ++gene_index) {
        size_t res_index = gene_index * 3;
        if (res_index >= frag3.size()) {
            res_index -= res_index - frag3.size() + 1;
        }
        size_t frag_index = floor(encoded_chromosome[gene_index] * frag3[res_index].size());
        vector<vector<double>> fragment = frag3[res_index][frag_index];
        for (size_t res_offset = 0; res_offset < fragment.size(); ++res_offset) {
            decoded_chromosome.push_back(fragment[res_offset][0]);
            decoded_chromosome.push_back(fragment[res_offset][1]);
            decoded_chromosome.push_back(fragment[res_offset][2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index + res_offset]; ++chi_idx)
                decoded_chromosome.push_back(0.0);
        }
    }

    return decoded_chromosome;
}

vector<double> frag9Decoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t gene_index = 0; gene_index < encoded_chromosome.size(); ++gene_index) {
        size_t res_index = gene_index * 9;
        if (res_index >= frag9.size()) {
            res_index -= res_index - frag9.size() + 1;
        }
        size_t frag_index = floor(encoded_chromosome[gene_index] * frag9[res_index].size());
        vector<vector<double>> fragment = frag9[res_index][frag_index];
        for (size_t res_offset = 0; res_offset < fragment.size(); ++res_offset) {
            decoded_chromosome.push_back(fragment[res_offset][0]);
            decoded_chromosome.push_back(fragment[res_offset][1]);
            decoded_chromosome.push_back(fragment[res_offset][2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index + res_offset]; ++chi_idx)
                decoded_chromosome.push_back(0.0);
        }
    }

    return decoded_chromosome;
}

vector<double> residueDecoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t gene_index = 0; gene_index < encoded_chromosome.size(); gene_index += 2) {
        decoded_chromosome.push_back(-180 + encoded_chromosome[gene_index] * 360);
        decoded_chromosome.push_back(-180 + encoded_chromosome[gene_index + 1] * 360);
        decoded_chromosome.push_back(180);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[gene_index / 2]; ++chi_idx)
            decoded_chromosome.push_back(0.0);
    }

    return decoded_chromosome;
}

vector<double> protein_angles;

vector<double> refineDecoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t res_index = 0; res_index < residue_indexes.size(); ++res_index) {
        double phi_delta = -10 + encoded_chromosome[res_index * 2] * 20;
        double psi_delta = -10 + encoded_chromosome[res_index * 2 + 1] * 20;
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index]] + phi_delta);
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index] + 1] + psi_delta);
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index] + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index]; ++chi_idx)
            decoded_chromosome.push_back(0.0);
    }

    return decoded_chromosome;
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

    core::pose::Pose native_structure;
    core::import_pose::pose_from_file(native_structure, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);

    core::pose::Pose native_structure2;
    core::import_pose::pose_from_file(native_structure2, "proteins/" + protein_name + "/pdb", false,
                                      core::import_pose::PDB_file);
    core::util::switch_to_residue_type_set(native_structure2, core::chemical::CENTROID);

    core::scoring::sasa::SasaCalc sasa_calc;

    core::scoring::dssp::Dssp dssp_native(native_structure);

    // getSS(protein_name);
    getTrueSS(dssp_native.get_dssp_secstruct());
    getFragments(protein_name, 3, frag3, frag3_prob);
    getFragments(protein_name, 9, frag9, frag9_prob);
    getContactMap(protein_name);

    evo_alg::FitnessFunction<double>::shared_ptr fit_centroid(RosettaCentroidEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_centroid =
        core::scoring::ScoreFunctionFactory::create_score_function("score3");
    max_energy = 200;

    evo_alg::FitnessFunction<double>::shared_ptr fit_fa(RosettaFaEnergyFunction::create(sequence));
    core::scoring::ScoreFunctionOP energy_score_fa =
        core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    max_energy_fa = energy_score_fa->score(faPose);

    cout << endl << "Protein: " << endl;
    cout << protein_desc << endl;
    cout << sequence << endl << endl;

    cout << "Sec. Struct: " << ss << endl << endl;

    evo_alg::Population<evo_alg::Individual<double>> pop;
    evo_alg::Individual<double> best_ind;
    vector<double> best_fit, mean_fit, diversity;

    size_t frag3_chromosome_size = ceil(residue_indexes.size() / 3.0);
    size_t frag9_chromosome_size = ceil(residue_indexes.size() / 9.0);
    size_t simple_chromosome_size = residue_indexes.size() * 2;

    size_t iteration_num = 1000;
    size_t pop_size = 250;
    double elite_x = 0.3;
    double mut_x = 0.3;
    double cr_x = 0.55;
    double diversity_threshold = 0.2;

    auto update_fn = [&](evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>>& config, size_t it) {
        if (it > iteration_num / 3) {
            double progress = (3.0 * it / (double)iteration_num - 1);
            double min_cross_pr = 1 - 3.0 / config.chromosome_size;
            config.elite_fraction =
                max((5.0 / config.pop_size), elite_x - (elite_x - (5.0 / config.pop_size)) * progress);
            config.diversity_threshold = max(0.0, diversity_threshold - (diversity_threshold - 0) * progress);
            config.elite_cross_pr = min(min_cross_pr, cr_x + (min_cross_pr - cr_x) * progress);
            config.mut_fraction = min(0.4, mut_x + (0.4 - mut_x) * progress);
        }
    };

    evo_alg::brkga::config_t<evo_alg::FitnessFunction<double>> brkga_config(
        iteration_num, pop_size, frag9_chromosome_size, elite_x, mut_x, cr_x, fit_centroid, frag9Decoder);

    brkga_config.diversity_threshold = diversity_threshold;
    brkga_config.log_step = 1;
    brkga_config.update_fn = update_fn;

    tie(best_ind, pop, best_fit, mean_fit, diversity) =
        evo_alg::brkga::run<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);

    protein_angles = best_ind.getChromosome();
    vector<double> chromosome_pre = best_ind.getChromosome();

    brkga_config.chromosome_size = simple_chromosome_size;
    brkga_config.decoder = refineDecoder;
    brkga_config.initial_pop = vector<vector<double>>(pop_size, vector<double>(simple_chromosome_size, 0.5));

    vector<double> best_fit_2, mean_fit_2, diversity_2;
    tie(best_ind, pop, best_fit_2, mean_fit_2, diversity_2) =
        evo_alg::brkga::run<evo_alg::Individual<double>, evo_alg::FitnessFunction<double>>(brkga_config);

    for (size_t index = 0; index < best_fit_2.size(); ++index) {
        best_fit.push_back(best_fit_2[index]);
        mean_fit.push_back(mean_fit_2[index]);
        diversity.push_back(diversity_2[index]);
    }

    vector<double> chromosome_post = best_ind.getChromosome();

    cout.precision(9);
    cout << fixed << endl;
    cout << "Pre-refine" << endl << endl;

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

    cout << "Best fitness: " << best_ind.getFitnessValue()[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << " "
         << energy_score_centroid->score(pose) << endl
         << endl;

    vector<double> scores_ = dynamic_pointer_cast<RosettaCentroidEnergyFunction>(fit_centroid)->getScores(pose);
    vector<double> scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);
    vector<double> scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);

    cout << "Energy: " << scores[0] << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    cout << "Energy: " << scores_[0] << " / " << max_energy << endl;
    cout << "SS: " << scores_[1] << " / " << scores_[2] << endl;
    cout << "CM: " << scores_[3] << " / " << scores_[4] << endl << endl;

    cout << "Energy: " << scores_native[0] << endl;
    cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
    cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    cout << endl
         << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    core::scoring::dssp::Dssp dssp0(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp0.get_dssp_secstruct() << endl << endl;

    if (pose_start - 1 > 0) {
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
    }

    cout << "Post-refine" << endl << endl;

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

    cout << "Best fitness: " << best_ind.getFitnessValue()[0] << endl;
    cout << "Centroid energy: " << energy_score_centroid->score(native_structure2) << " "
         << energy_score_centroid->score(pose) << endl
         << endl;

    scores_ = dynamic_pointer_cast<RosettaCentroidEnergyFunction>(fit_centroid)->getScores(pose);
    scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);
    scores_native = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(native_structure);

    cout << "Energy: " << scores[0] << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    cout << "Energy: " << scores_[0] << " / " << max_energy << endl;
    cout << "SS: " << scores_[1] << " / " << scores_[2] << endl;
    cout << "CM: " << scores_[3] << " / " << scores_[4] << endl << endl;

    cout << "Energy: " << scores_native[0] << endl;
    cout << "SS: " << scores_native[1] << " / " << scores_native[2] << endl;
    cout << "CM: " << scores_native[3] << " / " << scores_native[4] << endl << endl;

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    cout << endl
         << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    core::scoring::dssp::Dssp dssp1(faPose);

    cout << "Native: " << dssp_native.get_dssp_secstruct() << endl;
    cout << "Found:  " << dssp1.get_dssp_secstruct() << endl << endl;

    if (pose_start - 1 > 0) {
        core::pose::make_pose_from_sequence(
            faPose, sequence,
            *(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD)));
        for (size_t residue_idx = 1; residue_idx <= residue_num; ++residue_idx) {
            size_t chromosome_res_idx = residue_indexes[residue_idx - 1];
            faPose.set_phi(residue_idx, chromosome_post[chromosome_res_idx]);
            faPose.set_psi(residue_idx, chromosome_post[chromosome_res_idx + 1]);
            faPose.set_omega(residue_idx, chromosome_post[chromosome_res_idx + 2]);

            for (size_t chi_idx = 0; chi_idx < residue_chi_num[residue_idx - 1]; ++chi_idx)
                faPose.set_chi(chi_idx + 1, residue_idx, chromosome_post[chromosome_res_idx + 3 + chi_idx]);
        }
    }

    protocols::relax::FastRelax fast_relax(energy_score_fa);
    fast_relax.apply(faPose);

    scores = dynamic_pointer_cast<RosettaFaEnergyFunction>(fit_fa)->getScores(faPose);

    cout << endl
         << "Full Atom energy: " << energy_score_fa->score(native_structure) << " " << energy_score_fa->score(faPose)
         << endl
         << endl;

    cout << "Energy: " << scores[0] << endl;
    cout << "SS: " << scores[1] << " / " << scores[2] << endl;
    cout << "CM: " << scores[3] << " / " << scores[4] << endl << endl;

    if (pose_start - 1 > 0) {
        faPose.delete_residue_range_slow(1, pose_start - 1);
    }

    cout << endl
         << core::scoring::all_atom_rmsd(faPose, native_structure) << " "
         << core::scoring::CA_rmsd(faPose, native_structure) << endl
         << endl;

    cout << sasa_calc.calculate(native_structure) << endl;
    cout << sasa_calc.calculate(faPose) << endl;

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