#ifndef GUARD_H_PROTEIN_COMMONS
#define GUARD_H_PROTEIN_COMMONS

#include "rosetta.hpp"

#include <evo_alg/commons/utils.hpp>
#include <evo_alg/core.hpp>
#include <evo_alg/operators.hpp>

#include "utils.hpp"

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

using frag_list_t = vector<vector<vector<vector<double>>>>;

frag_list_t frag3, frag9;

vector<tuple<size_t, size_t, double>> cm;

size_t residue_num = 0;

double max_energy = 0;
double max_energy_fa = 0;

double ss_total = 0;
double cm_total = 0;
double energy_min = -100, energy_max = 10000;

double ss_w_coil = 1, ss_w_helix = 1, ss_w_sheet = 1;

double fix_angle(double angle) {
    return angle + (angle < -180 ? 360 : (angle > 180 ? -360 : 0));
}

void getFragments(string frag_path, size_t frag_size, frag_list_t& frag_list) {
    ifstream frag_in(frag_path);

    string _, line;
    size_t pos, n;
    size_t res_index = 0;
    while (getline(frag_in, line)) {
        istringstream iss(line);
        iss >> _ >> pos >> _ >> n;
        getline(frag_in, line);
        vector<vector<vector<double>>> frags;
        for (size_t index = 0; index < n; ++index) {
            vector<vector<double>> frag;
            for (size_t index2 = 0; index2 < frag_size; ++index2) {
                getline(frag_in, line);
                istringstream iss2(line);
                string res_ss;
                double phi, psi, omega;

                iss2 >> _ >> _ >> _ >> _;
                iss2 >> res_ss >> phi >> psi >> omega;

                frag.push_back({fix_angle(phi), fix_angle(psi), 180});
            }
            getline(frag_in, line);
            frags.push_back(frag);
        }
        frag_list.push_back(frags);

        ++res_index;
    }
}

vector<vector<double>> getPsipredSS(string ss_path) {
    ifstream ss_in(ss_path + "/ss_psipred");
    string line;

    string _;
    vector<vector<double>> ss_prob;
    double c_prob, h_prob, e_prob;
    while (ss_in >> _ >> _ >> _ >> c_prob >> h_prob >> e_prob) {
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
    return ss_prob;
}

vector<vector<double>> getRaptorxSS(string ss_path) {
    ifstream ss_in(ss_path + "/ss_raptorx");
    string line;

    string _;
    vector<vector<double>> ss_prob;
    double c_prob, h_prob, e_prob;
    while (ss_in >> _ >> _ >> _ >> h_prob >> e_prob >> c_prob) {
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
    return ss_prob;
}

vector<vector<double>> getSpider3SS(string ss_path) {
    ifstream ss_in(ss_path + "/ss_spider3");
    string line;

    string _;
    vector<vector<double>> ss_prob;
    double c_prob, h_prob, e_prob;
    while (ss_in >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> _ >> c_prob >> h_prob >> e_prob) {
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
    return ss_prob;
}

void getSS(string ss_path) {
    vector<vector<double>> psipred_prob = getPsipredSS(ss_path);
    vector<vector<double>> raptorx_prob = getRaptorxSS(ss_path);
    vector<vector<double>> spider3_prob = getSpider3SS(ss_path);

    for (size_t i = 0; i < psipred_prob.size(); ++i) {
        double c_prob, h_prob, e_prob;
        c_prob = psipred_prob[i][0];
        h_prob = psipred_prob[i][1];
        e_prob = psipred_prob[i][2];

        if (c_prob > h_prob && c_prob > e_prob) {
            ss += "L";
            ss_total += c_prob * ss_w_coil;
        } else if (h_prob > e_prob) {
            ss += "H";
            ss_total += h_prob * ss_w_helix;
        } else {
            ss += "E";
            ss_total += e_prob * ss_w_sheet;
        }

        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
}

void getTrueSS(string true_ss) {
    for (char c : true_ss) {
        double c_prob = 0, h_prob = 0, e_prob = 0;
        if (c == 'L') {
            ss_total += ss_w_coil;
            c_prob = 1;
        } else if (c == 'H') {
            ss_total += ss_w_helix;
            h_prob = 1;
        } else if (c == 'E') {
            ss_total += ss_w_sheet;
            e_prob = 1;
        }
        ss += c;
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
}

pair<double, double> ssScore(core::pose::Pose const& pose) {
    core::scoring::dssp::Dssp dssp(pose);

    string protein_ss = dssp.get_dssp_secstruct();

    double ss_score = 0, ss_total = 0;
    for (size_t ss_index = 0; ss_index < protein_ss.size(); ++ss_index) {
        double k;
        if (ss[ss_index] == 'L') {
            k = ss_w_coil;
            ss_total += ss_prob[ss_index][0] * k;
        }
        if (ss[ss_index] == 'H') {
            k = ss_w_helix;
            ss_total += ss_prob[ss_index][1] * k;
        }
        if (ss[ss_index] == 'E') {
            k = ss_w_sheet;
            ss_total += ss_prob[ss_index][2] * k;
        }

        if (protein_ss[ss_index] == 'L')
            ss_score += ss_prob[ss_index][0] * k;
        if (protein_ss[ss_index] == 'H')
            ss_score += ss_prob[ss_index][1] * k;
        if (protein_ss[ss_index] == 'E')
            ss_score += ss_prob[ss_index][2] * k;
    }

    return {ss_score, ss_total};
}

vector<vector<double>> getRaptorxCM(string cm_path) {
    ifstream cm_in(cm_path + "/con_raptorx");

    vector<vector<double>> contact_map(residue_num + 1, vector<double>(residue_num + 1));

    string line;
    while (getline(cm_in, line)) {
        istringstream iss(line);
        size_t i, j;
        double prob;
        iss >> i >> j >> prob;
        contact_map[i][j] = contact_map[j][i] = prob;
    }

    return contact_map;
}

vector<vector<double>> getMetaCM(string cm_path) {
    ifstream cm_in(cm_path + "/con_meta");

    vector<vector<double>> contact_map(residue_num + 1, vector<double>(residue_num + 1));

    string line, _;
    while (getline(cm_in, line)) {
        istringstream iss(line);
        size_t i, j;
        double prob;
        iss >> i >> j >> _ >> _ >> prob;
        contact_map[i][j] = contact_map[j][i] = prob;
    }

    return contact_map;
}

vector<vector<double>> getSpotCM(string cm_path) {
    ifstream cm_in(cm_path + "/con_spot");

    vector<vector<double>> contact_map(residue_num + 1, vector<double>(residue_num + 1));

    string line, _;
    while (getline(cm_in, line)) {
        istringstream iss(line);
        size_t i, j;
        double prob;
        iss >> i >> j >> _ >> _ >> prob;
        contact_map[i][j] = contact_map[j][i] = prob;
    }

    return contact_map;
}

vector<vector<double>> getGremlinCM(string cm_path) {
    ifstream cm_in(cm_path + "/con_gremlin");

    vector<vector<double>> contact_map(residue_num + 1, vector<double>(residue_num + 1));

    string line, _;
    while (getline(cm_in, line)) {
        istringstream iss(line);
        size_t i, j;
        double prob;
        iss >> i >> j >> _ >> _ >> _ >> _ >> prob;
        contact_map[i][j] = contact_map[j][i] = prob;
    }

    return contact_map;
}

void getContactMap(string cm_path) {
    ifstream cm_in(cm_path);

    vector<vector<double>> raptorx_cm = getRaptorxCM(cm_path);
    vector<vector<double>> meta_cm = getMetaCM(cm_path);
    vector<vector<double>> spot_cm = getSpotCM(cm_path);

    vector<vector<double>> contact_map(residue_num + 1, vector<double>(residue_num + 1));
    vector<tuple<double, int, int>> contact_pairs;
    for (size_t i = 1; i <= residue_num; ++i) {
        for (size_t j = i + 1; j <= residue_num; ++j) {
            contact_map[i][j] = raptorx_cm[i][j];
            contact_pairs.push_back({contact_map[i][j], i, j});
        }
    }

    sort(contact_pairs.rbegin(), contact_pairs.rend());

    for (size_t index = 0; index < contact_pairs.size(); ++index) {
        double prob;
        int i, j;
        tie(prob, i, j) = contact_pairs[index];
        if (cm.size() > residue_num)
            break;
        cm.push_back({i, j, prob});
        cm_total += prob;
    }
}

void getTrueCM(core::pose::Pose const& pose) {
    const size_t residue_num = pose.total_residue();
    vector<tuple<size_t, size_t, double>> contacts;
    for (size_t res_index_1 = 1; res_index_1 <= residue_num; ++res_index_1) {
        for (size_t res_index_2 = res_index_1 + 1; res_index_2 <= residue_num; ++res_index_2) {
            string atom = "CB";

            if (pose.residue(res_index_1).name1() == 'G' || pose.residue(res_index_2).name1() == 'G') {
                atom = "CA";
            }

            double dist = pose.residue(res_index_1).xyz(atom).distance(pose.residue(res_index_2).xyz(atom));
            if (evo_alg::utils::numericLowerEqual(dist, 8.0)) {
                contacts.push_back({res_index_1, res_index_2, 1.0});
            }
        }
    }
    shuffle(contacts.begin(), contacts.end(), evo_alg::utils::rng);
    for (size_t index = 0; index < contacts.size(); ++index) {
        cm.push_back(contacts[index]);
        cm_total += 1;
    }
}

pair<double, double> cmScore(core::pose::Pose const& pose) {
    double cm_score = 0, cm_total = 0;
    for (size_t cm_index = 0; cm_index < cm.size(); ++cm_index) {
        size_t i, j;
        double prob;
        tie(i, j, prob) = cm[cm_index];

        string atom = "CB";

        if (pose.residue(i).name1() == 'G' || pose.residue(j).name1() == 'G') {
            atom = "CA";
        }

        double dist = pose.residue(i).xyz(atom).distance(pose.residue(j).xyz(atom));
        if (evo_alg::utils::numericLowerEqual(dist, 8.0)) {
            cm_score += prob;
        } else {
            cm_score += prob / exp(dist - 8);
        }

        cm_total += prob;
    }

    return {cm_score, cm_total};
}

void setEnergyFunctionLimits(vector<tuple<vector<vector<long double>>, evo_alg::fitness::frontier_t>>& best_frontiers) {
    for (auto& best_frontier : best_frontiers) {
        evo_alg::fitness::frontier_t fitness = get<1>(best_frontier);
        for (auto& fit_values : fitness) {
            energy_min = min(energy_min, -fit_values[0]);
            energy_max = max(energy_max, -fit_values[0]);
        }
    }
}

void repackProtein(core::pose::Pose& pose, core::scoring::ScoreFunctionOP& scorefnx) {
    core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(pose);
    task->restrict_to_repacking();
    core::pack::pack_rotamers(pose, *scorefnx, task);
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

    vector<int8_t> getDirection() const override {
        return {1};
    }

    RosettaFaEnergyFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    evo_alg::fitness::FitnessValue normalize(evo_alg::fitness::FitnessValue const& fitness_value) override {
        return {0};
    }

    evo_alg::fitness::FitnessValue operator()(evo_alg::Genotype<double> const& genotype) override {
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

        double ss_score, ss_total;
        tie(ss_score, ss_total) = ssScore(protein_structure_);
        ss_score = ss_score / ss_total;

        double cm_score, cm_total;
        tie(cm_score, cm_total) = cmScore(protein_structure_);
        cm_score = cm_score / cm_total;

        double rosetta_score = 1 - energy_score_->score(protein_structure_) / max_energy_fa;

        double result = (rosetta_score + ss_score) / 2;

        return {result};
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
    RosettaFaEnergyFunction(size_t num_angles, core::pose::Pose& protein_structure)
        : FitnessFunction<double>({num_angles, {0, 1}}), protein_structure_{protein_structure},
          energy_score_{core::scoring::ScoreFunctionFactory::create_score_function("ref2015")} {};

    static constexpr size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

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

    vector<int8_t> getDirection() const override {
        return {1, 1, -1};
    }

    RosettaCentroidEnergyMOFunction* clone() const override {
        return create(protein_structure_.sequence());
    }

    evo_alg::fitness::FitnessValue normalize(evo_alg::fitness::FitnessValue const& fitness_value) override {
        double const ss_value = fitness_value[0] / ss_total;
        double const cm_value = fitness_value[1] / cm_total;
        double const energy_value = (-fitness_value[2] + energy_min) / (energy_max - energy_min);

        return {ss_value, cm_value, max(0.0, min(1.0, energy_value))};
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

        return {ss_score, cm_score, -rosetta_score};
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

#endif