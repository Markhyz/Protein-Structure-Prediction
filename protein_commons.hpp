#ifndef GUARD_H_PROTEIN_COMMONS
#define GUARD_H_PROTEIN_COMMONS

#include "rosetta.hpp"

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

using frag_list_t = vector<vector<vector<vector<double>>>>;

frag_list_t frag3, frag9;
vector<vector<double>> frag3_prob, frag9_prob;

vector<tuple<size_t, size_t, double>> cm;

double max_energy = 0;
double max_energy_fa = 0;

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

void getTrueSS(string true_ss) {
    for (char c : true_ss) {
        double c_prob = 0, h_prob = 0, e_prob = 0;
        if (c == 'L') {
            c_prob = 1;
        } else if (c == 'H') {
            h_prob = 1;
        } else if (c == 'E') {
            e_prob = 1;
        }
        ss += c;
        ss_prob.push_back({c_prob, h_prob, e_prob});
    }
}

double fix_angle(double angle) {
    return angle + (angle < -180 ? 360 : (angle > 180 ? -360 : 0));
}

void getFragments(string protein_name, size_t frag_size, frag_list_t& frag_list,
                  vector<vector<double>>& frag_prob_list) {
    ifstream frag_in("proteins/" + protein_name + "/frag" + to_string(frag_size));

    string _, line;
    size_t pos, n;
    size_t res_index = 0;
    while (getline(frag_in, line)) {
        istringstream iss(line);
        iss >> _ >> pos >> _ >> n;
        getline(frag_in, line);
        vector<vector<vector<double>>> frags;
        vector<double> frags_prob;
        for (size_t index = 0; index < n; ++index) {
            vector<vector<double>> frag;
            double frag_prob = evo_alg::utils::eps;
            for (size_t index2 = 0; index2 < frag_size; ++index2) {
                getline(frag_in, line);
                istringstream iss2(line);
                string res_ss;
                double phi, psi, omega;

                iss2 >> _ >> _ >> _ >> _;
                iss2 >> res_ss >> phi >> psi >> omega;

                if (res_ss == "L")
                    frag_prob += sqrt(ss_prob[res_index + index2][0]);
                if (res_ss == "H")
                    frag_prob += sqrt(ss_prob[res_index + index2][1]);
                if (res_ss == "E")
                    frag_prob += sqrt(ss_prob[res_index + index2][2]);

                frag.push_back({fix_angle(phi), fix_angle(psi), fix_angle(omega)});
            }
            getline(frag_in, line);
            frags.push_back(frag);
            frags_prob.push_back(frag_prob);
        }
        frag_list.push_back(frags);

        double total_prob = accumulate(frags_prob.begin(), frags_prob.end(), 0.0);
        for (double& prob : frags_prob)
            prob /= total_prob;
        assert(evo_alg::utils::numericEqual(accumulate(frags_prob.begin(), frags_prob.end(), 0.0), 1.0));
        frag_prob_list.push_back(frags_prob);

        ++res_index;
    }
}

void getContactMap(string protein_name) {
    ifstream cm_in("proteins/" + protein_name + "/con");

    size_t i, j;
    double lb, ub_, prob;

    vector<pair<size_t, size_t>> res_pairs;
    vector<pair<double, int>> contact_prob;

    size_t idx = 0;
    while (cm_in >> i >> j >> prob) {
        res_pairs.push_back({i, j});
        contact_prob.push_back({prob, idx++});
    }

    sort(contact_prob.rbegin(), contact_prob.rend());

    for (size_t index = 0; index < residue_indexes.size(); ++index) {
        i = res_pairs[contact_prob[index].second].first;
        j = res_pairs[contact_prob[index].second].second;
        cm.push_back({i, j, contact_prob[index].first});
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
    for (size_t index = 0; index < residue_indexes.size(); ++index) {
        cm.push_back(contacts[index]);
    }
}

pair<double, double> ssScore(core::pose::Pose const& pose) {
    core::scoring::dssp::Dssp dssp(pose);

    string protein_ss = dssp.get_dssp_secstruct();

    double ss_score = 0, ss_total = 0;
    for (size_t ss_index = 0; ss_index < protein_ss.size(); ++ss_index) {
        if (ss[ss_index] == 'L')
            continue;

        double og_prob;
        if (ss[ss_index] == 'H') {
            og_prob = ss_prob[ss_index][1];
        }
        if (ss[ss_index] == 'E') {
            og_prob = ss_prob[ss_index][2];
        }

        ss_total += og_prob;

        if (protein_ss[ss_index] == 'L')
            ss_score += ss_prob[ss_index][0];
        if (protein_ss[ss_index] == 'H')
            ss_score += ss_prob[ss_index][1];
        if (protein_ss[ss_index] == 'E')
            ss_score += ss_prob[ss_index][2];
    }

    return {ss_score, ss_total};
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

    RosettaFaEnergyFunction* clone() const override {
        return create(protein_structure_.sequence());
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

    size_t dimension_ = 1;
    core::pose::Pose protein_structure_;
    core::scoring::ScoreFunctionOP energy_score_;
};

#endif