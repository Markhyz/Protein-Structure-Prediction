#ifndef GUARD_H_BRKGA_COMMONS
#define GUARD_H_BRKGA_COMMONS

#include "protein_commons.hpp"

#include <evo_alg/algorithms/brkga.hpp>

vector<double> frag3Decoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome(angle_types.size());
    for (size_t gene_index = 0; gene_index < encoded_chromosome.size(); ++gene_index) {
        size_t res_index = gene_index * 3;
        if (res_index >= frag3.size()) {
            res_index -= res_index - frag3.size() + 1;
        }
        size_t frag_index = floor(encoded_chromosome[gene_index] * frag3[res_index].size());
        vector<vector<double>> fragment = frag3[res_index][frag_index];
        for (size_t res_offset = 0; res_offset < fragment.size(); ++res_offset) {
            size_t chromosome_residue_index = residue_indexes[res_index + res_offset];
            decoded_chromosome[chromosome_residue_index] = fragment[res_offset][0];
            decoded_chromosome[chromosome_residue_index + 1] = fragment[res_offset][1];
            decoded_chromosome[chromosome_residue_index + 2] = fragment[res_offset][2];
        }
    }

    return decoded_chromosome;
}

vector<double> frag9Decoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome(angle_types.size());
    for (size_t gene_index = 0; gene_index < encoded_chromosome.size(); ++gene_index) {
        size_t res_index = gene_index * 9;
        if (res_index >= frag9.size()) {
            res_index -= res_index - frag9.size() + 1;
        }
        size_t frag_index = floor(encoded_chromosome[gene_index] * frag9[res_index].size());
        vector<vector<double>> fragment = frag9[res_index][frag_index];
        for (size_t res_offset = 0; res_offset < fragment.size(); ++res_offset) {
            size_t chromosome_residue_index = residue_indexes[res_index + res_offset];
            decoded_chromosome[chromosome_residue_index] = fragment[res_offset][0];
            decoded_chromosome[chromosome_residue_index + 1] = fragment[res_offset][1];
            decoded_chromosome[chromosome_residue_index + 2] = fragment[res_offset][2];
        }
    }

    return decoded_chromosome;
}

vector<double> angleDecoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t res_index = 0; res_index < residue_indexes.size(); ++res_index) {
        decoded_chromosome.push_back(-180 + encoded_chromosome[res_index * 3] * 360);
        decoded_chromosome.push_back(-180 + encoded_chromosome[res_index * 3 + 1] * 360);
        decoded_chromosome.push_back(-180 + encoded_chromosome[res_index * 3 + 2] * 360);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index]; ++chi_idx)
            decoded_chromosome.push_back(0.0);
    }

    return decoded_chromosome;
}

vector<double> residueDecoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t res_index = 0; res_index < residue_indexes.size(); ++res_index) {
        double gene = encoded_chromosome[res_index];
        double phi = trunc(gene * 1e5) / 1e5;
        gene = (gene * 1e5 - trunc(gene * 1e5));
        double psi = trunc(gene * 1e5) / 1e5;
        gene = (gene * 1e5 - trunc(gene * 1e5));
        double omega = trunc(gene * 1e5) / 1e5;

        decoded_chromosome.push_back(-180 + phi * 360);
        decoded_chromosome.push_back(-180 + psi * 360);
        decoded_chromosome.push_back(-180 + omega * 360);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index]; ++chi_idx)
            decoded_chromosome.push_back(0.0);
    }

    return decoded_chromosome;
}

vector<double> protein_angles;

vector<double> refineDecoder(vector<double> const& encoded_chromosome) {
    vector<double> decoded_chromosome;
    for (size_t res_index = 0; res_index < residue_indexes.size(); ++res_index) {
        double phi_delta = -60 + encoded_chromosome[res_index * 2] * 120;
        double psi_delta = -60 + encoded_chromosome[res_index * 2 + 1] * 120;
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index]] + phi_delta);
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index] + 1] + psi_delta);
        decoded_chromosome.push_back(protein_angles[residue_indexes[res_index] + 2]);

        for (size_t chi_idx = 0; chi_idx < residue_chi_num[res_index]; ++chi_idx)
            decoded_chromosome.push_back(0.0);
    }

    return decoded_chromosome;
}

#endif