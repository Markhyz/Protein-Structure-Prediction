#include "protein_brkga.hpp"

#include <fstream>
#include <iostream>

using namespace std;

struct config_t {
    string name;
    string protein_fasta_path;
    string protein_pdb_path;
    string protein_frag3_path;
    string protein_frag9_path;
    string protein_ss_path;
    string protein_cm_path;
    string score_output_dir;
    string decoy_output_dir;
    string generation_output_dir;

    size_t iteration_num;
    size_t population_size;
    double elite_fraction = 0.5;
    double mutant_fraction = 0.2;
    double crossover_prob = 0.55;
    double diversity_threshold = 0.2;
    double diversity_enforcement = 0.5;

    size_t pose_start = 1;
};

string normalizeDirPath(string path) {
    return path.back() != '/' ? path + "/" : path;
}

void setConfigValue(config_t& config, map<string, string>& config_values, string value_name, bool required) {
    string value = config_values[value_name];
    if (value.empty()) {
        if (required) {
            cout << "error: missing required parameter \"" << value_name << "\"" << endl;

            exit(EXIT_FAILURE);
        }
    } else {
        if (value_name == "name") {
            config.name = value;

        } else if (value_name == "fasta_path") {
            config.protein_fasta_path = value;

        } else if (value_name == "pdb_path") {
            config.protein_pdb_path = value;
            
        } else if (value_name == "frag3_path") {
            config.protein_frag3_path = value;

        } else if (value_name == "frag9_path") {
            config.protein_frag9_path = value;

        } else if (value_name == "ss_path") {
            config.protein_ss_path = value;

        } else if (value_name == "cm_path") {
            config.protein_cm_path = value;

        } else if (value_name == "score_output_dir") {
            config.score_output_dir = normalizeDirPath(value);

        } else if (value_name == "decoy_output_dir") {
            config.decoy_output_dir = normalizeDirPath(value);

        } else if (value_name == "generation_output_dir") {
            config.generation_output_dir = normalizeDirPath(value);

        } else if (value_name == "iterations") {
            config.iteration_num = stoul(value);

        } else if (value_name == "population_size") {
            config.population_size = stoul(value);

        } else if (value_name == "elite_fraction") {
            config.elite_fraction = stod(value);

        } else if (value_name == "mutant_fraction") {
            config.mutant_fraction = stod(value);

        } else if (value_name == "crossover_prob") {
            config.crossover_prob = stod(value);

        } else if (value_name == "diversity_threshold") {
            config.diversity_threshold = stod(value);

        } else if (value_name == "diversity_enforcement") {
            config.diversity_enforcement = stod(value);
            
        } else if (value_name == "protein_offset") {
            config.pose_start = stoul(value);
        }
    }
}

void readConfigFile(string config_file_name, map<string, string>& config_values) {
    ifstream config_file(config_file_name);
    string config_line;
    while (getline(config_file, config_line)) {
        vector<string> tokens;
        string token;
        for (char c : config_line) {
            if (c == ' ') {
                if (!token.empty()) {
                    tokens.push_back(token);
                    token = "";
                }
            } else if (c == '=') {
                if (!token.empty()) {
                    tokens.push_back(token);
                }
                tokens.push_back("=");
                token = "";
            } else if (c == '#') {
                if (!token.empty()) {
                    tokens.push_back(token);
                }
                token = "";
                break;
            } else token += tolower(c);
        }
        if (!token.empty()) {
            tokens.push_back(token);
        }
        if (tokens.size() > 0) {
            if (tokens.size() != 3 || tokens[0] == "=" || tokens[1] != "=" || tokens[2] == "=") {
                cout << "error: invalid config file" << endl;

                exit(EXIT_FAILURE);
            }
            config_values[tokens[0]] = tokens[2];
        }
    }
}

void setConfig(string config_file_name, config_t& config) {
    map<string, string> config_values;

    readConfigFile(config_file_name, config_values);

    setConfigValue(config, config_values, "name", true);
    setConfigValue(config, config_values, "fasta_path", true);
    setConfigValue(config, config_values, "pdb_path", true);
    setConfigValue(config, config_values, "frag3_path", true);
    setConfigValue(config, config_values, "frag9_path", true);
    setConfigValue(config, config_values, "ss_path", true);
    setConfigValue(config, config_values, "cm_path", true);
    setConfigValue(config, config_values, "score_output_dir", true);
    setConfigValue(config, config_values, "decoy_output_dir", true);
    setConfigValue(config, config_values, "generation_output_dir", true);
    setConfigValue(config, config_values, "iterations", true);
    setConfigValue(config, config_values, "population_size", true);
    setConfigValue(config, config_values, "elite_fraction", false);
    setConfigValue(config, config_values, "mutant_fraction", false);
    setConfigValue(config, config_values, "crossover_prob", false);
    setConfigValue(config, config_values, "diversity_threshold", false);
    setConfigValue(config, config_values, "diversity_enforcement", false);
    setConfigValue(config, config_values, "protein_offset", false);
}