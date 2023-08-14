#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

using namespace std;

struct config_t {
    string protein_fasta_path;
    string protein_pdb_path;
    string protein_frag3_path;
    string protein_frag9_path;
    string protein_ss_path;
    string protein_cm_path;
    string score_output_dir;
    string decoy_output_dir;
    string algorithm_output_dir;

    size_t iteration_num;
    size_t population_size;

    double ss_w_coil;
    double ss_w_helix;
    double ss_w_sheet;

    string frag_type;

    size_t pose_start = 1;

    string output_level = "all";
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
        if (value_name == "fasta_path") {
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
        } else if (value_name == "algorithm_output_dir") {
            config.algorithm_output_dir = normalizeDirPath(value);
        } else if (value_name == "iterations") {
            config.iteration_num = stoul(value);
        } else if (value_name == "population_size") {
            config.population_size = stoul(value);
        } else if (value_name == "ss_coil_weight") {
            config.ss_w_coil = stod(value);
        } else if (value_name == "ss_helix_weight") {
            config.ss_w_helix = stod(value);
        } else if (value_name == "ss_sheet_weight") {
            config.ss_w_sheet = stod(value);
        } else if (value_name == "frag_type") {
            config.frag_type = value;
        } else if (value_name == "protein_offset") {
            config.pose_start = stoul(value);
        } else if (value_name == "output_level") {
            config.output_level = value;
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
            } else {
                token += c;
            }
        }
        if (!token.empty()) {
            tokens.push_back(token);
        }
        if (tokens.size() > 0) {
            if (tokens.size() != 3 || tokens[0] == "=" || tokens[1] != "=" || tokens[2] == "=") {
                cout << "error: invalid config file" << endl;

                exit(EXIT_FAILURE);
            }
            string key;
            for (char c : tokens[0])
                key += tolower(c);
            config_values[key] = tokens[2];
        }
    }
}

void setConfig(string config_file_name, config_t& config) {
    map<string, string> config_values;

    readConfigFile(config_file_name, config_values);

    setConfigValue(config, config_values, "fasta_path", true);
    setConfigValue(config, config_values, "pdb_path", true);
    setConfigValue(config, config_values, "frag3_path", true);
    setConfigValue(config, config_values, "frag9_path", true);
    setConfigValue(config, config_values, "ss_path", true);
    setConfigValue(config, config_values, "cm_path", true);
    setConfigValue(config, config_values, "score_output_dir", true);
    setConfigValue(config, config_values, "decoy_output_dir", true);
    setConfigValue(config, config_values, "algorithm_output_dir", true);
    setConfigValue(config, config_values, "iterations", true);
    setConfigValue(config, config_values, "population_size", true);
    setConfigValue(config, config_values, "ss_coil_weight", false);
    setConfigValue(config, config_values, "ss_helix_weight", false);
    setConfigValue(config, config_values, "ss_sheet_weight", false);
    setConfigValue(config, config_values, "frag_type", false);
    setConfigValue(config, config_values, "protein_offset", false);
    setConfigValue(config, config_values, "output_level", false);
}