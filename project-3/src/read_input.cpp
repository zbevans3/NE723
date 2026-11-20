#include "read_input.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

void Input::read(const std::string &filename, coreData &core, quadratureData &quadrature,
                 std::vector<materialData> &materials) {
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error("Cannot open input file: " + filename);
    }

    // read core block
    read_core(file, core);

    // allocate materials and read material block
    materials.assign(core.nmat, materialData{});
    read_material(file, core, materials);

    // read quadrature block
    read_quadrature(file, quadrature);
}

void Input::read_core(std::ifstream &file, coreData &core) {
    file.clear();
    file.seekg(0, std::ios::beg);

    if (!file) {
        throw std::runtime_error("Cannot rewind input file");
    }

    bool in_block = false;
    std::string line;

    while (std::getline(file, line)) {
        const auto comment = line.find('!');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream ss(line);
        std::string card;

        if (!(ss >> card)) {
            continue;
        }

        if (!in_block) {
            if (card == "[CORE]") {
                in_block = true;
            }
            continue;
        }

        // exit if another block is reached
        if (card.front() == '[') {
            break;
        }

        // read cards
        if (card == "xedge") {
            double vals;
            core.xedge.clear();
            while (ss >> vals) {
                core.xedge.push_back(vals);
            }
            core.nregions = int(core.xedge.size());
        } else if (card == "nx") {
            double vals;
            core.nx.clear();
            while (ss >> vals) {
                core.nx.push_back(vals);
            }
        } else if (card == "matid") {
            double vals;
            core.matid.clear();
            while (ss >> vals) {
                core.matid.push_back(vals);
            }
            for (size_t i = 0; i < core.matid.size(); i++) {
                if (core.matid[i] < 1) {
                    throw std::runtime_error("Error in material ids. Ids must be greater than 0.");
                }
            }
            core.nmat = *std::max_element(core.matid.begin(), core.matid.end());
        } else if (card == "maxinner") {
            ss >> core.maxin;
        } else if (card == "maxouter") {
            ss >> core.maxout;
        } else if (card == "epsk") {
            ss >> core.epsk;
        } else if (card == "epsflx") {
            ss >> core.epsflx;
        } else if (card == "epspow") {
            ss >> core.epspow;
        } else if (card == "search") {
            ss >> core.search;
        } else if (card == "scheme") {
            ss >> core.scheme;
        } else if (card == "acceleration") {
            ss >> core.acceleration;
        } else if (card == "ngroups") {
            ss >> core.ngroups;
        } else if (card == "kcrit") {
            ss >> core.kcrit;
        } else if (card == "bcleft") {
            ss >> core.bcleft;
            if (core.bcleft == "incoming") {
                ss >> core.psileft;
            }
        } else if (card == "bcright") {
            ss >> core.bcright;
            if (core.bcright == "incoming") {
                ss >> core.psiright;
            }
        } else if (card == "itdebug") {
            std::string tmpstr;
            ss >> tmpstr;
            std::transform(tmpstr.begin(), tmpstr.end(), tmpstr.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            if (tmpstr == "T" or tmpstr == "TRUE") {
                core.itdebug = true;
            } else {
                core.itdebug = false;
            }
        } else if (card == "otdebug") {
            std::string tmpstr;
            ss >> tmpstr;
            std::transform(tmpstr.begin(), tmpstr.end(), tmpstr.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            if (tmpstr == "T" or tmpstr == "TRUE") {
                core.otdebug = true;
            } else {
                core.otdebug = false;
            }
        } else {
            throw std::runtime_error("Unknown CORE card: " + card);
        }
    }

    if (!in_block) {
        throw std::runtime_error("Missing [CORE] block");
    }

    if (core.nmat <= 0) {
        throw std::runtime_error("[CORE] requires a positive nmat");
    }
    if (int(core.nx.size()) != core.nregions) {
        throw std::runtime_error("xedge and nx must have the same length.");
    }
    if (int(core.matid.size()) != core.nregions) {
        throw std::runtime_error("xedge and matid must have the same length.");
    }
}

void Input::read_material(std::ifstream &file, const coreData &core, std::vector<materialData> &materials) {
    file.clear();
    file.seekg(0, std::ios::beg);

    if (!file) {
        throw std::runtime_error("Cannot rewind input file");
    }

    bool in_block = false;
    std::vector<bool> seen(core.nmat, false);
    std::string line;
    int _ngroups = core.ngroups;
    int id = 0; // id of the material currently being read
    double xs = 0.0;

    while (std::getline(file, line)) {
        const auto comment = line.find('!');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream ss(line);
        std::string card;

        if (!(ss >> card)) {
            continue;
        }

        if (!in_block) {
            if (card == "[MATERIAL]") {
                in_block = true;
            }
            continue;
        }

        // exit if another block is reached
        if (card.front() == '[') {
            break;
        }

        // found material
        if (card == "material") {
            ss >> id;
            if (id < 1 || id > core.nmat) {
                throw std::runtime_error("Material ID out of range: " + std::to_string(id));
            }
            // check if duplicate
            if (seen[id - 1]) {
                throw std::runtime_error("Duplicate material ID: " + std::to_string(id));
            }
            seen[id - 1] = true;
            materials[id - 1].id = id;
            continue;
        }

        if (id == 0) {
            throw std::runtime_error("Card " + card + " comes before any material line");
        }
        materialData &mat = materials[id - 1];

        if (card == "total") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.xstot.push_back(xs);
            }
        } else if (card == "absorption") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.xsabs.push_back(xs);
            }
        } else if (card == "scattering") {
            {
                // full matrix: the next ngroups lines are the rows of the matrix,
                // row g = scattering from group g into each group gp (same layout as C5G7 tables)
                mat.xsscatmat.assign(_ngroups, std::vector<double>(_ngroups, 0.0));
                for (int g = 0; g < _ngroups; g++) {
                    std::string row;
                    std::getline(file, row);
                    std::istringstream sr(row);
                    for (int gp = 0; gp < _ngroups; gp++) {
                        sr >> mat.xsscatmat[g][gp];
                    }
                    // keep the within-group value in xsscat as well
                    mat.xsscat.push_back(mat.xsscatmat[g][g]);
                }
            }
        } else if (card == "fission") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.xsfiss.push_back(xs);
            }
        } else if (card == "nu") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.nu.push_back(xs);
            }
        } else if (card == "chi") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.xschi.push_back(xs);
            }
        } else if (card == "fixed_source") {
            for (int g = 0; g < _ngroups; g++) {
                ss >> xs;
                mat.q.push_back(xs);
            }
        } else {
            throw std::runtime_error("Unknown MATERIAL card: " + card);
        }

        // check if correct number of XS are provided for ngroups
        if (ss.fail() || ss >> xs) {
            throw std::runtime_error("Material " + std::to_string(id) + ": " + card + " needs exactly " +
                                     std::to_string(_ngroups) + " values");
        }
    }

    // throw error if materials never actually defined...
    if (!in_block) {
        throw std::runtime_error("Missing [MATERIAL] block");
    }

    // check for missing material
    for (int i = 0; i < core.nmat; ++i) {
        if (!seen[i]) {
            throw std::runtime_error("Missing material ID: " + std::to_string(i + 1));
        }
    }

    // fission, nu, chi and fixed_source are optional
    for (int i = 0; i < core.nmat; ++i) {
        materialData &mat = materials[i];
        if (mat.xsfiss.empty())
            mat.xsfiss.assign(_ngroups, 0.0);
        if (mat.nu.empty())
            mat.nu.assign(_ngroups, 0.0);
        if (mat.xschi.empty())
            mat.xschi.assign(_ngroups, 0.0);
        if (mat.q.empty())
            mat.q.assign(_ngroups, 0.0);

        if (int(mat.xstot.size()) != _ngroups || int(mat.xsabs.size()) != _ngroups ||
            int(mat.xsscat.size()) != _ngroups || int(mat.xsfiss.size()) != _ngroups ||
            int(mat.nu.size()) != _ngroups || int(mat.xschi.size()) != _ngroups || int(mat.q.size()) != _ngroups) {
            throw std::runtime_error("Material " + std::to_string(i + 1) +
                                     ": total, absorption and scattering are required");
        }
    }
}

void Input::read_quadrature(std::ifstream &file, quadratureData &quadrature) {
    file.clear();
    file.seekg(0, std::ios::beg);

    if (!file) {
        throw std::runtime_error("Cannot rewind input file");
    }

    bool in_block = false;
    std::string line;

    while (std::getline(file, line)) {
        const auto comment = line.find('!');
        if (comment != std::string::npos) {
            line.erase(comment);
        }

        std::istringstream ss(line);
        std::string card;

        if (!(ss >> card)) {
            continue;
        }

        if (!in_block) {
            if (card == "[QUADRATURE]") {
                in_block = true;
            }
            continue;
        }

        // exit if another block is reached
        if (card.front() == '[') {
            break;
        }

        // read cards
        if (card == "order") {
            ss >> quadrature.order;
            if (quadrature.order < 2) {
                throw std::runtime_error("Quadrature order must be greater than or equal to 2.");
            }
        } else if (card == "angle") {
            double vals;
            quadrature.angle.clear();
            while (ss >> vals) {
                quadrature.angle.push_back(vals);
            }
            if (int(quadrature.angle.size()) != quadrature.order) {
                throw std::runtime_error("Number of angles must be equal to quadrature order.");
            }
        } else if (card == "weight") {
            double vals;
            quadrature.weight.clear();
            while (ss >> vals) {
                quadrature.weight.push_back(vals);
            }
            if (int(quadrature.weight.size()) != quadrature.order) {
                throw std::runtime_error("Number of weights must be equal to quadrature order.");
            }
        } else {
            throw std::runtime_error("Unknown QUADRATURE card: " + card);
        }
    }

    if (!in_block) {
        throw std::runtime_error("Missing [QUADRATURE] block");
    }
}
