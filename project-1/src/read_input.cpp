#include "read_input.h"

#include <sstream>
#include <stdexcept>
#include <algorithm>

void Input::read(
    const std::string& filename,
    coreData& core,
    quadratureData& quadrature,
    std::vector<materialData>& materials)
{
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "Cannot open input file: " + filename
        );
    }

    // read core block
    read_core(file, core);

    // allocate materials and read material block
    materials.assign(core.nmat, materialData{});
    read_material(file, core, materials);

    // read quadrature block
    read_quadrature(file, quadrature);
}

void Input::read_core(std::ifstream& file, coreData& core)
{
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
            while(ss >> vals){
                core.xedge.push_back(vals);
            }
            core.nregions = int(core.xedge.size());
        } else if (card == "nx") {
            double vals;
            core.nx.clear();
            while(ss >> vals){
                core.nx.push_back(vals);
            }
        } else if (card == "matid") {
            double vals;
            core.matid.clear();
            while(ss >> vals){
                core.matid.push_back(vals);
            }
            for (size_t i = 0 ; i < core.matid.size() ; i++){
               if(core.matid[i] < 1) {
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
        } else if (card == "kcrit") {
            ss >> core.kcrit;            
        } else if (card == "bcleft") {
            ss >> core.bcleft;
            if (core.bcleft == "incoming"){
                ss >> core.psileft;
            }
        } else if (card == "bcright") {
            ss >> core.bcright;
            if (core.bcleft == "incoming"){
                ss >> core.psiright;
            }
        } else if (card == "itdebug") {
            std::string tmpstr;
            ss >> tmpstr;
            std::transform(tmpstr.begin(), tmpstr.end(), tmpstr.begin(), [](unsigned char c) {
                return std::toupper(c);
            });            
            if (tmpstr == "T" or tmpstr == "TRUE"){
                core.itdebug = true;
            } else {
                core.itdebug = false;
            }
        } else if (card == "otdebug") {
            std::string tmpstr;
            ss >> tmpstr;
            std::transform(tmpstr.begin(), tmpstr.end(), tmpstr.begin(), [](unsigned char c) {
                return std::toupper(c);
            });            
            if (tmpstr == "T" or tmpstr == "TRUE"){
                core.otdebug = true;
            } else {
                core.otdebug = false;
            }            
        } else {
            throw std::runtime_error(
                "Unknown CORE card: " + card
            );
        }
    }

    if (!in_block) {
        throw std::runtime_error("Missing [CORE] block");
    }

    if (core.nmat <= 0) {
        throw std::runtime_error(
            "[CORE] requires a positive nmat"
        );
    }
    if (int(core.nx.size()) != core.nregions) {
       throw std::runtime_error("xedge and nx must have the same length.");
    }
    if (int(core.matid.size()) != core.nregions) {
       throw std::runtime_error("xedge and nx must have the same length.");
    }
}

void Input::read_material(
    std::ifstream& file,
    const coreData& core,
    std::vector<materialData>& materials)
{
    file.clear();
    file.seekg(0, std::ios::beg);

    if (!file) {
        throw std::runtime_error("Cannot rewind input file");
    }

    bool in_block = false;
    std::vector<bool> seen(core.nmat, false);
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
            if (card == "[MATERIAL]") {
                in_block = true;
            }
            continue;
        }

        // exit if another block is reached        
        if (card.front() == '[') {
            break;
        }

        // throw error if wrong input is provided
        if (card != "material") {
            throw std::runtime_error(
                "Unknown MATERIAL card: " + card
            );
        }

        // read material into temporary struct
        materialData mat;

        if (!(ss >> mat.id
                 >> mat.xstot
                 >> mat.xsscat
                 >> mat.xsfiss
                 >> mat.nu
                 >> mat.q)) {
            throw std::runtime_error(
                "Expected: material <id> <xstot> <xsscat> "
                "<xsfiss> <nu> <q>"
            );
        }

        if (mat.id < 1 || mat.id > core.nmat) {
            throw std::runtime_error(
                "Material ID out of range: " +
                std::to_string(mat.id)
            );
        }

        // check if duplicate
        if (seen[mat.id - 1]) {
            throw std::runtime_error(
                "Duplicate material ID: " +
                std::to_string(mat.id)
            );
        }

        // fill actual structu
        materials[mat.id - 1] = mat;
        seen[mat.id - 1] = true;
    }

    // throw error if materials never actually defined...
    if (!in_block) {
        throw std::runtime_error("Missing [MATERIAL] block");
    }

    // check for missing material
    for (int i = 0; i < core.nmat; ++i) {
        if (!seen[i]) {
            throw std::runtime_error(
                "Missing material ID: " + std::to_string(i + 1)
            );
        }
    }
}

void Input::read_quadrature(std::ifstream& file, quadratureData& quadrature)
{
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
            if (quadrature.order < 2){
                throw std::runtime_error("Quadrature order must be greater than or equal to 2.");
            }
        } else if (card == "angle") {
            double vals;
            quadrature.angle.clear();
            while(ss >> vals){
                quadrature.angle.push_back(vals);
            }
            if(int(quadrature.angle.size()) != quadrature.order){
                throw std::runtime_error(
                    "Number of angles must be equal to quadrature order."
                );
            }
        } else if (card == "weight") {
            double vals;
            quadrature.weight.clear();
            while(ss >> vals){
                quadrature.weight.push_back(vals);
            }
            if(int(quadrature.weight.size()) != quadrature.order){
                throw std::runtime_error(
                    "Number of weights must be equal to quadrature order."
                );
            }
        } else {
            throw std::runtime_error(
                "Unknown QUADRATURE card: " + card
            );
        }
    }

    if (!in_block) {
        throw std::runtime_error("Missing [QUADRATURE] block");
    }
}