#include "read_input.h"

#include <sstream>
#include <stdexcept>

void Input::read(
    const std::string& filename,
    coreData& core,
    std::vector<materialData>& materials)
{
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "Cannot open input file: " + filename
        );
    }

    read_core(file, core);

    // CORE provides the size needed for MATERIAL.
    materials.assign(core.nmat, materialData{});

    read_material(file, core, materials);
}

void Input::read_core(std::ifstream& file, coreData& core)
{
    file.clear();
    file.seekg(0, std::ios::beg);

    if (!file) {
        throw std::runtime_error("Cannot rewind input file");
    }

    bool in_block = false;
    bool found_nmat = false;
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

        if (card.front() == '[') {
            break;
        }

        if (card == "nsize") {
            ss >> core.nsize;
        } else if (card == "isize") {
            ss >> core.isize;
        } else if (card == "maxin") {
            ss >> core.maxin;
        } else if (card == "maxout") {
            ss >> core.maxout;
        } else if (card == "hx") {
            ss >> core.hx;
        } else if (card == "epsk") {
            ss >> core.epsk;
        } else if (card == "epspow") {
            ss >> core.epspow;
        } else if (card == "nmat") {
            ss >> core.nmat;
            found_nmat = true;
        } else if (card == "search") {
            ss >> core.search;
        } else {
            throw std::runtime_error(
                "Unknown CORE card: " + card
            );
        }

        // Check the extraction performed by the matching branch.
        if (!ss) {
            throw std::runtime_error(
                "Invalid value for CORE card: " + card
            );
        }
    }

    if (!in_block) {
        throw std::runtime_error("Missing [CORE] block");
    }

    if (!found_nmat || core.nmat <= 0) {
        throw std::runtime_error(
            "[CORE] requires a positive nmat"
        );
    }

    // Add mapbc/mapmat parsing according to your map input format.
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

        if (card.front() == '[') {
            break;
        }

        if (card != "material") {
            throw std::runtime_error(
                "Unknown MATERIAL card: " + card
            );
        }

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

        if (seen[mat.id - 1]) {
            throw std::runtime_error(
                "Duplicate material ID: " +
                std::to_string(mat.id)
            );
        }

        materials[mat.id - 1] = mat;
        seen[mat.id - 1] = true;
    }

    if (!in_block) {
        throw std::runtime_error("Missing [MATERIAL] block");
    }

    for (int i = 0; i < core.nmat; ++i) {
        if (!seen[i]) {
            throw std::runtime_error(
                "Missing material ID: " + std::to_string(i + 1)
            );
        }
    }
}
