#ifndef READ_INPUT_H
#define READ_INPUT_H

#include "data.h"

#include <fstream>
#include <string>

class Input {
public:
    void read(
        const std::string& filename,
        coreData& core,
        std::vector<materialData>& materials
    );

private:
    void read_core(
        std::ifstream& file,
        coreData& core
    );

    void read_material(
        std::ifstream& file,
        const coreData& core,
        std::vector<materialData>& materials
    );
};

#endif
