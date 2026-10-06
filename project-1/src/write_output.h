#ifndef WRITE_OUTPUT_H
#define WRITE_OUTPUT_H

#include "data.h"

#include <fstream>
#include <sstream>

class Output {
public:
    void write_header(std::ofstream &outfile);

    void echo_input(
        std::ofstream &outputFile,
        coreData& core,
        quadratureData& quadrature,
        std::vector<materialData>& materials
    );

    void write_mesh(
        std::ofstream &outfile,
        coreData& core,
        meshData& mesh
    );
};

#endif
