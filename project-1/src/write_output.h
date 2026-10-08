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
        meshData& mesh,
        std::vector<materialData>& materials
    );

    void write_angular_flux(
        std::ofstream &outfile,
        meshData& mesh,
        solverData& solver
    );

    void write_angular_flux(
        std::ofstream& outfile,
        meshData& mesh,
        quadratureData& quadrature,
        solverData& solver
    );    

    void write_scalar_flux(
        std::ofstream& outfile,
        meshData& mesh,
        quadratureData& quadrature,
        solverData& solver
    );

    void write_current(
        std::ofstream& outfile,
        meshData& mesh,
        quadratureData& quadrature,
        solverData& solver
    );

    void write_residuals(
        std::ofstream& outfile,
        meshData& mesh,
        quadratureData& quadrature,
        solverData& solver
    ); 

    void write_inner_iteration(
        std::ofstream& outfile,
        coreData& core,        
        solverData& solver
    );

    void write_outer_iteration(
        std::ofstream& outfile,
        coreData& core,        
        solverData& solver
    );
};

#endif
