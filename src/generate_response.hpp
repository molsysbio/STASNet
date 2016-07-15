#include <ginac/ginac.h>                 // symbolic tool box for C++ 
#include <ginac/matrix.h>                // to obtain matrices which can contain
#include <vector>
#include "helper_types.hpp"
#include "modelstructure.hpp"

void generate_response( GiNaC::matrix &response,
            std::vector <GiNaC::symbol> &x,
            const ModelStructure &structure,
            const ExperimentalDesign &exp_design
            ) ;
