///////////////////////////////////////////////////////////////////////////////
//
// Generation of the Model numeric equations from the symbolic equations
// Copyright (C) 2013- Mathurin Dorel, Bertram Klinger, Nils Bluthgen
//
// Institute of Pathology and Institute for Theoretical Biology
// Charite - Universitätsmedizin Berlin - Chariteplatz 1, 10117 Berlin, Germany
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
///////////////////////////////////////////////////////////////////////////////

#include "generate_response.hpp"
#include "helper_types.hpp"
#include <boost/lexical_cast.hpp>
#include <iostream>                      // cin, cout, iostrstream...

extern bool debug;
extern bool verbosity;

void generate_response( GiNaC::matrix &response,
            std::vector <GiNaC::symbol> &x,
            const ModelStructure &structure,
            const ExperimentalDesign &exp_design) {
  //  (I) -Retrieve the symbolic "r" matrix, adjacency matrix "adm"
  const std::vector<std::string> names = structure.getNames();
  int_matrix r_idx;
  copy_matrix(structure.getRidx(), r_idx);
  int_matrix adm;
  copy_matrix(structure.getAdjacencyMatrix(), adm);
  const GiNaC::matrix r = structure.getR();
  x = structure.getSymbols();
  const size_t size=adm.shape()[0];

  // (II) -Include inhibitors as symbols-
  size_t c = x.size();
  if (debug) { std::cerr << "Including the inhibitors..." << std::endl; }
  int inh_idx[exp_design.inhib_nodes.size()];
  for (size_t i=0; i<exp_design.inhib_nodes.size() ;i++ ) {
    std::string tmpstr;
    if (names.size()==size) {
      tmpstr = "i" + names[exp_design.inhib_nodes[i]];
    } else {
      tmpstr = "i" + boost::lexical_cast<std::string>(i);
    }
    x.push_back(GiNaC::symbol(tmpstr));     // Put the inhibitor symbol in the vector
    inh_idx[i]=c++;                         // Remember the position of the inhibitor symbol
  }
  if (debug) { for(size_t ii=0; ii<x.size();ii++) { std::cout << " " << x[ii]; } std::cout << std::endl; }

  //  (III) -Create the global response matrix "R": R = -1 * inv(r)
  //std::cout << r << std::endl;
  
  if (debug) {
    std::cerr << "Inverting the response matrix..." << std::endl;
    std::cerr << "dim = " << size << "x" << size << std::endl;
  }
  GiNaC::matrix R=GiNaC::matrix(r.inverse()).mul_scalar(-1);
  //std::cout << R.expand() << std::endl; // Variation here due to different simplification by GiNaC


  //  (IV) -Create the actual response vector (of what was measured and perturbed including inhibitor effects)-
  if (debug) { std::cerr << "Creating of the actual response vector..." << std::endl; }
  response = GiNaC::matrix(exp_design.stimuli.shape()[0]*exp_design.measured_nodes.size(),1);  
  c=0;
  for (size_t i=0;i<exp_design.measured_nodes.size();i++){
    for (size_t j=0;j<exp_design.stimuli.shape()[0];j++){
          response(c,0)=0;

          // Effect of the inhibitor, record which links should be multiplied by the inhibitory term to build r tilde
          std::vector< std::pair<GiNaC::symbol, GiNaC::ex> > inhibited_links;
          for (size_t l=0;l<exp_design.inhibitor.shape()[1];l++){
            if (exp_design.inhibitor[j][l] == 1) { 
                int i_temp = exp_design.inhib_nodes[l];    
                for ( int m=0; m<size; m++) {
                    if (r_idx[m][i_temp] > -1) { 
                        inhibited_links.push_back(std::make_pair( x[r_idx[m][i_temp]], x[r_idx[m][i_temp]] *
                        (exp(-pow(pow(x[inh_idx[l]],2),0.5))) )); // Force the inhibitor action to be between 0 and 1
                        if (verbosity) {
                            std::cerr << "Substitute " << x[r_idx[m][i_temp]] << " by " << x[r_idx[m][i_temp]] *
                            (exp(-pow(pow(x[inh_idx[l]],2),0.5))) << std::endl << std::endl;
                        }
                    }
                }
            }
          }
          
          // Effect of the simuli
          for (size_t k=0; k<exp_design.stimuli.shape()[1]; k++) {
            if (exp_design.stimuli[j][k] == 1) { 
                int s_temp= exp_design.stim_nodes[k];    
                GiNaC::ex Rtmp=R(exp_design.measured_nodes[i],s_temp);
                // Substitute the inhibited links (r by r*i)
                for (size_t l=0 ; l < inhibited_links.size() ; l++) {
                    Rtmp = Rtmp.subs(inhibited_links[l].first == inhibited_links[l].second);
                }
                response(c,0)+=Rtmp;
            }
          }

          // Effect of the inhibition on the basal activity
          for (size_t k=0;k<exp_design.inhibitor.shape()[1];k++){
                if (exp_design.inhibitor[j][k] == 1) {
                    int i_temp = exp_design.inhib_nodes[k];
                    GiNaC::ex Rtmp2=0;
                    for (int l=0; l<size;l++){                                        
                        // Propagate the effect from the nodes downstream to the inhibited node to the measured node
                        if (exp_design.basal_activity[i_temp]==1) { 
                            if (r_idx[l][i_temp] > -1) {
                                // Add the dampening factors to the propagation
                                // The dampening of the directly inhibited links will only take place in case of feed-back
                                GiNaC::ex Rtmp3 = R(exp_design.measured_nodes[i],l);
                                for (size_t m=0 ; m < inhibited_links.size() ; m++) {
                                    Rtmp3 = Rtmp3.subs(inhibited_links[m].first == inhibited_links[m].second);
                                }
                                // Negative perturbation due to the inhibition
                                Rtmp2 += Rtmp3 * x[r_idx[l][i_temp]] *
                                (-pow(pow(x[inh_idx[k]],2),0.5)); // Force the inhibitor effect to be negative
                            }
                        }
                    }
                    // Substitute the inhibited links (r by r*i)
                    response(c,0) += Rtmp2;
                }
          }
      /*
          // Suppress the effect on all the nodes without basal activity that are not directly stimulated
          for (size_t m=0; m<size;m++) {
            bool remove_activity=false;
            if (exp_design.basal_activity[m]==0) { 
                // Select nodes without basal activity, but unselect those directly targeted by the stimulation (i.e. the receptor for the stimulation)
                remove_activity=true;
                for (size_t stim=0; stim<exp_design.stimuli.shape()[1]; stim++) {
                    if (exp_design.stimuli[j][stim] == 1 && adm[m][exp_design.stim_nodes[stim]]!=0) { 
                        remove_activity=false;
                    }
                }
                // Remove the transmission by all input nodes (it includes the feedbacks for the non activated receptors)
                if (remove_activity) {
                    for (int from=0; from<size;from++) {
                        if ((adm[m][from] != 0)) { 
                            response(c,0)=response(c,0).subs(x[r_idx[m][from]]==0);
                        }
                    }
                }

            }
          }
      */
          c++; // Next condition
    }
  }
  //  std::cout << response << std::endl;
}
