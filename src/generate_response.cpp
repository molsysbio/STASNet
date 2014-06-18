#include "generate_response.hpp"
#include "helper_types.hpp"
#include <boost/lexical_cast.hpp>
#include <iostream>                      // cin, cout, iostrstream...

void generate_response( GiNaC::matrix &response,
			std::vector <GiNaC::symbol> &x,
			const int_matrix &adm, 
			const ExperimentalDesign &exp_design,
			const std::vector<std::string> &names
			) {
  x.clear();
 //  (I) -Build symbolic "r" matrix from adjacency matrix "adm"-
  size_t size=adm.shape()[0];                  // because we need it all the time
  int r_idx[size][size];
  GiNaC::matrix r(size,size);               // local response matrix r

  size_t c=0;  
  for (size_t i=0; i<size; i++){
    for (size_t j=0; j<size; j++){
      std::string tmpstr;                
      if (i==j) r(j,i)=-1;                         // -1 on the diagonal of the local response is required for MRA
      if (adm[j][i]!=0){
	    if (names.size()==size) {
	      tmpstr= "r_" + names[j] + "_" + names[i];      // concatenates to "r_j_i" in tmpstr which can then be converted to str
	    } else {
	      tmpstr = "r_" + boost::lexical_cast<std::string>(j) + "_" + boost::lexical_cast<std::string>(i); 
	    }
            x.push_back(GiNaC::symbol(tmpstr));  // creates another symbol r with the position in the matrix (and adds one more entry to x)
  	    r(j,i)=x[c];                               // puts the symbol in the matrix position
	    r_idx[j][i]=c++;                           // to remember the position of the symbols
      }
      else r_idx[j][i]=-1;
    }
  }
  

  //  (II) -Include inhibitors as symbols-
  int inh_idx[exp_design.inhib_nodes.size()];
  for (size_t i=0; i<exp_design.inhib_nodes.size() ;i++ ) {
    std::string tmpstr;
    if (names.size()==size) {
      tmpstr = "i" + names[exp_design.inhib_nodes[i]];
    } else {
      tmpstr = "i" + boost::lexical_cast<std::string>(i);
    }
    x.push_back(GiNaC::symbol(tmpstr));
    inh_idx[i]=c++;
  }

  //  (III) -Create the global response matrix "R": R = -1 * inv(r)
  //std::cout << r << std::endl;
  GiNaC::matrix R=GiNaC::ex_to<GiNaC::matrix>(r.inverse()).mul_scalar(-1);
  //std::cout << R.expand() << std::endl; // Variation here due to different simplification by GiNaC


  //  (IV) -Create the actual response vector (of what was measured and perturbed including inhibitor effects)-

  response = GiNaC::matrix(exp_design.stimuli.shape()[0]*exp_design.measured_nodes.size(),1);  
  c=0;
  for (size_t i=0;i<exp_design.measured_nodes.size();i++){
    for (size_t j=0;j<exp_design.stimuli.shape()[0];j++){
          response(c,0)=0;
          
          for (size_t k=0; k<exp_design.stimuli.shape()[1]; k++) {   // Add stimuli
            if (exp_design.stimuli[j][k] == 1) { 
	            int s_temp= exp_design.stim_nodes[k];    
	            GiNaC::ex Rtmp=R(exp_design.measured_nodes[i],s_temp);
	            for (size_t l=0;l<exp_design.inhibitor.shape()[1];l++){    // Add the effect of the inhibitor
	                if (exp_design.inhibitor[j][l] == 1) { 
	                    int i_temp = exp_design.inhib_nodes[l];	   
	                    for ( int m=0; m<size; m++) {
		                    if (r_idx[m][i_temp] > -1) { 
		                        Rtmp = Rtmp.subs(x[r_idx[m][i_temp]]==x[r_idx[m][i_temp]] *
			                    (exp(-pow(pow(x[inh_idx[l]],2),0.5)))); // We force the inhibitor to have a negative effect
		                    }
	                    }
	                }
	            }
	            response(c,0)+=Rtmp;
	        }
          }

          for (size_t k=0;k<exp_design.inhibitor.shape()[1];k++){      // Basal inhibition
      	        if (exp_design.inhibitor[j][k] == 1) {
	                int i_temp = exp_design.inhib_nodes[k];
	                GiNaC::ex Rtmp2=0;
	                for (int l=0; l<size;l++){                                        
	                    if ((r_idx[l][i_temp] > -1) & (exp_design.basal_activity[i_temp]==1)) { 
	                        Rtmp2 += R(exp_design.measured_nodes[i],l)
		                    * x[r_idx[l][i_temp]] *
		                    (-pow(pow(x[inh_idx[k]],2),0.5));// !!! for multiple inhibitors acting in the same path this equation is incorrect!
	                    }
	                }
	                response(c,0) += Rtmp2;
	            }
          }
          
          // Basal nodes can only be modulated if they are directly stimulated.
          // Works only for receptors, not for downstream molecules

          for (size_t m=0; m<size;m++) {
	        bool remove_activity=false;
	        if (exp_design.basal_activity[m]==0) { // if node m has no basal activity
	          remove_activity=true;
	          for (size_t stim=0; stim<exp_design.stimuli.shape()[1]; stim++) {
	            if (exp_design.stimuli[j][stim] == 1) 
	              if (adm[m][exp_design.stim_nodes[stim]]!=0) { // if node m is stimulated
		            remove_activity=false;
	              }
	          }
	          if (remove_activity) {
	            for (int from=0; from<size;from++) {
	              if ((adm[m][from] != 0)) { 
		            response(c,0)=response(c,0).subs(x[r_idx[m][from]]==0);
	              }
	            }
	          }

	        }
          }

          c++;
    }
  }
}
