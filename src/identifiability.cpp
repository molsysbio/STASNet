#include "identifiability.hpp"


#include <boost/shared_ptr.hpp>
#include <boost/rational.hpp>
#include <iostream>  
#include <sstream>
#include <cassert>
#include <utility>
#include <fstream>

#include "mathtree.hpp"
#include "rref.hpp"


boost::shared_ptr<MathTree::parameter> parameterlist::getParameterForExpression(GiNaC::ex e) {
  // If expression e has allready been assined to a parameter, return that.
  for (iterator iter=begin(); iter!=end(); iter++) {
    if (iter->second==e) 
      return iter->first;
  }

  // otherwise generate a new parameter
  MathTree::parameter::Ptr par(new MathTree::parameter());
  par->set_parameter(boost::shared_ptr<double>(new double(1.0)));
  push_back(std::make_pair(par,e));
  return par;
}

/*
// Helper function to output a parameterlist does not work 

std::ostream &operator<<(std::ostream &os, parameterlist &a){
  for (int i=0;i<a.size();i++){
      os << a[i] << "\t";
    }
    os << std::endl;
  return os;
}
*/

// Converts an GiNaC expression into a MathTree object. 
// If reduce_products is true, products of symbols will be assigned to one parameter.
// Parameters are stored in parameterlist param
MathTree::math_item::Ptr put_into_mathtree_format (GiNaC::ex e, parameterlist &param, bool reduce_products) {
  
  MathTree::math_item::Ptr item;

  if (GiNaC::is_a<GiNaC::mul>(e)) {
    
    item = MathTree::mul::Ptr(new MathTree::mul());
    GiNaC::ex multiplier=1;
    
    for (size_t i=0; i<e.nops(); ++i) 
      {
	// This is the place where parameter reduction is happening:
	// All multiplies which are symbolic will be reduced to one...
	if (reduce_products) {
	  if (!GiNaC::is_a<GiNaC::symbol>(e.op(i))) {
	    boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));
	  } else {
	    multiplier*=e.op(i);
	  }
	} else {
	  boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));
	}
      }
    
    if (multiplier!=1) {
      MathTree::parameter::Ptr item2=param.getParameterForExpression(multiplier);
      boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(item2);
    }

  } else if (GiNaC::is_a<GiNaC::add>(e)) {

    item = MathTree::add::Ptr(new MathTree::add());
    for (size_t i=0; i<e.nops(); ++i) 
      boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

  } else if (GiNaC::is_a<GiNaC::power>(e)) {

    item = MathTree::pow::Ptr(new MathTree::pow());
    for (size_t i=0; i<e.nops(); ++i) 
      boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

  } else if (GiNaC::is_a<GiNaC::symbol>(e)) {

    item = param.getParameterForExpression(e);

  } else if (GiNaC::is_a<GiNaC::numeric>(e)) {

    item = MathTree::numeric::Ptr(new MathTree::numeric());
    boost::dynamic_pointer_cast<MathTree::numeric>(item)->set_numeric(GiNaC::ex_to<GiNaC::numeric>(e).to_double());

  } else if (GiNaC::is_a<GiNaC::function>(e)) {
    item = MathTree::pow::Ptr(new MathTree::exp());
    for (size_t i=0; i<e.nops(); ++i) 
      boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

  } else {

    std::cerr << "UNKNOWN TYPE: " << e;
    exit(-1);
  }

  return item;

}

template<class T> struct index_cmp {
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  { return arr[a] > arr[b]; }
  const T arr;
};

//resorts the parameter_dependency_matrix and the parameterlist
void sort ( int_matrix &old_pdmu,
	    parameterlist &param) {

  //get the number of variables
  size_t i=0,j=0,vars;  
  while(old_pdmu[i][j]!=-1){j++;}
  vars=j;
  
  // std::cerr << "number of variables: " << vars << std::endl;

  //get the sum of variables for each parameter
  int nr_of_vars[param.size()];  
  int indices[param.size()];
  for (size_t i=0; i<param.size();i++){
    nr_of_vars[i]=0;
    indices[i]=i;
    for(size_t j=0; j<vars;j++){
      nr_of_vars[i]+=std::abs(old_pdmu[i][j]);
    }
  }
  
  std::sort(indices,indices+param.size(),index_cmp<int*>(nr_of_vars) );
  
  
  //change the matrix and the parameter vector accordingly
  parameterlist tmpparam = param;
  int_matrix new_pdmu = old_pdmu; 
  for (size_t i=0;i<param.size();i++){
    tmpparam[i]=param[indices[i]];
    for (size_t j=0; j<vars;j++){
      new_pdmu[i][j]=old_pdmu[indices[i]][j];
    }
  }
  old_pdmu=new_pdmu;
  param=tmpparam;
  
}


// Performs identifiabiliy analysis on the GiNaC matrix input_matrix
void identifiability_analysis(	 equation_matrix &output_matrix, 
				 std::vector<MathTree::parameter::Ptr> &parameters,
				 double_matrix &parameter_dependency_matrix,
				 int_matrix &parameter_dependency_matrix_unreduced,
				 const symbolic_matrix &input_matrix, 
				 const std::vector<GiNaC::symbol> & vars) {


  // Gernerate new parameterisation
  output_matrix.resize(boost::extents[input_matrix.rows()][input_matrix.cols()]);
  parameterlist param;

  for (size_t i=0; i<input_matrix.rows(); i++) 
    for (size_t j=0; j<input_matrix.cols(); j++) 
      output_matrix[i][j]= put_into_mathtree_format(input_matrix(i,j).expand(),param);
  // Write parameter dependencies into matrix (this is the matrix which will 
  // be put into Row Echelon form.
  // Each row represnets one new parameter, the first colums are the old parameters, then the 
  // new ones.
  // It has form A|B where B is -unity matrix and A has the exponents of the products.

  parameter_dependency_matrix_unreduced.resize(boost::extents[param.size()][param.size()+vars.size()]);
  size_t x=0, y=0;
  for (parameterlist::iterator iter=param.begin(); iter!=param.end(); ++iter) {
    y=0;
    assert(GiNaC::is_a<GiNaC::mul>(iter->second) || GiNaC::is_a<GiNaC::symbol>(iter->second) );
    if (GiNaC::is_a<GiNaC::mul>(iter->second)) {
      GiNaC::mul m=GiNaC::ex_to<GiNaC::mul>(iter->second);
      for (std::vector<GiNaC::symbol>::const_iterator iter2=vars.begin(); iter2!=vars.end(); ++iter2) {
	parameter_dependency_matrix_unreduced[x][y++]=(int)round(m.degree(*iter2));
      }
    } else if (GiNaC::is_a<GiNaC::symbol>(iter->second)) {
      for (std::vector<GiNaC::symbol>::const_iterator iter2=vars.begin(); iter2!=vars.end(); ++iter2) {
	if (*iter2 == iter->second)
	  parameter_dependency_matrix_unreduced[x][y++] = 1;
	else 
	  parameter_dependency_matrix_unreduced[x][y++] = 0;
      }
    }
    parameter_dependency_matrix_unreduced[x][vars.size()+x]=-1;
    x++;
  }

  //sort the parameter dependency matrix
  sort (parameter_dependency_matrix_unreduced,param);

  parameter_dependency_matrix.resize(boost::extents[param.size()][param.size()+vars.size()]);
  //  bring into reduced row echelon form
  rational_matrix rational_matrix_for_rref;
  convert_int_to_rational_matrix(parameter_dependency_matrix_unreduced,
				 rational_matrix_for_rref);

  to_reduced_row_echelon_form(rational_matrix_for_rref);

  convert_rational_to_double_matrix(rational_matrix_for_rref, 
				    parameter_dependency_matrix);
  
    

  // Store parameters
  parameters.resize(param.size());
  for (size_t i=0; i<param.size(); i++){ 
    parameters[i]=param[i].first;
  }


}
 
