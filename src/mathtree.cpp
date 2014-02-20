#include "mathtree.hpp"
#include <cmath>
#include <iostream>

// These are helper classes to evaluate the model quickly
// They will be used to construct a tree that represents a formula

namespace MathTree {

  int parameter_counter=0;

  double math_item::eval() const { return -10000000.0; }
  void math_item::print(std::ostream &os) const {}
  
  
  parameter::parameter() : math_item() { parameter_index_=parameter_counter; parameter_counter++ ;}
  double parameter::eval() const { 
    return *num_; 
  }
  void parameter::print(std::ostream &os) const { os << "p" << parameter_index_; }
  void parameter::set_parameter(boost::shared_ptr<double> d) { num_=d; }
  void parameter::set_parameter(double d) { *num_=d; }
  void parameter::operator=(double d) {  *num_=d; }

  
  double numeric::eval() const { return num_; }
  void numeric::print(std::ostream &os) const { os << num_; }
  void numeric::set_numeric(double d) { num_=d; }
 


  void container::add_item(math_item::Ptr item) { subitems_.push_back(item); }
  void container::replace_subitem(math_item::Ptr what, 
				  math_item::Ptr replace_with) {
    for (std::vector<math_item::Ptr>::iterator iter=subitems_.begin();
	 iter!=subitems_.end(); ++iter) {
      if (boost::dynamic_pointer_cast<container>(*iter).get()==0) {
	if (*iter==what)
	  (*iter)=replace_with;
      } else {
	boost::dynamic_pointer_cast<container>(*iter)->replace_subitem(what, replace_with);
      }
    }
  }

  double add::eval() const { 
      double tmp=0.0; 
      for (std::vector<math_item::Ptr>::const_iterator iter=subitems_.begin(); iter!=subitems_.end(); ++iter) 
	tmp+=(*iter)->eval();
      return tmp;
    };
  void add::print(std::ostream &os) const { 
      os << " ( ";
      for (std::vector<math_item::Ptr>::const_iterator iter=subitems_.begin(); iter!=subitems_.end(); ++iter) 
	{  if (iter!=subitems_.begin()) os << " + "; 
	  (*iter)->print(os); }
      os << " ) ";
    }
 
 
  double mul::eval() const { 
      double tmp=1.0; 
      for (std::vector<math_item::Ptr>::const_iterator iter=subitems_.begin(); iter!=subitems_.end(); ++iter) 
	tmp*=(*iter)->eval();
      return tmp;
    }
  void mul::print(std::ostream &os) const { 
      os << " ( ";
      for (std::vector<math_item::Ptr>::const_iterator iter=subitems_.begin(); iter!=subitems_.end(); ++iter) 
	{ if (iter!=subitems_.begin())
	    os << " * "; 
	  (*iter)->print(os); }
      os << " ) ";

    }
 

  double pow::eval() const { 
      assert(subitems_.size()==2);
      return ::pow(subitems_[0]->eval(),subitems_[1]->eval());
    };
  void pow::print(std::ostream &os) const { 
      assert(subitems_.size()==2);
      os << *subitems_[0] << "^" << *subitems_[1];
  }


  double exp::eval() const { 
      assert(subitems_.size()==1);
      return ::exp(subitems_[0]->eval());
    };
  void exp::print(std::ostream &os) const { 
      assert(subitems_.size()==1);
      os << "exp(" << *subitems_[0] << ")";
  }
  
}

std::ostream  &operator<<(std::ostream &os, const MathTree::math_item &mi) { mi.print(os); return os;}
