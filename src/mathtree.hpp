///////////////////////////////////////////////////////////////////////////////
//
// Prototype of the MathTree class to represent the formulas of the model
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

#ifndef MATHTREE_HPP
#define MATHTREE_HPP 

#include <ostream>
#include <vector>
#include <boost/shared_ptr.hpp>


// These are helper classes to evaluate the model quickly
// They will be used to construct a tree that represents a formula

namespace MathTree {

  class math_item {
  public:
    virtual ~math_item() {}
    virtual double eval() const;
    virtual void print(std::ostream &os) const;
    typedef boost::shared_ptr<math_item> Ptr;
  };
  
  class parameter : public math_item{
  public:
    parameter() ;
    virtual double eval() const ;
    virtual void print(std::ostream &os) const ; 
    void set_parameter(boost::shared_ptr<double> d);
    void set_parameter(double d);
    void operator=(double d);
    
    typedef boost::shared_ptr<parameter> Ptr;
  private:
    int parameter_index_;
    boost::shared_ptr<double> num_;
  };
  
  
  class numeric : public math_item{
  public:
    numeric()  {}
    numeric(double d) : num_(d) {}
    virtual double eval() const;
    virtual void print(std::ostream &os) const;
    void set_numeric(double d);
    typedef boost::shared_ptr<numeric> Ptr;
  private:
    double num_;
  };


  class container : public math_item {
  public:
    void add_item(math_item::Ptr item) ;
    void add_item(math_item  * item) { add_item(math_item::Ptr(item)); }
    void replace_subitem(math_item::Ptr what, math_item::Ptr replace_with);
    typedef boost::shared_ptr<container> Ptr;
  protected:
    std::vector<math_item::Ptr> subitems_;

  };
  
  class add : public container {
  public:
    virtual double eval() const;
    virtual void print(std::ostream &os) const;
  };

 
  class mul : public container {
  public:
    virtual double eval() const;
    virtual void print(std::ostream &os) const;

  };  

  class pow : public container {
  public:
    virtual double eval() const;
    virtual void print(std::ostream &os) const;
  };  

  class exp : public container {
  public:
    virtual double eval() const;
    virtual void print(std::ostream &os) const;
  };  


};


std::ostream  &operator<<(std::ostream &os, const MathTree::math_item &mi);


#endif // MATHTREE_HPP
