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
