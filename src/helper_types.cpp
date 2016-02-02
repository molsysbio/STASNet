#include "helper_types.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/normal.hpp>
#include "model.hpp"
#include <fstream>

extern int verbosity;
extern bool debug;

// Performs a Latin Hypercube Sampling on ]0, 1[
// Each element of the return vector is a sample
std::vector< std::vector<double> > LHSampling (const int nb_samples, const int sample_size, const int decimals) {
    // Remember which "band" is filled for each dimension
    std::vector< std::vector<bool> > dimension_filled;
    dimension_filled.resize(sample_size);
    for (int j=0 ; j < sample_size ; j++) {
        dimension_filled[j] = std::vector<bool>(nb_samples, false);
    }
    // Will contain the samples
    std::vector< std::vector<double> > sampling;
    sampling.resize(nb_samples);
    for (int i=0 ; i < nb_samples ; i++) {
        sampling[i] = std::vector<double>(sample_size, 0);
    }

    int part = 0;
    double precision = pow(10, decimals);
    for (int i=0 ; i < nb_samples ; i++) {
        for (int j=0 ; j < sample_size ; j++) {
            part = rand() % nb_samples;
            // Change the "band" if the selected one has already been sampled
            // Select randomly if there is "enough" available choices or simply shift on the ring
            if (i < (double)(0.5 * nb_samples)) {
                while (dimension_filled[j][part]) {
                    part = rand() % nb_samples;
                }
            } else {
                while (dimension_filled[j][part]) {
                    part = ++part % nb_samples;
                }
            }
            dimension_filled[j][part] = true;
            sampling[i][j] = (double)(part + uniform_sampling(precision)) / nb_samples;
        }
    }

    return sampling;
}

// Perform LHS on a normal distribution
std::vector< std::vector<double> > normalLHS (const int nb_samples, const int sample_size, const int sd, const double decimals) {
    std::vector< std::vector<double> > sampling = LHSampling(nb_samples, sample_size, decimals);
    for (size_t i=0 ; i < sample_size ; i++) {
        for (size_t j=0 ; j < sample_size ; j++) {
            sampling[i][j] = boost::math::quantile(boost::math::normal(0, sd), sampling[i][j]);
        }
    }
    return sampling;
}

// Provides a random number in ]0, 1[
double uniform_sampling(double precision) {
    if (precision <= 1) {
      precision = 2;
    }
    return ((double)(rand() % ((int)precision-1)) + 1.0) / precision;
}

template<typename T> 
void convert_negatives_into_NaN(boost::multi_array< T, 2 > &a) {
  for (size_t i=0; i<a.shape()[0]; i++) {
    for (size_t j=0; j<a.shape()[1]; j++) {
      if (a[i][j]<0) { 
	    a[i][j]=std::numeric_limits<T>::quiet_NaN();
      }
    }
  }
}

// Helper function to output a double matrix
std::ostream &operator<<(std::ostream &os, const double_matrix &a) {

  for (size_t i=0; i<a.shape()[0]; i++) {
    for (size_t j=0; j<a.shape()[1]; j++) {
      os << a[i][j] << "\t";
    }
    os << std::endl;
  }
  return os;
}

// Helper function to output an int vector
std::ostream &operator<<(std::ostream &os, const std::vector<int> &a) {

  for (size_t i=0; i<a.size(); i++) {
    os << a[i] << std::endl;
    }
    os << std::endl;
  return os;
}
// Helper function to output an double vector
std::ostream &operator<<(std::ostream &os, const std::vector<double> &a) {

  for (size_t i=0; i<a.size(); i++) {
    os << a[i] << std::endl;
    }
    os << std::endl;
  return os;
}

// Helper function to output a int matrix
std::ostream &operator<<(std::ostream &os, const int_matrix &a) {

  for (size_t i=0; i<a.shape()[0]; i++) {
    for (size_t j=0; j<a.shape()[1]; j++) {
      os << a[i][j] << "\t";
    }
    os << std::endl;
  }
  return os;
}

// Helper function to read in parametermap
std::istream &operator>>(std::istream &is, std::map<std::string,double> &p) {
  while (is.good()) {
    std::string paramname;
    is >> paramname;
    if (is.good()) {
      std::string tmp;

      is >> tmp;
      if (is.good()) {
	try { p[paramname]=boost::lexical_cast<double>(tmp);
	}
	catch (boost::bad_lexical_cast &) { 
	  std::cerr << " problem reading param file. Exiting " << std::endl;
	  exit(-1);
	}
      }
    }
  }
  return is;
}

// Helper function to output an MathTree equation matrix
std::ostream &operator<<(std::ostream &os, const equation_matrix &a) {

  for (size_t i=0; i<a.shape()[0]; i++) {
    for (size_t j=0; j<a.shape()[1]; j++) {
      os << *a[i][j] << "\t";
    }
    os << std::endl;
  }
  return os;
}

void convert_int_to_rational_matrix(const int_matrix &i, rational_matrix &r) 
{
  r.resize(boost::extents[i.shape()[0]][i.shape()[1]]);
  for (size_t x=0; x< i.shape()[0]; ++x) {
    for (size_t y=0; y< i.shape()[1]; ++y) {
      r[x][y]=boost::rational<int>(i[x][y]);
    }
  }
}

void convert_rational_to_double_matrix(const rational_matrix &r, double_matrix &d) {
  d.resize(boost::extents[r.shape()[0]][r.shape()[1]]);
  for (size_t x=0; x< r.shape()[0]; ++x) {
    for (size_t y=0; y< r.shape()[1]; ++y) {
      d[x][y]=boost::rational_cast<double>(r[x][y]);
    }
  }
}

bool ExperimentalDesign::read_from_stream(std::istream &is) {
  std::string label;
  is >> label;
  if (label != "inhib_nodes") return false;
  if (!read_array(is,inhib_nodes)) { std::cerr << "problem reading inhib_nodes" << std::endl; return false; }

  is >> label;
  if (label != "stim_nodes") return false;
  if (!read_array(is,stim_nodes))  { std::cerr << "problem reading stim_nodes" << std::endl; return false; }

  is >> label;
  if (label != "measured_nodes") return false;
  if (!read_array(is,measured_nodes)) { std::cerr << "problem reading measured_nodes" << std::endl; return false; }

  is >> label;
  if (label != "stimuli") return false;
  if (!read_array(is,stimuli)) { std::cerr << "problem reading stimuli" << std::endl; return false; }

  is >> label;
  if (label != "inhibitor") return false;
  if (!read_array(is,inhibitor)) { std::cerr << "problem reading inhibitor" << std::endl; return false; }

  is >> label;
  if (label != "basal_activity") return false;
  if (!read_array(is,basal_activity)) { std::cerr << "problem reading basal_activity" << std::endl; return false; }

  return true;
  
}

Data::Data() {
    dataVector = new double[1];
    dataVectorComputed = false;
}

bool Data::read_from_stream(std::istream &is) {
  std::string label;
  is >> label;
  if (label != "unstim_data") return false;
  if (!read_array(is,unstim_data)) { std::cerr << "problem reading unstim_data" << std::endl; return false; }

  is >> label;
  if (label != "stim_data") return false;
  if (!read_array(is,stim_data)) { std::cerr << "problem reading stim_data" << std::endl; return false; }
  //  convert_negatives_into_NaN(stim_data);

  is >> label;
  if (label != "error") return false;
  if (!read_array(is,error))  { std::cerr << "problem reading error" << std::endl; return false; }
  convert_negatives_into_NaN(error);

  is >> label;
  if (label != "scale") return false;
  if (!read_array(is,scale))  { std::cerr << "problem reading scale" << std::endl; return false; }
  
  return true;
  
}

bool Data::data_consistent(const ExperimentalDesign &expdesign) const {
  return ( (unstim_data.shape()[0]==stim_data.shape()[0]) &&
	   (unstim_data.shape()[1]==stim_data.shape()[1]) &&
	   (error.shape()[1]==stim_data.shape()[1]) &&
	   (error.shape()[0]==stim_data.shape()[0]) &&
	   (scale.shape()[1]==stim_data.shape()[1]) &&
	   (scale.shape()[0]==stim_data.shape()[0]) &&
	   (expdesign.stim_nodes.size()==expdesign.stimuli.shape()[1]) &&
	   (stim_data.shape()[0]==expdesign.stimuli.shape()[0]) &&
	   (expdesign.inhib_nodes.size()==expdesign.inhibitor.shape()[1]) &&
	   (expdesign.measured_nodes.size()==unstim_data.shape()[1]) &&
	   (stim_data.shape()[0]==expdesign.inhibitor.shape()[0]) 
/*(basal_activity.size()==names.size())*/);
}

ExperimentalDesign &ExperimentalDesign::operator=(const ExperimentalDesign &exper) {
  stim_nodes=exper.stim_nodes;
  inhib_nodes=exper.inhib_nodes;
  measured_nodes=exper.measured_nodes;
  copy_matrix(exper.stimuli,stimuli);
  copy_matrix(exper.inhibitor,inhibitor);
  basal_activity=exper.basal_activity;
  return *this;
}

Data &Data::operator=(const Data &data) {
  copy_matrix(data.unstim_data,unstim_data);
  copy_matrix(data.stim_data,stim_data);
  copy_matrix(data.error,error);
  copy_matrix(data.scale,scale);
  delete[] dataVector;
  dataVector = new double[1];
  dataVectorComputed = false;
  computeDataVector();
  return *this;
}

std::ostream &operator<<(std::ostream &os, const Data &d) {
  return os;
}

long unsigned int getSeed()
{
  std::ifstream rand("/dev/urandom");
  char tmp[sizeof(long unsigned int)];
  rand.read(tmp,sizeof(long unsigned int));
  rand.close();
  int* number = reinterpret_cast<int*>(tmp);
  return (*number);
}

void Data::computeDataVector() {
    // define the measurement value to compare with simulated values divided by the error
    if (!dataVectorComputed && stim_data.shape()[0] == error.shape()[0] && stim_data.shape()[1] == error.shape()[1]) {
        size_t rows=stim_data.shape()[0], cols=stim_data.shape()[1];
        nb_measurements = rows * cols;

        if (dataVector != NULL) {
            //delete[] dataVector; // TODO Fix to avoid memory leaks
        }
        dataVector = new double[rows * cols];
        for (unsigned int i=0; i<cols;i++) { 
            for (unsigned int j=0; j<rows;j++) {
                if (std::isnan(error[j][i]) || std::isnan(stim_data[j][i])) {
                    dataVector[i*stim_data.shape()[0]+j]=0;
                } else {
                    dataVector[i*stim_data.shape()[0]+j]=stim_data[j][i]/error[j][i];
                }
            }
        }
        dataVectorComputed = true;
    }
}

DataSet::DataSet() {
}

DataSet::~DataSet() {
}

void DataSet::addData(Data &data, bool doDataVectorComputation) {
    datas_.push_back(data);
    //std::cout << "rows=" << data.unstim_data.shape()[0] << ", cols=" << data.unstim_data.shape()[1] << std::endl;

    // rbind_matrix generates an 'memory corruption' on the second call, could not figure out why
    // As a consequence, the rbind matrix must be provided in R (or an override of computeDataVector is necessary)
    /*
    rbind_matrix(data.error, error);
    std::cout << "Finished error !!" << std::endl;
    std::cout << "Go for scale !!" << std::endl;
    rbind_matrix(data.scale, scale);
    std::cout << "Go for unstim_data !!" << std::endl;
    rbind_matrix(data.unstim_data, unstim_data);
    std::cout << "Finished unstim_data !!" << std::endl;
    rbind_matrix(data.stim_data, stim_data);
    std::cout << "Finished stim_data !!" << std::endl;
    */

    /*
    boost::multi_array_types::index_gen indices; // To generate the views
    typedef boost::multi_array_types::index_range index_range;

    std::vector<size_t> extendlist(2);
    const size_t* ushape = unstim_data.shape();
    extendlist[0] = ushape[0]+data.unstim_data.shape()[0]; extendlist[1] = ushape[1];
    unstim_data.resize(extendlist);
    unstim_data[ indices[index_range(ushape[0]+1, unstim_data.shape()[0])][index_range(0, ushape[1])] ] = data.unstim_data;

    ushape = stim_data.shape();
    extendlist[0] = ushape[0]+data.stim_data.shape()[0]; extendlist[1] = ushape[1];
    stim_data.resize(extendlist);
    stim_data[ indices[index_range(ushape[0]+1, stim_data.shape()[0])][index_range(0, ushape[1])] ] = data.stim_data;

    ushape = error.shape();
    extendlist[0] = ushape[0]+data.error.shape()[0]; extendlist[1] = ushape[1];
    error.resize(extendlist);
    error[ indices[index_range(ushape[0]+1, error.shape()[0])][index_range(0, ushape[1])] ] = data.error;

    ushape = scale.shape();
    extendlist[0] = ushape[0]+data.scale.shape()[0]; extendlist[1] = ushape[1];
    scale.resize(extendlist);
    scale[ indices[index_range(ushape[0]+1, scale.shape()[0])][index_range(0, ushape[1])] ] = data.scale;
    */

    if (doDataVectorComputation) {
        computeDataVector();
    }
}

void DataSet::addDataFromMatrices(double_matrix unstim_data, double_matrix stim_data, double_matrix error, double_matrix scale, bool doDataVectorComputation) {
    Data data;
    data.setUnstimData(unstim_data);
    data.setStimData(stim_data);
    data.setError(error);
    data.setScale(scale);
    
    addData(data, doDataVectorComputation);
}

bool DataSet::data_consistent(const ExperimentalDesign &expdesign) const {
    if (verbosity > 7) { std::cout << "Calling DataSet data_consistent" << std::endl; }
    if (datas_.size() > 0) {
        return(datas_[0].data_consistent(expdesign));
    }
    return(false);
}

