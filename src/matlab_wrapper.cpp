
extern "C" {
#include <math.h>
#include "mex.h"

}

#include "fitmodel_in_CPP.hpp"
#include "helper_types.hpp"   

template<typename T>
double_matrix convert_matlab_to_double_matrix(T *tmp, unsigned M, unsigned N) {
  double_matrix mat(boost::extents[M][N]);
  for (unsigned i=0;i<N; i++) {
    for (unsigned j=0;j<M; j++) {
      mat[j][i]=*tmp++;
    }
  }
  return mat;
}

template<typename T>
int_matrix convert_matlab_to_int_matrix(T *tmp, unsigned M, unsigned N) {
  int_matrix mat(boost::extents[M][N]);
  for (unsigned i=0;i<N; i++) {
    for (unsigned j=0;j<M; j++) {
      mat[j][i]=*tmp++;
    }
  }
  return mat;
}



void cpp_matlab_call_wrapper(
			     int           nlhs,
			     mxArray       *plhs[],
			     int           nrhs,
			     const mxArray *prhs[] )
{
  if ((nrhs!=10) || (nlhs!=4)) {
    mexErrMsgTxt(
		 "call as follows:\n "
		 "[bestfit,bestresid,prediction,dependent] "
		 " = fitmodel_in_CPP(adm, unstim_data, stim_data, error, stim_nodes, inhib_nodes, stimuli, inhibitor, measured_nodes, scale) ");
      
  }
    
  if ( !mxIsLogical(prhs[0] ) )
    mexErrMsgTxt("first input must be logical");
    
  // Test arguments type 
  if ( !mxIsDouble(prhs[1]) )
    mexErrMsgTxt("second input must be double");
    
  if ( !mxIsDouble(prhs[2]) )
    mexErrMsgTxt("third input must be double");
    
  if ( !mxIsDouble(prhs[3]) )
    mexErrMsgTxt("forth input must be double");
    
  if ( !mxIsLogical(prhs[4] ) )
    mexErrMsgTxt("fifth input must be logical");

  if ( !mxIsLogical(prhs[5] ) )
    mexErrMsgTxt("sixth input must be logical");

  if ( !mxIsLogical(prhs[6] ) )
    mexErrMsgTxt("seventh input must be logical");

  if ( !mxIsLogical(prhs[7] ) )
    mexErrMsgTxt("eigth input must be logical");
 
  if ( !mxIsLogical(prhs[8] ) )
    mexErrMsgTxt("ninth input must be logical");

  if ( !mxIsDouble(prhs[9] ) )
    mexErrMsgTxt("tenth input must be double");
  
  // Read in input data in their appropriate format
  int_matrix adm=convert_matlab_to_int_matrix(mxGetLogicals(prhs[0]),mxGetM(prhs[0]),mxGetN(prhs[0]));

  Data data;
  data.unstim_data.resize(boost::extents[mxGetM(prhs[1])][mxGetN(prhs[1])]);
  data.unstim_data=convert_matlab_to_double_matrix(mxGetPr(prhs[1]),mxGetM(prhs[1]),mxGetN(prhs[1]));

  data.stim_data.resize(boost::extents[mxGetM(prhs[2])][mxGetN(prhs[2])]);
  data.stim_data=convert_matlab_to_double_matrix(mxGetPr(prhs[2]),mxGetM(prhs[2]),mxGetN(prhs[2]));
 
  data.error.resize(boost::extents[mxGetM(prhs[3])][mxGetN(prhs[3])]);
  data.error=convert_matlab_to_double_matrix(mxGetPr(prhs[3]),mxGetM(prhs[3]),mxGetN(prhs[3])); 
  
  int_matrix stim_nodestmp=convert_matlab_to_int_matrix(mxGetLogicals(prhs[4]),mxGetM(prhs[4]),mxGetN(prhs[4]));
  for(unsigned int i=0;i<stim_nodestmp.shape()[1];i++){if(stim_nodestmp[0][i]==1) {data.stim_nodes.push_back(i);}}
  
  int_matrix inhib_nodestmp=convert_matlab_to_int_matrix(mxGetLogicals(prhs[5]),mxGetM(prhs[5]),mxGetN(prhs[5]));
  for(unsigned int i=0;i<inhib_nodestmp.shape()[1];i++){if(inhib_nodestmp[0][i]==1) {data.inhib_nodes.push_back(i);}} 
  
  data.stimuli.resize(boost::extents[mxGetM(prhs[6])][mxGetN(prhs[6])]);
  data.stimuli=convert_matlab_to_int_matrix(mxGetLogicals(prhs[6]),mxGetM(prhs[6]),mxGetN(prhs[6]));
  
  data.inhibitor.resize(boost::extents[mxGetM(prhs[7])][mxGetN(prhs[7])]);
  data.inhibitor=convert_matlab_to_int_matrix(mxGetLogicals(prhs[7]),mxGetM(prhs[7]),mxGetN(prhs[7]));
  
  int_matrix measured_nodestmp=convert_matlab_to_int_matrix(mxGetLogicals(prhs[8]),mxGetM(prhs[8]),mxGetN(prhs[8]));
  for(unsigned int i=0;i<measured_nodestmp.shape()[1];i++){if(measured_nodestmp[0][i]==1) {data.measured_nodes.push_back(i);}}
  
  data.scale.resize(boost::extents[mxGetM(prhs[9])][mxGetN(prhs[9])]);
  data.scale=convert_matlab_to_double_matrix(mxGetPr(prhs[9]),mxGetM(prhs[9]),mxGetN(prhs[9]));
  
  //Sanity checks   
  if ( ! data.data_consistent() ) mexErrMsgTxt("data in shitty format");

  std::vector<double> bestfit;
  double bestresid;
  double_matrix prediction;
  std::vector<double> dependent;
    
  fitmodel(bestfit,
	   bestresid,
	   prediction,
	   dependent,
	   adm, 
	   data);

 double *tmp;
  int bestfit_rows=bestfit.size();
  plhs[0] = mxCreateDoubleMatrix(bestfit_rows,1,mxREAL);  
  tmp=mxGetPr(plhs[0]);
  for (unsigned i=0;i<mxGetM(plhs[0]); i++) {
    (*tmp++)=bestfit[i];
  }

  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  tmp = mxGetPr(plhs[1]);  
  *tmp=bestresid;
     
  int prediction_rows=prediction.shape()[0];
  int prediction_cols=prediction.shape()[1];
  plhs[2] = mxCreateDoubleMatrix(prediction_rows,prediction_cols,mxREAL);
  tmp=mxGetPr(plhs[2]);
  for (unsigned j=0;j<mxGetN(plhs[2]); j++) {
    for (unsigned i=0;i<mxGetM(plhs[2]); i++) {
      (*tmp++)=prediction[i][j];
    }
  }
  
  int dependent_rows=dependent.size();
  plhs[3]=mxCreateDoubleMatrix(dependent_rows,1,mxREAL);
  tmp=mxGetPr(plhs[3]);
  for(unsigned i=0;i<mxGetM(plhs[3]);i++){
      (*tmp++)=dependent[i];
    }
  
}




extern "C" {

  void mexFunction(
                   int           nlhs,
                   mxArray       *plhs[],
                   int           nrhs,
                   const mxArray *prhs[]
                   )
  {
    cpp_matlab_call_wrapper(nlhs, plhs, nrhs, prhs);
    return;
  }
  
}


