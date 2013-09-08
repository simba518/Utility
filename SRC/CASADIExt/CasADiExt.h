#ifndef _CASADIEXT_H_
#define _CASADIEXT_H_

#include <iostream>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <eigen3/Eigen/Dense>
#include <symbolic/sx/sx.hpp>
#include <symbolic/fx/sx_function.hpp>
#include <symbolic/matrix/matrix_tools.hpp>
#include <assertext.h>

/**
 * extension for the CasADi library.
 * 
 */
namespace LSW_CASADI{

  typedef std::vector<CasADi::SX> VSX;
  typedef std::vector<CasADi::SXMatrix> VMatSX;
  typedef std::vector<CasADi::SX>::iterator vsx_it;

  class CasADiExt{
	
  public:
	static void make_symbolic(vsx_it begin,vsx_it end,const std::string name){

	  int counter = 0;
	  for (vsx_it i = begin; i != end; ++i){
		(*i)=SX(name+std::string("_")+boost::lexical_cast<std::string>(counter));
		counter ++;
	  }
	}
	
	static VSX make_symbolic(const int len, const std::string name){
	  
	  VSX v(len);
	  make_symbolic(v.begin(), v.end(), name);
	  return v;
	}

	static void Sx2Mat(const VSX &in, const int num_sub_vec, VMatSX &out){

	  const int r = in.size()/num_sub_vec;
	  assert_eq ((int)in.size() , r*num_sub_vec);

	  out.resize(num_sub_vec);
	  for (int i = 0; i < num_sub_vec; ++i){
		out[i] = VSX(in.begin()+i*r,in.begin()+(i+1)*r);
	  }
	}

	static VMatSX Sx2Mat(const VSX &in, const int num_sub_vec){
	  
	  VMatSX out;
	  Sx2Mat (in, num_sub_vec, out);
	  return out;
	}

	static void Mat2Sx(const VMatSX &in, VSX &out){

	  int len = 0;
	  for (size_t i = 0; i < in.size(); ++i){
		len += in[i].size();
	  }

	  out.reserve(len);

	  for (size_t i = 0; i < in.size(); ++i){
		for (int r = 0; r < in[i].size1(); ++r){
		  for (int c = 0; c < in[i].size2(); ++c){
			out.push_back(in[i].elem(r,c));
		  }
		}
	  }
	}

	static VSX Mat2Sx(const VMatSX &in){

	  VSX out;
	  Mat2Sx (in, out);
	  return out;
	}

	static CasADi::SXMatrix convert(const VSX &in){
	  CasADi::SXMatrix out(in.size(),1);
	  for (size_t i = 0; i < in.size(); ++i){
		out(i) = in[i];
	  }
	  return out;
	}

	static VSX append(const VSX &v1, const VSX &v2){

	  VSX v = v1;
	  v.reserve(v1.size()+v2.size());
	  for (size_t i = 0; i < v2.size(); ++i){
		v.push_back(v2[i]);
	  }
	  return v;
	}

	static const Eigen::VectorXd &getResult (CasADi::SXFunction fun, 
											 const Eigen::VectorXd &x, 
											 Eigen::VectorXd &rlst){
	  fun.setInput(&x[0]);
	  fun.evaluate();
	  const CasADi::DMatrix &out = fun.output();
	  rlst.resize(out.size1()*out.size2());
	  assert_gt (rlst.size(),0);
	  out.get (&rlst[0], CasADi::DENSE);
	  return rlst;
	}

	static const Eigen::MatrixXd &getResult (CasADi::SXFunction fun, 
											 const Eigen::VectorXd &x, 
											 Eigen::MatrixXd &rlst){
	  fun.setInput(&x[0]);
	  fun.evaluate();
	  const CasADi::DMatrix &out = fun.output();
	  rlst.resize(out.size1(),out.size2());
	  for (int i = 0; i < out.size1(); ++i){
		for (int j = 0; j < out.size2(); ++j){
		  rlst(i,j) = out(i,j).toScalar();
		}
	  }
	  return rlst;
	}

  };  
}//end of namespace

#endif /*_CASADIEXT_H_*/
