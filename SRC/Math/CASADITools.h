#ifndef _CASADITOOLS_H_
#define _CASADITOOLS_H_

#include <iostream>
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <symbolic/sx/sx.hpp>
#include <symbolic/fx/sx_function.hpp>
#include <symbolic/matrix/matrix_tools.hpp>
#include <assertext.h>

namespace CASADI{

  typedef std::vector<CasADi::SX> VSX;
  typedef std::vector<CasADi::SXMatrix> VMatSX;
  typedef std::vector<CasADi::SX>::iterator vsx_it;

  inline void makeSymbolic(vsx_it begin,vsx_it end,const std::string name){
	int counter = 0;
	for (vsx_it i = begin; i != end; ++i){
	  (*i)=CasADi::SX(name+std::string("_")+boost::lexical_cast<std::string>(counter));
	  counter ++;
	}
  }
	
  inline VSX makeSymbolic(const int len, const std::string name){
	VSX v(len);
	makeSymbolic(v.begin(), v.end(), name);
	return v;
  }

  inline void convert(const VSX &in, const int num_sub_vec, VMatSX &out){
	const int r = in.size()/num_sub_vec;
	assert_eq ((int)in.size() , r*num_sub_vec);
	out.resize(num_sub_vec);
	for (int i = 0; i < num_sub_vec; ++i){
	  out[i] = VSX(in.begin()+i*r,in.begin()+(i+1)*r);
	}
  }
  inline VMatSX Sx2Mat(const VSX &in, const int num_sub_vec){
	VMatSX out;
	convert (in, num_sub_vec, out);
	return out;
  }

  inline void convert(const VMatSX &in, VSX &out){
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
  inline VSX Mat2Sx(const VMatSX &in){
	VSX out;
	convert (in, out);
	return out;
  }

  inline CasADi::SXMatrix convert(const VSX &in,const int cols=1,const bool columnMajor=true){

	assert_ge(cols,1);
	assert_ge(in.size()%cols,0);
	assert_ge(in.size(),1);
	CasADi::SXMatrix out(in.size()/cols,cols);
	int i = 0;
	if(columnMajor){
	  for (int c = 0; c < out.size2(); ++c){
		for (int r = 0; r < out.size1(); ++r){
		  out.elem(r,c) = in[i];
		  i ++;
		}
	  }
	}else{
	  for (int r = 0; r < out.size1(); ++r){
		for (int c = 0; c < out.size2(); ++c){
		  out.elem(r,c) = in[i];
		  i ++;
		}
	  }	  
	}
	return out;
  }

  inline VSX connect(const VSX &v1,const VSX &v2){
	VSX v = v1;
	v.reserve(v1.size()+v2.size());
	for (size_t i = 0; i < v2.size(); ++i){
	  v.push_back(v2[i]);
	}
	return v;
  }

  template<class VEC>
  inline CasADi::SXMatrix makeEyeMatrix(const VEC &v){
	CasADi::SXMatrix out(v.size(),v.size());
	for (size_t i = 0; i < v.size(); ++i){
	  out(i,i) = v[i];
	}
	return out;
  }

  template<class T>
  inline void evaluate(CasADi::SXFunction &fun,const Eigen::Matrix<T,-1,1> &x,Eigen::Matrix<T,-1,1> &rlst){
	fun.setInput(&x[0]);
	fun.evaluate();
	const CasADi::DMatrix &out = fun.output();
	rlst.resize(out.size1()*out.size2());
	assert_gt (rlst.size(),0);
	out.get (&rlst[0], CasADi::DENSE);//@bug when T==float?
  }

  template<class T> 
  inline void evaluate(CasADi::SXFunction &fun,const Eigen::Matrix<T,-1,1> &x,Eigen::Matrix<T,-1,-1> &rlst){
	fun.setInput(&x[0]);
	fun.evaluate();
	const CasADi::DMatrix &out = fun.output();
	rlst.resize(out.size1(),out.size2());
	for (int i = 0; i < out.size1(); ++i){
	  for (int j = 0; j < out.size2(); ++j){
		rlst(i,j) = out(i,j).toScalar();
	  }
	}
  }

  template<class T> 
  inline void convert(const Eigen::Matrix<T,-1,-1> &M,CasADi::SXMatrix &SM){

	SM.makeDense(M.rows(), M.cols(), 0.0f);
	for (int i = 0; i < M.rows(); ++i){
	  for (int j = 0; j < M.cols(); ++j){
		SM.elem(i,j) = M(i,j);
	  }
	}
  }

  template<class T> 
  inline void convert(const Eigen::SparseMatrix<T> &M, CasADi::SXMatrix &SM){

  	std::vector<int> cols;
  	std::vector<int> rows;
  	std::vector<T> values;
  	cols.reserve(M.nonZeros());
  	rows.reserve(M.nonZeros());
  	values.reserve(M.nonZeros());
  	for (int k=0; k<M.outerSize(); ++k) {
  	  for (typename Eigen::SparseMatrix<T>::InnerIterator it(M,k); it; ++it) {
  		cols.push_back(it.col());
  		rows.push_back(it.row());
  		values.push_back(it.value());
  	  }
  	}
  	SM = CasADi::DMatrix::sparse( rows, cols, values, M.rows(), M.cols() );
  }

  template<class T> 
  inline void convert(const CasADi::SXMatrix &SM, Eigen::SparseMatrix<T> &M){

  	std::vector<int> rows,cols;
  	SM.sparsity().getSparsity (rows, cols);
  	const std::vector<CasADi::SX> &data = SM.data ();

  	typedef Eigen::Triplet<T> E_Triplet;
  	std::vector<E_Triplet> SM_triplet;
  	const int nz = (int)data.size();
  	SM_triplet.reserve(nz);
  	for (int i = 0; i < nz; ++i){
  	  SM_triplet.push_back( E_Triplet(rows[i],cols[i],data[i].getValue()) );
  	}
  	M.resize (SM.size1(),SM.size2());
  	M.reserve (nz);
  	M.setFromTriplets( SM_triplet.begin(), SM_triplet.end() );
  }

  template<class T> 
  inline void convert(const std::vector<Eigen::Matrix<T,-1,-1> >&M,std::vector<CasADi::SXMatrix> &SM){
	SM.reserve(M.size());
	for (size_t i = 0; i < M.size(); ++i){
	  CasADi::SXMatrix s;
	  convert(M[i],s);
	  SM.push_back(s);
	}
  }

  template<class T> 
  inline void convert(const Eigen::Matrix<T,-1,1> &V, CasADi::SXMatrix &SV){
	SV.makeDense(V.size(),1, 0.0f);
	for (int i = 0; i < V.size(); ++i){
	  SV.elem(i) = V[i];
	}
  }

  template<class T>
  inline void convert(const std::vector<Eigen::Matrix<T,-1,1> >&V,std::vector<CasADi::SXMatrix> &SV){
	SV.reserve(V.size());
	for (size_t i = 0; i < V.size(); ++i){
	  CasADi::SXMatrix s;
	  convert(V[i],s);
	  SV.push_back(s);
	}
  }

  template<class T> 
  inline void convert(const CasADi::SXMatrix &SV, Eigen::Matrix<T,-1,1> &V){

	assert_eq(SV.size2(),1);
	V.resize(SV.size1());
	for (int i = 0; i < SV.size1(); ++i){
	  V[i] = SV.elem(i,0).getValue();
	}
  }

  template<class T> 
  inline void convert(const CasADi::SXMatrix &SM, Eigen::Matrix<T,-1,-1> &M){

	M.resize(SM.size1(),SM.size2());
	for (int i = 0; i < SM.size1(); ++i){
	  for (int j = 0; j < SM.size2(); ++j){
		M(i,j) = SM.elem(i,j).getValue();
		assert_eq_ext(M(i,j),M(i,j),"i="<<i<<",j="<<j<<",SM(i,j)="<<SM.elem(i,j));
	  }
	}
  }

  template<class T> 
  inline CasADi::SXMatrix convert(const Eigen::Matrix<T,-1,1> &V){
	CasADi::SXMatrix SV;
	SV.makeDense(V.size(),1,0.0f);
	for (int i = 0; i < V.size(); ++i){
	  SV.elem(i) = V[i];
	}
	return SV;
  }

  template<class T> 
  inline CasADi::SXMatrix convert(const Eigen::Matrix<T,-1,-1> &M){
	CasADi::SXMatrix SM;
	convert(M,SM);
	return SM;
  }

  template<class T> 
  inline Eigen::Matrix<T,-1,1> convert2Vec(const CasADi::SXMatrix &SV){
	assert_eq(SV.size2(),1);
	Eigen::Matrix<T,-1,1> V;
	convert(SV,V);
	return V;
  }

  template<class T> 
  inline Eigen::Matrix<T,-1,-1> convert(const CasADi::SXMatrix &SM){
	Eigen::Matrix<T,-1,-1> M(SM.size1(),SM.size2());
	convert(SM,M);
	return M;
  }
}

#endif /* _CASADITOOLS_H_ */
