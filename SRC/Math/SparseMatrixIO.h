#ifndef _SPARSEMATRIXIO_H_
#define _SPARSEMATRIXIO_H_

#include <vector>
#include <string>
using std::string;
using std::vector;

#include <log4cplus/logger.h>
#include <loggingmacrosext.h>
using namespace log4cplus;

#include <assertext.h>
#include <eigen3/Eigen/Sparse>
#include <VectorIO.h>
#include "TextVectorIO.h"
#include "BinaryVectorIO.h"

/**
 * read or save sparse matrix using triple formats:
 * rows    : number of rows
 * columns : number of columns
 * nz      : number of nonzeros, int type.
 * row_ids : index of rows, int type.
 * col_ids : index of columns, int type.
 * values  : data values, T or float or int type.
 */
namespace LSW_UTILITY{

  template <typename T> 
  bool read(Eigen::SparseMatrix<T>&S,const string f_name,IO_TYPE io_type=BINARY){

	// read from file
	bool succ = false;
	vector<int> row_ids, col_ids;
	vector<T> values;
	int rows, cols;
	if (io_type == TEXT){
	  succ = readText(rows,cols,row_ids,col_ids,values,f_name);
	}else{
	  succ = readBinary(rows,cols,row_ids,col_ids,values,f_name);
	}

	// initialize the sparse matrix
	if (succ){

	  assert_ge(rows,0);
	  assert_ge(cols,0);
	  assert_eq(row_ids.size(),values.size());
	  assert_eq(col_ids.size(),values.size());

	  vector<Eigen::Triplet<T> > triplets;
	  triplets.reserve(values.size());
	  for (int i = 0; i < (int)values.size(); ++i){
		assert_in(row_ids[i],0,rows-1);
		assert_in(col_ids[i],0,cols-1);
		triplets.push_back( Eigen::Triplet<T>(row_ids[i],col_ids[i],values[i]) );
	  }
	  S.resize(rows,cols);
	  S.reserve(values.size());
	  S.setFromTriplets(triplets.begin(), triplets.end());
	}
	return succ;
  }

  template <typename T> 
  bool save(const Eigen::SparseMatrix<T>&S,const string f_name,IO_TYPE io_type = BINARY){

	bool succ = false;
	// the rows, cols and values
	const int nz = S.nonZeros();
	vector<int> row_ids(nz), col_ids(nz);
	vector<T> values(nz);
	int i = 0;
	for (int k=0; k<S.outerSize(); ++k){ 
	  for (typename Eigen::SparseMatrix<T>::InnerIterator it(S,k); it; ++it){
		row_ids[i] = it.row();
		col_ids[i] = it.col();
		values[i] = it.value();
		i ++;
	  }
	}

	// save to file
	const int rows = S.rows();
	const int cols = S.cols();
	if (io_type == TEXT){
	  succ = saveText(rows,cols,row_ids,col_ids,values,f_name);
	}else{
	  succ = saveBinary(rows,cols,row_ids,col_ids,values,f_name);
	}
	return succ;
  }

  template <typename T> 
  bool saveBinary(const int rows,const int cols,
				  const vector<int> &row_ids,const vector<int> &col_ids,
				  const vector<T> &values,const string &f_name){

	// open file
	ofstream out_f;
	out_f.open(f_name.c_str(),ios_base::out|ios_base::binary);
	bool succ = out_f.is_open();

	// write data
	if (succ){
	  
	  const int nz = values.size();
	  out_f.write((char*)(&nz),sizeof(int));    // write nz
	  out_f.write((char*)(&rows),sizeof(int));  // write rows
	  out_f.write((char*)(&cols),sizeof(int));  // write cols
	  succ = !(out_f.fail());
	  if (succ && nz > 0){

		assert_eq ((int)row_ids.size(), nz);
		assert_eq ((int)col_ids.size(), nz);
		out_f.write((char*)(&row_ids[0]), sizeof(int)*(row_ids.size()));
		out_f.write((char*)(&col_ids[0]), sizeof(int)*(col_ids.size()));
		out_f.write((char*)(&values[0]), sizeof(T)*(values.size()));
		succ = !(out_f.fail());
	  }
	  out_f.close(); // close file
	  if (!succ){
		Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
		LOG4CPLUS_ERROR(logger, "failed to save the sparse matrix to: "<<f_name);
	  }
	}else{
	  Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
	  LOG4CPLUS_ERROR(logger, "failed to open file: "<<f_name);
	}
	return succ; 
  }

  template <typename T> 
  bool saveText(const int rows,const int cols,
				  const vector<int> &row_ids,const vector<int> &col_ids,
				  const vector<T> &values,const string &f_name){
	/// @todo undefined function
	Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
	const string msg=string("undefined function is called: ")+string(__FUNCTION__);
	LOG4CPLUS_WARN(logger, msg);
	return false;
  }

  template <typename T> 
  bool readBinary(int &rows,int &cols,
				  vector<int> &row_ids,vector<int> &col_ids,
				  vector<T> &values,const string &f_name){

	
	// open file
	ifstream in_f;
	in_f.open(f_name.c_str(),ios_base::in|ios_base::binary);
	bool succ = in_f.is_open();

	// read data
	if (succ){
	  
	  int nz = 0;
	  in_f.read((char*)(&nz),sizeof(int));    // read number of nonzeros
	  in_f.read((char*)(&rows),sizeof(int));  // read rows
	  in_f.read((char*)(&cols),sizeof(int));  // read cols
	  succ = !(in_f.fail());
	  
	  assert_ge(rows,0);
	  assert_ge(cols,0);
	  assert_ge(rows*cols,nz);
	  assert_ge(nz,0);
	  
	  row_ids.resize(nz);
	  col_ids.resize(nz);
	  values.resize(nz);
	  
	  if (succ && nz > 0){

		in_f.read((char*)(&row_ids[0]), sizeof(int)*(row_ids.size()));
		in_f.read((char*)(&col_ids[0]), sizeof(int)*(col_ids.size()));
		in_f.read((char*)(&values[0]), sizeof(T)*(values.size()));
		succ = !(in_f.fail());
	  }
	  in_f.close(); // close file
	  if (!succ){
		Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
		LOG4CPLUS_ERROR(logger, "failed to save the sparse matrix to: "<<f_name);
	  }
	}else{
	  Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
	  LOG4CPLUS_ERROR(logger, "failed to open file: "<<f_name);
	}
	return succ; 
  }

  template <typename T> 
  bool readText(int &rows,int &cols,
				  vector<int> &row_ids,vector<int> &col_ids,
				  vector<T> &values,const string &f_name){
	/// @todo undefined function
	Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("Utility"));
	const string msg=string("undefined function is called: ")+string(__FUNCTION__);
	LOG4CPLUS_WARN(logger, msg);
	return false;
  }
  
}

#endif /* _SPARSEMATRIXIO_H_ */
