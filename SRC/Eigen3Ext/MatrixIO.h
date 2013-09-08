#ifndef _MATRIXIO_H_
#define _MATRIXIO_H_

#include <eigen3/Eigen/Dense>
#include <IO.h>

namespace EIGEN3EXT{

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
	return UTILITY::loadVec(fname,data,io_type);
  }

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,-1,1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
	return UTILITY::loadVec(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeVec(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,-1,1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeVec(fname,data,io_type);
  }

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,-1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
    return UTILITY::loadMat(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,-1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeMat(fname,data,io_type,space);
  }

  template<class T>
  inline bool load(const std::string fname,std::vector<Eigen::Matrix<T,-1,1> > &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){

	bool succ = false;
	std::ifstream in;	
	if( openInputFile(in,fname,io_type) ){
	  int Total,subLen;
	  if(load(in,&subLen,1,io_type,fname) && load(in,&Total,1,io_type,fname) ){
		assert_gt(Total,0);
		assert_gt(subLen,0);
		data.resize(Total);
		for (size_t i = 0; i < data.size(); ++i){
		  data[i].resize(subLen);
		  succ = load(in,&(data[i][0]),subLen,io_type,fname);
		  if(!succ){
			break;
		  }
		}
	  }
	}
	in.close();
	return succ;
  }

  template<class T>
  inline bool write(const std::string fname,const std::vector<Eigen::Matrix<T,-1,1> > &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	
	bool succ = false;
	std::ofstream out;
	if( openOutputFile(out,fname,io_type) ){
	  const int Total = (int)data.size();
	  assert_gt(Total,0);
	  const int subLen = (int)data[0].size();
	  assert_gt(subLen,0);
	  if(write(out,&subLen,1,io_type,fname,space) && write(out,&Total,1,io_type,fname,space)){
		for (int i = 0; i < Total; ++i){
		  succ = write(out,&data[i][0],subLen,io_type,space,fname);
		  if(!succ)
			break;
		}
	  }
	}
	return succ;
  }

}

#endif /* _MATRIXIO_H_ */
