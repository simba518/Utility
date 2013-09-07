#ifndef _IO_H_
#define _IO_H_

#include <fstream>
#include <assertext.h>
#include <Log.h>

namespace UTILITY{

  template<class T>
  bool loadText(std::ifstream &in,T *data,const int len,const std::string fname=""){

	assert(data);
	assert_ge(len,0);
	bool succ = false;
	if (in.is_open()){
	  succ = true;
	  for (int i = 0; i < len; ++i){
		if ( in.eof() || !(in >> (*data[i])) ){
		  ERROR_LOG("failed to read data,"<<fname<<",i="<<i);
		  succ = false;
		  break;
		}
	  }
	}
	return succ;
  }

  template<class T>
  bool writeText(std::ofstream &out,const T *data,const int len,const std::string space="\t",const std::string fname=""){

	bool succ = false;
	if (out.is_open() && data != NULL){
	  succ = true;
	  for ( int i = 0; i < len ; i++){
		if ( !(out <<data[i]<<space) ) {
		  ERROR_LOG("failed to save data,"<<fname<<",i="<<i);
		  succ = false;
		  break;
		}
	  }
	}
	return succ;
  }

  template<class T>
  bool loadBinary(std::ifstream &in,T *data,const int len,const std::string fname=""){
		
	bool succ = true;
	const int nRead = in.read((char*)(data), sizeof(T)*len).gcount();
	if (nRead != len*(int)sizeof(T)){
	  ERROR_LOG("failed to read data,"<<fname);
	  succ = false;
	}
	return succ;
  }

  template<class T>
  bool writeBinary(std::ofstream &out,const T *data,const int len,const std::string fname=""){

	bool succ = true;
	out.write((char*)(data),sizeof(T)*len);
	if (out.fail()){
	  ERROR_LOG("failed to save data,"<<fname);
	  succ = false;
	}
	return succ;
  }

}

#endif /* _IO_H_ */
