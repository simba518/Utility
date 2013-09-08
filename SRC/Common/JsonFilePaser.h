#ifndef _JSONFILEPASER_H_
#define _JSONFILEPASER_H_

#include <stdlib.h>
#include <vector>
#include <set>
#include <fstream>
#include <boost/filesystem.hpp>
#include <eigen3/Eigen/Dense>
#include <Log.h>

namespace UTILITY{
  
  /**
   * @class JsonFilePaser base interface for reading the initfile.
   * 
   */
  class JsonFilePaser{
	
  public:
	JsonFilePaser(bool print = false){
	  _printData = print;
	}

	~JsonFilePaser(){
	  JsonFilePaser::close();
	}

	virtual bool open(const std::string filename){
	  _initFilename = filename;
	  JsonFilePaser::close();
	  _file.open(filename.c_str());
	  const bool succ = _file.is_open();
	  ERROR_LOG_COND("failed to open init file: "<<filename,succ);
	  return succ;
	}
	virtual void close(){ if(_file.is_open()) _file.close();}

	std::string getFileDir()const{
	  return _initFileDir;
	}
	std::string getFileName()const{
	  return _initFilename;
	}

	template<class T>
	bool read(const std::string eleName, T &value, const bool checkValidation = true){
	  if( !_file.is_open() ){
		return false;
	  }
	  bool succ = nodeIsFound(eleName);
	  if(succ){
		succ = this->actualRead(eleName,value);
	  }else{
		nodeNotFound(eleName);
	  }
	  return succ;
	}
	bool readFilePath(const std::string eleName, std::string &filePath,const bool checkFileExist = true){
	  if( !_file.is_open() ){
		return false;
	  }
	  bool succ = read(eleName,filePath);
	  if(succ){
		filePath = _initFileDir+filePath;
		if(checkFileExist)
		  succ = fileIsExist(eleName,filePath);
	  }
	  return succ;
	}
	bool readFilePathes(const std::string eleName, std::vector<std::string> &filePathes,const bool checkFileExist = true){
	  if( !_file.is_open() ){
		return false;
	  }
	  bool succ = read(eleName,filePathes);
	  if(succ){
		for (size_t i = 0; i < filePathes.size(); ++i){
		  filePathes[i] = _initFileDir+filePathes[i];
		  if(checkFileExist)
			succ &= fileIsExist(eleName,filePathes[i]);
		}
	  }
	  return succ;
	}
	/// replace the first "~" as the home directory
	std::string replaceHomeDir(const std::string &path,const bool checkPathExist = true){
	  std::string newPath = path;
	  if(path.size() > 0 && path[0] == '~')
		newPath = getenv("HOME") + newPath.erase(0,1);
	  if(checkPathExist)
		dirIsExist("",newPath);
	  return newPath;
	}
	
  protected:
	// print the name and value
	template<class T>
	void print(const std::string eleName, const T &value)const{
	  INFO_LOG_COND(eleName << "\t" << value,_printData);
	}
	void nodeNotFound(const std::string eleName)const{
	  WARN_LOG("the node '"<< eleName <<"' is not found!");
	}

	// checking
	bool fileIsExist(const std::string eleName,const std::string filePath)const{
	  const bool exist = (boost::filesystem::exists(filePath) && 
						  !(boost::filesystem::is_directory(filePath)));
	  WARN_LOG_COND("file '"<< filePath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
	  return exist;
	}
	bool dirIsExist(const std::string eleName,const std::string dirPath)const{
	  const bool exist = boost::filesystem::is_directory(dirPath);
	  WARN_LOG_COND("dir '"<< dirPath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
	  return exist;
	}
	virtual bool nodeIsFound(const std::string eleName) = 0;

	// perform the actual reading process
	template<class T>
	bool actualRead(const std::string eleName, std::set<T> &value){
	  value.clear();
	  std::vector<T> vec;
	  if ( this->actualRead(eleName, vec) ){
		for (size_t i = 0; i < vec.size(); ++i){
		  value.insert(vec[i]);
		}
		return true;
	  }
	  return false;
	}

	virtual bool actualRead(const std::string eleName, bool &value) = 0;
	virtual bool actualRead(const std::string eleName, int &value) = 0;
	virtual bool actualRead(const std::string eleName, float &value) = 0;
	virtual bool actualRead(const std::string eleName, double &value) = 0;
	virtual bool actualRead(const std::string eleName, std::string &value) = 0;

	virtual bool actualRead(const std::string eleName, std::vector<int> &value) = 0;
	virtual bool actualRead(const std::string eleName, std::vector<float> &value) = 0;
	virtual bool actualRead(const std::string eleName, std::vector<double> &value) = 0;
	virtual bool actualRead(const std::string eleName, std::vector<std::string> &value) = 0;

	// read in the eigen3 data types
	virtual bool actualRead(const std::string eleName, Eigen::VectorXd &value) = 0;
	virtual bool actualRead(const std::string eleName, Eigen::MatrixXd &value) = 0;
	virtual bool actualRead(const std::string eleName, std::vector<Eigen::VectorXd> &value) = 0;

  protected:
	std::ifstream _file;
	std::string _initFilename; //< filename of this initfile.
	std::string _initFileDir; //< directory of the initfile.
	bool _printData; // if true, function print() will work. false in default.
  };

}//end of namespace

#endif /* _JSONFILEPASER_H_ */

