#ifndef _JSONFILEPASER_H_
#define _JSONFILEPASER_H_

#include <stdlib.h>
#include <vector>
#include <set>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/filesystem.hpp>
#include <Log.h>
using namespace boost::property_tree;

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
	  close();
	}

	bool open(const std::string filename){

	  bool succ = false;
	  close();
	  _initFilename = filename;
	  _file.open(filename.c_str());
	  ERROR_LOG_COND("failed to open init file: "<<filename,_file.is_open());
	  if(_file.is_open()){
		try{
		  json_parser::read_json(_file,_jsonData);
		  succ = read("init_file_dir",_initFileDir);
		  _initFileDir = replaceHomeDir(_initFileDir);
		}catch(std::exception& ex){
		  ERROR_LOG("file: "<< filename << ex.what());
		}
	  }
	  return succ;
	}
	void close(){ if(_file.is_open()) _file.close();}

	std::string getFileDir()const{
	  return _initFileDir;
	}
	std::string getFileName()const{
	  return _initFilename;
	}

	template<class T>
	bool read(const std::string eleName, T &value, const bool checkValidation = true){
	  if( !_file.is_open() ){
		ERROR_LOG("json file is not open");
		return false;
	  }
	  bool succ = true;
	  try{
		actualRead(eleName,value);
	  }catch(std::exception& ex){
		succ = false;
		ERROR_LOG(ex.what());
	  }
	  return succ;
	}
	bool readFilePath(const std::string eleName, std::string &filePath,const bool checkFileExist = true){
	  bool succ = false;
	  if( read(eleName,filePath) ){
		filePath = getFileDir()+filePath;
		succ = fileIsExist(eleName,filePath,checkFileExist);
	  }
	  return succ;
	}
	bool readFilePath(const std::string eleName, std::vector<std::string> &filePathes,const bool checkFileExist = true){
	  bool succ = false;
	  if(read(eleName,filePathes)){
		for (size_t i = 0; i < filePathes.size(); ++i){
		  filePathes[i] = getFileDir()+filePathes[i];
		  succ &= fileIsExist(eleName,filePathes[i],checkFileExist);
		}
	  }
	  return succ;
	}
	std::string replaceHomeDir(const std::string &path,const bool checkPathExist = true){
	  /// replace the first "~" as the home directory
	  std::string newPath = path;
	  if(path.size() > 0 && path[0] == '~')
		newPath = getenv("HOME") + newPath.erase(0,1);
	  dirIsExist("",newPath,checkPathExist);
	  return newPath;
	}
	
  protected:
	bool fileIsExist(const std::string eleName,const std::string filePath,bool check=true)const{
	  bool exist = true;
	  if(check){
		exist = (boost::filesystem::exists(filePath) && 
				 !(boost::filesystem::is_directory(filePath)));
		WARN_LOG_COND("file '"<< filePath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
	  }
	  return exist;
	}
	bool dirIsExist(const std::string eleName,const std::string dirPath,bool check)const{
	  bool exist = true;
	  if(check){
		exist = boost::filesystem::is_directory(dirPath);
		WARN_LOG_COND("dir '"<< dirPath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
	  }
	  return exist;
	}

	template<class T>
	void actualRead(const std::string eleName, T &value){
	  value = _jsonData.get<T>(eleName);
	}

	template<class T>
	void actualRead(const std::string eleName, std::vector<T> &value){
	  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,_jsonData.get_child(eleName)){
		value.push_back(child.second.get<T>(""));
	  }
	}

	template<class T>
	void actualRead(const std::string eleName, std::set<T> &value){
	  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,_jsonData.get_child(eleName)){
		value.insert(child.second.get<T>(""));
	  }
	}

  private:
	std::ifstream _file;
	std::string _initFilename; 
	std::string _initFileDir;
	bool _printData;
	ptree _jsonData;
  };

}//end of namespace

#endif /* _JSONFILEPASER_H_ */

