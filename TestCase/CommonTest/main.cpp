#define LOG_WARN
#define LOG_ERROR
#define LOG_DEBUG
#define LOG_INFO
#define LOG_TRACE

#include <JsonFilePaser.h>
#include <Timer.h>
#include <IO.h>
#include <assertext.h>
#include <iostream>
#include <Log.h>
#include <MatrixIO.h>
#include <eigen3/Eigen/Dense>
using namespace UTILITY;
using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]){

  JsonFilePaser file;
  int a;
  std::vector<std::string> av;
  file.read("hh",av);

  return 0;
}


// #include <stdio.h>

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <locale>

// #include "boost/property_tree/ptree.hpp"
// #include "boost/foreach.hpp"
// #include "boost/property_tree/json_parser.hpp"
// #include "boost/property_tree/xml_parser.hpp"

// int main(int argc, char **argv){
//   try {

// 	std::ifstream jsonFile;
// 	jsonFile.open("/home/simba/Workspace/Utility/TestCase/TestData/json.in");

// 	boost::property_tree::ptree ptParse;
// 	boost::property_tree::json_parser::read_json(jsonFile, ptParse);
// 	const int num = ptParse.get<int>("intnumber");
// 	std::string strVal = ptParse.get<std::string>("filename");
// 	std::cout << "Num=" << std::dec << num << " Str=" << strVal << std::endl << std::endl;

// 	BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,ptParse.get_child("strArray")){
// 	  std::cout << child.second.get<int>("") << std::endl;
//     }

//   }catch (std::exception& ex){
// 	std::cout << "ERR:" << ex.what() << std::endl;
//   }

//   return 0;
// }
