#ifndef _DRAWCURVES_H_
#define _DRAWCURVES_H_

#include <string>
#include <sstream>
#include <fstream>
#include <assertext.h>
#include <Log.h>
using std::string;
using std::stringstream;

namespace UTILITY{

  /**
   * given a set of points (xi,yi), generate a python script to for drawing this
   * curve, and save it to a file.
   */
  class PythonScriptDraw2DCurves{
  
  public:
	template<class VECTOR>
	static bool write(const string fname,const VECTOR &y,const double dx=1.0,const double x0=0.0f){

	  VECTOR x(y.size());
	  for (int i = 0; i < x.size(); ++i){
		x[i] = x0+i*dx;
	  }
	  return write(fname,y,x);
	}

	template<class VECTOR>
	static bool write(const string fname,const VECTOR &y,const VECTOR &x){

	  assert_eq(x.size(),y.size());
	  if(y.size() <= 0){
		ERROR_LOG("the input data is zero");
		return false;
	  }

	  stringstream script;
	  script << head();

	  script << "x = [" << x[0];
	  for (int i = 1; i < x.size(); ++i){
		script << "," << x[i];
	  }	  
	  script << "];\n";

	  script << "y = [" << y[0];
	  for (int i = 1; i < y.size(); ++i){
		script << "," << y[i];
	  }
	  script << "];\n";

	  script << "plt.plot(x,y);\n";
	  script << end();

	  std::ofstream file;
	  file.open(fname.c_str());
	  ERROR_LOG_COND("failed to open file for writing: "<<fname,file.is_open());
	  file << script.str();
	  const bool succ = file.is_open();
	  file.close();
	  return succ;
	}

  protected:
	static string head(){
	  const string line1("\"\"\"Script generated by PythonScriptDraw2DCurves.\"\"\"\n");
	  const string line2("import numpy as np\n");
	  const string line3("import matplotlib.pyplot as plt\n");
	  return line1+line2+line3;
	}
	static string end(){
	  return string("plt.show()");
	}

  };
}
#endif /* _DRAWCURVES_H_ */
