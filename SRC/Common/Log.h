#ifndef _LOG_H_
#define _LOG_H_

#include <string>
#include <iostream>

namespace UTILITY{

  template<class T>
  void PRINT_MSG(const std::string &title,const T&event,bool cond=true,const std::string file="",const int line=0){
	if(!cond){
	  return ;
	}
	std::cout<<"*"<<title<<"\t: "<<event;
	if(file.size() > 1){
	  std::cout<<"\t at\t: "<<file<<":"<<line<<std::endl;
	}else{
	  std::cout<<std::endl;
	}
  }

#define PRINT_MSG_MICRO(title,event,cond)					\
  if(cond){													\
	std::cout<<"*"<<title<<"\t: "<<event<<std::endl;		\
  }

#define PRINT_MSG_MICRO_EXT(title,event,cond,file,line)		\
  if(cond){													\
	std::cout<<"*"<<title<<"\t: "<<event;					\
	std::cout<<"\t at\t: "<<file<<":"<<line<<std::endl;		\
  }

  /************************ERROR******************************/
#ifdef LOG_ERROR
#define ERROR_LOG_COND(event,cond) {PRINT_MSG_MICRO_EXT("ERROR",event,cond,__FILE__,__LINE__);}
#define ERROR_LOG(event) {PRINT_MSG_MICRO_EXT("ERROR",event,true,__FILE__,__LINE__);}
#else
#define ERROR_LOG_COND(event,cond)
#define ERROR_LOG(event)
#endif

  /************************WARN******************************/
#ifdef LOG_WARN
  template<class T> 
  inline void WARN_LOG(const T&event,bool cond=true){ PRINT_MSG("WARN",event,cond);}
#else
#define WARN_LOG(event,...)
#endif

  /************************TRACE******************************/
#ifdef LOG_TRACE
  class TRACE_CLASS{
  public:
	TRACE_CLASS(std::string funInfo):_funInfo(funInfo){
	  PRINT_MSG("ENETER_FUN: ",_funInfo);
	}
	~TRACE_CLASS(){
	  PRINT_MSG("OUT_FUN: ",_funInfo);
	}
  private:
	const std::string _funInfo;
  };
#define TRACE_FUN() TRACE_CLASS _m_trace_class(__PRETTY_FUNCTION__);
#else
#define TRACE_FUN()
#endif

  /************************INFO******************************/
#ifdef LOG_INFO
#define INFO_LOG(event) {PRINT_MSG("INFO",event);}
#define INFO_LOG_COND(event,cond) {PRINT_MSG_MICRO("INFO",event,cond);}
#else
#define INFO_LOG(event,)
#define INFO_LOG_COND(event,cond)
#endif

  /************************DEBUG******************************/
#ifdef LOG_DEBUG
  template<class T> 
  inline void DEBUG_LOG(const T&event,bool cond=true) {PRINT_MSG("DEBUG",event,cond);}
#else
#define DEBUG_LOG(event,...)
#endif
  
}

#endif /* _LOG_H_ */
