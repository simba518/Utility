#ifndef _LOG_H_
#define _LOG_H_

#include <iostream>

namespace UTILITY{

#define PRINT_MSG(TITLE , logEvent)						\
  std::cout<< TITLE << ": " << logEvent << std::endl;	

  /************************ERROR******************************/
#define ERROR_LOG(logEvent)					
#define ERROR_COND_LOG(logEvent, cond)		
#ifdef LOG_ERROR
#define ERROR_LOG(logEvent)						\
  PRINT_MSG("ERROR",logEvent);					
#define ERROR_COND_LOG(logEvent, cond)			\
  if (cond)										\
	{											\
	  ERROR_LOG(logEvent);						\
	}														
#endif

}

#endif /* _LOG_H_ */
