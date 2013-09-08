#ifndef _TIMER_H_
#define _TIMER_H_

#include <string>
#include <iomanip>
#include <Log.h>
#include <assertext.h>

namespace UTILITY{

#if (defined __unix__) || (defined __APPLE__)

#include <stdlib.h>
#include <sys/time.h>

  class PerformanceCounter
  {
  protected:

	long startCountSec,stopCountSec,startCountMicroSec,stopCountMicroSec;

  public:

	PerformanceCounter()
	{
	  // also, reset the starting counter
	  StartCounter();
	}

	void StartCounter(); // call this before your code block
	void StopCounter(); // call this after your code block

	// read elapsed time (units are seconds, accuracy is up to microseconds)
	double GetElapsedTime();

  };

  inline void PerformanceCounter::StartCounter()
  {
	struct timeval tv;

	gettimeofday(&tv,NULL);

	startCountSec = tv.tv_sec;
	startCountMicroSec = tv.tv_usec;

  }

  inline void PerformanceCounter::StopCounter()
  {
	struct timeval tv;

	gettimeofday(&tv,NULL);

	stopCountSec = tv.tv_sec;
	stopCountMicroSec = tv.tv_usec;
  }


  inline double PerformanceCounter::GetElapsedTime()
  {
	float elapsedTime = 1.0 * (stopCountSec-startCountSec) + 1E-6 * (stopCountMicroSec - startCountMicroSec);
	return elapsedTime;
  }

#endif

#ifdef WIN32

  /**************** WINDOWS COUNTER *******************/

#include <windows.h>

  class PerformanceCounter
  {
  protected:

	LARGE_INTEGER timerFrequency;
	LARGE_INTEGER startCount,stopCount;

  public:

	PerformanceCounter() 
	{
	  // reset the counter frequency
	  QueryPerformanceFrequency(&timerFrequency);
	  // also, reset the starting counter
	  StartCounter();
	}

	void StartCounter(); // call this before your code block
	void StopCounter(); // call this after your code block

	// read elapsed time (units are seconds, accuracy is up to microseconds)
	double GetElapsedTime();
  };

  inline void PerformanceCounter::StartCounter()
  {
	QueryPerformanceCounter(&startCount);
  }

  inline void PerformanceCounter::StopCounter()
  {
	QueryPerformanceCounter(&stopCount);
  }

  inline double PerformanceCounter::GetElapsedTime()
  {
	return ((double)(stopCount.QuadPart - startCount.QuadPart))
	  / ((double)timerFrequency.QuadPart);
  }

#endif

  class Timer{
	
  public:
	void start(){
	  _PerformanceCounter.StartCounter();
	}
	double stop(const std::string msg="time = ",const bool printOut=true,const int outPrecision=10){

	  assert_gt(outPrecision,0);
	  _PerformanceCounter.StopCounter();
	  const double elapsed_time = _PerformanceCounter.GetElapsedTime();
	  INFO_LOG_COND(msg<<std::setprecision(outPrecision)<<elapsed_time,printOut);
	  return elapsed_time;
	}
	
  private:
  	PerformanceCounter _PerformanceCounter;

  };
  
}//end of namespace

#endif /*_TIMER_H_*/
