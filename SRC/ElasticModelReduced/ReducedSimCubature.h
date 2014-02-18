#ifndef _REDUCEDSIMCUBATURE_H_
#define _REDUCEDSIMCUBATURE_H_

#include <boost/shared_ptr.hpp>

namespace SIMULATOR{
  
  /**
   * @class ReducedSimCubature the cubature operator for the reduced simulation,
   * which compute the weights and sampled tetrahedrons.
   * 
   */
  class ReducedSimCubature{
	
  public:
	ReducedSimCubature(){}
	
  protected:
	
  private:
	
  };
  
  typedef boost::shared_ptr<ReducedSimCubature> pReducedSimCubature;
  
}//end of namespace

#endif /*_REDUCEDSIMCUBATURE_H_*/
