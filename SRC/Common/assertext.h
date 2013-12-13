#ifndef _ASSERTEXT_H_
#define _ASSERTEXT_H_
/**
 * this is some extetions for the assert command.
 */

#include <string>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <assert.h>

#if defined(WIN32) || defined(NDEBUG)/* Not NDEBUG.  */

# define assert_eq(value_a,value_b)		
# define assert_ne(value_a,value_b)		
# define assert_ge(value_a,value_b)		
# define assert_gt(value_a,value_b)		
# define assert_le(value_a,value_b)		
# define assert_lt(value_a,value_b)		
# define assert_in(value_a,min,max)

# define assert_eq_ext(value_a,value_b,info)		
# define assert_ne_ext(value_a,value_b,info)		
# define assert_ge_ext(value_a,value_b,info)		
# define assert_gt_ext(value_a,value_b,info)		
# define assert_le_ext(value_a,value_b,info)		
# define assert_lt_ext(value_a,value_b,info)		
# define assert_in_ext(value_a,min,max,info)

#else

// define a new assert micro to throw expections instead of abort
#define my_assert( e, file, line ) ( throw std::runtime_error(std::string(file)+":"+boost::lexical_cast<std::string>(line)+": failed assertion "+e))

#define new_assert(e) ((void) ((e) ? 0 : my_assert (#e, __FILE__, __LINE__)))

// assert "value_a == value_b", if not, print the values and abort.
#define assert_eq(value_a, value_b)										\
  if (value_a != value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a == value_b);									\
	}																	\
  
#define assert_eq_ext(value_a, value_b, info)							\
  if (value_a != value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a == value_b);									\
	}


// assert "value_a != value_b", if not, print the values and abort.
#define assert_ne(value_a, value_b)										\
  if (value_a == value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a != value_b);									\
	}														

#define assert_ne_ext(value_a, value_b, info)							\
  if (value_a == value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a != value_b);									\
	}														

// assert "value_a >= value_b", if not, print the values and abort.
#define assert_ge(value_a, value_b)										\
  if (value_a < value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a >= value_b);									\
	}														

#define assert_ge_ext(value_a, value_b, info)							\
  if (value_a < value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a >= value_b);									\
	}														

// assert "value_a > value_b", if not, print the values and abort.
#define assert_gt(value_a, value_b)										\
  if (value_a <= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a > value_b);									\
	}														

#define assert_gt_ext(value_a, value_b, info)							\
  if (value_a <= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a > value_b);									\
	}														

// assert "value_a >= value_b", if not, print the values and abort.
#define assert_le(value_a, value_b)										\
  if (value_a > value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a <= value_b);									\
	}														

#define assert_le_ext(value_a, value_b, info)							\
  if (value_a > value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a <= value_b);									\
	}														

// assert "value_a > value_b", if not, print the values and abort.
#define assert_lt(value_a, value_b)										\
  if (value_a >= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  new_assert(value_a < value_b);									\
	}														

#define assert_lt_ext(value_a, value_b, info)							\
  if (value_a >= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a < value_b);									\
	}														

// assert "value_a in [min,max]", if not, print the values and abort.
#define assert_in(value_a, min, max)									\
  if ( !(value_a >= min && value_a <= max) )							\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(min) <<" = " << min << std::endl;			\
	  std::cout << __STRING(max) <<" = " << max << std::endl;			\
	  new_assert(value_a >= min && value_a <= max);						\
	}														

#define assert_in_ext(value_a, min, max, info)							\
  if ( !(value_a >= min && value_a <= max) )							\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(min) <<" = " << min << std::endl;			\
	  std::cout << __STRING(max) <<" = " << max << std::endl;			\
	  std::cout << "additional info: "<< info << std::endl;				\
	  new_assert(value_a >= min && value_a <= max);						\
	}														

#endif /* NDEBUG.  */

#endif /* _ASSERTEXT_H_ */
