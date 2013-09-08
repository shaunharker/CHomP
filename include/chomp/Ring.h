// Ring.h
// Shaun Harker
// 9/14/11

#ifndef CHOMP_RING_H
#define CHOMP_RING_H

#include <cstdlib>
#include <iostream>
#include <stdint.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include "Field.h"

namespace chomp {

/*********************
 *      Ring         *
 *********************/

template < class R >
bool invertible ( const R & element ) {
  return element . invertible ();
}

template <>
inline bool invertible < int64_t > ( const int64_t & element ) {
  if ( (element == 1) || (element == -1) ) return true;
  return false;
}

/// VARIOUS RING TYPES

class Long {
private:
  static const int64_t LongGone = ((int64_t)1 << 60);
  int64_t x;
public:
  Long ( void ) : x ( 0 ) {}
  Long ( int64_t x ) : x((int64_t)x) {}
  bool invertible ( void ) const {
    if ( x == 1 ) return true;
    if ( x == -1 ) return true;
    return false;
  }
  void integrity ( void ) {
    if ( (x > LongGone) || (-x > LongGone) ) {
      std::cout << "INTEGER OVERFLOW.\n";
      exit(1);
    }
  }
  
  Long operator - ( void ) const {
    Long result ( - x );
    result . integrity ();
    return result;
  }
  
  Long & operator += ( const Long & rhs ) {
    x += rhs . x;
    integrity ();
    return *this;
  }
  
  Long & operator -= ( const Long & rhs ) {
    x -= rhs . x;
    integrity ();
    return *this;
  }
  
  Long & operator *= ( const Long & rhs ) {
    x *= rhs . x;
    integrity ();
    return *this;
  }
  
  Long operator + ( const Long & rhs ) const {
    Long result ( x + rhs . x );
    result . integrity ();
    return result;
  }
  
  Long operator - ( const Long & rhs ) const {
    Long result ( x - rhs . x );
    result . integrity ();
    return result;
  }
  
  Long operator * ( const Long & rhs ) const {
    Long result ( x * rhs . x );
    result . integrity ();
    return result;
  }

  Long operator / ( const Long & rhs ) const {
    Long result ( x / rhs . x );
    result . integrity ();
    return result;
  }
  
  bool operator == ( const Long & rhs ) const {
    return x == rhs . x;
  }

  bool operator != ( const Long & rhs ) const {
    return x != rhs . x;
  }
  
  bool operator < ( const Long & rhs ) const {
    return x < rhs . x;
  }
  
  bool operator > ( const Long & rhs ) const {
    return x > rhs . x;
  }
  
  friend std::ostream & operator << ( std::ostream & outstream, const Long & rhs );
  
  friend class boost::serialization::access;
  template < class Archive >
  void serialize ( Archive & ar , const unsigned int version ) {
    ar & boost::serialization::make_nvp("e",x);
  }
};

inline std::ostream & operator << ( std::ostream & outstream, const Long & rhs ) {
  outstream << rhs . x;
  return outstream;
}

//typedef Long Ring;
#include "Field.h"


#ifndef RINGDEFINED
#define RINGDEFINED
//typedef Long Ring;
typedef Zp<3> Ring;
#endif
  
} // namespace chomp

#endif
