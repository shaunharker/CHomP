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

#include <gmpxx.h>

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
  
  Long balanced_value ( void ) const {
    return x;
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

class GMP_Integer {
private:
  mpz_class x;
public:
  GMP_Integer ( void ) : x ( 0 ) {}
  GMP_Integer ( int64_t x ) : x(x) {}
  GMP_Integer ( const mpz_class& x) : x(x) {}
  
  bool invertible ( void ) const {
    if ( x == 1 ) return true;
    if ( x == -1 ) return true;
    return false;
  }
  
  GMP_Integer balanced_value ( void ) const {
    return x;
  }

  GMP_Integer operator - ( void ) const {
    return mpz_class(-x);
  }
  
  GMP_Integer & operator += ( const GMP_Integer & rhs ) {
    x += rhs . x;
    return *this;
  }
  
  GMP_Integer & operator -= ( const GMP_Integer & rhs ) {
    x -= rhs . x;
    return *this;
  }
  
  GMP_Integer & operator *= ( const GMP_Integer & rhs ) {
    x *= rhs . x;
    return *this;
  }
  
  GMP_Integer operator + ( const GMP_Integer & rhs ) const {
    return mpz_class(x + rhs . x);
  }
  
  GMP_Integer operator - ( const GMP_Integer & rhs ) const {
    return mpz_class(x - rhs . x);
  }
  
  GMP_Integer operator * ( const GMP_Integer & rhs ) const {
    return mpz_class(x * rhs . x);
  }

  GMP_Integer operator / ( const GMP_Integer & rhs ) const {
    return mpz_class(x / rhs . x);
  }
  
  bool operator == ( const GMP_Integer & rhs ) const {
    return x == rhs . x;
  }

  bool operator != ( const GMP_Integer & rhs ) const {
    return x != rhs . x;
  }
  
  bool operator < ( const GMP_Integer & rhs ) const {
    return x < rhs . x;
  }
  
  bool operator > ( const GMP_Integer & rhs ) const {
    return x > rhs . x;
  }
  
  friend std::ostream & operator << ( std::ostream & outstream, const GMP_Integer & rhs );
  
  friend class boost::serialization::access;
  template < class Archive >
  void save(Archive & ar, const unsigned int version) const
  {
      const std::string mpz_string = x.get_str();
      ar & boost::serialization::make_nvp("mpz",mpz_string);
  }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
      std::string mpz_string;
      ar & boost::serialization::make_nvp("mpz",mpz_string);
      x.set_str(mpz_string,10);
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
};

inline std::ostream & operator << ( std::ostream & outstream, const GMP_Integer & rhs ) {
  outstream << rhs . x;
  return outstream;
}

//typedef Long Ring;
#include "Field.h"


#ifndef RINGDEFINED
#define RINGDEFINED
//typedef Long Ring;
typedef Zp<5> Ring;
#endif
  
} // namespace chomp

#endif
