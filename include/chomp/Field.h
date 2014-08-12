// Field.h
// Justin Bush and Shaun Harker
// 10/2/11

#ifndef CHOMP_FIELD_H
#define CHOMP_FIELD_H

#include <iostream>
#include <stdint.h>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

namespace chomp {
/**********************
 * FINITE FIELD CLASS *
 **********************/

typedef int64_t Prime;

template <Prime p>
class Zp {
private:
	int64_t value_;
public:
	Zp ( void );
	Zp ( int z );
  bool invertible ( void ) const {
    if ( value_ == 0 ) return false;
    return true;
  }
	int64_t & value ( void ) { return value_; };
	const int64_t & value ( void ) const { return value_; };
	int64_t balanced_value ( void ) const;
	Zp & operator += ( const Zp & rhs );
  Zp & operator *= ( const Zp & rhs );
  bool operator < ( const Zp & rhs ) const { return false; } // :(  Bezout in SmithNormalForm requires this
                                             // That needs to be fixed, then this can 
                                             // be removed.
  /// The serialization method.
  friend class boost::serialization::access;
  template < class Archive >
  void serialize ( Archive & ar , const unsigned int version ) {
    ar & boost::serialization::make_nvp("v", value_);
  }
};

template < Prime p >
std::ostream & operator << ( std::ostream & outstream, const Zp<p> & print_me ) {
	outstream << print_me . balanced_value ();
	return outstream;
}

template < Prime p >
int64_t Zp<p>::balanced_value ( void ) const {
	int64_t v = value_;
	/*
	// DEBUG
	if ( v > value_ || v < 0 ) {
		std::cout << "Zp coefficient unexpectedly out of range.\n";
		abort ();
	}
	// END DEBUG
	*/
	if ( v > p / 2 ) v -= p;
	return v;
}

// constructors

template < Prime p >
inline Zp<p>::Zp ( void ) {
	value_ = 0;
}

template < Prime p >
inline Zp<p>::Zp ( int z ) {
	if ( z < 0 ) {
		int64_t n = - z / p;
		z = z + p*(n+1);
	}
	value_ = z % p;
}

template < Prime p >
inline Zp<p> & Zp<p>::operator += ( const Zp & rhs ) {
	value_ += rhs . value ();
	value_ = value_ % p;
  return *this;
}

template < Prime p >
inline Zp<p> & Zp<p>::operator *= ( const Zp & rhs ) {
	value_ *= rhs . value ();
	value_ = value_ % p;
  return *this;
}

template < Prime p >
inline bool operator == (const Zp<p> & lhs, const Zp<p> & rhs) {
	if ( lhs . value () == rhs . value () ) return true;
	return false;
}

template < Prime p >
inline bool operator == (const Zp<p> & lhs, const int & rhs) {
	if ( lhs . value () == rhs ) return true;
	return false;
}

template < Prime p >
inline bool operator != (const Zp<p> & lhs, const Zp<p> & rhs) {
	if ( lhs . value () != rhs . value () ) return true;
	return false;
}


template < Prime p >
inline bool operator != (const Zp<p> & lhs, const int & rhs) {
	if ( lhs . value () != rhs ) return true;
	return false;
}

// additive and multiplicative inverse
template < Prime p >
inline Zp<p> operator - ( const Zp<p> & x) {
	return Zp<p> ( - x . value () );
}

template < Prime p >
inline Zp<p> inverse (const Zp<p> & x) {
	// First, invert rhs. We use repeated squaring.
	int64_t sqr = x . value (); 
	int64_t inv = 1;
	int64_t p0 = p - 2;
	while ( p0 ) {
		if ( p0 & 1 ) {
			inv *= sqr;
			inv = inv % p;
		}
		p0 >>= 1;
		sqr *= sqr;
		sqr = sqr % p;
	}
	return  Zp<p> ( inv );
}


// arithmetic operations
template < Prime p >
inline Zp<p> operator + (const Zp<p> & lhs, const Zp<p> & rhs) {
	return Zp<p> ( lhs . value () + rhs . value () );
}

template < Prime p >
inline Zp<p> operator - (const Zp<p> & lhs, const Zp<p> & rhs) {
	return Zp<p> ( lhs . value () - rhs . value () );
}

template < Prime p >
inline Zp<p> operator * (const Zp<p> & lhs, const Zp<p> & rhs) {
	return Zp<p> ( lhs . value () * rhs . value () );
  
}

template < Prime p >
inline Zp<p> operator / (const Zp<p> & lhs, const Zp<p> & rhs) {
	return lhs * inverse ( rhs );
}
  
} // namespace chomp


#endif
