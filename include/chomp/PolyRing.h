// PolyRing.h
// Justin Bush
// 10/2/11

#ifndef CHOMP_POLYRING_H
#define CHOMP_POLYRING_H

#include <iostream>
#include <vector>

#include "chomp/Field.h"

namespace chomp {
  
template < class Field >
class PolyRing {
private:
	std::vector < Field > coefficients_;
public:
	PolyRing ( void ) {}
	PolyRing ( int value ) { coefficients_ . push_back ( Field ( value ) ); }
	PolyRing ( const Field & value ) { coefficients_ . push_back ( value ); }
	void resize ( int d ) { coefficients_ . resize ( d, Field ( 0 ) ); }
	int degree ( void );
	int degree ( void ) const;
	Field & operator [] ( const int index ) { return coefficients_ [ index ]; }
	const Field & operator [] ( const int index ) const { return coefficients_ [ index ]; }
	PolyRing & operator += ( const PolyRing & rhs);
};

template < class Field >
std::ostream & operator << ( std::ostream & outstream, const PolyRing<Field> & print_me ) {
	if ( print_me . degree () == -1 ) {
		outstream << "0";
		return outstream;
	}
	if ( print_me . degree () == 0 ) {
		outstream << print_me [ 0 ];
		return outstream;
	}
	
	if ( print_me [ print_me . degree () ] == 1 ) {
		outstream << "x^" << print_me . degree ();
	} else {
		outstream << print_me [ print_me . degree () ] << "x^" << print_me . degree ();
	}
	
	for ( int i = (int) print_me . degree () - 1; i >= 0; -- i ) {
	  if ( print_me [ i ] == 0 ) continue;
	  if ( print_me [ i ] == 1 ) {
		  if ( i > 0 ) {
			outstream << " + " << "x^" << i;
		  } else {
			outstream << " + 1";
		  }
		} else {
			if ( i > 0 ) {
				outstream << " + " << print_me [ i ] << "x^" << i;
			} else {
				outstream << " + " << print_me [ i ];
			}
		}
	}
	return outstream;
}


template < class Field >
inline int PolyRing<Field>::degree ( void ) {
	int d = coefficients_ . size () - 1;
	while ( d >= 0 && ( coefficients_ [ d ] == Field ( 0 ) ) ) -- d;
	coefficients_ . resize ( d + 1 );
	return d;
}

template < class Field >
inline int PolyRing<Field>::degree ( void ) const {
	int d = coefficients_ . size () - 1;
	while ( d >= 0 && ( coefficients_ [ d ] == Field ( 0 ) ) ) -- d;
	return d;
}

template < class Field >
PolyRing<Field> & PolyRing<Field>::operator += ( const PolyRing & rhs ) {
	if ( rhs . degree () > degree () ) {
		coefficients_ . resize ( rhs . degree () + 1 );
	}
	for ( int i = 0; i <= rhs . degree (); ++ i ) {
		coefficients_ [ i ] += rhs [ i ];
	}
	return *this;
}


template < class Field >
PolyRing<Field> operator - (const PolyRing<Field> & x) {
	int degree = x . degree ();
	PolyRing<Field> result;
  result . resize ( degree + 1 );
	for ( int i = 0; i <= degree; ++ i ) {
		result [ i ] =  - x [ i ];
	}
	result . degree (); // resizes.
	return result;
}



template < class Field >
inline PolyRing<Field> operator + (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	int degree = std::max ( lhs . degree, rhs . degree );
	PolyRing<Field> result;
  result . resize ( degree + 1 );
	for ( int i = 0; i <= degree; ++ i ) {
    Field lvalue, rvalue;
    if ( lhs . degree () < i ) lvalue = Field ( 0 ); 
    else lvalue = lhs [ i ];
    if ( rhs . degree () < i ) rvalue = Field ( 0 ); 
    else rvalue = rhs [ i ];
		result [ i ] = lvalue + rvalue;
  }
	result . degree (); // resizes.
	return result;
}

template < class Field >
PolyRing<Field> operator - (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	int degree = std::max ( lhs . degree (), rhs . degree ());
	PolyRing<Field> result;
  result . resize ( degree + 1 );
	for ( int i = 0; i <= degree; ++ i ) {
    Field lvalue, rvalue;
    if ( lhs . degree () < i ) lvalue = Field ( 0 ); 
    else lvalue = lhs [ i ];
    if ( rhs . degree () < i ) rvalue = Field ( 0 ); 
    else rvalue = rhs [ i ];
		result [ i ] = lvalue - rvalue;
	}
	result . degree (); // resizes.
	return result;
}

template < class Field >
PolyRing<Field> operator * (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	int degree = lhs . degree () + rhs . degree ();
	PolyRing<Field> result;
  result . resize ( degree + 1 );
	for ( int i = 0; i <= lhs . degree (); ++ i ) {
		for ( int j = 0; j <= rhs . degree (); ++ j ) {
			result [ i + j ] += lhs [ i ] * rhs [ j ];
		}
	}
	return result;
}

template < class Field >
PolyRing<Field> operator * (const Field & lhs, const PolyRing<Field> & rhs) {
	PolyRing<Field> result = rhs;
	for ( int i = 0; i <= rhs . degree (); ++ i ) {
		result [ i ] = lhs * result [ i ];
	}
	return result;
}


template < class Field >
PolyRing<Field> operator / (const PolyRing<Field> & dividend, const PolyRing<Field> & divisor) {
  //std::cout << " Polyring division.\n";
  //std::cout << " dividend = " << dividend << ", has degree " << dividend . degree () << "\n";
  //std::cout << " divisor = " << divisor << ", has degree " << divisor . degree () << "\n";
  //NSLog ( @" polyring / \n" );
	PolyRing<Field> remainder = dividend;
	PolyRing<Field> quotient;
  int quotient_degree = dividend . degree () - divisor . degree ();
  if ( quotient_degree < 0 ) quotient_degree = -1;
  quotient . resize ( quotient_degree + 1 );
  Field x = divisor [ divisor . degree () ];
	while ( divisor . degree () <= remainder . degree () ) {
    //NSLog ( @" polyring / step \n" );
    //std::cout << " dividend = " << dividend << ", has degree " 
    //                            << dividend . degree () << "\n";
    //std::cout << " divisor = " << divisor << ", has degree " 
    //                           << divisor . degree () << "\n";
  
    //std::cout << "Remainder = " << remainder << "\n";
    Field a = remainder [ remainder . degree () ];
		Field q = a / x;
		int exponent = remainder . degree () - divisor . degree ();
		quotient [ exponent ] = q;
		for ( int i = 0; i <= divisor . degree (); ++ i ) {
			remainder [ exponent + i ] = remainder [ exponent + i ] - q * divisor [ i ];
		}
	}
  //std::cout << " division returning with quotient = " << quotient << "\n";
	return quotient;
}

template < class Field >
bool operator == (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	if ( lhs . degree () != rhs . degree () ) return false;
	for ( int i = 0; i <= lhs . degree (); ++ i ) {
		if ( lhs [ i ] != rhs [ i ] ) return false;
	}
	return true;
}

template < class Field >
bool operator != (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	if ( lhs . degree () != rhs . degree () ) return true;
	for ( int i = 0; i <= lhs . degree (); ++ i ) {
		if ( lhs [ i ] != rhs [ i ] ) return true;
	}
	return false;
}


template < class Field >
bool operator < (const PolyRing<Field> & lhs, const PolyRing<Field> & rhs) {
	return lhs . degree () < rhs . degree ();
}

} // namespace chomp

#endif