// Prism.h
// Shaun Harker
// 2/21/12 02/21/2012

#ifndef CHOMP_PRISM_H
#define CHOMP_PRISM_H

//#warning prism included
#define BOOST_UBLAS_NDEBUG

#include <iostream>
#include <vector>
#include <cmath>

#include "chomp/Rect.h"
#include "boost/serialization/serialization.hpp"

#include "boost/foreach.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas;

#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP
// REMEMBER to update "lu.hpp" header includes from boost-CVS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

#include "chomp/Real.h"

namespace chomp {

/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) {
 	using namespace boost::numeric::ublas;
 	typedef permutation_matrix<std::size_t> pmatrix;
 	// create a working copy of the input
 	matrix<T> A(input);
 	// create a permutation matrix for the LU-factorization
 	pmatrix pm(A.size1());
 	// perform LU-factorization
 	int res = lu_factorize(A,pm);
  if( res != 0 ) return false;
 	// create identity matrix of "inverse"
 	inverse.assign(ublas::identity_matrix<T>(A.size1()));
 	// backsubstitute to get the inverse
 	lu_substitute(A, pm, inverse);
 	return true;
}
#endif //INVERT_MATRIX_HPP

/*********
 * Prism *
 *********/

typedef ublas::matrix<Real, ublas::row_major> uMatrix;
typedef ublas::vector<Real> uVector;

class Prism {
public:
  // { Ax + c : -1 <= x <= 1 }
  int dim;
  uMatrix A; // edge vectors (half length)
  uVector c; // center
private:
  // member variables for intersection optimization
  mutable uMatrix Ainv;
  mutable bool Ainv_computed;
  mutable uVector D;
  mutable uVector d;
  mutable uVector a;
  mutable uVector b;
  mutable uVector x;
  mutable uVector u;
  mutable uVector v;

public:
  Prism ( void ) { Ainv_computed = false; dim = 0;}
  Prism ( int dim ) : dim ( dim ) {
    A . resize ( dim, dim );
    A = ublas::identity_matrix<Real> ( dim );
    c . resize ( dim );
    c = ublas::scalar_vector<Real> ( dim );
    
    Ainv . resize ( dim, dim );
    Ainv_computed = false;
    D . resize ( dim );
    d . resize ( dim );
    a . resize ( dim );
    b . resize ( dim );
    x . resize ( dim );
    u . resize ( dim );
    v . resize ( dim );

  }
  bool intersects ( const Rect & R ) const;
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & dim;
    ar & A;
    ar & c;
    Ainv . resize ( dim, dim );
    Ainv_computed = false;
    D . resize ( dim );
    d . resize ( dim );
    a . resize ( dim );
    b . resize ( dim );
    x . resize ( dim );
    u . resize ( dim );
    v . resize ( dim );
  }
  
};

std::ostream & operator << ( std::ostream & output_stream, const Prism & print_me );

inline bool Prism::intersects ( const Rect & r ) const {
  
  // This computes "weak intersection", the notion of
  // whether the prism and rectangle can be seperated with a hyperplane
  // which is either rect-aligned or prism-aligned (not an arbitrary hyperplane)
  // Compute "Ainv" if not already computed.
  if ( Ainv_computed == false ) {
    InvertMatrix ( A , Ainv ); 
    Ainv_computed = true;
  }
  
  // Now express the rectangle r in the form
  //      { x : -1 <= Dinv ( x - d )<= 1 }
  for ( int i = 0; i < dim; ++ i) {
    D ( i ) = (r . upper_bounds [ i ] - r . lower_bounds [ i ] ) / 2.0;
    d ( i ) = (r . upper_bounds [ i ] + r . lower_bounds [ i ] ) / 2.0;
  }
  
  // Note: A, D, are inverted compared to my notes.
  
  // Note prism may be expressed as { x : -1 <= Ainv(x-c) <= 1 }
  // a = (id + |Ainv|D) 1
  // b = (id + D^-1|A|)1

  for ( int i = 0; i < dim; ++ i) {
    a ( i ) = 1.0;
    b ( i ) = 1.0;
    for ( int j = 0; j < dim; ++ j) {
      a ( i ) += std::abs ( Ainv ( i, j ) ) * D ( j );
      b ( i ) += std::abs ( A ( i, j ) ) / D ( i );
    }
  }

  ublas::noalias(x) = c - d;

  // Check -a <= Ainv x <= a  && -b <= Dinv x <= b
  ublas::noalias(u) = ublas::prod ( Ainv, x );
  for ( int i = 0; i < dim; ++ i) {
    v ( i ) = x ( i ) / D ( i );
  }
  
  /*
  std::cout << "intersection check.\n prism = " << *this << "\n";
  std::cout << "rectangle = " << r << "\n";
  std::cout << "Ainv = " << Ainv << "\n";
  std::cout << "D = " << D << "\n";
  std::cout << "d = " << d << "\n";
  std::cout << "a = " << a << "\n";
  std::cout << "b = " << b << "\n";
  std::cout << "x = " << x << "\n";
  std::cout << "u = " << u << "\n";
  std::cout << "v = " << v << "\n";  
   */
  for ( int i = 0; i < dim; ++ i ) {
    if ( std::abs ( u ( i ) ) > a ( i ) ||
         std::abs ( v ( i ) ) > b ( i ) ) {
      //std::cout << "No intersection!\n";
      return false;
    }
  }
  //std::cout << "Intersection!\n";
  return true;
}

inline std::ostream & operator << ( std::ostream & output_stream, const Prism & print_me ) {
  output_stream << " A = " << print_me . A << ", c = " << print_me . c << "\n";
  return output_stream;
} 

} // namespace chomp

#endif
