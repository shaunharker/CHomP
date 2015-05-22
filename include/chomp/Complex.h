// Complex.h
// Shaun Harker
// 9/10/11


#ifndef CHOMP_COMPLEX_H
#define CHOMP_COMPLEX_H

#undef Complex

#include <vector>
#include <utility>

#include "boost/foreach.hpp"
#include "boost/unordered_set.hpp"

#include "chomp/Ring.h"
#include "chomp/Chain.h"
#include "chomp/SparseMatrix.h"

namespace chomp {
  
class Complex {
protected:
  int dimension_;
  std::vector<Index> sizes_;
public:
  // Boundary and Coboundary Of Cells -- User Supplied.
  virtual void boundary ( Chain * output, const Index input, int dim ) const = 0;
  virtual void coboundary ( Chain * output, const Index input, int dim ) const = 0;

  // Boundary and Coboundary of Chains.
  virtual void boundary ( Chain * output, const Chain & input) const {
    int d = input . dimension ();
    BOOST_FOREACH ( const Term & term, input () ) {
      (*output) += term . coef () * boundary ( term . index (), d );
    }
    output -> dimension () = d - 1;
    // WARNING: UNSIMPLIFIED, MAY BE LIKE TERMS (GOTCHA)
  }
  
  virtual void coboundary ( Chain * output, const Chain & input ) const {
    int d = input . dimension ();
    BOOST_FOREACH ( const Term & term, input () ) {
      (*output) += term . coef () * coboundary ( term . index (), d );
    }
    output -> dimension () = d + 1;
    // WARNING: UNSIMPLIFIED, MAY BE LIKE TERMS (GOTCHA)
  }
  
  // Convenience Forms
  Chain boundary ( const Index input, int dim ) const {
    Chain output; boundary ( &output, input, dim ); return output; 
  }
  Chain coboundary ( const Index input, int dim ) const {
    Chain output; coboundary ( &output, input, dim ); return output; 
  }
  Chain boundary ( const Chain & input ) const {
    Chain output; boundary ( &output, input ); return output; 
  }
  Chain coboundary ( const Chain & input ) const {
    Chain output; coboundary ( &output, input ); return output; 
  }

  // Constructor
  Complex ( void ) : dimension_ (0) {}
  virtual ~Complex (void) {}
  
  // dimension accessor
  int & dimension ( void ) { return dimension_; }
  const int & dimension ( void ) const { return dimension_; }
  
  // size
  Index size ( int d ) const { 
    if ( d < 0 || d > dimension_ ) return 0;
    if ( (size_t) d >= sizes_ . size () ) return 0;
    return sizes_ [ d ]; 
  }
  Index size ( void ) const { 
    Index result = 0;
    for ( int d = 0; d <= dimension (); ++ d ) result += size ( d );
    return result;
  }

};

typedef SparseMatrix<Ring> Matrix;

/// BoundaryMatrix
/// Creates a sparse matrix representation of the boundary map 
/// The columns are the boundaries of the d-dimensional cells.
/// output is returned in "output_matrix".
void BoundaryMatrix ( SparseMatrix < Ring > * output_matrix, 
                      const Complex & complex, 
                     const int d );

inline void BoundaryMatrix ( SparseMatrix < Ring > * output_matrix, 
                            const Complex & complex, 
                            const int d ) {
  uint64_t rows, columns;
  if ( d > 0 ) rows = complex . size ( d - 1 ); else rows = 0;
  if ( d <= complex . dimension () ) columns = complex . size ( d ); else columns = 0;
  output_matrix -> resize ( rows, columns );
  //std::cout << " d = " << d << "\n";
  //std::cout << " complex . dimension () = " << complex . dimension () << "\n";
  //std::cout << " rows = " << rows << " and columns = " << columns << "\n";
  if ( d >= 0 && d <= complex . dimension () )
  for ( Index i = 0; i < complex . size ( d ); ++ i ) {
		Chain bd = complex . boundary ( i, d );
    //std::cout << " boundary(" << i << ", " << d << ") = " << bd << "\n";
    BOOST_FOREACH ( const Term & t, bd () ) {
      //std::cout << "Write to " << t . index () << ", " << i << "\n";
      //std::cout << "value to write = " << t . coef () << "\n";
      output_matrix -> write ( t . index (), i, t . coef () ); 
    }
  } /* for */
} 

} // namespace chomp

///  CHOMP_COMPLEX(Cell) Macro
///  This provides several pre-written functions to a class 
///    derived from Complex. The implementor should write
///      "CHOMP_COMPLEX(Cell)" as a line in his class definition, where
///       Cell is the class the implementor has written specifically for
///       the complex.
///
///  The methods provided by this macro are:
///    (1) insertCell
///    (2) indexToCell
///    (3) cellToIndex
///   
///  In order to prepare the cell indexing, the following must be done:
///    ====> User calls "insertCell" on all cells. Repeats are OK.
///    Now the complex is ready for use by algorithms using index interface
///   Note that indexes are assigned contiguously, first come, first served.
///   Repeat inserts are ignored.
///   The indexToCell method assumes the index, dimension pair is valid.
///   The cellToIndex method does not assume the cell, dimension pair is valid. If
///     an invalid pair is given, the number of dim-dimensional cells is returned,
///     which the user can check is a valid index or not.
///
///  In order to write "bd" and "cbd" methods, the implementor will:
///    (A) Get the index of a cell as an argment
///    (B) Use "indexToCell" to convert it to a cell
///    (C) Apply a boundary/coboundary algorithm to the cell
///    (D) Use "cellToIndex" to convert the cells in the result to indexes
///  

#define CHOMP_COMPLEX(MyCell)                                     \
using Complex::boundary;                                          \
using Complex::coboundary;                                        \
typedef MyCell Cell;                                              \
std::vector < boost::unordered_map < Cell, Index > > CCindexes_;  \
std::vector < std::vector < Cell > > CCcells_;                    \
void insertCell ( const Cell & c, const int dim ) {               \
  if ( CCindexes_ . size () <= (size_t) dim )                     \
    CCindexes_ . resize ( dim + 1 );                              \
  if ( CCcells_ . size () <= (size_t) dim )                       \
    CCcells_ . resize ( dim + 1 );                                \
  if ( sizes_ . size () <= (size_t) dim )                         \
    sizes_ . resize ( dim + 1, 0 );                               \
  if ( dim > dimension () ) {                                     \
    dimension () = dim;                                           \
  }                                                               \
  if ( CCindexes_ [ dim ] . count ( c ) == 0 ) {                  \
    CCindexes_ [ dim ] [ c ] = CCcells_ [ dim ] . size ();        \
    CCcells_ [ dim ] . push_back ( c );                           \
    sizes_ [ dim ] = CCcells_ [ dim ] . size ();                  \
  }                                                               \
}                                                                 \
Cell indexToCell ( const Index i, const int dim ) const {         \
  return CCcells_ [ dim ] [ i ];                                  \
}                                                                 \
Index cellToIndex ( const Cell & cell, const int dim ) const {    \
  if ( (size_t) dim >= CCindexes_ . size () ) return 0;           \
  boost::unordered_map < Cell, Index >::const_iterator it =       \
   CCindexes_ [ dim ] . find ( cell );                            \
  if ( it == CCindexes_ [ dim ] . end () ) return size ( dim );   \
  return it -> second;                                            \
}                                                               


#endif

