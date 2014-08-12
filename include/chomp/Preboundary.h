// Preboundary.h
// Shaun Harker
// 9/21/11

#ifndef CHOMP_PREBOUNDARY_H
#define CHOMP_PREBOUNDARY_H
#include "chomp/Complex.h"
#include "chomp/Chain.h"
#include "chomp/Matrix.h"

namespace chomp {
  
/// Compute a preboundary of a chain in a complex
/// method: turn into an algebra problem, use Smith Normal Form.
inline Chain SmithPreboundary ( const Chain & input, 
                         const Complex & complex ) {
  // REQUIRED ENUMERATION
  int dim = input . dimension ();
  Matrix M;
  BoundaryMatrix ( &M, complex, dim + 1 );
  Matrix B ( complex . size ( dim ), 1 );
  Chain simplified_input = simplify ( input );
  BOOST_FOREACH ( const Term & t, simplified_input () ) {
    B . write ( t . index (), 0,
                t . coef () );
  }
  // Solve MX = B
  Matrix X = SmithSolve ( M, B );
  //print_matrix ( M*X);
  //print_matrix ( B);

  Chain return_value;
  return_value . dimension () = dim + 1;
  for ( Matrix::MatrixPosition entry = X . column_begin ( 0 );
        entry != X . end (); X . column_advance ( entry ) ) {
    return_value += Term (X . row ( entry ), 
                          X . read ( entry ) );
  }
  
  /*
   // DEBUG
  std::cout << "SmithPreboundary\n";
  std::cout << "  Input: " << input << "\n";
  std::cout << "  Output: " << return_value << "\n";
  Chain bd = complex . boundary ( return_value );
  std::cout << "  Difference: " << simplify ( input - bd ) << "\n";
   */
  return return_value;
}

} // namespace chomp

#endif
