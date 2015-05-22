// Matrix.h
// Shaun Harker
// 9/26/11

#ifndef CHOMP_MATRIX_H
#define CHOMP_MATRIX_H

#include <cstdlib>
#include "chomp/SparseMatrix.h"
#include "chomp/Chain.h"
#include "chomp/Complex.h"
#include "chomp/SmithNormalForm.h"

#undef Complex
#undef complex
namespace chomp {
  
Matrix chainsToMatrix ( const std::vector < Chain > & chains, 
                        const Complex & complex,
                        int dim );

Matrix chainsToMatrix ( const std::vector < std::pair < Chain, Ring > > & chains, 
                       const Complex & complex,
                       int dim );

Matrix SmithSolve ( const Matrix & A, const Matrix & B );

/// Definitions

inline Matrix chainsToMatrix ( const std::vector < Chain > & chains, 
                        const Complex & complex,
                        int dim ) {
  //std::cout << "chainsToMatrix.  d= " << dim << " (i,j) = (" << complex . size ( dim ) << ", " << chains . size () << ")\n";
  Matrix result ( complex . size ( dim ), chains . size () );
  for ( uint64_t chain_num = 0; 
        chain_num < chains . size (); 
        ++ chain_num ) {
    Chain chain = simplify ( chains [ chain_num ] );
    BOOST_FOREACH ( const Term & t, chain () ) {
      result . write ( t . index (), chain_num, t . coef () ); 
    }
  }
  return result;
}

inline Matrix chainsToMatrix ( const std::vector < std::pair < Chain, Ring > > & chains, 
                       const Complex & complex,
                       int dim ) {
    //std::cout << "chainsToMatrix.  d= " << dim << " (i,j) = (" << complex . size ( dim ) << ", " << chains . size () << ")\n";
  Matrix result ( complex . size ( dim ), chains . size () );
  for ( uint64_t chain_num = 0; 
        chain_num < chains . size (); 
        ++ chain_num ) {
    Chain chain = simplify ( chains [ chain_num ] . first );
    
    BOOST_FOREACH ( const Term & t, chain () ) {
      result . write ( t . index (), chain_num, t . coef () ); 
    }
  }
  return result;
}

/// Solve AX = B for X over PID modules using Smith Normal Form
inline Matrix SmithSolve ( const Matrix & A, const Matrix & B ) {
  Matrix U, Uinv, V, Vinv, D;
  SmithNormalForm ( &U, &Uinv, &V, &Vinv, &D, A );
  int n = A . number_of_rows ();
  int m = A . number_of_columns ();
  Matrix DT ( m, n );
  for (int i = 0; i < n && i < m; ++ i ) {
    if ( D . read (i, i) == Ring ( 0 ) ) continue; // avoid division by zero error
    DT . write (i, i, Ring ( 1 ) );
  }
  Matrix Y = DT * Uinv * B;
  for ( Matrix::size_type column = 0; 
       column < Y . number_of_columns (); 
       ++ column ) {
    for ( Matrix::MatrixPosition entry = Y . column_begin ( column ); 
         entry != Y . end (); 
         Y . column_advance ( entry ) ) {
      Matrix::size_type i = Y . row ( entry );
      Ring value = D . read ( i, i );
      if ( value == Ring ( 0 ) ) continue; // avoid division by zero error
                                           // DEBUG
      Ring newvalue = Y . read ( entry ) / value;
      if ( newvalue * value != Y . read ( entry ) ) {
      std::cout << "preboundary: INSOLUBLE MATRIX EQUATION! " 
        << Y . read ( entry ) 
        << " is not divisible by " 
        << value << ".\n"; 
        exit ( 1 );
      } 
      Y . write ( entry , newvalue );
    } /* for */
  } /* for */
  Matrix X = Vinv * Y;
  return X;
}

} // namespace chomp

#endif
