// Generators.h
// Shaun Harker
// 9/14/11

#ifndef CHOMP_GENERATORS_H
#define CHOMP_GENERATORS_H

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <iostream>

#include "chomp/Complex.h"
#include "chomp/Chain.h"
#include "chomp/SparseMatrix.h"
#include "chomp/SmithNormalForm.h"
#include "chomp/MorseComplex.h"

namespace chomp {
  
/*************************
 *    DECLARATIONS       *
 *************************/

/// Generators_t
/// Organized by dimension on the outermost vector,
/// each pair is a Chain x along with a coefficient k for which
/// kx is a boundary. The coefficient k is not a unit (not invertible).
/// It is only 0 if there are no other values for which kx is a boundary
/// but x is not.
typedef std::vector < std::vector < std::pair < Chain, Ring > > >  Generators_t;

Generators_t SmithGenerators ( const Complex & complex, int cutoff_dimension );
Generators_t MorseGenerators ( const Complex & complex, int cutoff_dimension );


/* Compute Homology Groups and also Homology generators */
/*
 
 First, calculate the smith form of each boundary matrix, along with the transformation matrices
 CAPD smithForm gives us 
 B_in = Q*B_out*Rinv and 
 s = number of "1"'s in the diagonal
 t = number of non-zero diagonal entries
 
 The trivial (boundary) generators are the left-most columns of Q
 The torsion (weak boundary) generators are the middle columns of Q.
 The betti number generators require a little work. The relevant information
 is involved in the Smith Normal Form decomposition of two consecutive boundary matrices.
 One also needs (R, Rinv) from the previous decomposition. Then the betti generators are given by
 a basis for the intersection of the spaces spanned by the right-most columns of R and the
 right-most columns of Q. 
 
 This is because:
 a) The right-most columns of R are a basis for the cycle subspace.
 b) The right-most columns of Q are a basis for the non-boundary subspace. 
 
 The matrix given by the bottom-most rows of Qinv (corresponding to the right-most columns of Q)
 is thus interesting because it gives a matrix which takes a chain, 'projects the boundary part out',
 and expresses it in the Q-basis.
 
 In particular the product M=(bottom-Qinv)*(right-R) is important; to obtain the Betti generators
 we must determine the rank r of this matrix (which is the betti number) and obtain exactly r
 linearly independent column vectors which span the range of M. We then multiply these columns
 by right-Q to get back into the ambient space.
 
 Proof:
 
 (It needs to be shown that any r linearly independent columns in this
 matrix has the same span as the entire matrix itself.)
 
 */

inline Generators_t SmithGenerators (const Complex & complex,
                                     int cutoff_dimension = -1 ) {
  if ( cutoff_dimension == -1 ) cutoff_dimension = complex . dimension ();
   if ( cutoff_dimension > complex . dimension () ) cutoff_dimension = complex . dimension ();
  // Prepare output.
  Generators_t return_value ( complex . dimension () + 1 );
	
  typedef SparseMatrix<Ring> Matrix;
	/* Loop through Chain Groups */
  Matrix first_Vinv;
  int first_t, second_s, second_t;
	/* The d_{0} boundary matrix is zero, so its SNF has no non-zero elements on its diagonal. */
	first_t = 0;
	for ( int d = 0; d <= cutoff_dimension; ++ d ) {
		/*******************************************************************
     * Compute the SNF of the new boundary matrix, d_{d} *             *
     *******************************************************************/
        
    Matrix second_U, second_Uinv, second_V, second_Vinv, smith_diagonal, boundary_matrix;
    BoundaryMatrix ( &boundary_matrix, complex, d + 1 );
    SmithNormalForm( &second_U, &second_Uinv, &second_V, &second_Vinv, &smith_diagonal, boundary_matrix);
    
    second_s = 0;
    for ( int i = 0; i < std::min(smith_diagonal . number_of_rows (), 
                                  smith_diagonal . number_of_columns ()); ++ i ) {
      if ( invertible ( smith_diagonal . read ( i, i ) ) ) ++ second_s; 
    }
    second_t = std::min(smith_diagonal . number_of_rows (), 
                        smith_diagonal . number_of_columns ());
    for ( int i = 0; i < std::min(smith_diagonal . number_of_rows (), 
                                  smith_diagonal . number_of_columns ()); ++ i ) {
      if ( smith_diagonal . read ( i, i ) == Ring(0) ) {
        second_t = i; 
        break;
      }
    }
    

    /* Here is the current method of calculation.
     We want to obtain the generators.
     The torsion generators can be read right off from the middle columns of U_2.
     To obtain the betti generators, here is a possible solution:
     1 Get the rightmost columns of U_2, call it A.
     2 Get the bottom rows of U^{-1}_2. Call it B. 
     3 Get the the rightmost columns of V^{-1}_1. Call it C.
     4 Form the product P = CAB.
     5 Find the smith normal form UDV of P. Claim: D is a projection.
     Now the betti generators are the non-zero columns of U corresponding to non-zero entries in D
     */
		/* Obtain the relevant sub-matrix from second_Qinv */
		Matrix A, B, C, U, Uinv, V, Vinv, D, P, G;
    Submatrix ( & A, 
               0, second_U . number_of_rows () - 1, 
               second_t, second_U . number_of_columns () - 1,
               second_U );
    Submatrix ( & B, 
               second_t, second_Uinv . number_of_rows () - 1, 
               0, second_Uinv . number_of_columns () - 1,
               second_Uinv );
    if ( d == 0 ) {
      Matrix::MatrixPosition n = complex . size ( 0 );
      C . resize ( n, n );
      for ( Matrix::MatrixPosition i = 0; i < n; ++ i ) C . write ( i, i, Ring ( 1 ) );
    } else {
      Submatrix ( & C, 
                 0, first_Vinv . number_of_rows () - 1, 
                 first_t, first_Vinv . number_of_columns () - 1,
                 first_Vinv );
    }

    P = A * (B * C);
    SmithNormalForm( &U, &Uinv, &V, &Vinv, &D, P);
    G = U * D;

		/* Prepare the output */
		const unsigned int betti_number = second_U . number_of_rows () - first_t - second_t;
		const unsigned int torsion_number = second_t - second_s;
    //const unsigned int trivial_number = second_s;
    //std::cout << " betti number = " << betti_number << ", torsion number = " << torsion_number << ", trivial number = " << trivial_number << "\n";
		std::vector < std::pair < Chain, Ring > > & generators = return_value [ d ];		
    generators . resize ( betti_number + torsion_number );
    /* Insert the betti generators */
    unsigned int betti_index = 0;
		for (Matrix::size_type column_number = 0; 
         column_number < G . number_of_columns (); ++ column_number ) {
      if ( G . column_size ( column_number ) == 0 ) continue;
			generators [ betti_index ] . second = 0;
			Chain & generator_chain = generators [ betti_index ] . first;
      generator_chain . dimension () = d;
      ++ betti_index;
			for (Matrix::MatrixPosition entry = G . column_begin ( column_number ); 
           entry != G . end (); G . column_advance ( entry ) ) {
				generator_chain += Term ( G . row ( entry ), 
                                  G . read ( entry ) );
      } /* for */
      //std::cout << "Betti Chain: " << generator_chain << "\n";
		}
    
		/* Insert the torsion generators */
		for ( unsigned int torsion_index = 0; torsion_index < torsion_number; ++ torsion_index ) {
      //std::cout << "torsion_index = " << torsion_index << " and second_s = " << second_s << "\n";
      //std::cout << "smith_diagonal = \n";
      //print_matrix ( smith_diagonal );
			generators [ betti_number + torsion_index ] . second = 
      smith_diagonal . read ( second_s + torsion_index, 
                             second_s + torsion_index);
			Chain & generator_chain = generators [ betti_number + torsion_index ] . first;
      generator_chain . dimension () = d;
      for (Matrix::MatrixPosition entry = second_U . column_begin ( second_s + torsion_index ); 
           entry != second_U . end (); second_U . column_advance ( entry ) ) {
        generator_chain += Term ( second_U . row ( entry ), 
                                  second_U . read ( entry ) );
      } /* for */
      //std::cout << "Torsion Chain: (order = " << generator_chain << "\n";
		}
    
		/* Store second_* into first_* */
    
		first_Vinv = second_Vinv;
		first_t = second_t;
		
	} /* for */
  //std::cout << "Returning from HGSNF\n";
	return return_value; 
} /* void SmithGenerators(...) */

  inline Generators_t MorseGenerators ( Complex & complex,
                                        int cutoff_dimension = -1 ) {
    if ( cutoff_dimension == -1 ) cutoff_dimension = complex . dimension ();
    if ( cutoff_dimension > complex . dimension () ) cutoff_dimension = complex . dimension ();
                                     
  Generators_t return_value;
  /* Perform Single Morse Reductions to complex */
  // TODO (possibly): make more than one. Use a tower.
  MorseComplex morse_complex ( complex );
  /* Get the Homology Generators on the Morse level */
  Generators_t deep_gen = SmithGenerators ( morse_complex, cutoff_dimension );
  /* Lift the Homology Generators to the top level */
  return_value . resize ( complex . dimension () + 1 );
  for (unsigned int d = 0; d < deep_gen . size (); ++ d) {
    return_value [ d ] . resize ( deep_gen [ d ] . size () );
    for (unsigned int gi = 0; gi < deep_gen [ d ] . size (); ++ gi) {
      Chain lifted = morse_complex . lift ( deep_gen [ d ] [ gi ] . first );
      return_value [ d ] [ gi ] = std::make_pair ( lifted, deep_gen [ d ] [ gi ] . second);
    } /* for */
  } /* for */
  return return_value;
} 

} // namespace chomp

#endif
