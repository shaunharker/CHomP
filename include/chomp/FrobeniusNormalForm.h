/// FrobeniusNormalForm.h
/// Shaun Harker 6/19/14

#ifndef CHOMP_FROBENIUSNORMALFORM_H
#define CHOMP_FROBENIUSNORMALFORM_H

#include <fstream>
#include <cstdlib>
#include <stdint.h>
#include "chomp/SparseMatrix.h"
#include "chomp/Ring.h"
#include "chomp/PolyRing.h"
#include "chomp/Algebra.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/thread.hpp>

namespace chomp {

/// FrobeniusNormalForm
/// Input: a square matrix A
/// Output: an array of polynomials P such that the diagonal
///         concatenation of their companion matrices is the
///         Frobenius Normal Form of A
/// Algorithm:  Storjohann, Arne. "An O(n^3) algorithm for the 
///             Frobenius normal form."  Proceedings of the 1998 
///             international symposium on Symbolic and algebraic 
///             computation. ACM, 1998. [Storjohann98]
///              + personal communication with Arne Storjohann
template < class R > 
std::vector<PolyRing<R> >
FrobeniusNormalForm ( SparseMatrix<R> const & A );

// Definitions.


namespace frobenius_detail {


/// ZigZagFormSubroutine
///    Perform algorithm of Lemma 1 in
///    [Storjohann98] to block [r,n)x[r,n)
template < class Matrix > void
ZigZagFormSubroutine ( Matrix * Z,
                       int64_t & r ) {


  //std::cout << "ZigZagFormSubroutine.\n";
  //print_matrix ( *Z );

  int64_t n = Z -> number_of_columns ();
  //std::cout << "Matrix is " << n << " by " << n << "\n";
  typedef typename Matrix::entry_type R;
  typedef typename SparseMatrix<R>::MatrixPosition MatrixPosition;

  MatrixPosition q;

  //    Stage 1.
  //std::cout << "ZigZagFormSubroutine. Stage 1.\n";
  int64_t j0 = r;
  int64_t j = r;
  while ( j < n ) {
    //std::cout << "ZigZagFormSubroutine. Stage 1. Top of loop. j = " << j << "\n";
    // Find a pivot row i \in [j+1, n)
    int64_t pivot_row;
    R pivot_value ( 0 );
    //std::cout << "ZigZagFormSubroutine. Stage 1. Search for pivot.\n";
    
    q = Z -> column_begin ( j );
    while ( q != Z -> end () ) {
      MatrixPosition p = q;
      Z -> column_advance ( q );
    
      pivot_row = Z -> row ( p );
      if ( pivot_row <= j ) continue;
      pivot_value = Z -> read ( p );
      if ( pivot_value != R(0) ) break;
    }
    if ( pivot_value == R(0) ) {
      break;
    }
    //std::cout << "ZigZagFormSubroutine. Stage 1. Swap pivot.\n";

    // Swap row i with row j+1 (and also the corresponding column op)
    Z -> similarity_operation ( R(0), R(1),
                                R(1), R(0),
                                pivot_row, j+1 );
    // Do a similarity operation so Z(j+1, j) == 1
    Z -> similarity_operation ( inverse(pivot_value) , j+1 );
    // Now we can pivot off (j+1,j) to clear the j column
    pivot_row = j+1;
    pivot_value = R(1);
    //std::cout << "ZigZagFormSubroutine. Stage 1. Elimination steps.\n";
    //std::cout << "ZigZagFormSubroutine. Pivot_value = " << Z -> read ( j+1, j ) << "\n";
    
    //for ( MatrixPosition p = Z -> column_begin ( j );
    //      p != Z -> end ();
    //      Z -> column_advance ( p ) ) {
    
    q = Z -> column_begin ( j );
    while ( q != Z -> end () ) {
      MatrixPosition p = q;
      Z -> column_advance ( q );

      //sane_matrix ( *Z );
      int64_t elim_row = Z -> row ( p );
      R val = Z -> read ( p );
      if ( elim_row == pivot_row ) continue;

      MatrixPosition q1 = p;
      Z -> column_advance ( q1 );
      //std::cout << "ZigZagFormSubroutine. TOP OF ELIMINATION LOOP. \n";
      //std::cout << "   p = " << p << " = (" << Z -> row ( p ) << ", " 
      //          << Z -> column ( p ) << ").\n";
      //std::cout << "   Column of interest = " << j << "\n"
      //          << "   Pivot Row = " << pivot_row << "\n"
      //          << "   Pivot_value = " << Z -> read ( pivot_row, j ) << "\n";
      //if ( q1 != Z -> end () ) {
      //   std::cout << "   Next in line before pivot: q1 = " << q1 << "\n";
      // } else {
      //   std::cout << "   Next in line before pivot: q1 =  END\n";
      // }
      // if ( j != Z -> column ( p ) ) {
      //   throw std::logic_error ( "ZigZagFormSubroutine. Error. Fell off column.\n");
      // }
      
      // std::cout << "   --->Pivoting. Eliminating val = " << val << " on elim row " << elim_row << " with pivot_row = " << pivot_row << "\n";
      Z -> similarity_operation ( R(1), 0,
                                   -val, R(1), 
                                   pivot_row, elim_row );
      // std::cout << "   p = " << p << " = (" << Z -> row ( p ) << ", " 
      //           << Z -> column ( p ) << ").\n";
      //sane_matrix ( *Z );
      /*
      MatrixPosition q2 = p;
      Z -> column_advance ( q2 );
      if ( q2 != Z -> end () ) {
        std::cout << "   Next in line after pivot: q2 = " << q2 << "\n";
      } else {
        std::cout << "   Next in line before pivot: q2 =  END\n";
      }
      if ( q != q2 ) {
        throw std::logic_error ( "ZigZagFormSubroutine. Error. Iteration structure damaged.\n");
      }
      */
      // if ( Z -> read ( pivot_row, j ) != R(1) ) {
      //   throw std::logic_error ( "ZigZagFormSubroutine. Error. Pivot element clobbered.\n");
      // }
      // if ( Z -> read ( elim_row, j ) != R(0) ) {
      //   throw std::logic_error ( "ZigZagFormSubroutine. Error. Pivot operation failed.\n");
      // }
    }
    ++ j;
  }

  if ( j == n ) {
    r = n; 
    return;
  }

  //std::cout << "Stage 1 complete. j0 = " << j0 << " and j = " << j << "\n";
  //print_matrix ( *Z );

  //    Stage 2.
  //std::cout << "ZigZagFormSubroutine. Stage 2.\n";
  for ( int64_t pivot_column = j-1; pivot_column >= j0; -- pivot_column ) {
    // Use column i to eliminate entries in columns [j+1,n) 
    //  (i.e. column i is all zeros except a 1 at row i+1)

    if ( Z -> read ( pivot_column+1, pivot_column ) != R(1) ) {
      throw std::logic_error ( "ZigZagFormSubroutine. Bad pivot choice. "
                               " (unexpected; please report a bug)\n" );
    }
    //for ( MatrixPosition p = Z -> row_begin ( i+1 );
    //      p != Z -> end ();
    //      Z -> row_advance ( p ) ) {

    q = Z -> row_begin ( pivot_column+1 );
    while ( q != Z -> end () ) {
      MatrixPosition p = q;
      Z -> row_advance ( q );


      int64_t elim_column = Z -> column ( p );
      if ( elim_column <= j ) continue;
      R val = Z -> read ( p );
      Z -> similarity_operation ( R(1), val, 
                                  R(0), R(1), 
                                  pivot_column, elim_column  );
      if ( Z -> read ( pivot_column+1, elim_column ) != 0 ) {
        throw std::logic_error ( "ZigZagFormSubroutine. Stage 2 elimination failed"
                                 " (unexpected; please report a bug)\n" );
      }
      //print_matrix ( *Z );
      //char c;
      //std::cin >> c;
    }
  }

  //std::cout << "Stage 2 complete.\n";
  //print_matrix ( *Z );
  //    Stage 3.
  //std::cout << "ZigZagFormSubroutine. Stage 3.\n";
  // Find a column pivot_column in [j+1,n) so that
  // Z(j0, pivot_column) != 0 
  R pivot_value ( 0 );
  int64_t pivot_column;
  //std::cout << " j = " << j << " and n = " << n << "\n";
  //for ( MatrixPosition p = Z -> row_begin ( j0 );
  //      p != Z -> end ();
  //      Z -> row_advance ( p ) ) {

  q = Z -> row_begin ( j0 );
  while ( q != Z -> end () ) {
    MatrixPosition p = q;
    Z -> row_advance ( q );
    
    pivot_column = Z -> column ( p );
    if ( pivot_column <= j ) continue;
    pivot_value = Z -> read ( p );
    if ( pivot_value != R(0) ) break;
  }

  if ( pivot_value != R(0) ) {
    //std::cout << "s3 pivot_value = " << pivot_value << "\n";
    // Perform similarity operation so pivot_column is j+1
    // and has pivot_value 1
    Z -> similarity_operation ( R(0), R(1),  
                                R(1), R(0),
                                pivot_column, j+1 );
    Z -> similarity_operation ( pivot_value, j+1 );
    // Now we can pivot off (j0,j+1) to clear the j0 row for j_elim > j+1
    pivot_column = j+1;
    // note: pivot_value = R(1);
    // Elimination
    //for ( MatrixPosition p = Z -> row_begin ( j0 );
    //  p != Z -> end ();
    //  Z -> row_advance ( p ) ) {
    q = Z -> row_begin ( j0 );
    while ( q != Z -> end () ) {
      MatrixPosition p = q;
      Z -> row_advance ( q );


      int64_t elim_column = Z -> column ( p );
      if ( elim_column <= j+1 ) continue;
      R val = Z -> read ( p );
      Z -> similarity_operation ( R(1), val,  
                                  R(0), R(1),
                                 pivot_column, elim_column );     
    }  
  }
  r = j+1;
  //std::cout << "Stage 3 complete.\n";
  //print_matrix ( *Z );
}
/// ZigZagForm
/// Input: a square matrix A
/// Output: a square matrix Z in the zig-zag form 
///         described in [Storjohann98]
template < class R > void
ZigZagForm (SparseMatrix<R> * Z,
            std::vector<int64_t> * blocks,
            SparseMatrix<R> const & A ) {
  //std::cout << "ZigZagForm.\n";

  int64_t n = A . number_of_rows ();
  int64_t m = A . number_of_columns ();
  typedef typename SparseMatrix<R>::MatrixPosition MatrixPosition;
  if ( n != m ) {
    throw std::domain_error ( "ZigZagForm: input must "
                              "be a square matrix.");
  }
  *Z = A;

  //print_matrix ( *Z );

  TransposeWrapper<R> ZT ( Z );

  int64_t r = 0;
  while ( r < n ) {
    //std::cout << "ZigZagForm. r = " << r << ", n = " << n << "\n";
    // Perform algorithm of Lemma 1 in
    // [Storjohann98] to block [r,n)x[r,n)
    //print_matrix ( *Z );
    ZigZagFormSubroutine ( Z, r );
    blocks -> push_back ( r );

    if ( r == n ) break;
    //print_matrix ( *Z );

    //std::cout << "Transpose step:\n";
    // Perform transpose of algorithm of Lemma 1 in
    // [Storjohann98] to block [r,n]x[r,n]
    ZigZagFormSubroutine ( &ZT, r );
    blocks -> push_back ( r );
  }

}

} //namespace frobenius_detail

template < class R > 
std::vector<PolyRing<R> >
FrobeniusNormalForm (SparseMatrix<R> const & A ) {
  using namespace frobenius_detail;

  //std::cout << "FrobeniusNormalForm\n";
  typedef chomp::SparseMatrix < chomp::PolyRing < chomp::Ring > > PolyMatrix;
  typedef chomp::PolyRing < chomp::Ring > Polynomial;

  std::vector<int64_t> blocks;
  SparseMatrix<R> Z;
  ZigZagForm ( &Z, &blocks, A );

  //std::cout << "ZigZagForm found.\n";
/*
e-mail from Arne Storjohann  --

The recipe for transforming from Zigzag to Frobenius is not elaborated
on in the paper, so here is some pseudocode.  We are starting with 
an n x n matrix in Zigzag form as shown in (1) in (*).  Then 1<=k<=n, 
where k is the number of diagonal companion blocks.  The first phase is
to parse the n x n Zigzag form into two arrays of length k.

Initialize an array C[1..k] of monic polynomials corresponding to
the diagonal (maybe transposed) companion blocks.  Don't forget to
negate the coefficients: see the picture of C_g on page 1 of the
paper.

*/
  //print_matrix ( Z );

  //std::cout << "Building monic polynomials C.\n";

  int64_t k = blocks . size ();

  //std::cout << "There were " << k << " blocks.\n";

  std::vector < Polynomial > C;

  for ( int64_t i = 0; i < k; ++ i ) {
    //std::cout << " Top of C loop. i = " << i << "\n";
    int64_t r0 = (i>0) ? blocks[i-1] : 0;
    int64_t r1 = blocks [i];
    Polynomial p;
    int64_t degree = r1 - r0; 
    p . resize ( degree + 1 );
    //std::cout << "Middle of C loop. r0 = " << r0 << " and r1 = " << r1 << "\n";
    if ( i % 2 == 0 ) {
      for ( int64_t d = 0; d < degree; ++ d ) {
        //std::cout << "Top of C d loop, d = " << d << "\n";
        p [ d ] = - Z . read ( r0 + d, r1 - 1 );
      }
    } else {
      for ( int64_t d = 0; d < degree; ++ d ) {
        //std::cout << "Top of C d loop, d = " << d << "\n";
        p [ d ] = - Z . read ( r1 - 1, r0 + d );
      }
    }
    p [ degree ] = R(1);
    C . push_back ( p );
    //std::cout << "Pushing poly " << p << " onto C.\n";
  }

/*  (email continued) Initialize an array B[1..k] of booleans corresponding to the
off-diagonal blocks in the Zigzag form.  More precisely,
B[1],B[2],...,B[k] are set to true/false depending on whether the
off-diagonal blocks contain a one/zero: B[1]=true iff the block to
the right of C[1] contains a one; B[2]=true iff the block below
C[2] contains a one; B[3]=true iff the block to the right of C[3]
contains a one, etc.  Note that if k is odd then B[k]=false by
default.   */

  //std::cout << "Building boolean array B.\n";

  std::vector<bool> B;
  for ( int64_t i = 0; i < k-1; ++ i ) {
    int64_t r0 = (i>0) ? blocks[i-1] : 0;
    int64_t r1 = blocks [i];
    if ( i % 2 == 0 ) {
      // Check to right
      B . push_back ( ( Z . read ( r0, r1 ) == R(1) ) );
    } else {
      // Check below
      B . push_back ( ( Z . read ( r1, r0 ) == R(1) ) );
    }
  }
  B . push_back ( false );
/* (email continued) Now execute the following code, where V[1..k] is a array
of type polynomial.

 if B[1] then V[1] <- 1 else V[1] <- 0 end if
 for j = 2 to k {
     V[j-1] <- Gcd(V[j-1], C[j], C[j-1])
     for i = j-2 by -1 to 1 { V[i] <- Gcd(V[i], V[i+1], C[i]) }
     m <- 1
     for i to j-1 {           <MISSING FIRST BOUND -- assuming 1 -- SRH>
         g <- Gcd(Product(m, V[i]), C[i])
         q <- Quotient(C[i], g)
         C[i] <- g
         if B[j] then V[i] <- Remainder(m, C[i]) else V[i] <- 0 end if
         m <- Product(m, q)
     }
     C[j] <- Product(m, C[j])
     if B[j] then V[j] <- m else V[j] <-  0 end if
  }
*/

  //std::cout << "Determining invariant factors.\n";

  std::vector<Polynomial> V ( k );
  //std::cout << "k = " << k << "\n";
  //std::cout << "V.size() = " << V . size () << "\n";
  //std::cout << "B.size() = " << B . size () << "\n";

  if ( k > 0 ) V[0] = B[0] ? Polynomial ( 1 ) : Polynomial ( 0 );
  for ( int64_t j = 1; j < k; ++ j ) {
    //std::cout << "Top of loop 1. j = " << j << "\n";
    V[j-1] = gcd ( V[j-1], gcd ( C[j], C[j-1] ) );
    V[j-1] =  V[j-1] / Polynomial ( V[j-1][V[j-1].degree()] ); // gcd not returning monic.
    for ( int64_t i = j-2; i >= 0; -- i ) {
      //std::cout << "Top of loop 2. i = " << i << "\n";
      V[i] = gcd ( V[i], gcd ( V[i+1], C[i] ) );
      V[i] =  V[i] / Polynomial ( V[i][V[i].degree()] ); // gcd not returning monic.
    }
    Polynomial m ( 1 );
    for ( int64_t i = 0; i < j; ++ i ) {
      //std::cout << "Top of loop 3. i = " << i << "\n";
      Polynomial g = gcd ( m * V[i], C[i] );
      g = g / Polynomial ( g[g.degree () ] ); // gcd not returning monic
      Polynomial q = C[i] / g; 
      C[i] = g;
      V[i] = B[j] ? m - (m / C[i]) * C[i] : Polynomial ( 0 );
      m = m * q;
    }
    C[j] = m * C[j];
    V[j] = B[j] ? m : Polynomial ( 0 );
  }

/* (email continued) After the code completes the array C contains the 
invariant factors of the Zigzag form (i.e., C[i] divides C[i+1] for 1<=i<k).  
Note that some of the C[i] might be equal to 1, in which case they correspond to
a companion block of dimension 0. */

  // Well I better return them then, shouldn't I? :)
  //std::cout << "Frobenius Form Returning.\n";
  return C;
}


}

#endif
