// SmithNormalForm.h
// Shaun Harker
// 9/14/11

#ifndef CHOMP_SMITHNORMALFORM_H
#define CHOMP_SMITHNORMALFORM_H

#include <cstdlib>
#include <stdint.h>

#include "chomp/SparseMatrix.h"

namespace chomp {

/// SmithNormalForm
/// Input: A
/// Produce matrices U, Uinv, V, Vinv, and D
/// such that U D V = A
/// U and V are Z-invertible with inverses Uinv and Vinv
/// D is diagonal, such that d_i | d_{i+1}
template < class R > void
SmithNormalForm (SparseMatrix<R> * U,
                 SparseMatrix<R> * Uinv,
                 SparseMatrix<R> * V,
                 SparseMatrix<R> * Vinv,
                 SparseMatrix<R> * D,
                 const SparseMatrix<R> & A);



/*****************************
 *        DEFINITIONS        *
 *****************************/

/* SMITH FORM CALCULATION */


template < class R >
R div ( const R y, const R x ) {
  /* For some reason, chip-makers decided to make (-2^63) / -1 generate an exception
   rather than returning 2^63 (which overflows into -2^63 ).
   This isn't comprehensible, since all other operations support overflow in the expected
   manner of treating the register as Z_{2^64} . But it means we have to check for it. 
   It is especially disconcerting, since multiplication does not have the same problem.
   However there is a simple workaround. */
  if ( x != R ( -1 ) ) {
    ////NSLog ( @" div 1\n" );
    return y / x;
  } else {
    ////NSLog ( @" div 2\n" );
    return y * x;
  }
}

template < class R >
bool divisable ( const R y, const R x ) {
  if ( x == R ( 0 ) ) return ( y == R ( 0 ) );
  ////NSLog ( @" divisable 1\n" );
  if ( (div(y, x) * x) == y ) return true;
  ////NSLog ( @" divisable 2\n" );
  return false;
}

/// Bezout
/// Given a and b, calculate a, b, and gcd such that
/// s * a + t * b = gcd
/// and gcd is the greatest common divisor of a and b
/// Also, if a = +/- gcd, then choose t = 0.
template < class R >
void Bezout (R * s_out,
             R * t_out,
             R * gcd_out,
             const R & a, 
             const R & b) {
  //std::cout << "Bezout\n";
  static R one ( 1 );
  static R zero ( 0 );
  static R neg_one ( -1 );
  // Relies upon Euclidean Division being available, and also comparison.
  // The user doesn't need to sort the input.
  bool reversed = false;
  if ( a < b ) reversed = true;
  // For the proof to follow, we assume without loss a >= b
  R x = std::min ( a, b ); // Think x = b for proof below
  R y = std::max ( a, b ); // Think y = a for proof below
                           // For extended euclidean algorithm
  R s0 = one; R s1 = zero;
  R t0 = zero; R t1 = one;
  while ( x != zero ) {
    //std::cout << "top of bezout\n";
    //std::cout << " y = " << y << " and x = " << x << "\n";
    //std::cout << "y - 1 = " << y - 1 << "\n";
    R q = div ( y, x );
    R r = y - x * q;
    R s = s0 - q * s1;
    R t = t0 - q * t1;
    //std::cout << "Claim: " << r << " = " << s << " * " << a << " + " << t << " * " << b << "\n";
    /* r = s * a + t * b
     Proof:
     By induction. First, k = 0 case.
     Since s_0 = 1 and t_0 = -q_0, we get
     r_0 = a - q_0 * b, which is indeed the case.
     Now the k = 1 case.
     s_1 = 0 - q_1 * (1)
     t_1 = 1 - q_1 * (-q_0)
     s_1 * a + t_1 * b = 
     -q_1 * a + b + q_0 * q_1 * b =
     b - q_1 ( a - q_0 * b ) =
     b - q_1 * r_0 = r_1.
     Now the inductive step for k > 1. 
     s_k * a + t_k * b =
     (s_{k-2} - q_k * s_{k-1}) * a + (t_{k-2} - q_k * t_{k-1}) * b =
     (s_{k-2} * a + t_{k-2} * b ) - q_k * (s_{k-1} * a + t_{k-1} * b ) =
     r_{k-2} - q_k * r_{k-1} = r_k.
     */
    s0 = s1; s1 = s;
    t0 = t1; t1 = t;
    y = x;
    x = r;
  } /* while */
  // Set output
  // The Bezout coefficients s and t are the second to last ones
  // to be calculated (the last ones give 0 = s*a + t*b)
  if ( not reversed ) {
    * s_out = s0;
    * t_out = t0;
  } else {
    * t_out = s0;
    * s_out = t0;
  }
  * gcd_out = y;
  // TODO: generalize this to all unit multiples, somehow, for general Rs
  // TODO: should this really be here? perhaps the caller should worry about this
  if ( * gcd_out == a ) {
    * s_out = one; * t_out = zero;
  }
  if ( * gcd_out == -a ) {
    * s_out = neg_one; * t_out = zero;
  }
}



#ifdef USE_GMP
#include <gmpxx.h>
template < >
void Bezout < mpz_class > ( mpz_class * s_out,
                           mpz_class * t_out,
                           mpz_class * gcd_out,
                           const mpz_class & a, 
                           const mpz_class & b) {
  mpz_gcdext ( gcd_out -> get_mpz_t (), 
              s_out -> get_mpz_t (), 
              t_out -> get_mpz_t (), 
              a . get_mpz_t (), 
              b . get_mpz_t () );
  // TODO: should this really be here? perhaps the caller should worry about this
  static mpz_class one ( 1 );
  static mpz_class zero ( 0 );
  static mpz_class neg_one ( -1 );
  if ( * gcd_out == a ) {
    * s_out = one; * t_out = zero;
  }
  if ( * gcd_out == -a ) {
    * s_out = neg_one; * t_out = zero;
  }
}
#endif

/**********************
 *      PIVOTING      *
 **********************/

/// RowPivot
/// Using the pivot element (i, j), clear the jth column using row operations
template < class R >
void RowPivot (SparseMatrix<R> * U,
               SparseMatrix<R> * Uinv,
               SparseMatrix<R> * D,
               const int i,
               const int j) {
  
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::Index Index;
  // Main Loop
  Index pivot_index = D -> find ( i, j );
  Index index = D -> column_begin ( j );
  while ( index != D -> end () ) { 
    // Perform the pivot operation. Use (i, j) to eliminate element.
    R s, t, g, x, y;
    R a = D -> read ( pivot_index );
    R b = D -> read ( index ); // Remains valid.
    int k = D -> row ( index );
    D -> column_advance ( index );
    
    //std::cout << " a = " << a << " b = " << b << " k = " << k << " and j = " << j << "\n";
    //std::cout << " (column advanced) index = " << index << "\n";
    
    // Prevent self-elimination.
    if ( i == k ) continue;
    // Determine necessary row operations with Euclidean Algorithm
    Bezout ( &s, &t, &g, a, b );
    x = div ( a, g );
    y = div ( b, g );
    
    
    
    
    
    // DEBUG
    /*
     std::cout << " Value at pivot = " << a << "\n";
     std::cout << " Elimination value = " << b << "\n";
     std::cout << " Elimination index = " << k << "\n";
     std::cout << " Bezout Formula: " << s << " * " << a << " + " << t << " * " << b << " = " << g << "\n";
     */
    
    /* Explanation of the row operations:
     Apply the 2x2 matrix from the left:
     M=  [  s  t  ]     Minv = [ x  -t ]
     [ -y  x  ]            [ y   s ]
     This means the following: Let I and K represent the ith and kth rows, respectively.
     We do the following:  I' <- s I + t K
     K' <- -y I + x K.
     See that this makes a <- s * a + t * b = g
     and b <- (- b * a + a * b) / g = 0, as desired.
     Also, we have to update U and U_inv. 
     U is to be updated by multiplying on the right by M_inv.
     This involves column operations. Let I and K be the ith and kth columns of U, respectively.
     We set I' <- x * I + y * K
     K' <- -t * I + s * K
     Uinv is to be updated by multiplying on the left by M. This is the same as what we did to D.
     */
    
    D -> row_operation (s, t,
                        -y, x,
                        i, k);
    
    
    Uinv -> row_operation (s, t,
                           -y, x,
                           i, k);
    
    
    U -> column_operation (x, y,
                           -t, s,
                           i, k);
    //std::cout << "Result of elimination step:\n";
    //print_matrix ( *D );
  } /* boost_foreach */
}

/// ColumnPivot
/// Using the pivot element (i, j), clear the ith row using column operations
template < class R >
void ColumnPivot (SparseMatrix<R> * V,
                  SparseMatrix<R> * Vinv,
                  SparseMatrix<R> * D,
                  const int i,
                  const int j) {
  
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::Index Index;
  // Main Loop
  Index pivot_index = D -> find ( i, j );
  Index index = D -> row_begin ( i );
  
  //std::cout << "DEBUG: ColumnPivot. pivot_index = " << pivot_index << " and index = " << index << "\n";
  //std::cout << "DEBUG: D -> end () = " << D -> end () << "\n";
  
  while ( index != D -> end () ) { 
    //std::cout << "Inside while loop...\n";
    // Perform the pivot operation. Use (i, j) to eliminate element.
    R s, t, g, x, y;
    R a = D -> read ( pivot_index );
    R b = D -> read ( index ); // Remains valid.
    int k = D -> column ( index );
    //std::cout << " pivot_index = " << pivot_index << ", index = " << index << "\n";
    // Prevent self-elimination.
    D -> row_advance ( index );
    if ( j == k ) continue;
    //std::cout << " Eliminating (" << i << ", " << k << "; " << b << ") with (" << i << ", " << j << "; " << a << ")\n";
    // Determine necessary row operations with Euclidean Algorithm
    Bezout ( &s, &t, &g, a, b );
    x = div ( a, g );
    y = div ( b, g );
    
    //std::cout << " Bezout Formula: " << s << " * " << a << " + " << t << " * " << b << " = " << g << "\n";
    
    /* Explanation of the row operations:
     Apply the 2x2 matrix from the right:
     M=  [  s  -y  ]     Minv = [  x  y ]
     [  t   x  ]            [ -t  s ]
     This means the following: Let J and K represent the jth and kth rows, respectively.
     We do the following:  J' <-  s J + t K
     K' <-  -y J + x K.
     See that this makes a <- s * a + t * b = g
     and b <- (- b * a + a * b) / g = 0, as desired.
     Also, we have to update U and U_inv. 
     U is to be updated by multiplying on the left by M_inv.
     This involves column operations. Let J and K be the jth and kth columns of U, respectively.
     We set J' <-  x * J + y * K
     K' <-  -t * J + s * K
     Uinv is to be updated by multiplying on the left by M. This is the same as what we did to D.
     */
    
    D -> column_operation (s, t,
                           -y, x,
                           j, k);
    
    Vinv -> column_operation (s, t,
                              -y, x,
                              j, k);
    
    V -> row_operation (x, y,
                        -t, s,
                        j, k);
    
    
    //std::cout << "Result of elimination step:\n";
    //print_matrix ( *D );
  } /* while */
  
}

/// DoublePivot
/// Using the pivot element (i, j), clear the ith row 
/// and jth column using row and column operations
template < class R >
void SmithPivot (SparseMatrix<R> * U,
                 SparseMatrix<R> * Uinv,
                 SparseMatrix<R> * V,
                 SparseMatrix<R> * Vinv,
                 SparseMatrix<R> * D,
                 const int i,
                 const int j) {
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::Index Index;
  // Obtain the (i, j)th element of D
  Index pivot_index = D -> find ( i, j );
  R pivot_value = D -> read ( pivot_index ); 
  
  //DEBUG
  if ( D -> read (i, j) == R ( 0 ) ) {
    std::cout << "Pivot cannot be zero!\n";
    exit ( 1 );
  }
  
  
  R old_pivot_value ( 0 );
  //std::cout << " **** SMITH PIVOT (" << i << ", " <<
  // j << ") value = " << pivot_value << " **** \n";
  //print_matrix ( * D );
  // We assume pivot_value != 0, so this while loop will run at least once:
  while ( pivot_value != old_pivot_value ) {
    old_pivot_value = pivot_value;
    
    //sane_matrix ( *D );
    //std::cout << "***** Calling ColumnPivot on (" << i << ", " << j << ")\n";
    ColumnPivot ( V, Vinv, D, i, j);
    //print_matrix ( *D );
    
    //sane_matrix ( *D );
    //std::cout << "***** Calling RowPivot on (" << i << ", " << j << ")\n";
    RowPivot ( U, Uinv, D, i, j );
    //print_matrix ( *D );
    
    pivot_value = D -> read ( pivot_index ); 
    // if pivot_value has not changed, then no new work has been produced.
  }
  //std::cout << " **** SMITH PIVOT COMPLETE **** \n";
}

#ifdef SNF_DEBUG
uint64_t number_of_pivots = 0;
#endif

/// SmithNormalForm
/// Input: A
/// Produce matrices U, Uinv, V, Vinv, and D
/// such that U D V = A
/// U and V are Z-invertible with inverses Uinv and Vinv
/// D is diagonal, such that d_i | d_{i+1}
/// WARNING: output matrices MUST be empty.
template < class R >
void SmithNormalForm (SparseMatrix<R> * U,
                      SparseMatrix<R> * Uinv,
                      SparseMatrix<R> * V,
                      SparseMatrix<R> * Vinv,
                      SparseMatrix<R> * D,
                      const SparseMatrix<R> & A) {
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::Index Index;
  typedef typename Matrix::size_type size_type;
  // We copy A into D.
  *D = A;
  
  //NSLog(@"SNF1\n" );

  //std::cout << "SmithNormalForm Calculation.\n";
  //std::cout << "   this is " << A . number_of_rows () 
  // << " x " << A . number_of_columns () << "\n";
  //print_matrix ( *D );
  //sane_matrix ( *D );
  // Resize outputs, set to identity
  
  U -> resize ( A . number_of_rows (), A . number_of_rows () );
  for ( int i = 0; i < A . number_of_rows (); ++ i ) 
    U -> write ( i, i, R ( 1 ), true ); 
  *Uinv = *U;
    
  V -> resize ( A . number_of_columns (), A . number_of_columns () );
  for ( int i = 0; i < A . number_of_columns (); ++ i ) 
    V -> write ( i, i, R ( 1 ), true ); 
  *Vinv = *V;
  
  // The algorithm proceeds by selecting a sequence of pivots (t, j_t)
  int t = 0;
#ifdef SNF_DEBUG
  std::cout << "**** ELIMINATION STAGE ****\n";
#endif
  
  //NSLog(@"SNF2\n" );

  for ( int j = 0; j < D -> number_of_columns (); ++ j ) {
#ifdef SNF_DEBUG
    std::cout << "Working on column " << j << 
    " out of " << D -> number_of_columns () << ".\n";
#endif

    if ( D -> column_size ( j ) == 0 ) continue;
    
    // We want to find a good pivot row from the jth column
    // The simplest criterion is to pick the one with the smallest size
    size_type pivot_row = 0; // we find this
    size_type best_size = D -> number_of_columns (); // maximum size row could be
    Index index = D -> column_begin ( j );
    while ( index != D -> end () ) {
      size_type which_row = D -> row ( index );
      size_type row_size = D -> row_size ( which_row );
      //std::cout << "Looking at index " << index << 
      //", we see that row " << which_row << " has size " << row_size << "\n";
      if ( row_size <= best_size ) {
        best_size = row_size;
        pivot_row = which_row;
      } /* if */
      D -> column_advance ( index );
    } /* while */
    // At this point we have found a pivot. We'll need to swap rows.
    
    //std::cout << " Just before swapping:\n";
    //print_matrix ( *D );
    //sane_matrix ( *D );
    
    //std::cout << "Swapping rows. (" << t << ", " << pivot_row << ")...\n";
    D -> swap_rows ( t, pivot_row );
    U -> swap_columns ( t, pivot_row );
    Uinv -> swap_rows ( t, pivot_row );
    
    
    //print_matrix ( *D );
    //sane_matrix ( *D );
    // Now the pivot_row is really t
    
    // Now we use pivot off from the pivot choice,
    // zeroing out all elements in its row and column except itself.
    // At the end of this procedure, all remaining entries have
    // row number greater than t and column number greater than j
    //std::cout << "Performing the pivot step at (" << t << ", " << j << "):\n";
    SmithPivot ( U, Uinv, V, Vinv, D, t, j );
    //print_matrix ( *D );
    //std::cout << "Swapping columns.\n";
    // Now we swap columns, making this the t-th column rather than the jth column
    // (note that j >= t )
    D -> swap_columns ( t, j );
    V -> swap_rows ( t, j );
    Vinv -> swap_columns ( t, j );
    
    // Finally, increment t.
    ++ t;
  } /* for */
  // At this stage, we are almost done. What remains is to ensure the
  // divisability d_i | d_{i+1}.
  // The following trick accomplishes this: one by one, we 
  // copy the value of the diagonals (except the first) to the position
  // immediately to the left
  
#ifdef SNF_DEBUG
  std::cout << "*** DIAGONAL DIVISABILITY STAGE ****\n";
#endif
  //NSLog(@"SNF3\n" );

  //print_matrix ( *D );
  
  int r = std::min ( D -> number_of_rows (), D -> number_of_columns () );
  bool pass_again = true;
  while ( pass_again ) {
    //NSLog(@"SNF - PASS\n" );
    pass_again = false;
    for ( int i = 0; i < r - 1; ++ i ) {
      ////NSLog(@"SNF - PASS -- i = %d / r = %d\n", i, r );
#ifdef SNF_DEBUG
      std::cout << " i = " << i << " and r = " << r << "\n";
      std::cout << " size of almost diagonal matrix D: " << D -> size () << "\n";
      std::cout << " consider " << D -> read ( i, i ) << " and " << D -> read ( i + 1, i + 1 ) << "\n";
      if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) 
        std::cout << "we consider them divisible.\n";
      else std::cout << "we don't consider them divisible\n";
#endif
      //print_matrix ( *D );
      
      ///// DEBUG ////
      //std::cout << "Check:" << D -> read ( i + 1, i + 1 ) << " is divisible by " <<
      //D -> read ( i, i ) << "\n";
      //if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) {
      //  std::cout << "YES.\n";
      //};
      /////////////////
      
      ////NSLog(@"SNF - PASS -- A\n" );

      ////// DEBUG //////
      //R temp1 = D -> read ( i + 1, i + 1 );
      //R temp2 = D -> read ( i, i );
      //std::cout << "temp1 = " << temp1 << "\n";
      //std::cout << "temp2 = " << temp2 << "\n";
      ////NSLog( @"SNF - PASS - A2\n");
            
      if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) continue;
      
      ////NSLog(@"SNF - PASS -- B\n" );

      pass_again = true;
      //////// DEBUG /////////
      //std::cout << "NO.\n";
      ////////////////////////
      
      D -> write ( i + 1, i,  D -> read ( i + 1, i + 1 ) );
      
      //print_matrix ( *D );
      
      // This is a column operation; we need to update V and Vinv
      /* The matrix in question is
       M =  [ 1 0 ]    Minv = [  1 0 ]
       [ 1 1 ]           [ -1 1 ]
       We should multiply V on the left by Minv and Vinv on the right by M.
       */
      ////NSLog(@"SNF - PASS -- C\n" );

      Vinv -> column_operation (R ( 1 ), R ( 1 ),
                                R ( 0 ), R ( 1 ),
                                i, i + 1 );
      
      ////NSLog(@"SNF - PASS -- D\n" );

      V -> row_operation (R ( 1 ), R ( 0 ),
                          R ( -1 ), R ( 1 ),
                          i, i + 1 );
      
      ////NSLog(@"SNF - PASS -- E\n" );

      SmithPivot ( U, Uinv, V, Vinv, D, i, i );
      
      ////NSLog(@"SNF - PASS -- F\n" );

      ////// DEBUG ///////
      // std::cout << "Check again:" << D -> read ( i + 1, i + 1 ) << " is divisible by " <<
      // D -> read ( i, i ) << "\n";    
      // if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) {
      //   std::cout << "YES.\n";
      // };
      ////////////////////
      
      //print_matrix ( *D );
      
    }
  }
  
  //NSLog(@"SNF4\n" );
  /*
   Not legitimate unless you also fix the transformation matrices 
   for ( int i = 0; i < r; ++ i ) {
   D -> write ( i, i, abs ( D -> read ( i, i ) ) );
   }
   */
  // DEBUG.
  SparseMatrix<R> UDV = (*U) * (*D) * (*V);
  
  //NSLog(@"SNF4.5\n" );

  //NSLog ( @"Size of UDV = (%d, %d)", UDV . number_of_rows (), UDV . number_of_columns () );
  //NSLog ( @"Size of A = (%d, %d)", A . number_of_rows (), A . number_of_columns () );
  
  for ( int i = 0; i < UDV . number_of_rows (); ++ i ) {
    for ( int j = 0; j < UDV . number_of_columns (); ++ j ) {
      ////NSLog(@"about to read\n" );
      if ( UDV . read (i, j) != A . read (i, j) ) {
        //NSLog(@"problem at (%d, %d)\n", i, j );
        std::cout << "Smith Normal Form failure.\n";
        print_matrix ( A );
        //print_matrix ( UDV );
        return;
        //exit ( 1 );
      }
      ////NSLog(@"did read\n" );
    }
  }
  //NSLog(@"SNF5\n" );
#ifdef SNF_DEBUG
  std::cout << "Done!\n";
  //std::cout << "Total number of pivot moves: " << number_of_pivots << "\n";
#endif
}
  
} // namespace chomp

#endif
