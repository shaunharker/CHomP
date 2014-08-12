// SmithNormalForm.h
// Shaun Harker
// 9/14/11

#ifndef CHOMP_SMITHNORMALFORM_H
#define CHOMP_SMITHNORMALFORM_H

#include <fstream>
#include <cstdlib>
#include <stdint.h>

#include "chomp/SparseMatrix.h"
#include "chomp/Algebra.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/thread.hpp>

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
                 const SparseMatrix<R> & A );



/*****************************
 *        DEFINITIONS        *
 *****************************/

/* SMITH FORM CALCULATION */



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
  typedef typename Matrix::MatrixPosition MatrixPosition;
  // Main Loop
  MatrixPosition pivot_index = D -> find ( i, j );
  MatrixPosition index = D -> column_begin ( j );
  while ( index != D -> end () ) { 

    boost::this_thread::interruption_point ();

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
    
    
    #ifdef SNF_DEBUG
    std::cout << "Bezout: (s, t, x, y, a, b ) = (" << s 
      << ", " << t << ", " << x << ", " << y << ", " << a << ", " << b << ")\n";

     std::cout << " Value at pivot = " << a << "\n";
     std::cout << " Elimination value = " << b << "\n";
     std::cout << " Elimination index = " << k << "\n";
     std::cout << " Bezout Formula: " << s << " * " << a << " + " << t << " * " << b << " = " << g << "\n";
     std::cout << s*a + t*b << "\n";
    #endif
    
    /* Explanation of the row operations:
     Apply the 2x2 matrix from the left:
     M :=  [  s  t  ]     
           [ -y  x  ]            
     This means the following: Let I and K represent the ith and kth rows, respectively.
     We do the following:  
         I' <-  s I + t K
         K' <- -y I + x K.
     See that this makes a <- s * a + t * b = g
     and b <- (- b * a + a * b) / g = 0, as desired.
     Also, we have to update U and U_inv. 
     U is to be updated by multiplying on the right by M_inv.
     Minv := [ x  -t ]
             [ y   s ]
     This involves column operations. Let I and K be the ith and kth columns of U, respectively.
     We set 
         I' <-  x * I + y * K
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
  typedef typename Matrix::MatrixPosition MatrixPosition;
  // Main Loop
  MatrixPosition pivot_index = D -> find ( i, j );
  MatrixPosition index = D -> row_begin ( i );
  
  //std::cout << "DEBUG: ColumnPivot. pivot_index = " << pivot_index << " and index = " << index << "\n";
  //std::cout << "DEBUG: D -> end () = " << D -> end () << "\n";
  
  while ( index != D -> end () ) { 

    boost::this_thread::interruption_point ();

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
    #ifdef SNF_DEBUG
    std::cout << " Eliminating (" << i << ", " << k << "; " << b << ") with (" << i << ", " << j << "; " << a << ")\n";
    #endif
    // Determine necessary row operations with Euclidean Algorithm
    Bezout ( &s, &t, &g, a, b );
    x = div ( a, g );
    y = div ( b, g );
    
    #ifdef SNF_DEBUG
    std::cout << "Bezout: (s, t, x, y, a, b ) = (" << s 
      << ", " << t << ", " << x << ", " << y << ", " << a << ", " << b << ")\n";

     std::cout << " Value at pivot = " << a << "\n";
     std::cout << " Elimination value = " << b << "\n";
     std::cout << " Elimination index = " << k << "\n";
     std::cout << " Bezout Formula: " << s << " * " << a << " + " << t << " * " << b << " = " << g << "\n";
     std::cout << s*a + t*b << "\n";
    #endif    
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

    #ifdef SNF_DEBUG
     std::cout << ".\n";
    #endif

     D -> column_operation (s, t,
       -y, x,
       j, k);
    #ifdef SNF_DEBUG

     std::cout << ".\n";
    #endif

     Vinv -> column_operation (s, t,
      -y, x,
      j, k);
    #ifdef SNF_DEBUG

     std::cout << "col_op.\n";
    #endif

     V -> row_operation (x, y,
      -t, s,
      j, k);
    #ifdef SNF_DEBUG
     std::cout << "column row_op.\n";
    #endif
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
  typedef typename Matrix::MatrixPosition MatrixPosition;
  // Obtain the (i, j)th element of D
  MatrixPosition pivot_index = D -> find ( i, j );

  //DEBUG
  if ( pivot_index == D -> end () ) {
    std::cout << "Unable to find pivot " << i << ", " << j << "\n";
    std::cout << "Hence pivot_value = " << D -> read ( pivot_index );
    std::cout << "THROWING BAD SNF CALCULATION\n";
    throw 42; // bad snf calculation
  }
  //END DEBUG

  //DEBUG
  if ( D -> read (i, j) == R ( 0 ) ) {
    std::cout << "Pivot cannot be zero!\n";
    exit ( 1 );
  }
  //END DEBUG
  
  #ifdef SNF_DEBUG
  R pivot_value = D -> read ( pivot_index ); 
  std::cout << " **** SMITH PIVOT (" << i << ", " <<
   j << ") value = " << pivot_value << " **** \n";
  print_matrix ( * D );
  #endif
  // We assume pivot_value != 0, so this while loop will run at least once:
  while ( D -> row_size ( i ) > 1 ||
          D -> column_size ( j ) > 1 ) { 
    
    //sane_matrix ( *D );
    
    #ifdef SNF_DEBUG
    std::cout << "***** Calling ColumnPivot on (" << i << ", " << j << ")\n";
    #endif
    
    ColumnPivot ( V, Vinv, D, i, j);

    #ifdef SNF_DEBUG
    print_matrix ( *D );
    #endif

    //sane_matrix ( *D );
    #ifdef SNF_DEBUG
    std::cout << "***** Calling RowPivot on (" << i << ", " << j << ")\n";
    #endif

    RowPivot ( U, Uinv, D, i, j );

    #ifdef SNF_DEBUG
    print_matrix ( *D );
    #endif

  }
  #ifdef SNF_DEBUG
  std::cout << " **** SMITH PIVOT COMPLETE **** \n";
  #endif
}

#ifdef SNF_DEBUG
uint64_t number_of_pivots = 0;
#endif

#define CHECK_SNF_PRODUCT { SparseMatrix<R> UDV = (*U) * (*D) * (*V); \
for ( int i = 0; i < UDV . number_of_rows (); ++ i ) { \
    for ( int j = 0; j < UDV . number_of_columns (); ++ j ) { \
      if ( UDV . read (i, j) != A . read (i, j) ) abort (); \
    } \
  } }

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
                      const SparseMatrix<R> & A ) {
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::MatrixPosition MatrixPosition;
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
  
  #ifdef SNF_DEBUG
  CHECK_SNF_PRODUCT
  #endif

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
    if ( D -> read ( t, t ) != R ( 0 ) ) {
      pivot_row = t;
    } else {
    size_type best_size = D -> number_of_columns (); // maximum size row could be
    MatrixPosition index = D -> column_begin ( j );
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
    }
    // At this point we have found a pivot. We'll need to swap rows.
    
    #ifdef SNF_DEBUG
    std::cout << "Pivot Row selected: t = " << t << ", pivot_row = " << pivot_row << "\n";
    std::cout << " Just before swapping:\n";
    print_matrix ( *D );
    //sane_matrix ( *D );
    std::cout << "Swapping rows. (" << t << ", " << pivot_row << ")\n";
    #endif

    D -> swap_rows ( t, pivot_row );
    U -> swap_columns ( t, pivot_row );
    Uinv -> swap_rows ( t, pivot_row );
    
    #ifdef SNF_DEBUG
    CHECK_SNF_PRODUCT
    #endif

    #ifdef SNF_DEBUG
    print_matrix ( *D );
    //sane_matrix ( *D );
    #endif
    // Now the pivot_row is really t
    
    // Now we use pivot off from the pivot choice,
    // zeroing out all elements in its row and column except itself.
    // At the end of this procedure, all remaining entries have
    // row number greater than t and column number greater than j
    #ifdef SNF_DEBUG
    std::cout << "Performing the pivot step at (" << t << ", " << j << "):\n";
    #endif
    SmithPivot ( U, Uinv, V, Vinv, D, t, j );

    #ifdef SNF_DEBUG
    CHECK_SNF_PRODUCT
    #endif

    #ifdef SNF_DEBUG
    print_matrix ( *D );
    std::cout << "Swapping columns. (" << t << ", " << j << ")\n";
    #endif
    // Now we swap columns, making this the t-th column rather than the jth column
    // (note that j >= t )
    D -> swap_columns ( t, j );
    V -> swap_rows ( t, j );
    Vinv -> swap_columns ( t, j );
    
    #ifdef SNF_DEBUG
    CHECK_SNF_PRODUCT
    #endif
    
    #ifdef SNF_DEBUG
    print_matrix ( *D );
    //sane_matrix ( *D );
    #endif

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
    #ifdef SNF_DEBUG
   std::cout << " TOP OF OUTER DIAGONAL DIVISABILITY STAGE LOOP \n ";
#endif
    pass_again = false;
    for ( int i = 0; i < r - 1; ++ i ) {
      ////NSLog(@"SNF - PASS -- i = %d / r = %d\n", i, r );
#ifdef SNF_DEBUG
      std::cout << " TOP OF INNER DIAGONAL DIVISABILITY STAGE LOOP \n ";
      print_matrix ( *D );
      std::cout << " i = " << i << " and r = " << r << "\n";
      std::cout << " size of D: " << D -> size () << "\n";
      std::cout << " We look at (" << i << ", " << i << ") and (" << i+1 << ", " << i+1 << ")\n";
      std::cout << " We consider " << D -> read ( i, i ) << " and " << D -> read ( i + 1, i + 1 ) << "\n";
      if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) 
        std::cout << D -> read ( i + 1, i + 1 ) << " is  divisable by " << D -> read ( i, i )  << " \n";
      else std::cout << D -> read ( i + 1, i + 1 ) << " is NOT divisable by " << D -> read ( i, i )  << " \n";
#endif
      
      if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) continue;
      
      pass_again = true;
      
      D -> write ( i + 1, i,  D -> read ( i + 1, i + 1 ) );

      // This is a column operation; we need to update V and Vinv
      /* The matrix in question is
       M =  [ 1 0 ]    Minv = [  1 0 ]
       [ 1 1 ]           [ -1 1 ]
       We should multiply V on the left by Minv and Vinv on the right by M.
       */

      Vinv -> column_operation (R ( 1 ), R ( 1 ),
                                R ( 0 ), R ( 1 ),
                                i, i + 1 );
      V -> row_operation (R ( 1 ), R ( 0 ),
                          R ( -1 ), R ( 1 ),
                          i, i + 1 );
      
      #ifdef SNF_DEBUG
      std::cout << "Applied an operation:\n";
      print_matrix ( *D );
      #endif

      #ifdef SNF_DEBUG
      CHECK_SNF_PRODUCT
      #endif
    
      SmithPivot ( U, Uinv, V, Vinv, D, i, i );
    
      #ifdef SNF_DEBUG
      CHECK_SNF_PRODUCT
      #endif
      ////// DEBUG ///////
      #ifdef SNF_DEBUG
            print_matrix ( *D );

       std::cout << "Check again:" << D -> read ( i + 1, i + 1 ) << " is divisible by " <<
       D -> read ( i, i ) << "\n";    
       if ( divisable ( D -> read ( i + 1, i + 1 ), D -> read ( i, i ) ) ) {
         std::cout << "YES.\n";
       } else {
        std::cout << "NO!!!!!\n";
        abort ();
       }
      print_matrix ( *D );

      #endif
      ////////////////////
      
      
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

  //NSLog(@"SNF4.5\n" );

  //NSLog ( @"Size of UDV = (%d, %d)", UDV . number_of_rows (), UDV . number_of_columns () );
  //NSLog ( @"Size of A = (%d, %d)", A . number_of_rows (), A . number_of_columns () );
  
  //#ifdef SNF_DEBUG
  { SparseMatrix<R> UDV = (*U) * (*D) * (*V);
  for ( int i = 0; i < UDV . number_of_rows (); ++ i ) {
    for ( int j = 0; j < UDV . number_of_columns (); ++ j ) {
      ////NSLog(@"about to read\n" );
      if ( UDV . read (i, j) != A . read (i, j) ) {
        //NSLog(@"problem at (%d, %d)\n", i, j );
        std::cout << "Smith Normal Form failure.\n";
        print_matrix ( A );
        print_matrix ( UDV );
        {
           std::cout << "SNF FAILURE SAVE\n";
           std::ofstream ofs("SNFFailure.txt");
          assert(ofs.good());
          boost::archive::text_oarchive oa(ofs);
          oa << A;
        }
        abort ();
        throw 42;
        return;
      }
    }
  } }
  //#endif

#ifdef SNF_DEBUG
  std::cout << "Done!\n";
  //std::cout << "Total number of pivot moves: " << number_of_pivots << "\n";
#endif
}
  
} // namespace chomp

#endif
