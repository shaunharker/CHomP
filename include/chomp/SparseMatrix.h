/*
 * SparseMatrix.h
 * chomp-rutgers
 * 4/27/11
 * Shaun Harker
 *
 * This header contains declarations for Sparse Matrices.
 * And definitions, as of 9/14/11
 */ 

#ifndef CHOMP_SPARSEMATRIX_H
#define CHOMP_SPARSEMATRIX_H

#include <iomanip>
#include <cstdlib>
#include <vector>
#include <deque>
#include <utility>
#include <algorithm>

#include "boost/foreach.hpp"
#include "boost/functional/hash.hpp"
#include "boost/unordered_map.hpp"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/deque.hpp"
#include "boost/serialization/unordered_map.hpp"
#include <boost/serialization/nvp.hpp>

#define HASH_SWITCH 10

namespace chomp {

//int number_of_pivots = 0;

// Sparse Matrix Data Structure

typedef std::pair < int64_t, int64_t > Position;

// Forward Declarations
template < class R > struct Element;
template < class R > class SparseMatrix;

// struct Element
template < class R >
struct Element {
  typedef int64_t size_type;
  Position position;
  R value;
  size_type left;
  size_type right;
  size_type up;
  size_type down;
  
  Element ( void ) {}
  
  Element ( Position position, R value, size_type left, size_type right, size_type up, size_type down) :
  position(position), value(value), left(left), right(right), up(up), down(down) {}
  
  template < class T > Element ( const T & copy_me ) : 
  position( copy_me . position), value( copy_me . value), left( copy_me . left), 
  right( copy_me . right), up( copy_me . up), down( copy_me . down) {}

  /// The serialization method.
  friend class boost::serialization::access;  
  template < class Archive >
  void serialize ( Archive & ar , const unsigned int version ) {
    ar & BOOST_SERIALIZATION_NVP(position);
    ar & BOOST_SERIALIZATION_NVP(value);
    ar & BOOST_SERIALIZATION_NVP(left);
    ar & BOOST_SERIALIZATION_NVP(right);
    ar & BOOST_SERIALIZATION_NVP(up);
    ar & BOOST_SERIALIZATION_NVP(down);
      //std::cout << "Serializing Element\n";
  }
};

// friends of SparseMatrix
// Multiply
template < class R >
SparseMatrix<R> operator * (const SparseMatrix<R> & rhs,
                               const SparseMatrix<R> & lhs);

// Add
template < class R >
SparseMatrix<R> operator + (const SparseMatrix<R> & rhs,
                               const SparseMatrix<R> & lhs);

// Subtract
template < class R >
SparseMatrix<R> operator - (const SparseMatrix<R> & rhs,
                               const SparseMatrix<R> & lhs);

// forward declare
template < class Matrix >
void print_matrix ( const Matrix & print_me );
    
// class Sparse Matrix
template < class R >
class SparseMatrix {
public:
  // typedefs
  typedef R value_type;
  typedef R entry_type;
  typedef int64_t MatrixPosition;
  typedef int64_t size_type;
private:
  

public:
  
  /// erase 
  ///   Set the entry corresponding to MatrixPosition "index"
  ///   to zero
  void erase ( const MatrixPosition index );

  /// find
  ///    Find the MatrixPosition corresponding to (i,j)
  ///    If there is no such MatrixPosition (i.e. A(i,j) == 0 )
  ///    then return end ()
  MatrixPosition find ( const size_type i, const size_type j ) const;
  
  /// read 
  ///     Return the value at the (i,j) position
  R read ( const size_type i, const size_type j ) const;

  /// read
  ///     Return the value at the position 
  ///     represented by MatrixPosition "index"
  R read ( MatrixPosition index ) const;

  /// write ( i, j, value, insert = false)
  ///   Set (i,j) entry to "value".
  ///     Note: If it is known in advance that (i,j) holds zero before the call, 
  ///           then "insert" may be set to true (it is false by default)
  ///           which results in a speedup.
  MatrixPosition 
  write ( const size_type i, const size_type j, 
          const R value, bool insert = false ); 
  
  /// write ( MatrixPosition index, const R value )
  ///   Copy "value" to the position represented by "index"
  MatrixPosition 
  write ( MatrixPosition index, const R value );
  
  // begin and end
  MatrixPosition row_begin ( const size_type i ) const;
  MatrixPosition column_begin ( const size_type j ) const;
  MatrixPosition end ( void ) const;
  
  // traversal. row_advance advances within a row, column_advance within a column
  void row_advance ( MatrixPosition & index ) const;
  void column_advance ( MatrixPosition & index ) const;
  
  // row and column (to learn position)
  size_type row ( const MatrixPosition index ) const;
  size_type column ( const MatrixPosition index ) const;
  
  // add to some position
  void add ( size_type i, size_type j, const R value );
  
  /// SparseMatrix
  ///   Default constructor 
  SparseMatrix ( void );

  /// SparseMatrix ( i, j )
  ///   Construct matrix with shape i x j
  SparseMatrix ( size_type i, size_type j );

  /// SparseMatrix
  ///    Copy constructor. Can convert types.
  template < class T > SparseMatrix ( const T & copy_me );
  
  /// resize
  ///    Change the shape of the matrix
  void resize ( size_type i, size_type j );

  /// size
  ///   Report number of non-zero elements in
  ///   the sparse matrix
  size_type size ( void ) const; // sparsity size
  
  /// row_operation
  ///   Multiply row i by a           
  void row_operation (const R a, size_type i );

  /// column_operation
  ///   Multiply column j by a
  void column_operation (const R a, size_type j );

  /// similarity_operation
  ///   Multiply row i by a and divide column i by a
  void similarity_operation (const R a, size_type i);

  /// row_operation
  ///   Multiply  on the left by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a b ]
  ///               [ c d ]             
  void row_operation (const R a, const R b,
                      const R c, const R d,
                      size_type i, size_type j );

  /// column_operation
  ///   Multiply  on the right by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a c ]
  ///               [ b d ]  
  void column_operation (const R a, const R b,
                         const R c, const R d,
                         size_type i, size_type j );

  /// similarity_operation
  ///   Multiply  on the left by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a b ]
  ///               [ c d ] 
  ///   and multiply on the right by the inverse of M
  void similarity_operation (const R a, const R b,
                             const R c, const R d,
                             size_type i, size_type j );

  /// swap_rows
  ///    Swap row i and row j
  void swap_rows ( const size_type i, const size_type j );

  /// swap_columns
  ///    Swap column i and column j
  void swap_columns ( const size_type i, const size_type j );

  /// number_of_rows
  ///    Return number of rows in matrix
  size_type number_of_rows ( void ) const;

  /// number_of_columns
  ///    Return number of columns in matrix
  size_type number_of_columns ( void ) const;

  /// row_size
  ///     Return number of non-zero elements in row i
  size_type row_size ( const size_type i ) const;

  /// columns_size
  ///     Return number of non-zero elements in column j  
  size_type column_size ( const size_type j ) const;
  

public: // not friends with different templated versions, weirdly
  // data to store the Sparse Matrix
  std::vector < Element<R> > data_;
  // Garbage Collection structure
  std::deque < MatrixPosition > garbage_;
  // data to assist in O(1) random access times
  boost::unordered_map < Position, MatrixPosition, boost::hash< Position > > access_;
  typedef boost::unordered_map < Position, MatrixPosition, boost::hash< Position > >::const_iterator access_iterator;
  // data to store the beginning of the rows and columns
  std::vector < MatrixPosition > row_begin_;
  std::vector < MatrixPosition > column_begin_;
  // data to store amount of non-zero elements per row and column
  std::vector < size_type > row_sizes_;
  std::vector < size_type > column_sizes_;
  // data to handle quick permutations of rows and columns
  std::vector < size_type > row_names_;
  std::vector < size_type > column_names_;
  // data to handle quick linear algebra
  std::vector < R > cache_A;
  std::vector < size_type > cache_A_TS;
  std::deque < size_type > cache_A_S;
  std::vector < R > cache_B;
  std::vector < size_type > cache_B_TS;
  std::deque < size_type > cache_B_S;
  size_type timestamp;
  // technicals
  MatrixPosition new_index ( void );
  void delete_index ( const MatrixPosition index );
  void hash_check ( const MatrixPosition index );
  void increment_row_size ( const size_type i );
  void increment_column_size ( const size_type j );
  void decrement_row_size ( const size_type i );
  void decrement_column_size ( const size_type j );

  // friends
  friend SparseMatrix<R> operator * <> ( const SparseMatrix<R> &, const SparseMatrix<R> & );
  friend SparseMatrix<R> operator + <> ( const SparseMatrix<R> &, const SparseMatrix<R> & );
  friend SparseMatrix<R> operator - <> ( const SparseMatrix<R> &, const SparseMatrix<R> & );
  
  /// The serialization method.
  friend class boost::serialization::access;  
  template < class Archive >
  void serialize ( Archive & ar , const unsigned int version ) {
    // TODO: More efficient to split into save/load routines, saving has unnecessary operations
      //std::cout << "Serializing SparseMatrix\n";
    std::vector < Element<R> > data_copy = data_;
    size_type number_of_rows = row_sizes_ . size ();
    size_type number_of_columns = column_sizes_ . size ();
    ar & boost::serialization::make_nvp("DATA",data_copy);
    ar & boost::serialization::make_nvp("NUMROWS",number_of_rows);
    ar & boost::serialization::make_nvp("NUMCOLS",number_of_columns);
    resize ( number_of_rows, number_of_columns );
    BOOST_FOREACH ( const Element<R> & e, data_copy ) {
      write ( e . position . first, e . position . second, e . value );
    }    
  }

/// sanityCheck
///   debugging code    
  void sanityCheck ( void ) const {
    std::cout << "sanityCheck.\n";
    std::cout << "data_.size() = " << data_ . size () << "\n";
    BOOST_FOREACH ( const Element<R> & e, data_ ) {
        std::cout << "(" << e . position . first << ", " << e . position . second << "; " << e  . value << ") ";
    }
    std::cout << "\n";
    std::cout << "access_.size() = " << access_ . size () << "\n";
    typedef std::pair < Position, MatrixPosition > access_value_t;
    BOOST_FOREACH ( const access_value_t & v, access_ ) {
        std::cout << "(" << v . first . first << ", " << v. first . second << "; " << v  . second << ") ";
    }
    std::cout << "\n";
    std::cout << "row_sizes_.size() = " << row_sizes_ . size () << "\n";
    for ( size_type i = 0; i < row_sizes_ . size (); ++ i ) {
        std::cout << "(" << i << ", " << row_sizes_ [ i ] << "\n";
    }
    std::cout << "\n";
    std::cout << "column_sizes_.size() = " << column_sizes_ . size () << "\n";
    for ( size_type i = 0; i < column_sizes_ . size (); ++ i ) {
        std::cout << "(" << i << ", " << column_sizes_ [ i ] << "\n";
    }
    std::cout << "row_names_.size() = " << row_names_ . size () << "\n";
    std::cout << "column_sizes_.size() = " << column_names_ . size () << "\n";
  }

};


template < class Matrix >
void print_matrix ( const Matrix & print_me ) {
  typedef int64_t size_type;
  size_type I = print_me . number_of_rows ();
  size_type J = print_me . number_of_columns ();
  
  std::cout << " Matrix is " << I << " x " << J << "\n";
  for ( size_type i = 0; i < I; ++ i ) {
    std::cout << "[";
    for ( size_type j = 0; j < J; ++ j ) {
      std::cout << std::setw(2) << print_me . read ( i, j ) << " ";
    } /* for */
    std::cout << "]\n";
  } /* for */
  //char c;
  //std::cin >> c;  
}

template < class Matrix >
void sane_matrix ( const Matrix & print_me ) {
  typedef typename Matrix::size_type size_type;
  //std::cout <<  "MATRIX SANITY CHECK\n";
  size_type I = print_me . number_of_rows ();
  size_type J = print_me . number_of_columns ();
  //std::cout << "Checking rows.\n";
  for ( size_type i = 0; i < I; ++ i ) {
    int64_t ind = print_me . row_begin ( i );
    int64_t count = 0;
    while ( ind != print_me . end () ) {
      //std::cout << "(" << print_me . row ( ind ) << ", " << print_me . column ( ind ) << "; " << print_me . read ( ind ) << " | " << ind << ") \n";
      if ( print_me . row ( ind ) != i ) {
        std::cout << "Row iterator fault.\n";
        throw std::logic_error ( "Row iterator fault.\n" );
      }
      print_me . row_advance ( ind );
      ++ count;
    }
    if ( count != print_me . row_size ( i ) ) {
      std::cout << "Discrepancy on row i = " << i << "\n";
      throw std::logic_error ( "Row size fault.\n" );
    }
  } /* for */
  //std::cout << "Checking columns.\n";

  for ( size_type j = 0; j < J; ++ j ) {
    int64_t ind = print_me . column_begin ( j );
    int64_t count = 0;
    while ( ind != print_me . end () ) {
      //std::cout << "(" << print_me . row ( ind ) << ", " << print_me . column ( ind ) << "; " << print_me . read ( ind ) << " | " << ind << ") \n";
      if ( print_me . column ( ind ) != j ) {
        std::cout << "Column iterator fault.\n";
        throw std::logic_error ( "Column iterator fault.\n" );
      }
      print_me . column_advance ( ind );
      ++ count;
    }
    if ( count != print_me . column_size ( j ) ) {
      std::cout << "Discrepancy on column j = " << j << "\n";
      throw std::logic_error ( "Column size fault.\n" );
    }  
  } /* for */
  //char c;
  //std::cin >> c;
  
  //std::cout << "\n";
}


/// Submatrix
/// Input: A, top, bottom, left, right
/// Produce matrix B = A(top:bottom, left:right)
template < class R > void
Submatrix (SparseMatrix<R> * B, 
           const typename SparseMatrix<R>::size_type top,
           const typename SparseMatrix<R>::size_type bottom,
           const typename SparseMatrix<R>::size_type left,
           const typename SparseMatrix<R>::size_type right,
           const SparseMatrix<R> & A);


/*********************************************************************************
 *             DEFINITIONS                   *
 *********************************************/

namespace SparseMatrix_detail {
  template < class R >
  void mul ( R * output, const R & lhs, const R & rhs ) {
    *output = lhs * rhs;
  }
  
  template < class R >
  void addmul ( R * output, const R & lhs, const R & rhs ) {
    *output += lhs * rhs;
  }
  
#ifdef USE_GMP
  template < >
  void mul < mpz_class > ( mpz_class * output, const mpz_class & lhs, const mpz_class & rhs ) {
    mpz_mul ( output -> get_mpz_t (), lhs . get_mpz_t (), rhs . get_mpz_t () );
  }
  
  template < >
  void addmul < mpz_class > ( mpz_class * output, const mpz_class & lhs, const mpz_class & rhs ) {
    mpz_addmul ( output -> get_mpz_t (), lhs . get_mpz_t (), rhs . get_mpz_t () );
  }
#endif
}
// Sparse Matrix Algorithm Definitions

/// submatrix (NOT IMPLEMENTED)
/// Copy the submatrix with rows [i0, i1) 
///                  and columns [j0, j1) into a new Sparse Matrix.
template < class R >
SparseMatrix<R> submatrix ( int64_t i0, int64_t i1, 
                            int64_t j0, int64_t j1, 
                            const SparseMatrix<R> & A ) {
  SparseMatrix<R> result;
  std::cout << "submatrix NOT IMPLEMENTED\n";
  return result;
}


/// operator * (SparseMatrix overload)
template < class R >
SparseMatrix<R> operator * (const SparseMatrix<R> & lhs,
                               const SparseMatrix<R> & rhs) {
  //std::cout << " operator * \n";
  typedef SparseMatrix<R> Matrix;
  typedef typename Matrix::MatrixPosition MatrixPosition;
  typedef typename Matrix::size_type size_type;
  
  int64_t I = lhs . number_of_rows ();
  int64_t K = lhs . number_of_columns ();
  int64_t J = rhs . number_of_columns ();
  Matrix result ( I, J );
  for ( int64_t k = 0; k < K; ++ k ) {
    MatrixPosition left_column_index = lhs . column_begin ( k );
    MatrixPosition right_row_index = rhs . row_begin ( k );
    
    if ( lhs . column_size ( k ) < rhs . row_size ( k ) ) {
      left_column_index = lhs . column_begin ( k );
      while ( left_column_index != lhs . end () ) {
        size_type roww = lhs . row ( left_column_index );
        R left_value = lhs . read ( left_column_index );
        right_row_index = rhs . row_begin ( k );
        while ( right_row_index != rhs . end () ) {   
          //std::cout << "M(" << roww << ", " << rhs . column ( right_row_index )
          //<< " += " << left_value * rhs . read ( right_row_index ) << "\n";
          result . add ( roww, rhs . column ( right_row_index ), 
                        left_value * rhs . read ( right_row_index ) );
          rhs . row_advance ( right_row_index );
        } /* while */
        lhs . column_advance ( left_column_index );
      } /* while */
    } else {
      right_row_index = rhs . row_begin ( k );
      while ( right_row_index != rhs . end () ) {
        size_type col = rhs . column ( right_row_index );
        R right_value = rhs . read ( right_row_index );
        left_column_index = lhs . column_begin ( k );
        while ( left_column_index != lhs . end () ) {
          //std::cout << "M(" << lhs . row ( left_column_index ) << ", " << col
          //<< " += " << lhs . read ( left_column_index ) * right_value << "\n";
          result . add ( lhs . row ( left_column_index ), col,
                        lhs . read ( left_column_index ) * right_value);
          lhs . column_advance ( left_column_index );
        } /* while */
        rhs . row_advance ( right_row_index );
      } /* while */
    } /* if-else */
  } /* for */
  return result;
}

/// operator + (SparseMatrix overload)
template < class R >
SparseMatrix<R> operator + (const SparseMatrix<R> & lhs,
                               const SparseMatrix<R> & rhs) {
  SparseMatrix<R> result;
  std::cout << "operator + NOT IMPLEMENTED\n";
  return result;
}

/// operator - (SparseMatrix overload)
template < class R >
SparseMatrix<R> operator - (const SparseMatrix<R> & lhs,
                               const SparseMatrix<R> & rhs) {
  SparseMatrix<R> result;
  std::cout << "operator - NOT IMPLEMENTED\n";
  return result;
}

/* SparseMatrix<> Method Definitions */

template < class R >
SparseMatrix<R>::
SparseMatrix ( void ) {
  resize ( 0, 0 );
}

template < class R >
SparseMatrix<R>::
SparseMatrix ( int64_t i, int64_t j ) {
  resize ( i, j );
}

template < class R >
template < class T >
SparseMatrix<R>::
SparseMatrix ( const T & copy_me ) :
  data_ ( copy_me . data_ . begin (), copy_me . data_ . end () ),
  garbage_ ( copy_me . garbage_ ), 
  access_ ( copy_me . access_ ),
  row_begin_ ( copy_me . row_begin_ ),
  column_begin_ ( copy_me . column_begin_ ),
  row_sizes_ ( copy_me . row_sizes_ ),
  column_sizes_ ( copy_me . column_sizes_ ),
  row_names_ ( copy_me . row_names_ ),
  column_names_ ( copy_me . column_names_ ) {
  resize ( copy_me . number_of_rows (), copy_me . number_of_columns () );
}

template < class R >
void SparseMatrix<R>::
resize ( int64_t i, int64_t j ) {
  row_begin_ . resize ( i, end () );
  column_begin_ . resize ( j, end () );
  row_sizes_ . resize ( i, 0 );
  column_sizes_ . resize ( j, 0 );
  // row_names and column_names aren't used yet TODO
  cache_A . resize ( std::max (i, j) );
  cache_A_TS . resize ( std::max (i, j), 0 );
  cache_B . resize ( std::max (i, j) );
  cache_B_TS . resize ( std::max (i, j), 0 );
  timestamp = 0;
  // WARNING, no support for shrinking at the moment.
}

template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::
new_index ( void ) {
  MatrixPosition result;
  // First we reclaim garbage space
  if ( not garbage_ . empty () ) {
    result = garbage_ . back ();
    garbage_ . pop_back ();
    return result;
  }
  // No garbage space is left to reclaim.
  // We make use of the automatically resizing "data_"
  result = data_ . size ();
  data_ . push_back ( Element < R > () );
  // This should be fine, as std::vector is smart enough to do length doubling
  // or other tricks
  return result;
}

template < class R >
void SparseMatrix<R>::
delete_index ( const MatrixPosition index ) {
  garbage_ . push_back ( index );
}

template < class R >
void SparseMatrix<R>::
hash_check ( const MatrixPosition index ) {
  // TODO: Probably this could save some time with extra checks
  size_type i = row ( index );
  size_type j = column ( index );
  if ( std::min ( row_sizes_ [ i ], column_sizes_ [ j ] ) <= HASH_SWITCH ) {
    access_ . erase ( Position ( i, j ) );
  } else {
    access_ . insert ( std::make_pair ( Position (i, j ), index ) );
  } /* if-else */
} /* SparseMatrix<R>::hash_check */

// TODO: put some hysteresis in the hash switching, 
//       to prevent adversarial conditions
template < class R >
void SparseMatrix<R>::
increment_row_size ( const size_type i ) {
  ++ row_sizes_ [ i ];
  if ( row_sizes_ [ i ] == 1 + HASH_SWITCH ) {
    MatrixPosition index = row_begin ( i );
    while ( index != end () ) {
      hash_check ( index );
      row_advance ( index );
    } /* while */
  } /* if */
} /* SparseMatrix<R>::increment_row_size */

template < class R >
void SparseMatrix<R>::
increment_column_size ( const size_type j ) {
  ++ column_sizes_ [ j ];
  if ( column_sizes_ [ j ] == 1 + HASH_SWITCH ) {
    MatrixPosition index = column_begin ( j );
    while ( index != end () ) {
      hash_check ( index );
      column_advance ( index );
    } /* while */
  } /* if */
} /* SparseMatrix<R>::increment_row_size */

template < class R >
void SparseMatrix<R>::
decrement_row_size ( const size_type i ) {
  -- row_sizes_ [ i ];
  if ( row_sizes_ [ i ] == HASH_SWITCH ) {
    MatrixPosition index = row_begin ( i );
    while ( index != end () ) {
      hash_check ( index );
      row_advance ( index );
    } /* while */
  } /* if */
} /* SparseMatrix<R>::decrement_row_size */

template < class R >
void SparseMatrix<R>::decrement_column_size ( const size_type j ) {
  -- column_sizes_ [ j ];
  if ( column_sizes_ [ j ] == HASH_SWITCH ) {
    MatrixPosition index = column_begin ( j );
    while ( index != end () ) {
      hash_check ( index );
      column_advance ( index );
    } /* while */
  } /* if */
} /* SparseMatrix<R>::decrement_row_size */

// erase, find
template < class R >
void SparseMatrix<R>::erase ( const MatrixPosition index ) {
  // Repair the links
  //std::cout << "ERASING " << index << "\n";
  Element < R > element = data_ [ index ];
  if ( element . left != end () ) data_ [ element . left ] . right = element . right;
  if ( element . right != end () ) data_ [ element . right ] . left = element . left;
  if ( element . up != end () ) data_ [ element . up ] . down = element . down;
  if ( element . down != end () ) data_ [ element . down ] . up = element . up;
  
  // Update the row_begin_ and column_begin_ vectors if needed
  if ( element . left == end () ) row_begin_ [ element . position . first ] = element . right;
  if ( element . up == end () ) column_begin_ [ element . position . second ] = element . down;
  
  // Empty value from Random Access Hash Table if needed
  // Check if need to include in Random Access Hash Table
  size_type i = row ( index );
  size_type j = column ( index );
  size_type r_size = row_size ( i );
  size_type c_size = column_size ( j );
  if ( std::min ( r_size, c_size ) > HASH_SWITCH ) {
    access_ . erase ( Position (i, j) ) ;
  } /* if */
  
  // Free up the index
  delete_index ( index );
  
  // Update the row and column sizes
  decrement_row_size ( i );
  decrement_column_size ( j );
}

template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::find ( const int64_t i, const int64_t j ) const {
  //std::cout << "Looking...\n";
  //static int find_count = 0;
  // Check the size of row i and the size of column j
  size_type r_size = row_size ( i );
  size_type c_size = column_size ( j );
  if ( std::min ( r_size, c_size ) > HASH_SWITCH ) {
    access_iterator it = access_ . find ( Position ( i, j ) );
    if ( it == access_ . end () ) return end ();
    //std::cout << "hash found: " << find_count ++ << "\n";
    return it -> second;
  }
  if ( r_size < c_size ) {
    // Search through row for (i, j)
    //std::cout << "row search\n";
    MatrixPosition index = row_begin ( i );
    while ( index != end () ) {
      if ( row ( index ) == i && column ( index ) == j ) {
        //std::cout << "found: " << find_count ++ << "\n";
        return index;
      } /* if */
      row_advance ( index );
    } /* while */
    //std::cout << "unfound: " << find_count ++ << "\n";
    return end ();
  } else {
    // Search through column for (i, j)
    //std::cout << "column search\n";
    MatrixPosition index = column_begin ( j );
    while ( index != end () ) {
      //std::cout << "Magically, index " << index << " is connected downward to " << data_ [ index ] . down << "\n";
      if ( row ( index ) == i && column ( index ) == j ) {
        //std::cout << "found: " << find_count ++ << "\n";
        return index;
      } /* if */
      column_advance ( index );
    } /* while */
    //std::cout << "unfound: " << find_count ++ << "\n";
    return end ();
  } /* if-else */
} /* SparseMatrix<R>::find */

// read and write

template < class R >
R SparseMatrix<R>::read ( const int64_t i, const int64_t j ) const {
  MatrixPosition index = find ( i, j );
  return read ( index );
}

template < class R >
R SparseMatrix<R>::read ( const MatrixPosition index ) const {
  if ( index == end () ) return R ( 0 );
  return data_ [ index ] . value;
}

template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::write ( const int64_t i, const int64_t j, const R value, bool insert ) {
 //////////////DEBUG BEGIN//////////////
 //std::cout << "Insert of (" << i << ", " << j << ")\n";
  /*
  static long count = 0;
  ++ count;
  if ( count % 100000 == 0) {
    std::cout << count << ".\n";
    std::cout << "size = " << size () << "\n";
  }
  if ( i >= number_of_rows () || j >= number_of_columns () ) {
    std::cout << "(i,j) = (" << i << ", " << j << ")\n";
    std::cout << "(nr,nc) = (" << number_of_rows () << ", " << number_of_columns () << ")\n";

    std::cout << "out of bounds!\n";
    exit ( 1 );
  }
  //std::cout << "Write Called.\n";
   //print_matrix ( *this );
   if ( insert ) {
     //std::cout << "Optimized insert of (" << i << ", " << j << ")\n";
     if ( find (i, j) != end () ) {
      std::cout << "Caller mistaken, the element already exists.\n";
       exit ( 1 );
     }
   }
  */
  ////////////////DEBUG END//////////////////
  MatrixPosition index;
  if ( not insert ) {
    index = find ( i, j );
    if ( index != -1 ) {
      //std::cout << " (" << i << ", " << j << ") was found!\n";
      return write ( index, value ); // deletes it if value is zero, notice
    } /* if */
    //std::cout << " (" << i << ", " << j << ") was not found!\n";
  } /* if */
  
  if ( value == R ( 0 ) ) return end (); // if value is zero, quit
                                            // Allocate new index
  index = new_index ();

  //////////////DEBUG BEGIN//////////////
  /*
  std::cout << "Allocated a new index " << index << "!\n";
  // Send it to beginning of ith row and jth column
  std::cout << "INSERTING index = " << index << " at (" << i << ", " << j << ") with value " << value << "\n";
   */
  //////////////DEBUG END////////////////

  Element < R > element ( Position ( i, j ), value, end (), row_begin ( i ), end (), column_begin ( j ) );

  data_ [ index ] = element;

  if ( row_begin ( i ) != end () ) data_ [ row_begin ( i ) ] . left = index;

  if ( column_begin ( j ) != end () ) data_ [ column_begin ( j ) ] . up = index;
  

  //////////////DEBUG BEGIN//////////////
  /*
   std::cout << "  row linking  " << data_ [ index ] . left << " " << index << " " << 
   data_ [ index ] . right << " = " << row_begin ( i ) << "\n";
   std::cout << "  col linking  " << data_ [ index ] . up << " " << index << " " << 
   data_ [ index ] . down << " = " << column_begin ( j ) << "\n";
   */
  //////////////DEBUG END////////////////

  row_begin_ [ i ] = index;
  column_begin_ [ j ] = index;

  // Check if need to include in Random Access Hash Table
  size_type r_size = row_size ( i ) + 1;

  size_type c_size = column_size ( j ) + 1;

  if ( std::min ( r_size, c_size ) > HASH_SWITCH ) {
    access_ [ Position (i, j) ] = index;
  } /* if */
  
  // Update the row and column sizes
  increment_row_size ( i );

  increment_column_size ( j );
  
  return index;
}

template < class R >
typename SparseMatrix<R>::MatrixPosition SparseMatrix<R>::write ( MatrixPosition index, const R value ) {
  if ( value == R ( 0 ) ) {
    erase ( index );
    return end ();
  } else {
    data_ [ index ] . value = value;
    return index;
  } /* if-else */
}  

// begin and end
template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::row_begin ( const int64_t i ) const {
  return row_begin_ [ i ];
}

template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::column_begin ( const int64_t j ) const {
  return column_begin_ [ j ];
}

template < class R >
typename SparseMatrix<R>::MatrixPosition 
SparseMatrix<R>::end ( void ) const {
  return -1;
}

// traversal
template < class R >
void SparseMatrix<R>::row_advance ( MatrixPosition & index ) const {
  //std::cout << "row_advance ( " << index << " ) = ";
  index = data_ [ index ] . right;
  //std::cout << index << "\n";
}

template < class R >
void SparseMatrix<R>::column_advance ( MatrixPosition & index ) const {
  //std::cout << "column_advance ( " << index << " ) = ";
  index = data_ [ index ] . down;
  // std::cout << index << "\n";
}

// row, column
template < class R >
typename SparseMatrix<R>::size_type SparseMatrix<R>::row ( const MatrixPosition index ) const {
  return data_ [ index ] . position . first;
}

template < class R >
typename SparseMatrix<R>::size_type SparseMatrix<R>::column ( const MatrixPosition index ) const {
  return data_ [ index ] . position . second;
}

// add to some position
template < class R >
void SparseMatrix<R>::add ( int64_t i, int64_t j, const R value ) {
  MatrixPosition index = find ( i, j );
  if ( index == end () ) {
    write ( i, j, value, true );
    return;
  } /* if */
  data_ [ index ] . value += value;
  // If we zeroed out the entry, delete it
  if ( data_ [ index ] . value == R ( 0 ) ) {
    erase ( index );
  }
  return;
}


template < class R >
typename SparseMatrix<R>::size_type SparseMatrix<R>::number_of_rows ( void ) const {
  return row_sizes_ . size ();
}

template < class R >
typename SparseMatrix<R>::size_type SparseMatrix<R>::number_of_columns ( void ) const {
  return column_sizes_ . size ();
}

template < class R >
typename SparseMatrix<R>::size_type SparseMatrix<R>::size ( void ) const {
  int64_t result = 0;
  BOOST_FOREACH ( size_type summand, row_sizes_ ) {
    result += summand;
  } /* boost_foreach */
  return result;
}

/***********************************
 *    ROW AND COLUMN OPERATIONS    *
 ***********************************/
         
template < class R > 
void SparseMatrix<R>::
row_operation (const R a, size_type i ) {
  for ( MatrixPosition p = row_begin ( i );
        p != end ();
        row_advance ( p ) ) {
    R value = read ( p );
    value = value * a;
    write ( p, value );
  }
}

template < class R > 
void SparseMatrix<R>::
column_operation (const R a, size_type j ) {
  for ( MatrixPosition p = column_begin ( j );
        p != end ();
        column_advance ( p ) ) {
    R value = read ( p );
    value = value * a;
    write ( p, value );
  }
}

template < class R > 
void SparseMatrix<R>::
similarity_operation (const R a, size_type i) {
  row_operation ( a, i );
  column_operation ( inverse(a), i );
}


template < class R > 
void SparseMatrix<R>::row_operation ( const R a, const R b,
                                        const R c, const R d,
                                        size_type i, size_type j ) {
  //++ number_of_pivots;
  // Increment timestamp
  ++ timestamp;
  // Note: we expect the stacks are clear. 
  
  // TODO: if timestamp has rolled over (as if) then clear cache
  /* Stage I. Copy a times the contents of row i into cache 1,
   Copy d times the contents of row j into cache 2.
   Mark the timestampts. Do not write to stack. */
  MatrixPosition i_index = row_begin ( i );
  MatrixPosition j_index = row_begin ( j );
  while ( i_index != end () ) {
    size_type col = column ( i_index );
#ifdef USE_GMP
    SparseMatrix_detail::mul ( & cache_A [ col ], a, read ( i_index ) );
#else
    cache_A [ col ] = a * read ( i_index );
#endif
    cache_A_TS [ col ] = timestamp;
    row_advance ( i_index );
  }
  
  while ( j_index != end () ) {
    size_type col = column ( j_index );
#ifdef USE_GMP
    SparseMatrix_detail::mul ( & cache_B [ col ], d, read ( j_index ) );
#else
    cache_B [ col ] = d * read ( j_index );
#endif
    cache_B_TS [ col ] = timestamp;
    row_advance ( j_index );
  }
  
  /* Stage II. Copy c * the contents of row i into cache 2,
   Copy b * the contents of row j into cache 1. 
   This time, do not report to the stack if the timestamp matches,
   and do not overwrite if the timestamp matches. If the timestamp doesn't match
   then overwrite and push to the stack */
  
  i_index = row_begin ( i );
  while ( i_index != end () ) {
    size_type col = column ( i_index );
    if ( cache_B_TS [ col ] == timestamp ) {
#ifdef USE_GMP
      SparseMatrix_detail::addmul ( & cache_B [ col ], c, read ( i_index ) );
#else
      cache_B [ col ] += c * read ( i_index );
#endif
    } else {
#ifdef USE_GMP
      SparseMatrix_detail::mul ( & cache_B [ col ], c, read ( i_index ) );
#else
      cache_B [ col ] = c * read ( i_index );
#endif
      cache_B_S . push_back ( col );
    }
    row_advance ( i_index );
  }
  
  j_index = row_begin ( j );
  while ( j_index != end () ) {
    size_type col = column ( j_index );
    if ( cache_A_TS [ col ] == timestamp ) {
#ifdef USE_GMP
      SparseMatrix_detail::addmul ( & cache_A [ col ], b, read ( j_index ) );
#else
      cache_A [ col ] += b * read ( j_index );
#endif
    } else {
#ifdef USE_GMP
      SparseMatrix_detail::mul ( & cache_A [ col ], b, read ( j_index ) );
#else
      cache_A [ col ] = b * read ( j_index );
#endif
      cache_A_S . push_back ( col );
    }
    row_advance ( j_index );
  }
  
  /* Stage III. Replace the entries in row i and row j that are already there with the new values */
  
  i_index = row_begin ( i );
  while ( i_index != end () ) {
    size_type col = column ( i_index );
    write ( i_index, cache_A [ col ] );
    row_advance ( i_index );
  }
  
  j_index = row_begin ( j );
  while ( j_index != end () ) {
    size_type col = column ( j_index );
    write ( j_index, cache_B [ col ] );
    row_advance ( j_index );
  }
  
  /* Stage IV. Unload the stacks and insert the new elements */
  while ( not cache_A_S . empty () ) {
    size_type col = cache_A_S . back ();
    cache_A_S . pop_back ();
    write ( i, col, cache_A [ col ], true );
  }
  
  while ( not cache_B_S . empty () ) {
    size_type col = cache_B_S . back ();
    cache_B_S . pop_back ();
    write ( j, col, cache_B [ col ], true );
  }
  
}

//   i -> ai+bj   [i j]  [a  c]
//   j -> ci+dj          [b  d]
template < class R > 
void SparseMatrix<R>::
column_operation (const R a, const R b,
                  const R c, const R d,
                  size_type i, size_type j ) {
  //++ number_of_pivots;
  // Increment timestamp
  ++ timestamp;
  // Note: we expect the stacks are clear. 
  
  //std::cout << "column operation.\n";
  //std::cout << a << " " << b << " " << c << " " << d << " " << i << " " << j << "\n";
  
  // TODO: if timestamp has rolled over (as if) then clear cache
  /* Stage I. Copy a times the contents of column i into cache 1,
   Copy d times the contents of column j into cache 2.
   Mark the timestampts. Do not write to stack. */
  MatrixPosition i_index = column_begin ( i );
  while ( i_index != end () ) {
    size_type roww = row ( i_index );
#ifdef USE_GMP
    SparseMatrix_detail::mul ( &cache_A [ roww ], a, read ( i_index ) );
#else
    cache_A [ roww ] = a * read ( i_index );
#endif
    cache_A_TS [ roww ] = timestamp;
    //std::cout << i_index << "\n";
    column_advance ( i_index );
  }
  //std::cout << "*****\n";
  
  MatrixPosition j_index = column_begin ( j );
  while ( j_index != end () ) {
    size_type roww = row ( j_index );
#ifdef USE_GMP
    SparseMatrix_detail::mul ( &cache_B [ roww ], d, read ( j_index ) );
#else
    cache_B [ roww ] = d * read ( j_index );
#endif
    cache_B_TS [ roww ] = timestamp;
    column_advance ( j_index );
  }
  
  /* Stage II. Copy c * the contents of column i into cache 2,
   Copy b * the contents of column j into cache 1. 
   This time, do not report to the stack if the timestamp matches,
   and do not overwrite if the timestamp matches. If the timestamp doesn't match
   then overwrite and push to the stack */
  
  i_index = column_begin ( i );
  while ( i_index != end () ) {
    size_type roww = row ( i_index );
    if ( cache_B_TS [ roww ] == timestamp ) {
#ifdef USE_GMP
      SparseMatrix_detail::addmul ( &cache_B [ roww ], c, read ( i_index ) );
#else
      cache_B [ roww ] += c * read ( i_index );
#endif
    } else {
#ifdef USE_GMP
      SparseMatrix_detail::mul ( &cache_B [ roww ], c, read ( i_index ) );
#else
      cache_B [ roww ] = c * read ( i_index );
#endif
      cache_B_S . push_back ( roww );
    }
    column_advance ( i_index );
  }
  
  j_index = column_begin ( j );
  while ( j_index != end () ) {
    size_type roww = row ( j_index );
    if ( cache_A_TS [ roww ] == timestamp ) {
#ifdef USE_GMP
      SparseMatrix_detail::addmul ( &cache_A [ roww ], b, read ( j_index ) );
#else
      cache_A [ roww ] += b * read ( j_index );
#endif
    } else {
#ifdef USE_GMP
      SparseMatrix_detail::mul ( &cache_A [ roww ], b, read ( j_index ) );
#else
      cache_A [ roww ] = b * read ( j_index );
#endif
      cache_A_S . push_back ( roww );
    }
    column_advance ( j_index );
  }
  
  /* Stage III. Replace the entries in column i and column j that are already there with the new values */
  
  i_index = column_begin ( i );
  while ( i_index != end () ) {
    size_type roww = row ( i_index );
    //std::cout << "III. (" << row ( i_index ) << ", " << column ( i_index ) << ") = " << cache_A [ roww ] << "\n";
    write ( i_index, cache_A [ roww ] );
    column_advance ( i_index );
  }
  
  j_index = column_begin ( j );
  while ( j_index != end () ) {
    size_type roww = row ( j_index );
    //std::cout << "III. (" << row ( j_index ) << ", " << column ( j_index ) << ") = " << cache_B [ roww ] << "\n";
    write ( j_index, cache_B [ roww ] );
    column_advance ( j_index );
  }
  
  /* Stage IV. Unload the stacks and insert the new elements */
  while ( not cache_A_S . empty () ) {
    size_type roww = cache_A_S . back ();
    cache_A_S . pop_back ();
    //std::cout << "IV. (" << roww << ", " << i << ") = " << cache_A [ roww ] << "\n";
    write ( roww, i, cache_A [ roww ], true );
  }
  
  while ( not cache_B_S . empty () ) {
    size_type roww = cache_B_S . back ();
    cache_B_S . pop_back ();
    //std::cout << "IV. (" << roww << ", " << j << ") = " << cache_B [ roww ] << "\n";
    write ( roww, j, cache_B [ roww ], true );
  }
}

template < class R >
void SparseMatrix<R>::
similarity_operation ( const R a, const R b,
                       const R c, const R d,
                       size_type i, size_type j ) {
  R det = a * d - b * c;
  row_operation ( a, b,
                  c, d,
                  i, j);
  column_operation ( d/det, -c/det,
                     -b/det, a/det, 
                     i, j );
}

template < class R >
void SparseMatrix<R>::
swap_rows ( const int64_t i, const int64_t j ) {
  // TODO: make more efficient
  row_operation (R ( 0 ), R ( 1 ),
                 R ( 1 ), R ( 0 ),
                 i, j );
}

template < class R >
void SparseMatrix<R>::
swap_columns ( const int64_t i, const int64_t j ) {
  // TODO: make more efficient
  column_operation (R ( 0 ), R ( 1 ),
                    R ( 1 ), R ( 0 ),
                    i, j );
}

template < class R >
typename SparseMatrix<R>::size_type 
SparseMatrix<R>::row_size ( const size_type i ) const {
  return row_sizes_ [ i ];
}

template < class R >
typename SparseMatrix<R>::size_type 
SparseMatrix<R>::column_size ( const size_type j ) const {
  return column_sizes_ [ j ];
}



/// Submatrix
/// Input: A, top, bottom, left, right
/// Produce matrix B = A(top:bottom, left:right)

template < class R >
void Submatrix (SparseMatrix<R> * B, 
                const typename SparseMatrix<R>::size_type top,
                const typename SparseMatrix<R>::size_type bottom,
                const typename SparseMatrix<R>::size_type left,
                const typename SparseMatrix<R>::size_type right,
                const SparseMatrix<R> & A) {
  typedef typename SparseMatrix<R>::size_type size_type;
  typedef typename SparseMatrix<R>::MatrixPosition MatrixPosition;
  B -> resize ( bottom - top + 1, right - left + 1 );
  for ( size_type i = top; i <= bottom; ++ i ) { 
    // warning: if top - bottom >> right - left and matrix is VERY sparse? then we should have make columns the outer loop
    MatrixPosition entry = A . row_begin ( i );
    while ( entry != A . end () ) {
      size_type j = A . column ( entry );
      if (  left <= j && j <= right ) B -> write ( i - top, j - left, A . read ( entry ) );
      A . row_advance ( entry );
    } // while
  } // for
}


template < class R >
class TransposeWrapper {
private:
  SparseMatrix<R> * m_;
public:
  typedef R value_type;
  typedef R entry_type;
  typedef int64_t MatrixPosition;
  typedef int64_t size_type;
  TransposeWrapper ( SparseMatrix<R> * m ) : m_(m) {}

  /// erase 
  ///   Set the entry corresponding to MatrixPosition "index"
  ///   to zero
  void erase ( const MatrixPosition index ) {
    m_ -> erase ( index );
  }

  /// find
  ///    Find the MatrixPosition corresponding to (i,j)
  ///    If there is no such MatrixPosition (i.e. A(i,j) == 0 )
  ///    then return end ()
  MatrixPosition find ( const size_type i, const size_type j ) const {
    return m_ -> find ( j, i );
  }
  
  /// read 
  ///     Return the value at the (i,j) position
  R read ( const size_type i, const size_type j ) const {
    return m_ -> read ( j, i );
  }

  /// read
  ///     Return the value at the position 
  ///     represented by MatrixPosition "index"
  R read ( MatrixPosition index ) const {
    return m_ -> read ( index );
  }

  /// write ( i, j, value, insert = false)
  ///   Set (i,j) entry to "value".
  ///     Note: If it is known in advance that (i,j) holds zero before the call, 
  ///           then "insert" may be set to true (it is false by default)
  ///           which results in a speedup.
  MatrixPosition 
  write ( const size_type i, const size_type j, 
          const R value, bool insert = false ) {  
    return m_ -> write ( j, i, value, insert );
  }

  /// write ( MatrixPosition index, const R value )
  ///   Copy "value" to the position represented by "index"
  MatrixPosition 
  write ( MatrixPosition index, const R value ) {
    return m_ -> write ( index, value );
  }
  
  // begin and end
  MatrixPosition row_begin ( const size_type i ) const {
    return m_ -> column_begin ( i );
  }

  MatrixPosition column_begin ( const size_type j ) const {
    return m_ -> row_begin ( j );
  }

  MatrixPosition end ( void ) const {
    return m_ -> end ();
  } 

  // traversal. row_advance advances within a row, column_advance within a column
  void row_advance ( MatrixPosition & index ) const {
    m_ -> column_advance ( index );
  }

  void column_advance ( MatrixPosition & index ) const {
    m_ -> row_advance ( index );
  }
  
  // row and column (to learn position)
  size_type row ( const MatrixPosition index ) const {
    return m_ -> column ( index );
  }

  size_type column ( const MatrixPosition index ) const {
    return m_ -> row ( index );
  }
  
  // add to some position
  void add ( size_type i, size_type j, const R value ) {
    m_ -> add ( j, i, value );
  }
  
  /// resize
  ///    Change the shape of the matrix
  void resize ( size_type i, size_type j ) {
    m_ -> resize ( j, i );
  }

  /// size
  ///   Report number of non-zero elements in
  ///   the sparse matrix
  size_type size ( void ) const {
    return m_ -> size ();
  }
  
  /// row_operation
  ///   Multiply row i by a           
  void row_operation (const R a, size_type i ) {
    m_ -> column_operation ( a, i );
  }

  /// column_operation
  ///   Multiply column j by a
  void column_operation (const R a, size_type j ) {
    m_ -> row_operation ( a, j );
  }

  /// similarity_operation
  ///   Multiply row i by a and divide column i by a
  void similarity_operation (const R a, size_type i) {
    m_ -> similarity_operation ( inverse(a), i );
  }

  /// row_operation
  ///   Multiply  on the left by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a b ]
  ///               [ c d ]             
  void row_operation (const R a, const R b,
                      const R c, const R d,
                      size_type i, size_type j ) {
    m_ -> column_operation ( a,b,c,d,i,j );
  }

  /// column_operation
  ///   Multiply  on the right by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a c ]
  ///               [ b d ]  
  void column_operation (const R a, const R b,
                         const R c, const R d,
                         size_type i, size_type j ) {
    m_ -> row_operation ( a,b,c,d,i,j );
  }

  /// similarity_operation
  ///   Multiply  on the left by the matrix M
  ///      where M_kk =  1  for  k \notin {i, j}
  ///      and the {i,j}x{i,j} submatrix M_s is
  ///         M_s = [ a b ]
  ///               [ c d ] 
  ///   and multiply on the right by the inverse of M
  void similarity_operation (const R a, const R b,
                             const R c, const R d,
                             size_type i, size_type j ) {
    R det = a*d - b*c;
    m_ -> similarity_operation ( d/det, -c/det, -b/det, a/det, i, j);
  }

  /// swap_rows
  ///    Swap row i and row j
  void swap_rows ( const size_type i, const size_type j ) {
    m_ -> swap_columns ( i, j );
  }

  /// swap_columns
  ///    Swap column i and column j
  void swap_columns ( const size_type i, const size_type j ) {
    m_ -> swap_rows ( i, j );
  }

  /// number_of_rows
  ///    Return number of rows in matrix
  size_type number_of_rows ( void ) const {
    return m_ -> number_of_columns ();
  }

  /// number_of_columns
  ///    Return number of columns in matrix
  size_type number_of_columns ( void ) const {
    return m_ -> number_of_rows ();
  }

  /// row_size
  ///     Return number of non-zero elements in row i
  size_type row_size ( const size_type i ) const {
    return m_ -> column_size ( i );
  }

  /// columns_size
  ///     Return number of non-zero elements in column j  
  size_type column_size ( const size_type j ) const {
    return m_ -> row_size ( j );
  }

};
} // namespace chomp

#endif
