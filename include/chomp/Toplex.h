// Toplex.h
// Shaun Harker
// 9/16/11

#ifndef CHOMP_TOPLEX_H
#define CHOMP_TOPLEX_H

#include <cstdlib>
#include <stdint.h>

#include <vector>
#include <stack>
#include <boost/unordered_set.hpp>

#include "chomp/Rect.h"
#include "chomp/Prism.h"
#include "chomp/RelativePair.h"
#include "chomp/CubicalComplex.h"
#include "chomp/ToplexDetail.h"
#include "chomp/BitmapSubcomplex.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace chomp {

/**********
 * Toplex *
 **********/

class Toplex {
public:
  /* typedefs */
  typedef GridElement value_type;
  typedef uint32_t size_type;
  typedef Toplex_const_iterator iterator;
  typedef iterator const_iterator;
  typedef GridElement Top_Cell; //compatibility
  typedef Rect Geometric_Description; // compatibility
  /* Basic Container Interface */
  void erase ( iterator erase_me );
  void clear ( void );
  iterator find ( const GridElement & find_me ) const;
  iterator begin ( void ) const;
  iterator end ( void ) const;
  size_type size ( void ) const;
  
  /// tree_size
  size_type tree_size ( void ) const;
  
  /// dimension
  int dimension ( void ) const;
  
  /// bounds
  const Rect & bounds ( void ) const;

  /// geometry
  Rect geometry ( const GridElement & GridElement ) const;
  Rect geometry ( const const_iterator & cell_iterator ) const;
  
  /// prefix ( GridElement )
  ///   Return a vector with the prefix string of tree moves necessary to 
  ///   navigate to its leaf in the Toplex tree structure
  
  std::vector < unsigned char > prefix ( const GridElement & element ) const;

  /// children ( GridElement)
  ///    Return an array of children GridElements,
  ///    or else empty if the GridElement is not further subdivided.
  template < class InsertIterator > InsertIterator 
  children ( InsertIterator ii, const GridElement & element ) const;

  /// leaves ( GridElement )
  template < class InsertIterator, class Container > InsertIterator
  leaves ( InsertIterator ii, 
           const Container & elements ) const;
  
  /// umbrella ( std::vector < GridElement > & elements )
  ///    Return the set of all GridElements whose every descendent is in "elements"
  
  template < class InsertIterator, class Container > InsertIterator
  umbrella ( InsertIterator ii, const Container & elements ) const;

  /// cover   (rect version)
  template < class InsertIterator > InsertIterator
  cover ( InsertIterator & ii, const Rect & r ) const;

  /// cover   (prism version)
  template < class InsertIterator > InsertIterator
  cover ( InsertIterator & ii, const Prism & p ) const;

  /// cover   (vector version)  FOR UNIONS
  template < class InsertIterator, class T > void
  cover ( InsertIterator & ii, const std::vector < T > & V ) const;
  
  /// cover   (pair version)    FOR INTERSECTIONS
  template < class InsertIterator, class S, class T > void
  cover ( InsertIterator & ii, const std::pair < S, T > & V ) const;
  
  /// coarse cover   (whenever node containment, report parent, not children)
  template < class InsertIterator > InsertIterator
  coarseCover ( InsertIterator ii, const Rect & geometric_region ) const;
  
  /// subdivide
  template < class InsertIterator > InsertIterator
  subdivide ( InsertIterator ii, GridElement divide_me );

  template < class InsertIterator > InsertIterator
  subdivide ( InsertIterator ii, iterator divide_me );
  
  template < class InsertIterator, class Container > InsertIterator
  subdivide ( InsertIterator ii, const Container & subset_to_divide );

  template < class InsertIterator > InsertIterator
  subdivide ( InsertIterator ii );
  
  void subdivide ( void );
  
  /// depth
  int getDepth ( const GridElement & ge ) const;
  
  template < class Container >
  int getDepth ( const Container & subset ) const;
  
  /// coarsen
  template < class Container > void 
  coarsen ( const Container & coarsen_to );
  
  /// GridElementToCubes
  void GridElementToCubes ( std::vector<std::vector < uint32_t > > * cubes, 
                        const GridElement e, int d ) const;
  
  /// relativeComplex
  template < class Container > void 
  relativeComplex ( RelativePair * pair,
                    const Container & XGridElements,
                    const Container & AGridElements,
                    int depth) const;
  // Construction
  Toplex ( void );
  Toplex ( const Rect & outer_bounds_of_toplex );
  void initialize ( const Rect & outer_bounds_of_toplex );
  void initialize ( const Rect & outer_bounds_of_toplex, const std::vector < bool > & periodic );

  ~Toplex ( void );
  
private:
  const_iterator begin_;
  const_iterator end_;
  size_type size_;
  size_type tree_size_;
  std::vector < iterator > find_;
  Node * root_;
  Rect bounds_;
  int dimension_;
  std::vector < bool > periodic_;
  friend class boost::serialization::access;
public:
  std::vector < bool > periodic ( void ) const { return periodic_; }
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & begin_;
    ar & end_;
    ar & size_;
    ar & tree_size_;
    ar & find_;
    ar & root_;
    ar & bounds_;
    ar & dimension_;
    ar & periodic_;
  }
  // file operations
  void save ( const char * filename ) const {
    std::ofstream ofs(filename);
    assert(ofs.good());
    boost::archive::text_oarchive oa(ofs);
    oa << * this;
    ofs . close ();
  }
  
  void load ( const char * filename ) {
    std::ifstream ifs(filename);
    if ( not ifs . good () ) {
      std::cout << "Could not load " << filename << "\n";
      exit ( 1 );
    }
    boost::archive::text_iarchive ia(ifs);
    ia >> * this;    
    ifs . close ();
  }
 
};

inline void Toplex::erase ( iterator erase_me ) {
  /* Update begin_ if necessary */
  /* TODO: if this is being used to kill an entire subtree, not just a leaf, then
   this will not update begin_ properly. */
  /* BUG: size not updated properly, because when we move onto parents
   it keeps decreasing the size erroneously */
  /* debug */
  if ( erase_me . node () -> left_ != NULL ||
      erase_me . node () -> right_ != NULL )
    std::cout << "Erasing a non-leaf node\n";
  /* end debug */
  if ( erase_me == begin_ ) ++ begin_;
  /* Disconnect from parent */
  Node * node = erase_me . node_;
  while ( 1 ) {
    Node * parent = node -> parent_;
    if ( parent == NULL ) {
      /* We are emptying the toplex -- the root is gone! */
      root_ = NULL;
      break;
    }
    /* This is not the root: a parent exists. */
    if ( parent -> left_ == node ) {
      /* This is a left-child; disconnect accordingly. */
      parent -> left_ = NULL;
    } else {
      /* This is a right-child; disconnect accordingly. */
      parent -> right_ = NULL;
    } /* if-else */
    if ( parent -> left_ != parent -> right_ ) break;
    /* We will erase this node and move on to erase its parent */
    find_ [ node -> contents_ ] = end_;
    delete node;
    node = parent;
  } /* while */
  /* Update find_ */
  find_ [ node -> contents_ ] = end_;
  /* Update size_ */
  -- size_;
  /* Deallocate the node */
  delete node;
} /* Adaptive_Cubical::Toplex::erase */

inline void Toplex::clear ( void ) {
  size_ = 0;
  tree_size_ = 0;
  find_ . clear ();
  begin_ = end_;
  if ( root_ != NULL ) delete root_;
} /* Adaptive_Cubical::Toplex::clear */

inline Toplex::iterator Toplex::find ( const GridElement & find_me ) const {
  return find_ [ find_me ];
} /* Adaptive_Cubical::Toplex::find */

inline Toplex::iterator Toplex::begin ( void ) const {
  return begin_;
} /* Adaptive_Cubical::Toplex::begin */

inline Toplex::iterator Toplex::end ( void ) const {
  return end_;
} /* Adaptive_Cubical::Toplex::end */

inline Toplex::size_type Toplex::size ( void ) const {
  return size_;
} /* Adaptive_Cubical::Toplex::size */

inline Toplex::size_type Toplex::tree_size ( void ) const {
  return tree_size_;
} /* Adaptive_Cubical::Toplex::tree_size */

inline int Toplex::dimension ( void ) const {
  return dimension_;
} /* Adaptive_Cubical::Toplex::dimension */

inline const Rect & Toplex::bounds ( void ) const {
  return bounds_;
} /* Adaptive_Cubical::Toplex::bounds */

inline std::vector < unsigned char > Toplex::prefix ( const GridElement & cell ) const {
  std::vector < unsigned char > result, reversed;
  iterator cell_iterator = find ( cell );
  Node * node_ptr = cell_iterator . node_;
  while ( node_ptr != root_ ) {
    Node * parent = node_ptr -> parent_;
    if ( parent -> left_ == node_ptr ) {
      /* This is a left-child */
      reversed . push_back ( 0 );
    } else {
      /* This is a right-child */
      reversed . push_back ( 1 );
    } /* if-else */
    node_ptr = parent;
  } /* while */
  /* Now reverse the order */
  for ( int i = reversed . size () - 1; i >= 0; -- i ) {
    result . push_back ( reversed [ i ] );
  }
  // WEIRD: the algorithm reverse_copy doesn't work.
  // I get this bizarre seg-fault, and gdb tells me we've have been hijacked
  // by some boost MPL crap (other than STL implementation). Apparently, it doesn't work right.
  // std::reverse_copy ( reversed . begin (), reversed . end (), result . begin () );
  return result;
} /* Adaptive_Cubical::Toplex::prefix */


template < class InsertIterator > InsertIterator
Toplex::children ( InsertIterator ii, 
                  const GridElement & element ) const {
  iterator cell_iterator = find ( element );
  Node * ptr = cell_iterator . node_;
  if ( ptr -> left_ != NULL ) * ii ++ = ptr -> left_ -> contents_;
  if ( ptr -> right_ != NULL ) * ii ++ = ptr -> right_ -> contents_;
  return ii;
}

template < class InsertIterator, class Container > InsertIterator
Toplex::leaves ( InsertIterator ii, 
                 const Container & elements ) const {
  std::stack < Node * > nodes;
  BOOST_FOREACH ( GridElement element, elements ) {
    iterator cell_iterator = find ( element );
    Node * node_ptr = cell_iterator . node_;
    nodes . push ( node_ptr );
  }
  while ( not nodes . empty () ) {
    Node * ptr = nodes . top ();
    nodes . pop ();
    if ( ptr -> left_ == NULL && ptr -> right_ == NULL ) {
      // It's a leaf.
      * ii ++ = ptr -> contents_;
    } else {
      // It's not a leaf.
      if ( ptr -> left_ != NULL ) {
        nodes . push ( ptr -> left_ );
      }
      if ( ptr -> right_ != NULL ) {
        nodes . push ( ptr -> right_ );
      }
    }
  }
  return ii;
}

#if 0
// oops, this isn't what i wanted, but maybe it could be useful later
template < class InsertIterator > void 
Toplex::children ( InsertIterator & ii, 
                   const GridElement & element ) const {
  iterator cell_iterator = find ( element );
  Node * node_ptr = cell_iterator . node_;
  std::stack < Node * > nodes;
  nodes . push ( node_ptr );
  while ( not nodes . empty () ) {
    Node * ptr = nodes . top ();
    nodes . pop ();
    if ( ptr -> dimension_ == 0 ) {
      * ii ++ = ptr -> contents_;
    } else {
      if ( ptr -> left_ != NULL ) {
        nodes . push ( ptr -> left_ );
      }
      if ( ptr -> right_ != NULL ) {
        nodes . push ( ptr -> right_ );
      }
    }
  }
}
#endif

template < class InsertIterator, class Container > InsertIterator
Toplex::umbrella ( InsertIterator ii, 
                   const Container & elements ) const {
  std::vector<GridElement> result ( elements . begin (), elements . end () );
  boost::unordered_set < GridElement > umb;
  int N = result . size ();
  for ( int i = 0; i < N; ++ i ) {
    GridElement element = result [ i ];
    iterator cell_iterator = find ( element );
    Node * parent_node = find ( element ) . node_ -> parent_;
    if ( parent_node == NULL ) continue;
    GridElement parent =  parent_node -> contents_;
    if ( umb . insert ( parent ) . second ) {
      result . push_back ( parent );
      ++ N;
    }
  }
  for ( int i = 0; i < N; ++ i ) * ii ++ = result [ i ];
  return ii;
}


inline Rect Toplex::geometry ( const const_iterator & cell_iterator ) const {
  Rect return_value ( dimension_, Real ( 0 ) );
  //std::cout << "geometry of " << cell_iterator . node_ << " (" << cell_iterator . node_ -> contents_ << ")\n";
  //std::cout << "root = " << root_ << "\n";
  /* Climb the tree */
  Node * node_ptr = cell_iterator . node_;
  while ( node_ptr != root_ ) {
    //std::cout << "visiting " << node_ptr << " with parent " << node_ptr -> parent_ << "\n";
    Node * parent = node_ptr -> parent_;
    int division_dimension = parent -> dimension_;
    if ( parent -> left_ == node_ptr ) {
      /* This is a left-child */
      return_value . upper_bounds [ division_dimension ] += Real ( 1 );
    } else {
      /* This is a right-child */
      return_value . lower_bounds [ division_dimension ] += Real ( 1 );
    } /* if-else */
    return_value . lower_bounds [ division_dimension ] /= Real ( 2 );
    return_value . upper_bounds [ division_dimension ] /= Real ( 2 );
    node_ptr = parent;
  } /* while */
  for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
    /* Produce convex combinations */
    return_value . lower_bounds [ dimension_index ] = return_value . lower_bounds [ dimension_index ] * bounds_ . upper_bounds [ dimension_index ] +
    ( Real ( 1 ) - return_value . lower_bounds [ dimension_index ] ) * bounds_ . lower_bounds [ dimension_index ];
    return_value . upper_bounds [ dimension_index ] = return_value . upper_bounds [ dimension_index ] * bounds_ . lower_bounds [ dimension_index ] +
    ( Real ( 1 ) - return_value . upper_bounds [ dimension_index ] ) * bounds_ . upper_bounds [ dimension_index ];
    //DEBUG
    if ( return_value . lower_bounds [ dimension_index ] > return_value . lower_bounds [ dimension_index ] ) {
      std::cout << "Toplex::geometry ERROR: constructed invalid region.\n";
      exit(1);
    }
  } /* for */
  //std::cout << "returning.\n";
  return return_value;
} /* Adaptive_Cubical::Toplex::geometry */

inline Rect Toplex::geometry ( const GridElement & cell  ) const {
  return geometry ( find ( cell ) );
} /* Adaptive_Cubical::Toplex::geometry */


// TODO:
// I have three cover routines below.
// The first two should be made into one, with a template parameter
// The third I don't know what to do with yet. Will the coarse cover idea pan out?

  template < class InsertIterator >
  inline InsertIterator
  Toplex::cover ( InsertIterator & ii, const Rect & geometric_region ) const {
    //std::cout << "Rect version of Cover\n";
    //std::cout << "Covering " << geometric_region << "\n";
    // Deal with periodicity
    
    boost::unordered_set < GridElement > redundancy_check;
    
    std::vector < double > width ( dimension_ );
    for ( int d = 0; d < dimension_; ++ d ) {
      width [ d ] = bounds_ . upper_bounds [ d ] - bounds_ . lower_bounds [ d ];
    }
    
    std::stack < Rect > work_stack;
    Rect R = geometric_region;
    for ( int d = 0; d < dimension_; ++ d ) {
      if ( periodic_ [ d ] == false ) continue;
      if ( R . upper_bounds [ d ] > bounds_ . upper_bounds [ d ] ) {
        R . lower_bounds [ d ] -= width [ d ];
        R . upper_bounds [ d ] -= width [ d ];
      }
      if ( R . upper_bounds [d] - R . lower_bounds [ d ] > width [ d ] )
        R . upper_bounds [ d ] = R . lower_bounds [ d ] + width [ d ];
    }
    
    long periodic_long = 0;
    for ( int d = 0; d < dimension_; ++ d ) {
      if ( periodic_ [ d ] ) periodic_long += (1 << d );
    }
    
    // loop through all 2^l periodic images, avoiding repeats
    std::set < long > periodic_images;
    long hypercube = 2 << dimension_;
    for ( long k = 0; k < hypercube; ++ k ) {
      if ( periodic_images . count ( k & periodic_long ) ) continue;
      periodic_images . insert ( k & periodic_long );
      Rect r = R;
      for ( int d = 0; d < dimension_; ++ d ) {
        if ( periodic_ [ d ] == false ) continue;
        if ( k & ( 1 << d ) ) {
          r . lower_bounds [ d ] += width [ d ];
          r . upper_bounds [ d ] += width [ d ];
        }
      }
      work_stack . push ( r );
      //std::cout << "Pushed " << r << "\n";
    }
    
    //std::cout << "ready to cover pushed things\n";
    /* Use a stack, not a queue, and do depth first search.
     The advantage of this is that we can maintain the geometry during our Euler Tour.
     We can maintain our geometry without any roundoff error if we use the standard box
     [0,1]^d. To avoid having to translate to real coordinates at each leaf, we instead
     convert the input to these standard coordinates, which we put into integers. */
    
    while ( not work_stack . empty () ) {
      //std::cout << "Top of cover loop. Size of work stack = " << work_stack . size () << "\n";
      Rect GR = work_stack . top ();
      work_stack . pop ();
      //std::cout << "Trying to cover " << GR << "\n";
      // Step 1. Convert input to standard coordinates. 
      Rect region ( dimension_ );
      static std::vector<uint64_t> LB ( dimension_);
      static std::vector<uint64_t> UB ( dimension_);
#define INTPHASEWIDTH (((uint64_t)1) << 60)
      static Real bignum ( INTPHASEWIDTH );
      bool out_of_bounds = false;
      for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
        region . lower_bounds [ dimension_index ] = 
        (GR . lower_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]) /
        (bounds_ . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]);
        region . upper_bounds [ dimension_index ] = 
        (GR . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]) /
        (bounds_ . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]);
        
        if (region . upper_bounds [ dimension_index ] < Real ( 0 ) ||
            region . lower_bounds [ dimension_index ] > Real ( 1 ) )  {
          out_of_bounds = true;
          break;
        }
        //std::cout << "dim " << dimension_index << ": [" << region . lower_bounds [ dimension_index ] << ", " << region . upper_bounds [ dimension_index ] << "]\n";
        if ( region . lower_bounds [ dimension_index ] < Real ( 0 ) ) 
          region . lower_bounds [ dimension_index ] = Real ( 0 );
        if ( region . lower_bounds [ dimension_index ] > Real ( 1 ) ) 
          region . lower_bounds [ dimension_index ] = Real ( 1 );
        LB [ dimension_index ] = (uint64_t) ( bignum * region . lower_bounds [ dimension_index ] );
        if ( region . upper_bounds [ dimension_index ] < Real ( 0 ) ) 
          region . upper_bounds [ dimension_index ] = Real ( 0 );
        if ( region . upper_bounds [ dimension_index ] > Real ( 1 ) ) 
          region . upper_bounds [ dimension_index ] = Real ( 1 );
        UB [ dimension_index ] = (uint64_t) ( bignum * region . upper_bounds [ dimension_index ] );
      }
      if ( out_of_bounds ) continue;
      // Step 2. Perform DFS on the Toplex tree, recursing whenever we have intersection,
      //         (or adding leaf to output when we have leaf intersection)
      static std::vector<uint64_t> NLB ( dimension_);
      static std::vector<uint64_t> NUB ( dimension_);
      for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
        //if ( LB [ dimension_index ] < (1 << 20) ) LB [ dimension_index ] = 0;
        //if ( LB [ dimension_index ] >= (1 << 20) ) LB [ dimension_index ] -= (1 << 20);
        //if ( UB [ dimension_index ] < (INTPHASEWIDTH - (1 << 20)) ) UB [ dimension_index ] += (1 << 20);
        //if ( UB [ dimension_index ] >= (INTPHASEWIDTH - (1 << 20)) ) UB [ dimension_index ] = INTPHASEWIDTH;
        NLB [ dimension_index ] = 0;
        NUB [ dimension_index ] = INTPHASEWIDTH;
      }
      //std::cout << "C\n";
      
      /* Strategy. 
       We will take the Euler Tour using a 4-state machine.
       There are Four states.
       0 = Just Descended. Check for an intersection.
       1 = Descend to the left
       2 = Descend to right 
       3 = Rise.
       */
      
      Node * N = root_;
      char state = 0; 
      
      //std::cout << "D\n";
      
      while ( 1 ) {
        //std::cout << "Entering Loop, state = " << (int) state << "\n";
        //std::cout << " N = " << N << "\n";
        if ( state == 0 ) {
          // If we have descended here, then we should check for intersection.
          bool intersect_flag = true;
          for ( int d = 0; d < dimension_; ++ d ) {
            if ( LB[d] > NUB[d] || UB[d] < NLB [d] ) {  // INTERSECTION CHECK
              intersect_flag = false;
              break;
            }
          }
          
          if ( intersect_flag ) {
            //std::cout << "Detected intersection.\n";
            // Check if its a leaf.
            if ( N -> left_ == NULL ) {
              if ( N -> right_ == NULL ) {
                // Here's what we are looking for.
                if ( redundancy_check . count ( N -> contents_ ) == 0 ) {
                  * ii ++ = N -> contents_; // OUTPUT
                                            //std::cout << "output " << N -> contents_ << "\n";
                            redundancy_check . insert ( N -> contents_ ); 
                                            }
                                          //std::cout << "cover -- " << N -> contents_ << "\n";
                                          // Issue the order to rise.
                                          //std::cout << "Issue rise.\n";
                state = 3;
              } else {
                // Issue the order to descend to the right.
                //std::cout << "Issue descend right.\n";
                state = 2;
              } 
            } else {
              // Issue the order to descend to the left.   
              //std::cout << "Issue descend left.\n";
              state = 1;
            }
          } else {
            // No intersection, issue order to rise.
            //std::cout << "No intersection. \n";
            //std::cout << "Issue Rise.\n";
            state = 3;
          } // intersection check complete
        } // state 0
        
        if ( state == 1 ) {
          // We have been ordered to descend to the left.
          //std::cout << "Descend left.\n";
          int div_dim = N -> dimension_;
          NUB[div_dim] -= ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
          N = N -> left_;
          state = 0;
          continue;
        } // state 1
        
        if ( state == 2 ) {
          // We have been ordered to descend to the right.
          //std::cout << "Descend right.\n";
          int & div_dim = N -> dimension_;
          NLB[div_dim] += ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
          N = N -> right_;
          state = 0;
          continue;
        } // state 2
        
        if ( state == 3 ) {
          // We have been ordered to rise.
          //std::cout << "Rise.\n";
          Node * P = N -> parent_;
          // Can't rise if root.
          if ( P == NULL ) break; // algorithm complete
          int & div_dim = P -> dimension_;
          if ( P -> left_ == N ) {
            // This is a left child.
            //std::cout << "We are rising from left.\n";
            NUB[div_dim] += NUB[div_dim]-NLB[div_dim];
            // If we rise from the left child, we order parent to go right.
            state = 2;
          } else {
            // This is the right child.
            //std::cout << "We are rising from right.\n";
            NLB[div_dim] -= NUB[div_dim]-NLB[div_dim];
            // If we rise from the right child, we order parent to rise.
            state = 3;
          }
          N = P;
        } // state 3
        
      } // while loop
    }
    return ii;
  } // cover

template < class InsertIterator > inline InsertIterator
Toplex::cover ( InsertIterator & ii, const Prism & P ) const {
  //std::cout << "Prism version of Cover\n";
  static Rect G ( dimension_ );

  /* Use a stack, not a queue, and do depth first search.
   The advantage of this is that we can maintain the geometry during our Euler Tour.
   We can maintain our geometry without any roundoff error if we use the standard box
   [0,1]^d. To avoid having to translate to real coordinates at each leaf, we instead
   convert the input to these standard coordinates, which we put into integers. */
  
  // Step 1. Convert input to standard coordinates. 
  
#define INTPHASEWIDTH (((uint64_t)1) << 60)
  
  // Step 2. Perform DFS on the Toplex tree, recursing whenever we have intersection,
  //         (or adding leaf to output when we have leaf intersection)
  static std::vector<uint64_t> NLB ( dimension_);
  static std::vector<uint64_t> NUB ( dimension_);
  for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
    NLB [ dimension_index ] = 0;
    NUB [ dimension_index ] = INTPHASEWIDTH;
  }
  //std::cout << "Cover\n";
  
  /* Strategy. 
   We will take the Euler Tour using a 4-state machine.
   There are Four states.
   0 = Just Descended. Check for an intersection.
   1 = Descend to the left
   2 = Descend to right 
   3 = Rise.
   */
  
  Node * N = root_;
  char state = 0; 
  
  //std::cout << "Above main loop.\n";

  while ( 1 ) {
    //std::cout << "Entering Loop, state = " << (int) state << "\n";
    //std::cout << " N = " << N << "\n";
    if ( state == 0 ) {
      // If we have descended here, then we should check for intersection.
      
      // INTERSECTION CHECK
      for ( int d = 0; d < dimension_; ++ d ) {
        G . lower_bounds [ d ] = bounds () . lower_bounds [ d ] +
        (bounds () . upper_bounds [ d ] - bounds () . lower_bounds [ d ] ) * ( (Real) NLB [ d ] / (Real) INTPHASEWIDTH );
        G . upper_bounds [ d ] = bounds () . lower_bounds [ d ] +
        (bounds () . upper_bounds [ d ] - bounds () . lower_bounds [ d ] ) * ( (Real) NUB [ d ] / (Real) INTPHASEWIDTH );
      }
      
      //std::cout << "checking intersection:\n";
      if ( P . intersects ( G ) ) {
        //std::cout << "Detected intersection.\n";
        // Check if its a leaf.
        if ( N -> left_ == NULL ) {
          if ( N -> right_ == NULL ) {
            // Here's what we are looking for.
            * ii ++ = N -> contents_; // OUTPUT
            //std::cout << "cover -- " << N -> contents_ << "\n";
            // Issue the order to rise.
            //std::cout << "Issue rise.\n";
            state = 3;
          } else {
            // Issue the order to descend to the right.
            //std::cout << "Issue descend right.\n";
            state = 2;
          } 
        } else {
          // Issue the order to descend to the left.   
          //std::cout << "Issue descend left.\n";
          state = 1;
        }
      } else {
        // No intersection, issue order to rise.
        //std::cout << "No intersection. \n";
        //std::cout << "Issue Rise.\n";
        state = 3;
      } // intersection check complete
    } // state 0
    
    if ( state == 1 ) {
      // We have been ordered to descend to the left.
      //std::cout << "Descend left.\n";
      int div_dim = N -> dimension_;
      NUB[div_dim] -= ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
      N = N -> left_;
      state = 0;
      continue;
    } // state 1
    
    if ( state == 2 ) {
      // We have been ordered to descend to the right.
      //std::cout << "Descend right.\n";
      int & div_dim = N -> dimension_;
      NLB[div_dim] += ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
      N = N -> right_;
      state = 0;
      continue;
    } // state 2
    
    if ( state == 3 ) {
      // We have been ordered to rise.
      //std::cout << "Rise.\n";
      Node * P = N -> parent_;
      // Can't rise if root.
      if ( P == NULL ) break; // algorithm complete
      int & div_dim = P -> dimension_;
      if ( P -> left_ == N ) {
        // This is a left child.
        //std::cout << "We are rising from left.\n";
        NUB[div_dim] += NUB[div_dim]-NLB[div_dim];
        // If we rise from the left child, we order parent to go right.
        state = 2;
      } else {
        // This is the right child.
        //std::cout << "We are rising from right.\n";
        NLB[div_dim] -= NUB[div_dim]-NLB[div_dim];
        // If we rise from the right child, we order parent to rise.
        state = 3;
      }
      N = P;
    } // state 3
    
  } // while loop
  return ii;
} // cover

  
  // UNION version of cover
template < class InsertIterator, class T >
inline void Toplex::cover ( InsertIterator & ii, const std::vector < T > & V ) const {
  BOOST_FOREACH ( const T & geo, V ) {
    cover ( ii, geo );
  }
}

  // INTERSECTION (pair) for cover
  template < class InsertIterator, class S, class T >
  inline void Toplex::cover ( InsertIterator & ii, const std::pair < S, T > & V ) const {
    // Cover V . first, store in "firstcover"
    boost::unordered_set < GridElement > firstcover;
    std::insert_iterator < boost::unordered_set < GridElement > > fcii ( firstcover, firstcover . begin () );
    cover ( fcii, V . first );
    // Cover V . second, store in "secondcover"
    boost::unordered_set < GridElement > secondcover;
    std::insert_iterator < boost::unordered_set < GridElement > > scii ( secondcover, secondcover . begin () );
    cover ( scii, V . second );
    // Without loss, let firstcover be no smaller than secondcover.
    if ( firstcover . size () < secondcover. size () ) std::swap ( firstcover, secondcover );
    // Compute intersection by checking if each element in smaller cover is in larger cover
    BOOST_FOREACH ( const GridElement & ge, secondcover ) {
      if ( firstcover . count ( ge ) ) {
        // Use "ii" to insert "ge" into output.
        * ii ++ = ge;
      }
    }
  }
  
template < class InsertIterator >
inline InsertIterator
Toplex::coarseCover ( InsertIterator ii, const Rect & geometric_region ) const {
  
  /* Use a stack, not a queue, and do depth first search.
   The advantage of this is that we can maintain the geometry during our Euler Tour.
   We can maintain our geometry without any roundoff error if we use the standard box
   [0,1]^d. To avoid having to translate to real coordinates at each leaf, we instead
   convert the input to these standard coordinates, which we put into integers. */
  
  // Step 1. Convert input to standard coordinates. 
  Rect region ( dimension_ );
  static std::vector<uint64_t> LB ( dimension_);
  static std::vector<uint64_t> UB ( dimension_);
#define INTPHASEWIDTH (((uint64_t)1) << 60)
  static Real bignum ( INTPHASEWIDTH );
  for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
    region . lower_bounds [ dimension_index ] = 
    (geometric_region . lower_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]) /
    (bounds_ . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]);
    region . upper_bounds [ dimension_index ] = 
    (geometric_region . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]) /
    (bounds_ . upper_bounds [ dimension_index ] - bounds_ . lower_bounds [ dimension_index ]);
    
    if ( region . upper_bounds [ dimension_index ] < Real ( 0 ) ) return ii;
    if ( region . lower_bounds [ dimension_index ] > Real ( 1 ) ) return ii;
    
    if ( region . lower_bounds [ dimension_index ] < Real ( 0 ) ) 
      region . lower_bounds [ dimension_index ] = Real ( 0 );
    if ( region . lower_bounds [ dimension_index ] > Real ( 1 ) ) 
      region . lower_bounds [ dimension_index ] = Real ( 1 );
    LB [ dimension_index ] = (uint64_t) ( bignum * region . lower_bounds [ dimension_index ] );
    if ( region . upper_bounds [ dimension_index ] < Real ( 0 ) ) 
      region . upper_bounds [ dimension_index ] = Real ( 0 );
    if ( region . upper_bounds [ dimension_index ] > Real ( 1 ) ) 
      region . upper_bounds [ dimension_index ] = Real ( 1 );
    UB [ dimension_index ] = (uint64_t) ( bignum * region . upper_bounds [ dimension_index ] );
  }
  
  // Step 2. Perform DFS on the Toplex tree, recursing whenever we have intersection,
  //         (or adding leaf to output when we have leaf intersection)
  static std::vector<uint64_t> NLB ( dimension_);
  static std::vector<uint64_t> NUB ( dimension_);
  for ( int dimension_index = 0; dimension_index < dimension_; ++ dimension_index ) {
    if ( LB [ dimension_index ] < (1 << 20) ) LB [ dimension_index ] = 0;
    if ( LB [ dimension_index ] >= (1 << 20) ) LB [ dimension_index ] -= (1 << 20);
    if ( UB [ dimension_index ] < (INTPHASEWIDTH - (1 << 20)) ) UB [ dimension_index ] += (1 << 20);
    if ( UB [ dimension_index ] >= (INTPHASEWIDTH - (1 << 20)) ) UB [ dimension_index ] = INTPHASEWIDTH;
    NLB [ dimension_index ] = 0;
    NUB [ dimension_index ] = INTPHASEWIDTH;
  }
  //std::cout << "C\n";
  
  /* Strategy. 
   We will take the Euler Tour using a 4-state machine.
   There are Four states.
   0 = Just Descended. Check for an intersection.
   1 = Descend to the left
   2 = Descend to right 
   3 = Rise.
   */
  
  Node * N = root_;
  char state = 0; 
  
  //std::cout << "D\n";
  
  while ( 1 ) {
    //std::cout << "Entering Loop, state = " << (int) state << "\n";
    //std::cout << " N = " << N << "\n";
    if ( state == 0 ) {
      // If we have descended here, then we should check for intersection.
      bool intersect_flag = true;
      bool contain_flag = true;
      for ( int d = 0; d < dimension_; ++ d ) {
        if ( LB[d] > NUB[d] || UB[d] < NLB [d] ) {  // INTERSECTION CHECK
          intersect_flag = false;
          break;
        }
        if ( LB[d] > NLB[d] || UB[d] < NUB [d] ) {  // CONTAINMENT CHECK
          contain_flag = false;
        }
      }
      
      if ( intersect_flag && contain_flag ) {
        // Here's what we are looking for.
        // std::cout << "Detected containment.\n";
        * ii ++ = N -> contents_; // OUTPUT
                                  // std::cout << "coarsecover -- " << N -> contents_ << "\n";
                                  // Issue the order to rise.
                                  // std::cout << "Issue rise.\n";
        /*
        for ( int d = 0; d < dimension_; ++ d ) {
          std::cout << "[" << LB[d] << ", " << UB[d] << "] x ";
        }
        std::cout << "\n";
        for ( int d = 0; d < dimension_; ++ d ) {
          std::cout << "[" << NLB[d] << ", " << NUB[d] << "] x ";
        }
        std::cout << "\n";
        */
        
        state = 3;
      } else if ( intersect_flag ) {
        //std::cout << "Detected intersection.\n";
        // Check if its a leaf.
        if ( N -> left_ == NULL ) {
          if ( N -> right_ == NULL ) {
            // Here's what we are looking for.
            * ii ++ = N -> contents_; // OUTPUT
                                      //std::cout << "cover -- " << N -> contents_ << "\n";
                                      // Issue the order to rise.
                                      //std::cout << "Issue rise.\n";
            state = 3;
          } else {
            // Issue the order to descend to the right.
            //std::cout << "Issue descend right.\n";
            state = 2;
          } 
        } else {
          // Issue the order to descend to the left.   
          //std::cout << "Issue descend left.\n";
          state = 1;
        }
      } else {
        // No intersection, issue order to rise.
        //std::cout << "No intersection. \n";
        //std::cout << "Issue Rise.\n";
        state = 3;
      } // intersection check complete
    } // state 0
    
    if ( state == 1 ) {
      // We have been ordered to descend to the left.
      //std::cout << "Descend left.\n";
      int div_dim = N -> dimension_;
      NUB[div_dim] -= ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
      N = N -> left_;
      state = 0;
      continue;
    } // state 1
    
    if ( state == 2 ) {
      // We have been ordered to descend to the right.
      //std::cout << "Descend right.\n";
      int & div_dim = N -> dimension_;
      NLB[div_dim] += ( (NUB[div_dim]-NLB[div_dim]) >> 1 );
      N = N -> right_;
      state = 0;
      continue;
    } // state 2
    
    if ( state == 3 ) {
      // We have been ordered to rise.
      //std::cout << "Rise.\n";
      Node * P = N -> parent_;
      // Can't rise if root.
      if ( P == NULL ) break; // algorithm complete
      int & div_dim = P -> dimension_;
      if ( P -> left_ == N ) {
        // This is a left child.
        //std::cout << "We are rising from left.\n";
        NUB[div_dim] += NUB[div_dim]-NLB[div_dim];
        // If we rise from the left child, we order parent to go right.
        state = 2;
      } else {
        // This is the right child.
        //std::cout << "We are rising from right.\n";
        NLB[div_dim] -= NUB[div_dim]-NLB[div_dim];
        // If we rise from the right child, we order parent to rise.
        state = 3;
      }
      N = P;
    } // state 3
    
  } // while loop
  return ii;
} // coarseCover

template < class InsertIterator > inline InsertIterator
Toplex::subdivide ( InsertIterator ii, GridElement divide_me ) {
  return subdivide ( ii, find ( divide_me ) );
}

template < class InsertIterator >
inline InsertIterator
Toplex::subdivide ( InsertIterator ii, iterator cell_to_divide ) {
  cell_to_divide . node_ -> left_ = new Node;
  cell_to_divide . node_ -> right_ = new Node;
  cell_to_divide . node_ -> left_ -> parent_ = cell_to_divide . node_;
  cell_to_divide . node_ -> right_ -> parent_ = cell_to_divide . node_;
  int new_dim = (cell_to_divide . node_ -> dimension_ + 1 ) % dimension ();
  cell_to_divide . node_ -> left_ -> dimension_ = new_dim;
  cell_to_divide . node_ -> right_ -> dimension_ = new_dim;
  
  /* Update begin_, size_, tree_size_, find_ and initialize new nodes.*/
  if ( begin_ == cell_to_divide . node_ )
    begin_ . node_ = cell_to_divide . node_ -> left_;
  ++ size_;
  * ii ++ = cell_to_divide . node_ -> left_ -> contents_ = tree_size_ ++;
  * ii ++ = cell_to_divide . node_ -> right_ -> contents_ = tree_size_ ++;
  find_ . push_back ( const_iterator ( cell_to_divide . node_ -> left_ ) );
  find_ . push_back ( const_iterator ( cell_to_divide . node_ -> right_ ) );

#if 0
  //std::cout << "Subdivide: called on cell " << *cell_to_divide << "\n";
  std::deque < std::pair < const_iterator, int > > work_deque;
  work_deque . push_back ( std::pair < const_iterator, int >
                          (cell_to_divide,
                           cell_to_divide . node () -> dimension_ ) );
  while ( not work_deque . empty () ) {
    std::pair < const_iterator, int >  work_pair = work_deque . front ();
    work_deque . pop_front ();
    if ( work_pair . second < dimension_ ) { // NON-CYCLIC WORK TODO
      work_pair . first . node_ -> dimension_ = work_pair . second;
      /* We must subdivide further */
      work_pair . first . node_ -> left_ = new Node;
      work_pair . first . node_ -> right_ = new Node;
      /* Update begin_, size_, tree_size_, find_ and initialize new nodes.*/
      if ( begin_ == work_pair . first . node_ )
        begin_ . node_ = work_pair . first . node_ -> left_;
      ++ size_;
      work_pair . first . node_ -> left_ -> contents_ = tree_size_ ++;
      work_pair . first . node_ -> left_ -> parent_ = work_pair . first . node_;
      find_ . push_back ( const_iterator ( work_pair . first . node_ -> left_ ) );
      work_pair . first . node_ -> right_ -> contents_ = tree_size_ ++;
      work_pair . first . node_ -> right_ -> parent_ = work_pair . first . node_;
      find_ . push_back ( const_iterator ( work_pair . first . node_ -> right_ ) );
      /* Push the children onto the work_deque */
      work_deque . push_back ( std::pair < const_iterator, int >
                              (const_iterator ( work_pair . first . node_ -> left_ ),
                               work_pair . second + 1 ) );
      work_deque . push_back ( std::pair < const_iterator, int >
                              (const_iterator ( work_pair . first . node_ -> right_ ),
                               work_pair . second + 1 ) );
    } else {
      work_pair . first . node_ -> dimension_ = 0;
      //std::cout << "subdivide: inserting " << work_pair . first . node_ -> contents_ << "\n";
      * ii ++ = work_pair . first . node_ -> contents_;
    } /* if-else */
  } /* while */
#endif
  return ii;
} /* Adaptive_Cubical::Toplex::subdivide */

template < class InsertIterator, class Container >
inline InsertIterator
Toplex::subdivide ( InsertIterator ii, const Container & subset_to_divide ) {
  //std::cout << " ENTER \n";
  BOOST_FOREACH ( GridElement cell, subset_to_divide ) {
    //std::cout << "subdividing " << cell << "\n";
    ii = subdivide ( ii, find ( cell ) );
  }
  //std::cout << " EXIT \n";

  return ii;
}

template < class InsertIterator >
inline InsertIterator Toplex::subdivide ( InsertIterator ii ) {
  // TODO: CHECK THIS
  std::vector < GridElement > all;
  std::insert_iterator<std::vector<GridElement> > all_ii ( all, all . end () );
  cover ( all_ii, bounds () );
  return subdivide ( ii, all );
}

inline void Toplex::subdivide ( void ) {
  std::vector < GridElement > dummy;
  std::insert_iterator<std::vector<GridElement> > ii ( dummy, dummy . end () );
  subdivide ( ii );
}

/// depth
inline int Toplex::getDepth ( const GridElement & ge ) const {
  std::vector < unsigned char > p = prefix ( ge );
  return p . size ();// / dimension ();
}

template < class Container >
int Toplex::getDepth ( const Container & subset ) const {
  int depth = 0;
  BOOST_FOREACH ( const GridElement & ge, subset ) {
    int ge_depth = getDepth ( ge );
    if ( ge_depth > depth ) depth = ge_depth;
  }
  return depth;
}

/////////////////////// COARSEN /////////////////////////////
inline void branch ( std::vector < Node * > & nodes, Node * node ) {
  nodes . push_back ( node );
  // std::cout << "  pushing " << node << "  " << node -> contents_ << "\n";
  if ( node -> left_ != NULL ) branch ( nodes, node -> left_ );
  if ( node -> right_ != NULL ) branch ( nodes, node -> right_ );
}

template < class CellContainer >
inline void Toplex::coarsen ( const CellContainer & coarsen_to ) {
  std::vector < Node * > nodes;
  // Produce a list "nodes" of all descendant nodes
  BOOST_FOREACH ( GridElement cell, coarsen_to ) {
    Node * node = find ( cell ) . node ();
    //std::cout << "Coarsen to: " << node << "  " << node -> contents_ << "\n";
    // If this will be a new leaf and was not before, increment size
    if ( node -> left_ != NULL || node -> right_ != NULL ) {
      ++ size_;
      //std::cout << "  Will be a new leaf.\n";
    } 
    if ( node -> left_ != NULL ) branch ( nodes, node -> left_ );
    if ( node -> right_ != NULL ) branch ( nodes, node -> right_ );  
    node -> left_ = node -> right_ = NULL;
  }
  // Remove the top cells from the find structure
  BOOST_FOREACH ( Node * node, nodes ) {
    //std::cout <<  "Removing node " << node << "  " << node -> contents_ << "\n";
    find_ [ node -> contents_ ] = end ();
    // If we are deleting what used to be a leaf, decrement size
    if ( node -> left_ == NULL && node -> right_ == NULL ) -- size_;
    // Single delete (not recursive, so forget children first)
    node -> left_ = node -> right_ = NULL;
    delete node;
  }
  // Recompute begin_
  Node * n = root_;
  while ( n -> left_ != NULL || n -> right_ != NULL ) {
    while ( n -> left_ != NULL ) n = n -> left_;
    if ( n -> right_ != NULL ) n = n -> right_;
  }
  begin_ = iterator ( n );
}

inline void Toplex::GridElementToCubes ( std::vector<std::vector < uint32_t > > * cubes, 
                     const GridElement e, int depth ) const {
  // Obtain the prefix
  //std::cout << "GEtoCubes: " << geometry ( e ) << "\n";
  int D = dimension ();
  std::vector < unsigned char > p = prefix ( e );
  int GridElement_depth = p . size (); // == getDepth ( e );
  //std::cout << "  gedepth = " << GridElement_depth << ", from " << p . size () << "\n";
  if ( GridElement_depth > depth ) GridElement_depth = depth; //effectively truncates the prefix
  // Determine width
  typedef std::vector < uint32_t > Cube;
  Cube cube ( D, 0 );
  int pos = 0;
  
  // NON-CYCLIC OLD CODE
  //for ( int d = 0; d < GridElement_depth; ++ d ) {
  //  for ( int dim = 0; dim < D; ++ dim ) {
  //    cube [ dim ] <<= 1;
  //    cube [ dim ] |= (uint32_t) p [ pos ++ ];
  //  }
  //}

  int dim = 0;
  for ( int d = 0; d < GridElement_depth; ++ d ) {
    if ( dim == D ) dim = 0;
    cube [ dim ] <<= 1;
    cube [ dim ] |= (uint32_t) p [ pos ++ ];
    ++ dim;
  }
  // make the cubes
  if ( GridElement_depth == depth ) {
    cubes -> push_back ( cube );
    //std::cout << "    simple case.\n";
    return;
  }

  // We must make more than one output cube; 
  // the user has requested a greater depth than
  // the toplex provides.
  //std::cout << "    hard case.\n";
  
  std::vector < Cube > work_stack, split_stack;
  work_stack . push_back ( cube );
  for ( int dim = GridElement_depth; dim < depth; ++ dim ) {
    std::vector < Cube > split_stack;
    BOOST_FOREACH ( Cube c, work_stack ) {
      c [ dim % D ] <<= 1;
      split_stack . push_back ( c );
      c [ dim % D ] |= 1;
      split_stack . push_back ( c );
    }
    std::swap ( work_stack, split_stack );
  }
  BOOST_FOREACH ( const Cube & outcube, work_stack ) {
    cubes -> push_back ( outcube );
  }
  
}

template < class Container >
inline void Toplex::relativeComplex ( RelativePair * pair,
                                      const Container & XGridElements,
                                      const Container & AGridElements,
                                      int depth ) const {
  // DEBUG
  // check A \subset X
  /*
  boost::unordered_set < GridElement > XS, AS;
  BOOST_FOREACH ( GridElement x, XGridElements ) XS . insert ( x );
  BOOST_FOREACH ( GridElement a, AGridElements ) {
    AS . insert ( a );
    if ( XS . count ( a ) == 0 ) {
      std::cout << "couldn't find " << a << "\n";
      exit ( 1 );
    }
  }
  */
  
  // Produce the full complex.
  CubicalComplex * full_complex = new CubicalComplex;
  CubicalComplex & X = *full_complex;
  int D = dimension ();

  //std::cout << "relativeComplex.\n";
  //std::cout << "XGridElements:\n";
  
  // Make set of cubes, and learn bounds of the cubes.
  typedef std::vector < uint32_t > Cube;
  Cube mincube ( D, -1 );
  Cube maxcube ( D, 0 );
  Rect newbounds ( D );
  for ( int d = 0; d < D; ++ d ) {
		newbounds . lower_bounds [ d ] = bounds () . upper_bounds [ d ];
		newbounds . upper_bounds [ d ] = bounds () . lower_bounds [ d ];
  }
  BOOST_FOREACH ( GridElement e, XGridElements ) {
    Rect geo = geometry ( e );
    for ( int d = 0; d < D; ++ d ) {
    if ( newbounds . lower_bounds [ d ] > geo . lower_bounds [ d ] )
    	newbounds . lower_bounds [ d ] = geo . lower_bounds [ d ];
    if ( newbounds . upper_bounds [ d ] < geo . upper_bounds [ d ] )
    	newbounds . upper_bounds [ d ] = geo . upper_bounds [ d ];
    	}
    std::vector < Cube > cubes;
    GridElementToCubes ( &cubes, e, depth );
    BOOST_FOREACH ( Cube & cube, cubes ) {
      for ( int d = 0; d < D; ++ d ) {
      	if ( mincube [ d ] > cube [ d ] ) mincube [ d ] = cube [ d ];
      	if ( maxcube [ d ] < cube [ d ] ) maxcube [ d ] = cube [ d ];
      }
    }
  }
  
  std::vector < uint32_t > dimension_sizes ( D, 1 );
  std::vector < bool > is_periodic = periodic_;
  for ( int d = 0; d < D; ++ d ) {
    dimension_sizes [ d ] = maxcube [ d ] - mincube [ d ] + 1;
    if ( newbounds . lower_bounds [ d ] > bounds () . lower_bounds [ d ] )
    	is_periodic [ d ] = false;
    if ( newbounds . upper_bounds [ d ] < bounds () . upper_bounds [ d ] )
    	is_periodic [ d ] = false;
  }
  
  X . bounds () = newbounds;
  X . initialize ( dimension_sizes, is_periodic );
  
  BOOST_FOREACH ( GridElement e, XGridElements ) {
    //std::cout << e << "\n";
    //std::cout << "GEOMETRY = " << geometry ( e ) << "\n"; 
    typedef std::vector < uint32_t > Cube;
    std::vector < Cube > cubes;
    GridElementToCubes ( &cubes, e, depth );
    BOOST_FOREACH ( Cube & cube, cubes ) {
      //std::cout << "XCube = ";
      //for ( int d = 0; d < dimension (); ++ d ) std::cout << cube [ d ] << " ";
      //std::cout << "\n";
      Cube offset ( D );
      for ( int d = 0; d < D; ++ d ) offset [ d ] = cube [ d ] - mincube [ d ];
      X . addFullCube ( offset );
      //std::cout << " CC-Geometry = " << X . geometryOfCube ( cube ) << "\n";
    }
  }
  X . finalize ();
  
  // Produce the relative complex
  BitmapSubcomplex * pair_complex = new BitmapSubcomplex ( X, true );
  BitmapSubcomplex * rel_complex = new BitmapSubcomplex ( X, false );
  BitmapSubcomplex & XA = * pair_complex;
  BitmapSubcomplex & A = * rel_complex;
  //std::cout << "AGridElements:\n";
  BOOST_FOREACH ( GridElement e, AGridElements ) {
    //std::cout << e << "\n";
    typedef std::vector < uint32_t > Cube;
    std::vector < Cube > cubes;
    GridElementToCubes ( &cubes, e, depth );
    BOOST_FOREACH ( Cube & cube, cubes ) {
      //std::cout << "ACube = ";
      //for ( int d = 0; d < dimension (); ++ d ) std::cout << cube [ d ] << " ";
      //std::cout << "\n";
      Cube offset ( D );
      for ( int d = 0; d < D; ++ d ) offset [ d ] = cube [ d ] - mincube [ d ];
      
      std::vector < std::vector < Index > > cells = 
      X . fullCubeIndexes ( offset );
      for ( int d = 0; d <= dimension (); ++ d ) {
        BOOST_FOREACH ( Index cell, cells [ d ] ) {
          XA . erase ( cell, d );
          A . insert ( cell, d );
        }
      }
    }
  }
  XA . finalize ();
  A . finalize ();
  pair -> initialize ( &X, &XA, &A ); 
}

inline void Toplex::initialize ( const Rect & outer_bounds_of_toplex ) {
  if ( root_ != NULL ) clear ();
  dimension_ = outer_bounds_of_toplex . lower_bounds . size ();
  periodic_ . resize ( dimension_, false );
  bounds_ = outer_bounds_of_toplex;
  root_ = new Node;
  tree_size_ = 1;
  size_ = 1;
  begin_ = const_iterator ( root_ );
  find_ . push_back ( begin_ );
} /* Adaptive_Cubical::Toplex::initialize */

inline void Toplex::initialize ( const Rect & outer_bounds_of_toplex,
                                const std::vector < bool > & periodic ) {
  initialize ( outer_bounds_of_toplex );
  periodic_ = periodic;
}

inline Toplex::Toplex ( void ) {
  end_ = const_iterator ( NULL );
  begin_ = end_;
  size_ = 0;
  tree_size_ = 0;
  root_ = NULL;
  dimension_ = 0;
} /* Adaptive_Cubical::Toplex::Toplex */

inline Toplex::Toplex ( const Rect & outer_bounds_of_toplex ) {
  end_ = const_iterator ( NULL );
  begin_ = end_;
  size_ = 0;
  tree_size_ = 0;
  root_ = NULL;
  dimension_ = 0;
  initialize ( outer_bounds_of_toplex );
} /* Adaptive_Cubical::Toplex::Toplex */

inline Toplex::~Toplex ( void ) {
  if ( root_ != NULL ) delete root_;
} /* Adaptive_Cubical::Toplex::~Toplex */

} // namespace chomp

#endif
