// CoreductionDecomposer.h
// Shaun Harker
// 9/13/11


/// TODO:
/// Inefficient mistake: I didn't bother ordering queens
/// according to their discrete morse functions.

#ifndef CHOMP_COREDUCTIONDECOMPOSER_H
#define CHOMP_COREDUCTIONDECOMPOSER_H 

#include <iostream>
#include <cstdlib>
#include <vector>
#include <stack>

#include "chomp/Decomposer.h"
#include "chomp/Complex.h"
#include "chomp/Chain.h"

namespace chomp {
  
/**************************
 *  CoreductionDecomposer *
 **************************/

class CoreductionDecomposer : public Decomposer {
private:
  std::vector < uint64_t > queen_end_; 
  std::vector < uint64_t > ace_end_;
  std::vector < uint64_t > king_rbegin_;
  std::vector < std::vector < uint64_t > > toOriginalIndexing_;
  std::vector < std::vector < uint64_t > > permute;

public:

  int type ( Index j, int d ) const { 
    Index i = permute [ d ] [ j ];
    if ( i < queen_end_ [ d ] ) return QUEEN;
    if ( i < ace_end_ [ d ] ) return ACE;
    return KING;
  }
  
  Index mate ( Index j, int d ) const { 
    int t = type (j, d);
    Index i = permute [ d ] [ j ];
    if ( t == QUEEN ) {
     return toOriginalIndexing_ [ d + 1 ] [ king_rbegin_ [ d + 1 ] - i ]; 
    }
    if ( t == ACE ) {
      return j;
    }
    //if ( t == KING) 
    return toOriginalIndexing_ [ d - 1 ] [ king_rbegin_ [ d ] - i ];
  }
  
  bool compare ( Index lhs, Index rhs, int d ) const {
    return permute[d][lhs] < permute[d][rhs];
  }
  
  void decompose ( Complex & complex );
  
  CoreductionDecomposer ( Complex & complex ) : Decomposer(complex) {
    decompose ( complex );
  }
  
  virtual ~CoreductionDecomposer ( void ) {}
};

/****************
 *  Definition  *
 ****************/

/*
 Algorithm Description.
 This is a greedy algorithm which finds all "coreduction pairs"
 and marks them as K-Q pairs. When no pairs are left to find, an arbitrary
 cell with no boundary as chosen as an ace. Then process repeats until all
 cells have been marked as either Ace, King, or Queen.
 
 Implementation details:
 The computational burden is to efficiently be aware of coreduction pairs.
 A coreduction pair is easily identified by knowing the set of cells
 with boundaries dK = uQ for some unit u and cell Q. 
 NOTE: If ever we discover a cell K with boundary dK = nQ, where n is NOT
 a unit, then we immediately know that K and Q should be aces.
 The method used is to keep track of the number of boundary elements of each
 cell. If a cell is labeled A, K, or Q, then it is no longer considered in this
 count, so whenever a cell is labeled, we must update the counts of the coboundary.
 
 After all cells have been processed and the decomposition is known, it must be
 stored in some fashion. This is done by invalidating the indexing of the complex
 and replacing it with one where the cells are stored in the ordering
 QAK, QAK, QAK, with the Queens and Kings corresponding to each other in reverse
 order. The rationale is the following:
 1) To do Morse Boundary, we need to find maximal queens. By simply ordering by
    indexing, this can be accomplished.
 2) By storing Queens and Kings in reverse order, we know
    indexToCell ( i, d )   MATCHES WITH  indexToCell ( size ( d + 1 ) - i - 1, d+1)
 3) We have to provide some extra information. The most convenient
    solution, given kings are iterated in reverse compared to associated queens,
    is to give the indices one past the end of the queens (which is the start
    of the aces), one past the end of the aces, (which is the end of the kings,
    iterating backwards) and the final king -- (rbegin). 
 
 Technical questions:
 1) Can a cell be placed into the king_candidate stack twice?
    Answer: No.
 2) Can a cell be placed into the ace_candidate stack twice?
    Answer: No.
 3) Can a cell be in the king_candidate stack and no longer be 
    a candidate for a king?
    Answer. Yes, but only in one way: its queen could be claimed by
    another cell. In fact, we may as well make this the only way a cell 
    gets into the ace candidate stack, apart from the initial sweep.
*/
inline void CoreductionDecomposer::decompose ( Complex & complex ) {
  /************************
   * INITIALIZE ALGORITHM *
   ************************/
  int D = complex . dimension ();

  // Prepare output variables
  queen_end_ . resize ( D + 1 );
  ace_end_ . resize ( D + 1 );
  king_rbegin_ . resize ( D + 1 );
  
  std::vector < std::vector < Index > > aces ( D + 1 );
  std::vector < std::vector < Index > > kings ( D + 2 );
  std::vector < std::vector < Index > > queens ( D + 1 );
  
  // Declare a structure to hold the boundary counts.
  std::vector < std::vector < int > > 
  bd_count ( D + 2 );
  
  // Initialize the boundary count structure. 
  for ( int d = 0; d <= D; ++ d ) {
    bd_count [ d ] . resize ( complex . size ( d ) );
    for ( Index i = 0; i < complex . size ( d ); ++ i ) {
      Chain bd = complex . boundary ( i, d );
      bd_count [ d ] [ i ] = bd () . size ();
    }
  }

  // Declare structures to keep track of cells 
  // which might become kings or aces.
  std::queue < std::pair < Index, int > > king_candidates;
  std::queue < std::pair < Index, int > > ace_candidates;
  
  // Initialize those structures. 
  for ( int d = 0; d <= D; ++ d ) {
    //std::cout << "Dim " << d << " has " << complex . size ( d ) << " cells.\n";
    for ( Index i = 0; i < complex . size ( d ); ++ i ) {
      std::pair < Index, int > candidate = std::make_pair ( i, d );
      if ( bd_count [ d ] [ i ] == 1 ) 
        king_candidates . push ( candidate );
      if ( bd_count [ d ] [ i ] == 0 ) 
        ace_candidates . push ( candidate );
    }
  }
  
  // Declare structures to keep track of which cells are paired
  std::vector < std::vector < bool > > paired;
  paired . resize ( D + 1 );
  for ( int d = 0; d <= D; ++ d ) {
    paired [ d ] . resize ( complex . size ( d ), false );
  }
  
  
  // A macro to clean up the code: 
  // This is the procedure to lower the boundary counts 
  // of all cells in the coboundary by one.
#define CBDPROCESS(i,d)                                     \
{                                                           \
  paired [ d ] [ i ] = true;                                \
  Chain cbd = complex . coboundary ( i, d );                \
  BOOST_FOREACH ( const Term & t, cbd () ) {                \
    Index newcount = -- bd_count [ d + 1 ] [ t . index () ];\
    if ( newcount == 1 )                                    \
      king_candidates . push ( std::make_pair (t . index (), d + 1) ); \
    if ( newcount == 0 )                                               \
      ace_candidates . push ( std::make_pair (t . index (), d + 1) );  \
  }                                                                    \
}                                                       

  /************* 
   * MAIN LOOP *
   *************/
  while ( (not ace_candidates . empty ()) || (not king_candidates . empty ()) ) {
    
    /********************** 
     * KING-QUEEN ROUTINE *
     **********************/
    //std::cout << "There are " << king_candidates . size () << " pairs\n";
    while ( not king_candidates . empty () ) {
      std::pair < Index, int > candidate = king_candidates . front ();
      king_candidates . pop ();

      Index ki = candidate . first;
      int kd = candidate . second;
      // Determine the connection
      Chain bd = complex . boundary ( ki, kd );
      Index qi;
      int qd = kd - 1;
      Ring connection ( 0 ); // If queen is already paired, this remains 0
      BOOST_FOREACH ( const Term & t, bd () ) {
        if ( paired [ qd ] [ t . index () ] == false ) {
          qi = t . index ();
          connection += t . coef ();
        }
      }
      // If the connection is a unit, mate them:
      if ( invertible ( connection ) ) {
        queens [ qd ] . push_back ( qi );
        kings [ kd ] . push_back ( ki );
        CBDPROCESS ( ki, kd );
        CBDPROCESS ( qi, qd );
        //std::cout << "pushing (" <<qi<<", "<<qd<< ") -- ("<<ki<< ", "<<kd<<")\n";
      } 
    }
    
    /*************** 
     * ACE ROUTINE *
     ***************/
    
    while ( king_candidates . empty () ) {
      if ( ace_candidates . empty () ) break;
      std::pair < Index, int > candidate = ace_candidates . front ();
      ace_candidates . pop ();
      Index ai = candidate . first;
      int ad = candidate . second;
      if ( not paired [ ad ] [ ai ] ) {
        CBDPROCESS ( ai, ad );
        aces [ ad ] . push_back ( ai );
        //std::cout << "pushing ace (" << ai << ", " << ad << ")\n";

      }
    }
  } /* while */
  
  /*********************** 
   * MORSE INDEX ROUTINE *
   ***********************/

  permute . resize ( D + 1 );
  for ( int d = 0; d <= D; ++ d ) {
    permute [ d ] . resize ( complex . size ( d ) );
  }
  for ( int d = 0; d <= D; ++ d ) {
    std::vector < Index > & A = aces [ d ];
    std::vector < Index > & K = kings [ d + 1 ];
    std::vector < Index > & Q = queens [ d ];

    uint64_t ai = Q . size ();
    uint64_t ki = complex . size ( d + 1 ) - K . size ();
    uint64_t qi = Q . size ();

    for ( uint64_t i = 0; i < Q . size (); ++ i ) {
      permute [ d ] [ Q [ i ] ] = -- qi;
      permute [ d + 1 ] [ K [ i ] ] = ki ++;
    }
    for ( uint64_t i = 0; i < A . size (); ++ i ) {
      permute [ d ] [ A [ i ] ] = ai ++;
    }
    
    queen_end_[d] = Q . size ();
    ace_end_[d] = A . size () + Q . size ();
    king_rbegin_[d] = complex . size ( d ) - 1;
  }

  toOriginalIndexing_ . resize ( D + 1 );
  for ( int d = 0; d <= D; ++ d ) {
    Index N = complex . size ( d );
    toOriginalIndexing_ [ d ] . resize ( N );
    for ( Index i = 0; i < N; ++ i ) {
      toOriginalIndexing_ [ d ] [ permute [ d ] [ i ] ] = i;
    }
  }
  
  // DEBUG
#if 0
  for ( int d = 0; d <= D; ++ d ) {
    for ( Index i = 0; i < complex . size ( d ); ++ i ) {
      if ( type ( i, d ) == QUEEN ) {
        Index j = mate ( i, d );
        Chain bd = complex . boundary ( j, d + 1 );
        if ( type ( j, d + 1 ) != KING ) {
          std::cout << "CD: mate issue (1).\n";
          exit ( 1 );
        }
        bool flag = false;
        BOOST_FOREACH ( const Term & t, bd () ) {
          if ( permute [ d ] [ t . index () ] < permute [ d ] [ i ] ) {
            std::cout << "COREDUCTION DECOMPOSER: Possible acyclicity violation.\n";
            exit ( 1 );
          }
          if ( t . index () == i ) flag = true;
        }
        if ( not flag ) {
          std::cout << "QUEEN is not paired with mate!\n";
          exit ( 1 );
        }
      }
      if ( type ( i, d ) == KING ) {
        Index j = mate ( i, d );
        if ( type ( j, d - 1 ) != QUEEN ) {
          std::cout << "CD: mate issue (2).\n";
          exit ( 1 );
        }
      }
    }
  }
#endif
  // ENDDEBUG
  
  //complex . reindex ( permute );

}

} // namespace chomp
#endif

