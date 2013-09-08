// MorseComplex.h
// Shaun Harker
// 9/10/11

#ifndef CHOMP_MORSECOMPLEX_H
#define CHOMP_MORSECOMPLEX_H

#include <cstdlib>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include "boost/unordered_set.hpp"
#include "boost/unordered_map.hpp"

#include "chomp/Chain.h"
#include "chomp/Complex.h"
#include "chomp/Decomposer.h"
#include "chomp/CoreductionDecomposer.h"

namespace chomp {
  
class MorseComplex : public Complex {

public:      
  MorseComplex ( void ) {}  
  MorseComplex ( Complex & base );
  virtual ~MorseComplex ( void );

  // Boundary and Coboundary
  virtual void boundary ( Chain * output, const Index input, int dim ) const;
  virtual void coboundary ( Chain * output, const Index input, int dim ) const;

  // Morse Theory
  template < class Decomposer >
  void initialize ( Complex & base );
  void include ( Chain * output, const Chain & input ) const;
  void project ( Chain * output, const Chain & input ) const;
  void flow ( Chain * canonical, Chain * gamma, const Chain & input ) const;
  void lift ( Chain * output, const Chain & input ) const;
  void lower ( Chain * output, const Chain & input ) const;
  void colift ( Chain * output, const Chain & input ) const;
  void colower ( Chain * output, const Chain & input ) const;
  
  // Convenience:

  Chain include ( const Chain & input ) const {
    Chain output; include ( &output, input ); return output;
  }
  Chain project ( const Chain & input ) const {
    Chain output; project ( &output, input ); return output;
  }  
  Chain lift ( const Chain & input ) const {
    Chain output; lift ( &output, input ); return output;
  }
  Chain lower ( const Chain & input ) const {
    Chain output; lower ( &output, input ); return output;
  }
  Chain colift ( const Chain & input ) const {
    Chain output; colift ( &output, input ); return output;
  }
  Chain colower ( const Chain & input ) const {
    Chain output; colower ( &output, input ); return output;
  }
  // Cell type = Index
  CHOMP_COMPLEX(Index)
//private:
  // Base Complex
  Complex * base_;
  Complex & base () { return *base_; }
  const Complex & base () const { return *base_; }
  // Decomposition
  Decomposer * decomposer_; 
  int type ( Index i, int d ) const { 
    return decomposer_ -> type ( i, d );
  }
  Index mate ( Index i, int d ) const { 
    return decomposer_ -> mate ( i, d );
  }
  // Boundary/Coboundary storage
  mutable std::vector < boost::unordered_map < Index, Chain > > 
  boundary_cache_;
  
  mutable std::vector < boost::unordered_map < Index, Chain > > 
  coboundary_cache_;

};

/***************
 * Definitions *
 ***************/

/**************************
 * Morse Complex Creation *
 **************************/
inline MorseComplex::MorseComplex ( Complex & base ) {
  initialize < CoreductionDecomposer > ( base );
}

inline MorseComplex::~MorseComplex ( void ) {
  if ( decomposer_ != NULL ) delete decomposer_;
}

template < class Decomposer >
void MorseComplex::initialize ( Complex & basearg ) {
  // Set the base complex
  base_ = &basearg;
    
  // Call the decomposition algorithm
  decomposer_ = new Decomposer ( base () );

  // Produce the complex
  //startInserting ();
  for ( int d = 0; d <= base () . dimension (); ++ d ) {
    Cell hoist = base () . size ( d );
    for ( Cell i = 0; i < hoist; ++ i ) {
      if ( type ( i, d ) == Decomposer::ACE ) insertCell ( i, d ); 
    }
  }
  //finishedInserting ();
  boundary_cache_ . resize ( dimension () + 1 );
  coboundary_cache_ . resize ( dimension () + 1 );
  
}

/************************
 * Morse Theory Algebra *
 ************************/

inline void MorseComplex::include ( Chain * output, 
                                    const Chain & input ) const {
  //std::cout << "morsecomplex INCLUDE " << input . dimension () << "\n";
  int dim = input . dimension ();
  BOOST_FOREACH ( Term term, input () ) {
    term . index () = indexToCell ( term . index (), dim );
    (*output) += term;
  }
  output -> dimension () = dim;
}

inline void MorseComplex::project ( Chain * output, 
                                   const Chain & input ) const {
  //std::cout << "PROJECT " << input . dimension () << "\n";
  int dim = input . dimension ();
  BOOST_FOREACH ( Term term, input () ) {
    if ( type ( term . index (), dim ) == Decomposer::ACE ) {
      term . index () = cellToIndex ( term . index (), dim );
      (*output) += term;
    }
  }
  output -> dimension () = dim;
}

namespace flow_detail {
  class CellCompare {
    Decomposer * decomposer_;
    int d_;
  public:
    CellCompare ( Decomposer * decomposer, int d ) :
    decomposer_(decomposer), d_(d) {}
    bool operator () ( const Index & lhs, const Index & rhs ) const {
      return decomposer_ -> compare ( rhs, lhs, d_ ); //reversal intentional
    }
  };
}

inline void MorseComplex::flow ( Chain * canonical, 
                                 Chain * gamma, 
                                 const Chain & input ) const {
  //std::cout << "FLOW " << input . dimension () << "\n";

  // Set the dimensions on the output chains
  int D = input . dimension ();
  canonical -> dimension () = D;
  gamma -> dimension () = D + 1;
  
  // data structures
  boost::unordered_map < Cell, Ring > work;
  boost::unordered_set < Cell > queens;

  flow_detail::CellCompare comparator ( decomposer_, D );
  std::priority_queue<Cell, std::vector<Cell>, flow_detail::CellCompare > priority ( comparator );


  // A MACRO TO ADD CHAINS TO THE "work" CHAIN
#define ADDTOWORKCHAIN(chain)                               \
BOOST_FOREACH ( const Term & t, chain () ) {                \
  if ( work . find ( t . index () ) == work . end () )      \
    work [ t . index () ] = Ring ( 0 );                     \
  if ( type ( t . index (), D ) == Decomposer::QUEEN        \
    && queens . count ( t . index () ) == 0           ) {   \
    queens . insert ( t . index () );                       \
    priority . push ( t . index () );                       \
  }                                                         \
  work [ t . index () ] += t . coef ();                     \
  if ( work [ t . index () ] == Ring ( 0 ) ) {              \
    work . erase ( t . index () );                          \
  }                                                         \
}



  // Copy input into data structures
  ADDTOWORKCHAIN(input)

  // Gamma algorithm
  while ( not priority . empty () ) {
    Cell queen = priority . top ();
    priority . pop ();
    if ( work . count ( queen ) == 0 ) continue;
    Ring queen_coef = work [ queen ];
    Cell king = mate ( queen, D );
    Chain king_bd = base () . boundary ( king, D + 1 );
    // Get coefficient unit := <Q, dK>
    Ring unit (0);
    BOOST_FOREACH ( const Term & king_bd_term, king_bd () ) {
      if ( king_bd_term . index () == queen ) {
        unit += king_bd_term . coef ();
      }
    } // for each king boundary
    
    // debug
    if ( not invertible ( unit ) ) {
      // problem.
      std::cout << "Flow problem.\n";
      std::cout << "King: " << king << "\n";
      std::cout << "Queen: " << queen << "\n";
      std::cout << "unit: " << unit << "\n";
      std::cout << "bd(king) = " << king_bd << "\n";
      exit ( 1 );
    } //else { std::cout << unit << " "; }
    
    Ring factor = - queen_coef  / unit;
    (*gamma) += Term ( king, factor );
    king_bd *= factor;
    ADDTOWORKCHAIN ( king_bd )
  } // while queens are in the boundary
  
  // Write canonical chain from "work"
  BOOST_FOREACH ( const Term & term, work ) (*canonical) += term;
  
}

inline void MorseComplex::lift ( Chain * output, const Chain & input ) const {
  //std::cout << "morsecomplex LIFT " << input . dimension () << "\n";
  Chain & answer = *output;
  // Include the chain into the original complex
  answer = include ( input );
  // Take the boundary of the chain
  Chain answer_bd = base () . boundary ( answer );
  // Apply gradient flow to the boundary chain
  Chain canonical; Chain gamma;
  flow ( &canonical, &gamma, answer_bd );
  // Return the sum of the original input and "gamma"
  answer += gamma;
}


inline void MorseComplex::lower ( Chain * output, const Chain & input ) const {
  //std::cout << "LOWER " << input . dimension () << "\n";
  Chain canonical;
  Chain gamma;
  flow ( &canonical, &gamma, input );
  project ( output, canonical );
}

inline void MorseComplex::colift ( Chain * output, const Chain & input ) const {
}


inline void MorseComplex::colower ( Chain * output, const Chain & input ) const {
}

/// DEBUG methods

inline void MorseSanity ( const MorseComplex & c ) {
  std::cout << "MorseSanity Check.\n";
  const Complex & base = c . base ();
  int D = base . dimension ();
  for ( int d = 0; d <= D; ++ d ) {
    Index N = base . size ( d );
    std::cout << "Dimension " << d << " has " << N << " cells.\n";
    for ( Index i = 0; i < N; ++ i ) {
      std::cout << "Index " << i << "\n";
      if ( c.type ( i, d ) == Decomposer::QUEEN ) {
        std::cout << "Queen.\n";
        if ( c.type ( c.mate ( i, d ), d + 1 ) != Decomposer::KING ) {
          std::cout << "problem with morse decomposition.\n";
          std::cout << "Queen(i,d) = (" << i << ", " << d << ")\n"; 
          std::cout << "Mate(i,d) = (" << c.mate(i,d) << ", " << d + 1 << ")\n"; 
          std::cout << c.type ( c.mate ( i, d ), d + 1 ) << "\n";
        }
      }
      if ( c.type ( i, d ) == Decomposer::KING ) {
        std::cout << "King.\n";
        if ( c.type ( c.mate ( i, d ), d - 1 ) != Decomposer::QUEEN ) {
          std::cout << "problem with morse decomposition.\n";
          std::cout << "King(i,d) = (" << i << ", " << d << ")\n"; 
          std::cout << "Mate(i,d) = (" << c.mate(i,d) << ", " << d - 1 << ")\n"; 
          std::cout << c.type ( c.mate ( i, d ), d - 1 ) << "\n";

        }
      }
    }
  }
}
 

/*******************
 * Complex Methods *
 *******************/
/// Boundary
inline void MorseComplex::boundary ( Chain * output, 
                                     const Index input, 
                                     int dim ) const {
  if ( boundary_cache_ [ dim ] . count ( input ) == 0 ) {
    //debug
    //std::cout << "debug1: " << base () . boundary ( indexToCell ( input, dim ), dim ) << "\n";
    boundary_cache_ [ dim ] [ input ] =
    lower ( base () . boundary ( indexToCell ( input, dim ), dim ) );
  }
  *output = boundary_cache_ [ dim ] [ input ];
}

/// Coboundary
inline void MorseComplex::coboundary ( Chain * output, 
                                       const Index input,
                                       int dim ) const {
  std::cout << "MorseComplex::coboundary Error, coboundary not " << 
               "working until colower and colift written\n";
  exit ( 1 );
  if ( coboundary_cache_ [ dim ] . count ( input ) == 0 ) {
    coboundary_cache_ [ dim ] [ input ] = 
    colower ( base () . coboundary ( indexToCell ( input, dim ), dim ) );
  }
  *output = coboundary_cache_ [ dim ] [ input ];
}

} // namespace chomp

#endif
