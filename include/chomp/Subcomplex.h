// Subcomplex.h
// Shaun Harker
// 9/21/11

#ifndef CHOMP_SUBCOMPLEX_H
#define CHOMP_SUBCOMPLEX_H

#include "chomp/Complex.h"
#include "chomp/Chain.h"

namespace chomp {

/// Subcomplex.
/// Currently this isn't a true complex
/// as it doesn't support "enumeration"
/// cellToIndex and indexToCell are trivialized --
/// they need to switch their behavior if one
/// enumerates; this isn't coded yet.


class CellSubset {
 public:
  virtual bool operator () ( const Index i, const int d ) = 0;
};

class Subcomplex : public Complex {
 private:
  Complex * supercomplex;
  mutable CellSubset * cells_;
  bool query ( const Index i, int d ) const 
  { return (*cells_)(i,d); } 
public:
  Subcomplex ( Complex * supercomplex, CellSubset * cells ) 
  : supercomplex_(supercomplex), cells_(cells) {};
  virtual ~Subcomplex ( void ) {}; 
  
  // Inclusion and Projection
  void include ( Chain * output, const Chain & input ) const;
  void project ( Chain * output, const Chain & input ) const;

  // Boundary and Coboundary
  virtual void boundary ( Chain * output, 
                          const Index input, int dim ) const;
  virtual void coboundary ( Chain * output, 
                            const Index input, int dim ) const;
  
  //
  Chain include ( const Chain & input ) const {
    Chain output; include ( &output, input ); return output;
  }
  Chain project ( const Chain & input ) const {
    Chain output; project ( &output, input ); return output;
  }  
  
  // cellToIndex and indexToCell
  // TODO: allow enumeration.
  typedef Index Cell;
  Index cellToIndex ( const Index i ) { return i; }
  Index indexToCell ( const Index i ) { return i; }

};

/********************
 *   DEFINITIONS    *
 ********************/

inline void Subcomplex::include ( Chain * output, 
                                   const Chain & input ) const {
  //std::cout << "INCLUDE " << input . dimension () << "\n";
  int dim = input . dimension ();
  BOOST_FOREACH ( Term term, input () ) {
    term . index () = indexToCell ( term . index (), dim );
    (*output) += term;
  }
  output -> dimension () = dim;
}

inline void Subcomplex::project ( Chain * output, 
                                   const Chain & input ) const {
  //std::cout << "PROJECT " << input . dimension () << "\n";
  int dim = input . dimension ();
  BOOST_FOREACH ( Term term, input () ) {
    if ( query ( type ( term . index (), dim ) ) ) {
      term . index () = cellToIndex ( term . index (), dim );
      (*output) += term;
    }
  }
  output -> dimension () = dim;
}

inline void Subcomplex::boundary 
( Chain * output, const Index input, int dim ) const {
  Chain input_chain;
  input_chain += Term ( input, Ring ( 1 ) );
  input_chain . dimension () = dim;
  boundary ( output, input_chain );
}

inline void Subcomplex::coboundary 
( Chain * output, const Index input, int dim ) const {
  Chain input_chain;
  input_chain += Term ( input, Ring ( 1 ) );
  input_chain . dimension () = dim;
  coboundary ( output, input_chain );
}

inline void Subcomplex::boundary 
( Chain * output, const Chain & input ) const {
  project ( output, 
            supercomplex_ -> boundary ( include ( input_chain ) ) );
}

inline void Subcomplex::coboundary 
( Chain * output, const Chain & input ) const {
  project ( output, 
            supercomplex_ -> coboundary ( include ( input_chain ) ) );

}
  
} // namespace chomp

#endif
