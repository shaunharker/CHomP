// BitmapSubcomplex.h
// Shaun Harker
// 9/16/11

#ifndef CHOMP_BITMAPSUBCOMPLEX_H 
#define CHOMP_BITMAPSUBCOMPLEX_H

#include "chomp/Complex.h"

#include "boost/foreach.hpp"

namespace chomp {
  
class BitmapSubcomplex : public Complex {
public:
  CHOMP_COMPLEX(Index)
  virtual void boundary ( Chain * output, const Index input, int dim ) const;
  virtual void coboundary ( Chain * output, const Index input, int dim ) const;
  const Complex & full ( void ) const;
  void include ( Chain * output, const Chain & input ) const;
  Chain include ( const Chain & input ) const;  
  void project ( Chain * output, const Chain & input ) const;  
  Chain project ( const Chain & input ) const;  
  void initialize ( Complex & base, bool full = true );
  void erase ( const Cell cell, int d );
  void insert ( const Cell cell, int d );
  void finalize ( void );
  BitmapSubcomplex ( Complex & base );
  BitmapSubcomplex ( Complex & b, bool flag );
  const Complex & base ( void ) const;
private:
  Complex * base_;
  std::vector<std::vector<bool> > pattern_;
};

/***************
 * DEFINITIONS *
 ***************/

inline void BitmapSubcomplex::boundary (Chain * output, 
                                       const Index input, 
                                       int dim ) const {
  Cell c = indexToCell ( input, dim );
  Chain bd = base () . boundary ( c, dim );
  project ( output, bd );
}

inline void BitmapSubcomplex::coboundary (Chain * output, 
                                         const Index input, 
                                         int dim ) const {
  Cell c = indexToCell ( input, dim );
  Chain cbd = base () . coboundary ( c, dim );
  project ( output, cbd );
}

inline const Complex & BitmapSubcomplex::full ( void ) const {
  return base ();
}

inline void BitmapSubcomplex::include ( Chain * output, 
                                 const Chain & input ) const {
  int dim = input . dimension ();
  output -> dimension () = dim;
  BOOST_FOREACH ( Term term, input () ) {
    term . index () = indexToCell ( term . index (), dim );
    (*output) += term;
  }
}

inline Chain BitmapSubcomplex::include ( const Chain & input ) const {
  Chain output; include( &output, input); return output;
}

inline void BitmapSubcomplex::project (Chain * output, 
                                      const Chain & input ) const {
  int d = input . dimension ();
  output -> dimension () = d;
  BOOST_FOREACH ( const Term & t, input () ) {
    //std::cout << "Bitsubcomplex, t .index() " << t . index () << "\n";
    //std::cout << "                pattern_ [ " << d << "] . size () = " 
    //   << pattern_ [ d ] . size () << "\n";
    if ( pattern_ [ d ] [ t . index () ] ) {
      Index i = cellToIndex ( t . index (), d );
      (*output) += Term ( i, t . coef () );
    }
  }
}

inline Chain BitmapSubcomplex::project ( const Chain & input ) const {
  Chain output; project( &output, input); return output;
}

inline void BitmapSubcomplex::initialize ( Complex & b, bool bit ) {
  //std::cout << "BitmapSubcomplex: initializing " << this << 
  // " from " << &b << "\n";
  base_ = &b;
  dimension () = base () . dimension ();
  pattern_ . resize ( dimension () + 1 );
  for ( int d = 0; d <= dimension (); ++ d ) {
    pattern_ [ d ] . resize ( base () . size ( d ), bit );
  }
}

inline void BitmapSubcomplex::erase ( const Cell cell, int d ) {
  pattern_ [ d ] [ cell ] = false;
}

inline void BitmapSubcomplex::insert ( const Cell cell, int d ) {
  pattern_ [ d ] [ cell ] = true;
}

inline void BitmapSubcomplex::finalize ( void ) {
  //std::cout << "BitmapSubcomplex: finalizing " << this << "\n";
  //startInserting ();
  for ( int d = 0; d <= dimension (); ++ d ) {
    Index hoist = base () . size ( d );
    for ( Index i = 0; i < hoist; ++ i ) {
      if ( pattern_ [ d ] [ i ] ) insertCell ( i, d );
    }
  }
  //finishedInserting ();
}

inline BitmapSubcomplex::BitmapSubcomplex ( Complex & b ) {
  initialize ( b, true );
}

inline BitmapSubcomplex::BitmapSubcomplex ( Complex & b, bool flag ) {
  initialize ( b, flag );
}

inline const Complex & BitmapSubcomplex::base ( void ) const {
  return *base_;
}

} // namespace chomp

#endif
