// Closure.h
// Shaun Harker
// 9/21/11

#ifndef CHOMP_CLOSURE_H
#define CHOMP_CLOSURE_H

#include "chomp/Complex.h"
#include "chomp/Chain.h"

namespace chomp {
  
  template < class VectorOfSets > void 
  closure ( VectorOfSets & indices,
            Complex const& complex ) {
    int D = complex . dimension ();
    indices . resize ( D + 1 );
    for ( int d = D; d > 0; -- d ) {
      for ( Index i : indices [ d ] ) {
        Chain bd = complex . boundary ( i, d );
        for ( const Term & t : bd () ) {
          indices [ d - 1 ] . insert ( t . index () );
        }
      }
    }
  }
  
  
  template < class VectorOfSets > void 
  star ( VectorOfSets & indices,
         Complex const& complex ) {
    int D = complex . dimension ();
    indices . resize ( D + 1 );
    for ( int d = 0; d < D; ++ d ) {
      for ( Index i : indices [ d ] ) {
        Chain cbd = complex . coboundary ( i, d );
        for ( const Term & t : cbd () ) {
          indices [ d + 1 ] . insert ( t . index () );
        }
      }
    }
  }
  
} // namespace chomp

#endif
