// Closure.h
// Shaun Harker
// 9/21/11

#ifndef CHOMP_CLOSURE_H
#define CHOMP_CLOSURE_H

#include "chomp/Complex.h"
#include "chomp/Chain.h"
#include <vector>
#include "boost/unordered_set.hpp"

namespace chomp {
  
  inline void closure
  ( std::vector < boost::unordered_set < Index > > & indices,
   const Complex & complex ) {
    int D = complex . dimension ();
    indices . resize ( D + 1 );
    for ( int d = D; d > 0; -- d ) {
      BOOST_FOREACH ( Index i, indices [ d ] ) {
        Chain bd = complex . boundary ( i, d );
        BOOST_FOREACH ( const Term & t, bd () ) {
          indices [ d - 1 ] . insert ( t . index () );
        }
      }
    }
  }
  
  inline void star
  ( std::vector < boost::unordered_set < Index > > & indices,
   const Complex & complex ) {
    int D = complex . dimension ();
    indices . resize ( D + 1 );
    for ( int d = 0; d < D; ++ d ) {
      BOOST_FOREACH ( Index i, indices [ d ] ) {
        Chain cbd = complex . coboundary ( i, d );
        BOOST_FOREACH ( const Term & t, cbd () ) {
          indices [ d + 1 ] . insert ( t . index () );
        }
      }
    }
  }
  
} // namespace chomp

#endif
