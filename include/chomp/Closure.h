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
( std::vector < boost::unordered_set < Index > > & cells,
  const Complex & complex ) {
  int D = complex . dimension ();
  for ( int d = D; d > 0; -- d ) {
    BOOST_FOREACH ( Index i, cells [ d ] ) {
      Chain bd = complex . boundary ( i, d );
      BOOST_FOREACH ( const Term & t, bd () ) {
        cells [ d - 1 ] . insert ( t . index () );
      }
    }
  }
}
  
} // namespace chomp

#endif
