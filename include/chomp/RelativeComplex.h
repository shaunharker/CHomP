// RelativeComplex.h

// Shaun Harker
// 7/08/2013

#ifndef CHOMP_RELATIVECOMPLEX_H
#define CHOMP_RELATIVECOMPLEX_H

#include <iostream>
#include <cstdlib>
#include <stdint.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
#include <queue>


#include "boost/unordered_set.hpp"
#include "boost/unordered_map.hpp"

#include "chomp/Complex.h"
#include "chomp/Chain.h"
#include "chomp/Rect.h"

namespace chomp {


/********************************
 *      RELATIVE COMPLEX        *
 ********************************/
class RelativeComplex : public Complex {
public:
  // Cell type: uint64_t
  CHOMP_COMPLEX(uint64_t)
	
  /*******************************
   *      COMPLEX INTERFACE      *
   *******************************/
  /// boundary. See Complex.h 
  virtual void boundary ( Chain * output, const Index input, int dim ) const;
  /// coboundary. See Complex.h
  virtual void coboundary ( Chain * output, const Index input, int dim ) const;

  /**********************************
   * SPECIFIC TO RELATIVE COMPLEX   *
   **********************************/
	
private:
  /* Relative Complex */
  Complex * full;

	public:
	
	RelativeComplex ( Complex * X, const std::vector < boost::unordered_set <uint64_t> > & A );

	void include ( Chain * output, const Chain & input ) const;
	void project ( Chain * output, const Chain & input ) const;

	Chain include ( const Chain & input ) const { Chain output; include ( &output, input ); return output; }
	Chain project ( const Chain & input ) const { Chain output; project ( &output, input ); return output; }

};

/*******************************
 *        DEFINITIONS          *
 *******************************/

inline RelativeComplex::RelativeComplex ( Complex * X, 
 	                                        const std::vector < boost::unordered_set <uint64_t> > & A ) : full ( X )  {
  for ( int d = 0; d <= X -> dimension (); ++ d ) {
    bool insert_all = false;
    if ( d >= A . size () ) insert_all = true;
 		for ( uint64_t i = 0; i < X -> size ( d ); ++ i ) {
 			if ( insert_all || (A [ d ] . count ( i ) == 0) ) insertCell ( i, d );
 		}
 	}
}

inline void RelativeComplex::boundary ( Chain * output, const Index input, int dim ) const {
	Chain input_chain;
	input_chain . dimension () = dim;
	input_chain += Term ( input, Ring ( 1 ) );
	Chain included;
	include ( &included, input_chain);
	Chain bd = full -> boundary ( included );
	project ( output, bd );
} /* RelativeComplex::boundary */

inline void RelativeComplex::coboundary ( Chain * output, const Index input, int dim ) const {
	Chain input_chain;
	input_chain . dimension () = dim;
	input_chain += Term ( input, Ring ( 1 ) );
	Chain included;
	include ( &included, input_chain);
	Chain cbd = full -> coboundary ( included );
	project ( output, cbd );
} /* RelativeComplex::coboundary */
  
inline void RelativeComplex::include (Chain * output, 
                                      const Chain & input ) const {
    int D = output -> dimension () = input . dimension ();
    BOOST_FOREACH ( const Term & t, input () ) {
    	* output += Term ( indexToCell ( t . index (), D ), t . coef () ); 
    }
    *output = simplify ( *output ); 
  }

inline void RelativeComplex::project (Chain * output, 
                                      const Chain & input ) const {
    int D = output -> dimension () = input . dimension ();
    BOOST_FOREACH ( const Term & t, input () ) {
    	Term s ( cellToIndex ( t . index (), D ), t . coef () );
    	if ( s == size ( D ) ) continue;
    	* output += s; 
    }
    *output = simplify ( *output ); 
  }
} // namespace chomp

#endif
