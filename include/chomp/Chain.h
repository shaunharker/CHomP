//
//  Chain.h
//  Shaun Harker 
//  9/10/11.
//

#ifndef CHOMP_CHAIN_H
#define CHOMP_CHAIN_H

#include <iostream>
#include <vector>
#include <cstdlib>

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"

#include "chomp/Ring.h"

namespace chomp {

/*********************
 *      Index        *
 *********************/
typedef uint64_t Index;

/*********************
 *      Term         *
 *********************/
class Term {
  Index index_;
  Ring coef_;
public:
  Index & index ( void ) { return index_; }
  const Index & index ( void ) const { return index_; }
  Ring & coef ( void ) { return coef_; }
  const Ring & coef ( void ) const { return coef_; }
  Term ( void ) {}
  Term ( const Index & index_, const Ring & coef_) 
  : index_(index_), coef_(coef_) {}
  Term ( const std::pair < const Index, Ring > & pair ) : index_(pair.first), coef_(pair.second) {}
};

inline bool operator==(const Term &x, const Term &y)
{
  return (x.index() == y.index()) && (x.coef() == y.coef());
}

inline bool operator!=(const Term &x, const Term &y)
{
  return !(x == y);
}

inline std::ostream & operator << ( std::ostream & outstream, const Term & print_me ) {
  outstream << print_me . coef () << "[" << print_me . index () << "]";
  return outstream;
}


/*********************
 *      Chain        *
 *********************/
template < class C >
class CustomChain {
private:
  typedef C Container;
  typedef uint64_t size_type;
  int dimension_;
  Container container_; 
public:
  int & dimension ( void ) { return dimension_; }
  const int & dimension ( void ) const { return dimension_; }
  CustomChain ( void ) {}
  //CustomChain ( size_type d ) : dimension_ ( d ) {}  SERIOUS GOTCHA! Was hijacking mistakes.
  Container & operator () ( void ) { return container_; }
  const Container & operator () ( void ) const { return container_; }
};



template < class C >
inline CustomChain<C> & operator += ( CustomChain<C> & lhs, 
                                     const Term & rhs ) {
  std::cout << "Sorry, this chain needs to have += implemented.\n";
  //std::cout << "Specialize += for CustomChain<C>, with C =" << typeid ( C ) << "\n";
  exit ( 1 );
  return lhs;
}

template < class C >
inline std::ostream & operator << ( std::ostream & outstream, const CustomChain<C> & print_me ) {
  bool first = true;
  BOOST_FOREACH ( const Term & t, print_me () ) {
    if ( not first ) {
      outstream << " + ";
    } else {
      first = false;
    }
    outstream << t;
  }
  if ( first ) outstream << "0";
  return outstream;
}

/*******************************
 *      Vector-Based Chain     *
 *******************************/

typedef CustomChain < std::vector < Term > > Chain;

template <>
inline CustomChain< std::vector < Term > > & 
operator += ( CustomChain< std::vector < Term > > & lhs, 
             const Term & rhs ) {
  lhs () . push_back ( rhs );
  return lhs;
}

/*******************************
 *      Chain Arithmetic       *
 *******************************/
template < class C >
inline CustomChain<C> & operator -= ( CustomChain<C> & lhs, 
                                      Term rhs ) {
  rhs . coef () = - rhs . coef ();
  return (lhs += rhs);
}

template < class L, class R >
inline CustomChain<L> & operator += ( CustomChain<L> & lhs, 
                                      const CustomChain<R> & rhs ) {
  BOOST_FOREACH ( const Term & term, rhs () ) {
    lhs += term;
  }
  return lhs;
}

template < class L, class R >
inline CustomChain<L> & operator -= ( CustomChain<L> & lhs, 
                                      const CustomChain<R> & rhs ) {
  BOOST_FOREACH ( Term term, rhs () ) {
    term . coef () = - term . coef ();
    lhs += term;
  }
  return lhs;
}

template < class C >
inline CustomChain<C> & operator *= ( CustomChain<C> & lhs, 
                                      const Ring & rhs ) {
  BOOST_FOREACH ( Term & term, lhs () ) {
    term . coef () *= rhs;
  }
  return lhs;
}


template < class C >
inline CustomChain<C> operator * ( const CustomChain<C> & lhs, 
                                   const Ring & rhs ) {
  CustomChain<C> result;
  result . dimension () = lhs . dimension ();
  BOOST_FOREACH ( Term term, lhs ()) {
    term . coef () *= rhs;
    result += term;
  }
  return result;
}

template < class C >
inline CustomChain<C> operator * ( const Ring & lhs, 
                                   const CustomChain<C> & rhs ) {
  CustomChain<C> result;
  result . dimension () = rhs . dimension ();
  BOOST_FOREACH ( Term term, rhs ()) {
    term . coef () *= lhs;
    result += term;
  }
  return result;
}

// These ones really don't seem like they should be used,
// except perhaps as convenience functions.
inline Chain operator + ( const Chain & lhs, const Chain & rhs ) {
  Chain result;
  result . dimension () = lhs . dimension ();
  BOOST_FOREACH ( const Term & term, lhs () ) {
    result += term;
  }
  BOOST_FOREACH ( const Term & term, rhs () ) {
    result += term;
  }
  return result;
}

inline Chain operator - ( const Chain & lhs, const Chain & rhs ) {
  Chain result;
  result . dimension () = lhs . dimension ();
  BOOST_FOREACH ( const Term & term, lhs ()) {
    result += term;
  }
  BOOST_FOREACH ( Term term, rhs ()) {
    term . coef () = - term . coef ();
    result += term;
  }
  return result;
}

inline Chain simplify ( const Chain & chain ) {
  boost::unordered_map < Index, Ring > map_chain;
  BOOST_FOREACH ( Term term, chain ()) {
    if ( map_chain . count ( term . index () ) ) 
      map_chain [ term . index () ] += term . coef ();
    else
      map_chain [ term . index () ] = term . coef ();
  }
  Chain result;
  typedef std::pair<Index, Ring> mapTerm;
  BOOST_FOREACH ( mapTerm term, map_chain ) {
    if ( term . second != Ring ( 0 ) )
      result += Term ( term . first, term . second );
  }
  result . dimension () = chain . dimension ();
  return result;
}
  
} // namespace chomp
#endif
