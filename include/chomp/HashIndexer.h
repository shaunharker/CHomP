// HashIndexer.h
// Shaun Harker
// 9/13/11

#ifndef CHOMP_HASHINDEXER_H
#define CHOMP_HASHINDEXER_H

#include <vector>
#include "boost/unordered_map.hpp"

namespace chomp {
  
template < class K, class M >
class HashIndexer {
public:
  typedef K Key;
  typedef M Mapped;
  void clear ( void );
  void initialize ( const std::vector < Key > & data );
  uint32_t size ( void ) const;
  const Key & key ( Mapped i ) const;
  //Mapped & rank ( const Key & k );
  const Mapped rank ( const Key & k ) const;
  void reindex ( const std::vector < Mapped > & permute );
private:
  boost::unordered_map < Key, Mapped > data_;
  std::vector < Key > keys_;
};

template < class K, class M > void
HashIndexer<K,M>::clear ( void ) {
  data_ . clear ();
  keys_ . clear ();
}
template < class K, class M > void
HashIndexer<K,M>::initialize ( const std::vector < Key > & keys ) {
  clear ();
  keys_ = keys;
  for ( Mapped i = 0; i < keys . size (); ++ i ) {
    data_ [ keys [ i ] ] = i;
    //std::cout << "ASSIGNED data_ [ " << keys [ i ] << " ] = " << i << "\n";
  }
}

template < class K, class M > uint32_t
HashIndexer<K,M>::size ( void ) const {
  return data_ . size ();
}

template < class K, class M > const K &
HashIndexer<K,M>::key ( Mapped i ) const {
  return keys_ [ i ];
}

  //template < class K, class M > M &
  //HashIndexer<K,M>::rank ( const Key & k ) {
  //return data_ [ k ];
  //}

template < class K, class M > const M
HashIndexer<K,M>::rank ( const Key & k ) const {
  typename boost::unordered_map < Key, Mapped >::const_iterator it = data_ . find ( k );
  if ( it == data_ . end () ) {
    std::cout << "keys_ . size () = " << keys_ . size () << "\n";
    return keys_ . size ();
  }
  if ( it -> second  >= keys_ . size () ) std::cout << "PROBLEM in hashindexer rank\n";
  return it -> second;
}

template < class K, class M > void 
HashIndexer<K,M>::reindex ( const std::vector < Mapped > & permute ) {
  std::vector < Key > new_keys ( size () );
  for ( uint32_t i = 0; i < size (); ++ i ) {
    new_keys [ permute [ i ] ] = keys_ [ i ];
  }
  initialize ( new_keys );
}

} // namespace chomp

#endif
