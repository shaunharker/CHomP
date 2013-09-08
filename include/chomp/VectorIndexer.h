// VectorIndexer.h
// Shaun Harker
// 9/15/11

#ifndef CHOMP_VECTORINDEXER_H
#define CHOMP_VECTORINDEXER_H

#include <vector>

namespace chomp {

// both K and M are unsigned integer types here
template < class K, class M >
class VectorIndexer {
public:
  typedef K Key;
  typedef M Mapped;
  void clear ( void );
  void initialize ( const std::vector < Key > & data );
  uint32_t size ( void ) const;
  const Key & key ( Mapped i ) const;
  Mapped & rank ( const Key & k );
  const Mapped & rank ( const Key & k ) const;
  void reindex ( const std::vector < Mapped > & permute );
private:
  std::vector < Key > keys_;
  std::vector < Mapped > mapped_;
};

template < class K, class M > void
VectorIndexer<K,M>::clear ( void ) {
  mapped_ . clear ();
  keys_ . clear ();
}
template < class K, class M > void
VectorIndexer<K,M>::initialize ( const std::vector < Key > & keys ) {
  clear ();
  keys_ = keys;
  for ( Mapped i = 0; i < keys . size (); ++ i ) {
    if ( mapped_ . size () <= keys [ i ] ) 
      mapped_ . resize ( keys[i] + 1);
    mapped_ [ keys [ i ] ] = i;
  }
}

template < class K, class M > uint32_t
VectorIndexer<K,M>::size ( void ) const {
  return keys_ . size ();
}

template < class K, class M > const K &
VectorIndexer<K,M>::key ( Mapped i ) const {
  return keys_ [ i ];
}

template < class K, class M > M &
VectorIndexer<K,M>::rank ( const Key & k ) {
  return mapped_ [ k ];
}

template < class K, class M > const M &
VectorIndexer<K,M>::rank ( const Key & k ) const {
  return mapped_ [ k ];
}

template < class K, class M > void 
VectorIndexer<K,M>::reindex ( const std::vector < Mapped > & permute ) {
  std::vector < Key > new_keys ( size () );
  for ( int i = 0; i < size (); ++ i ) {
    new_keys [ permute [ i ] ] = keys_ [ i ];
  }
  initialize ( new_keys );
}
} // namespace chomp
  
#endif
