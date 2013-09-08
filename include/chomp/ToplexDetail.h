/// ToplexDetail.h
/// Shaun Harker
/// 9/21/11

#ifndef CHOMP_TOPLEXDETAIL_H
#define CHOMP_TOPLEXDETAIL_H

#include "boost/serialization/serialization.hpp"

namespace chomp {

typedef uint32_t GridElement;

/***************
 * struct Node *
 ***************/
struct Node {
  Node * left_;
  Node * right_;
  Node * parent_;
  GridElement contents_;
  int dimension_;
  Node ( void );
  ~Node ( void );
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & left_;
    ar & right_;
    ar & parent_;
    ar & contents_;
    ar & dimension_;
  }

};

inline Node::Node ( void ) : left_ ( NULL ), right_ ( NULL ), parent_ ( NULL ), contents_ ( 0 ), dimension_ ( 0 ) {
} /* Adaptive_Cubical::Node::Node */

inline Node::~Node ( void ) {
  if ( left_ != NULL ) delete left_;
  if ( right_ != NULL ) delete right_;
} /* Adaptive_Cubical::Node::~Node */


/********************************
 * class Toplex_const_iterator  *
 ********************************/
class Toplex_const_iterator {
public:
  /* Iterator typedefs */
  typedef std::forward_iterator_tag iterator_category;
  typedef GridElement value_type;
  typedef GridElement * pointer;
  typedef std::ptrdiff_t difference_type;
  typedef const GridElement & reference;
  
  Toplex_const_iterator ( void );
  Toplex_const_iterator ( Node * node );
  Toplex_const_iterator & operator ++ ( void );
  bool operator != ( const Toplex_const_iterator & right_hand_side ) const;
  bool operator == ( const Toplex_const_iterator & right_hand_side ) const;
  bool operator < ( const Toplex_const_iterator & right_hand_side ) const;
  const GridElement & operator * ( void ) const; 
  Node * node ( void );
private:
  friend class Toplex;
  Node * node_;
  friend class boost::serialization::access;
public:
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & node_;
  }

};

inline Toplex_const_iterator::Toplex_const_iterator ( void ) {
} /* Adaptive_Cubical::Toplex_const_iterator::Toplex_const_iterator */

inline Toplex_const_iterator::Toplex_const_iterator ( Node * node ) : node_ ( node ) {
} /* Adaptive_Cubical::Toplex_const_iterator::Toplex_const_iterator */

inline Toplex_const_iterator & Toplex_const_iterator::operator ++ ( void ) {
  /* Go up until either root is reached, or node is a left child of a parent possessing a right child as well */
  while (  node_ -> parent_ != NULL && ( node_ -> parent_ -> left_ != node_ || node_ -> parent_ -> right_ == NULL ) ) node_ = node_ -> parent_;
  /* If root has been reached, return NULL */
  if ( node_ -> parent_ == NULL ) {
    node_ = NULL;
    return *this;
  } /* if */
  /* Move to right child of parent */
  node_ = node_ -> parent_ -> right_;
  /* Descend to a leaf, taking left route whenever possible */
  while ( 1 ) {
    if ( node_ -> left_ != NULL ) {
      node_ = node_ -> left_;
    } else if ( node_ -> right_ != NULL ) {
      node_ = node_ -> right_;
    } else break;
  } /* while */
  return *this;
} /* Adaptive_Cubical::Toplex_const_iterator::operator ++ */

inline bool Toplex_const_iterator::operator != ( const Toplex_const_iterator & right_hand_side ) const {
  return node_ != right_hand_side . node_;
} /* Adaptive_Cubical::Toplex_const_iterator::operator != */

inline bool Toplex_const_iterator::operator == ( const Toplex_const_iterator & right_hand_side ) const {
  return node_ == right_hand_side . node_;
} /* Adaptive_Cubical::Toplex_const_iterator::operator == */

inline bool Toplex_const_iterator::operator < ( const Toplex_const_iterator & right_hand_side ) const {
  return node_ < right_hand_side . node_;
} /* Adaptive_Cubical::Toplex_const_iterator::operator < */

inline const GridElement & Toplex_const_iterator::operator * ( void ) const {
  return node_ -> contents_;
} /* Adaptive_Cubical::Toplex_const_iterator::operator * */

inline Node * Toplex_const_iterator::node ( void ) {
  return node_;
} /* Adaptive_Cubical::Toplex_const_iterator::node */

} // namespace chomp



#endif
