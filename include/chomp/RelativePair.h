// RelativePair.h
// Shaun Harker
// 9/16/11

#ifndef CHOMP_RELATIVEPAIRCOMPLEX_H
#define CHOMP_RELATIVEPAIRCOMPLEX_H

#include "chomp/Complex.h"
#include "chomp/CubicalComplex.h"
#include "chomp/BitmapSubcomplex.h"

namespace chomp {

class RelativePair {
private:
  CubicalComplex * base_;
  BitmapSubcomplex * pair_;
  BitmapSubcomplex * relative_;
public:
  RelativePair ( void ) : base_ (NULL), pair_ (NULL), relative_ (NULL) {}
  RelativePair ( CubicalComplex * base_, BitmapSubcomplex * pair_, BitmapSubcomplex * relative_ )
  : base_(base_), pair_(pair_), relative_(relative_) {}
  
  void initialize ( CubicalComplex * b, BitmapSubcomplex * p, BitmapSubcomplex * r ) {
    base_ = b;
    pair_ = p;
    relative_ = r;
  }
  ~RelativePair ( void ) { 
    if ( base_ != NULL ) delete base_;
    if ( pair_ != NULL ) delete pair_;
    if ( relative_ != NULL ) delete relative_;
  }
  CubicalComplex & base ( void ) { return *base_; }
  BitmapSubcomplex & pair  ( void ) { return *pair_; }
  BitmapSubcomplex & relative ( void ) { return *relative_; }

  const CubicalComplex & base ( void ) const { return *base_; }
  const BitmapSubcomplex & pair  ( void ) const { return *pair_; }
  const BitmapSubcomplex & relative ( void ) const { return *relative_; }

};
} // namespace chomp


#endif
