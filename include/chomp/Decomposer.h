// Decomposer.h
// Shaun Harker
// 9/15/11

#ifndef CHOMP_DECOMPOSER_H
#define CHOMP_DECOMPOSER_H

#include "chomp/Complex.h"


namespace chomp {
  
/**************
 * Decomposer *
 **************/

class Decomposer {
protected:
  Complex * complex_;
public:
  static const int ACE = 0;
  static const int KING = 1;
  static const int QUEEN = 2;

  virtual int type ( Index x, int d ) const = 0;
  virtual Index mate ( Index x, int d ) const = 0;
  virtual bool compare ( Index lhs, Index rhs, int d ) const = 0;

  virtual void decompose ( Complex & complex ) = 0;
  Decomposer ( Complex & complex ) {
    complex_ = &complex;
  }
  virtual ~Decomposer ( void ) {}
};

} // namespace chomp

#endif
