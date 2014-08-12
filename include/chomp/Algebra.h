/// Algebra.h  6-19-14  Shaun Harker
#ifndef CHOMP_ALGEBRA_H
#define CHOMP_ALGEBRA_H

namespace chomp {

template < class R >
R div ( const R y, const R x ) {
  /* For some reason, chip-makers decided to make (-2^63) / -1 generate an exception
   rather than returning 2^63 (which overflows into -2^63 ).
   This isn't comprehensible, since all other operations support overflow in the expected
   manner of treating the register as Z_{2^64} . But it means we have to check for it. 
   It is especially disconcerting, since multiplication does not have the same problem.
   However there is a simple workaround. */
  if ( x != R ( -1 ) ) {
    ////NSLog ( @" div 1\n" );
    return y / x;
  } else {
    ////NSLog ( @" div 2\n" );
    return y * x;
  }
}

template < class R >
bool divisable ( const R y, const R x ) {
  if ( x == R ( 0 ) ) return ( y == R ( 0 ) );
  ////NSLog ( @" divisable 1\n" );
  if ( (div(y, x) * x) == y ) return true;
  ////NSLog ( @" divisable 2\n" );
  return false;
}

/// Bezout
/// Given a and b, calculate s, t, and gcd such that
/// s * a + t * b = gcd
/// and gcd is the greatest common divisor of a and b
/// Also, if a = +/- gcd, then choose t = 0.
template < class R >
void Bezout (R * s_out,
             R * t_out,
             R * gcd_out,
             const R & a, 
             const R & b) {
  //std::cout << "Bezout\n";
  static R one ( 1 );
  static R zero ( 0 );
  static R neg_one ( -1 );
  // Relies upon Euclidean Division being available, and also comparison.
  // The user doesn't need to sort the input.
  bool reversed = false;
  if ( a < b ) reversed = true;
  // For the proof to follow, we assume without loss a >= b

        //std::cout << " a = " << a << " and b = " << b << "\n";

  R x, y;
  if ( reversed ) {
    x = a;
    y = b;
  } else {
    x = b;
    y = a;
  }

  //std::cout << " y = " << y << " and x = " << x << "\n";

                           // For extended euclidean algorithm
  R s0 = one; R s1 = zero;
  R t0 = zero; R t1 = one;
  while ( x != zero ) {
    //std::cout << "top of bezout\n";
    //std::cout << " y = " << y << " and x = " << x << "\n";
    R q = div ( y, x );
    //std::cout << " q = y/x = " << q << "\n";
    R r = y - x * q;
    R s = s0 - q * s1;
    R t = t0 - q * t1;
    //std::cout << "Claim: " << r << " = " << s << " * " << a << " + " << t << " * " << b << "\n";
    /* r = s * a + t * b
     Proof:
     By induction. First, k = 0 case.
     Since s_0 = 1 and t_0 = -q_0, we get
     r_0 = a - q_0 * b, which is indeed the case.
     Now the k = 1 case.
     s_1 = 0 - q_1 * (1)
     t_1 = 1 - q_1 * (-q_0)
     s_1 * a + t_1 * b = 
     -q_1 * a + b + q_0 * q_1 * b =
     b - q_1 ( a - q_0 * b ) =
     b - q_1 * r_0 = r_1.
     Now the inductive step for k > 1. 
     s_k * a + t_k * b =
     (s_{k-2} - q_k * s_{k-1}) * a + (t_{k-2} - q_k * t_{k-1}) * b =
     (s_{k-2} * a + t_{k-2} * b ) - q_k * (s_{k-1} * a + t_{k-1} * b ) =
     r_{k-2} - q_k * r_{k-1} = r_k.
     */
    s0 = s1; s1 = s;
    t0 = t1; t1 = t;
    y = x;
    x = r;
  } /* while */
  // Set output
  // The Bezout coefficients s and t are the second to last ones
  // to be calculated (the last ones give 0 = s*a + t*b)
  if ( not reversed ) {
    * s_out = s0;
    * t_out = t0;
  } else {
    * t_out = s0;
    * s_out = t0;
  }
  * gcd_out = y;
  // TODO: generalize this to all unit multiples, somehow, for general Rs
  // TODO: should this really be here? perhaps the caller should worry about this
  if ( * gcd_out == a ) {
    * s_out = one; * t_out = zero;
  }
  if ( * gcd_out == -a ) {
    * s_out = neg_one; * t_out = zero;
  }
}

/// gcd 
template < class R >
R gcd ( const R & a, 
        const R & b) {
  R s, t, gcd;
  Bezout ( &s, &t, &gcd, a, b );
  return gcd;
}

#ifdef USE_GMP
#include <gmpxx.h>
template < >
void Bezout < mpz_class > ( mpz_class * s_out,
                           mpz_class * t_out,
                           mpz_class * gcd_out,
                           const mpz_class & a, 
                           const mpz_class & b) {
  mpz_gcdext ( gcd_out -> get_mpz_t (), 
              s_out -> get_mpz_t (), 
              t_out -> get_mpz_t (), 
              a . get_mpz_t (), 
              b . get_mpz_t () );
  // TODO: should this really be here? perhaps the caller should worry about this
  static mpz_class one ( 1 );
  static mpz_class zero ( 0 );
  static mpz_class neg_one ( -1 );
  if ( * gcd_out == a ) {
    * s_out = one; * t_out = zero;
  }
  if ( * gcd_out == -a ) {
    * s_out = neg_one; * t_out = zero;
  }
}
#endif


}

#endif
