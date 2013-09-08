#include <iostream>

#define RINGDEFINED
#include "chomp/Ring.h"
using namespace chomp;
typedef Long Ring;

#include "chomp/MatrixComplex.h"
#include "chomp/Generators.h"
#include "chomp/MorseComplex.h"

int main ( int argc, char * argv [] ) {
	
	MatrixComplex complex;
	typedef MatrixComplex::Cell Cell;
	if ( argc == 1 ) {
    std::cout << "Give a filename.\n";
    abort ();
	} else {
	  complex.loadFromFile ( argv [ 1 ] );
	}
  
  
  std::cout << "Betti Numbers: ";
  Generators_t gen = MorseGenerators ( complex );
  for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
    int betti = 0;
    for ( unsigned int gen_index = 0; gen_index < gen [ dim ] . size (); ++gen_index ) {
      if ( gen[dim][gen_index].second == 0 ) ++ betti;
    }
    std::cout << betti << " ";
  }
  std::cout << "\n";
  
  
  if ( argc > 2 ) {
    
    std::cout << "----------------------------------\n";
    std::cout << "  (Index,dim)     <=>    Cell     \n";
    std::cout << "----------------------------------\n";
    for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
      for ( Index i = 0; i < complex . size ( dim ); ++ i ) {
        Cell s = complex . indexToCell ( i, dim );
        std::cout << " (" << i << ", " << dim << ") <=> " << s << "\n";
      }
    }
    
    std::cout << "----------------\n";
    std::cout << "| bd and cbd   | \n";
    std::cout << "----------------\n";
    for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
      for ( Index i = 0; i < complex . size ( dim ); ++ i ) {
        std::cout << "bd(" << i << ", " << dim << ") = " << complex.boundary (i, dim) << "\n";
        std::cout << "cbd(" << i << ", " << dim << ") = " << complex.coboundary (i, dim) << "\n";
      }
    }
    
    std::cout << "----------------------\n";
    std::cout << "| complex property   | \n";
    std::cout << "----------------------\n";
    for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
      for ( Index i = 0; i < complex . size ( dim ); ++ i ) {
        std::cout << "bd(bd(" << i << ", " << dim << ")) = "
        << simplify(complex.boundary(complex.boundary (i, dim))) << "\n";
        std::cout << "cbd(cbd(" << i << ", " << dim << ")) = "
        << simplify(complex.coboundary(complex.coboundary (i, dim))) << "\n";
      }
    }
    
    
    std::cout << "------------------------------------\n";
    std::cout << "| Generators via SmithGenerators   | \n";
    std::cout << "------------------------------------\n";
    {
      Generators_t gen = SmithGenerators ( complex );
      
      for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
        std::cout << "H_" << dim << ":\n";
        for ( unsigned int gen_index = 0; gen_index < gen [ dim ] . size (); ++gen_index ) {
          std::cout << "  (" << gen[dim][gen_index].second << ") : " << gen [ dim ] [gen_index ] . first << "\n";
          std::cout << "bd(cycle) = "<< simplify ( complex . boundary ( gen [ dim ] [gen_index ] . first ) ) << "\n";
        }
      }
    }
    
    
    std::cout << "------------------------------------\n";
    std::cout << "| Generators via MorseGenerators   | \n";
    std::cout << "------------------------------------\n";
    {
      Generators_t gen = MorseGenerators ( complex );
      
      for ( int dim = 0; dim <= complex . dimension (); ++ dim ) {
        std::cout << "H_" << dim << ":\n";
        for ( unsigned int gen_index = 0; gen_index < gen [ dim ] . size (); ++gen_index ) {
          std::cout << "  (" << gen[dim][gen_index].second << ") : " << gen [ dim ] [gen_index ] . first << "\n";
          std::cout << "bd(cycle) = "<< simplify ( complex . boundary ( gen [ dim ] [gen_index ] . first ) ) << "\n";
        }
      }
    }
    
    std::cout << "---------------------\n";
    std::cout << "|   Morse complex   |\n";
    std::cout << "---------------------\n";
    
    MorseComplex morse ( complex );
    //MorseSanity ( morse );

    for ( int dim = 0; dim <= morse . dimension (); ++ dim ) {
      for ( Index i = 0; i < morse . size ( dim ); ++ i ) {
        std::cout << "bd(" << i << ", " << dim << ") = \n ";
        std::cout << morse.boundary (i, dim) << "\n";
        //std::cout << "cbd(" << i << ", " << dim << ") = " << morse.coboundary (i, dim) << "\n";
      }
    }
    for ( int dim = 0; dim <= morse . dimension (); ++ dim ) {
      std::cout << "There are " << morse . size ( dim ) << " "
      << dim <<"-dimensional critical cells.\n";
    }
  }
  return 0;
}
