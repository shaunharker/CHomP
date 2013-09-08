// MatrixComplex.h
// Shaun Harker
// 3/21/2013

#ifndef CHOMP_MATRIXCOMPLEX_H
#define CHOMP_MATRIXCOMPLEX_H

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
 *      MATRIX COMPLEX      *
 ********************************/
class MatrixComplex : public Complex {
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
   * SPECIFIC TO MATRIX COMPLEX *
   **********************************/
	
private:
  /* Matrix Complex */
	
  // boundary data
  std::vector < boost::unordered_map<Index, Chain > > boundaries_;

	// coboundary data
  std::vector < boost::unordered_map<Index, Chain > > coboundaries_;
	
	public:
	/// loadFromFile.
	/// dimension of boundary map, followed by triples
  /// example: d_2 takes C_3 to C_2
	void loadFromFile ( const char * FileName );
	
};

/*******************************
 *        DEFINITIONS          *
 *******************************/
inline void MatrixComplex::loadFromFile ( const char * FileName) {
	//char *ptr;
	std::ifstream input_file ( FileName ); 
	if ( not input_file . good () ) {
		std::cout << "MatrixComplex::loadFromFile. Fatal Error. " 
		<< FileName << " not found.\n";
    exit ( 1 );
  } /* if */
	
  int current_dimension = 0;
  boundaries_.resize(current_dimension + 2);
  coboundaries_.resize(current_dimension + 2);
  
	while ( not input_file . eof () ) {
		std::string line;
		getline( input_file, line );
		std::istringstream is( line );
		if ( line . length () == 0 ) continue;
    std::vector<int> linedata;
		while ( is . good () ) {
			int v;
			is >> v;
      linedata . push_back ( v );
		}
    if ( linedata . size () == 0 ) continue;
    if ( linedata . size () == 1 ) {
      current_dimension = linedata[0];
      //std::cout << "Dimension " << current_dimension << "\n";
      if ( current_dimension >= boundaries_ . size () ) {
        boundaries_.resize(current_dimension + 2);
        coboundaries_.resize(current_dimension + 2);
      }
    }
    if ( linedata . size () == 3 ) {
      //std::cout << "linedata.size() == " << linedata.size() << "\n";
      insertCell ( linedata[0], current_dimension );
      insertCell ( linedata[1], current_dimension + 1 );
      Index bdindex = cellToIndex(linedata[0],current_dimension);
      Index cbdindex = cellToIndex(linedata[1],current_dimension+1);
      
      boundaries_[current_dimension+1][cbdindex] += Term(bdindex, linedata[2]);
      coboundaries_[current_dimension][bdindex] += Term(cbdindex, linedata[2]);
      boundaries_[current_dimension+1][cbdindex].dimension () = current_dimension;
      coboundaries_[current_dimension][bdindex].dimension () = current_dimension + 1;
      
      
    }
    
	}
	input_file . close ();
	
	return;
}


inline void MatrixComplex::boundary ( Chain * output, const Index input, int dim ) const {
  if ( boundaries_[dim].find(input) == boundaries_[dim].end() ) return;
  *output = boundaries_[dim].find(input)->second;
} /* MatrixComplex::boundary */

inline void MatrixComplex::coboundary ( Chain * output, const Index input, int dim ) const {
  if ( coboundaries_[dim].find(input) == coboundaries_[dim].end() ) return;
	*output = coboundaries_[dim].find(input)->second;
} /* MatrixComplex::coboundary */
  
} // namespace chomp


#endif

