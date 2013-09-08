// SimplicialComplex.h
// Amit Patel
// 12/2/11

#ifndef CHOMP_SIMPLICIALCOMPLEX_H
#define CHOMP_SIMPLICIALCOMPLEX_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <iterator>
#include <sstream>
#include <queue>


#include "boost/unordered_set.hpp"

#include "chomp/Complex.h"
#include "chomp/Chain.h"
#include "chomp/Rect.h"

namespace chomp {

/************************
 *      SIMPLEX					*
 ************************/
class Simplex {
private:
	std::vector<int> simplex_;
	
public:
	Simplex ( void ) {};
	Simplex (const std::vector<int> & simplex ) : simplex_ ( simplex ) {}
	friend std::size_t hash_value(Simplex const &);
	
	const std::vector < int > & simplex ( void ) const { return simplex_; }
	std::vector < int > & simplex ( void ) { return simplex_; }
	
	unsigned int dimension ( void ) const { return simplex_ . size () - 1; }
	int operator [] ( int i ) const { return simplex_ [ i ]; }
	//Simplex & operator = ( const Simplex & rhs );
	bool operator == ( const Simplex & rhs ) const;
	friend std::ostream & operator << ( std::ostream & outstream,
									   const Simplex & simplex );
};

std::ostream & operator << ( std::ostream & outstream, const Simplex & s ) {
	int size = s . simplex () . size ();
	outstream << "[";
	for ( int i = 0; i < size - 1; ++ i ) {
		outstream << s . simplex () [ i ] << " ";
	}
	if ( s . simplex () . size () > 0 ) outstream << s . simplex () [ size - 1 ]; 
	outstream << "]";
	return outstream;
}

/*******************************
 *        DEFINITIONS          *
*******************************/

/*
Simplex & Simplex::operator = ( const Simplex & rhs ) {
	simplex . clear ();
	simplex = rhs . simplex;
	return *this;
}
 */

bool Simplex::operator == ( const Simplex & rhs ) const {
	if ( dimension () != rhs . dimension () ) return false;
	for ( unsigned int i = 0; i < simplex_ . size (); ++ i ) {
		if ( simplex_ [ i ] != rhs . simplex_ [ i ] ) return false;
	}
	return true;
}

/*
inline Simplex::Simplex (const int vertex_list[]) :
simplex (vertex_list, vertex_list + sizeof(vertex_list) / sizeof(int)) 
{
}
*/

inline unsigned int fvn_hash( unsigned int hash_me ) {
	// Fowler Noll Vo hash function.
	unsigned int hash = 2166136261u;
	hash ^= (hash_me >> 24);
	hash *= 16777619;
	hash ^= ((hash << 8) >> 24);
	hash *= 16777619;
	hash ^= ((hash_me << 16) >> 24);
	hash *= 16777619;
	hash ^= ((hash_me << 24) >> 24);
	hash *= 16777619;
	return hash;
}

inline std::size_t hash_value(Simplex const & simplex ) {
	//boost::hash<int> hasher;
	
	std::size_t seed = 0;
	int d = simplex . dimension ();
	for ( int i = 0; i <= d; ++ i ) {
		seed += fvn_hash ( simplex [ i ] );
		//boost::hash_combine ( seed, simplex [ i ] );
		//std::cout << "i = " << i << ", seed = " << seed << ", simplex [" << i << "] = " << simplex [ i ] << "\n";
	}
	return seed;
}


/********************************
 *      SIMPLICIAL COMPLEX      *
 ********************************/
class SimplicialComplex : public Complex {
public:
  // Cell type: Simplex
  CHOMP_COMPLEX(Simplex)
	
  /*******************************
   *      COMPLEX INTERFACE      *
   *******************************/
  /// boundary. See Complex.h 
  virtual void boundary ( Chain * output, const Index input, int dim ) const;
  /// coboundary. See Complex.h
  virtual void coboundary ( Chain * output, const Index input, int dim ) const;

  /**********************************
   * SPECIFIC TO SIMPLICIAL COMPLEX *
   **********************************/
	
private:
  /* Simplicial Complex */
	boost::unordered_set<Simplex> processed_;
	
	// coboundary data
	std::vector < boost::unordered_map<Index, Chain > > coboundaries_;
	
	public:
	/// loadFromFile.
	/// the file is a list of simplicies where an n-simplex is a list of integers 
	/// "int_1 int_2 ...  int_{n+1} \n" each representing a vertex of the complex,
	/// Assumption: The integers describing each simplex is listed in increasing value
	/// The function reads the file and adds each simplex as well as all its faces
	/// No one line can have more than 512 characters
	void loadFromFile ( const char * FileName );
	
	void generateCoboundaryData ( void );
	
	

	
};

/*******************************
 *        DEFINITIONS          *
 *******************************/
inline void SimplicialComplex::loadFromFile ( const char * FileName) {
	//char *ptr;
	std::ifstream input_file ( FileName ); 
	if ( not input_file . good () ) {
		std::cout << "SimplicialComplex::loadFromFile. Fatal Error. " 
		<< FileName << " not found.\n";
    exit ( 1 );
  } /* if */
	//int index = 0;
	
	//Simplex s;
	std::vector< std::vector<int> > max_simplices;
	while ( not input_file . eof () ) {
		std::string line;
		getline( input_file, line );
		std::istringstream is( line );
		if ( line . length () == 0 ) continue;
		max_simplices . push_back ( std::vector < int > () );
		while ( is . good () ) {
			int v;
			is >> v;
			max_simplices . back () . push_back ( v );
		}
	}
	input_file . close ();
	
  //DEBUG
	// for ( unsigned int i = 0; i < max_simplices . size (); ++ i ) {
	//	for ( unsigned int vi = 0; vi < max_simplices [ i ] . size (); ++ vi ) {
	//		std::cout << max_simplices [ i ] [ vi ] << " ";
	//	}
	//	std::cout << "\n";
	//}
	
	boost::unordered_set < Simplex > processed_;
	std::queue < Simplex > work_q;
	BOOST_FOREACH ( const std::vector < int > & maximal, max_simplices ) {
		processed_.insert(maximal);
		if (maximal.size() == 1) continue;
		work_q.push(maximal);
		//std::cout << maximal << "\n";
		while ( not work_q.empty() ) {
			//std::cout << "Front of the queue : " << work_q.front() << "\n";
			for ( unsigned int i = 0; i <= work_q.front().dimension(); i ++ ) {
				Simplex work_simplex;
				for (unsigned int j = 0; j <= work_q.front().dimension(); j ++ ) {
					if (j != i)
						work_simplex.simplex().push_back(work_q.front()[j]);
				}
				//std::cout << "   Inserting simplex : " << work_simplex << "\n";
				processed_.insert(work_simplex);
				if (work_simplex.dimension() > 0) work_q.push(work_simplex);
			}
			//std::cout << "Popping simplex : " << work_q.front() << "\n";
			work_q.pop();
		}
	}
	
	/*
	BOOST_FOREACH( const Simplex &s, processed_ ) {
		std::cout << s << "\n";
	}
	*/
	
	//startInserting ();
	BOOST_FOREACH( const Simplex &s, processed_ ) {
		//std::cout << "Inserting cell " << s << "\n";
		insertCell (s, s.dimension() );
	}
	//finishedInserting ();
	generateCoboundaryData ();
	return;
}


inline void SimplicialComplex::boundary ( Chain * output, const Index input, int dim ) const {
	if (dim == 0) {
		output -> dimension () = -1;
		return;
	}
	Ring positive = Ring (1); 
	Ring negative = - positive;
	
	bool sign = true;

	Simplex s = indexToCell(input, dim);
	//std::cout << "bd(" << input << ", " << dim << ") = ";
	for ( unsigned int i = 0; i <= s.dimension(); i++ ) {
		Simplex work_simplex;
		for ( unsigned int j = 0; j <= s.dimension(); j++ ) {
			if ( j != i) work_simplex.simplex().push_back( s[j] );
		}
		//std::cout << "Simplex " << work_simplex << "has index " << cellToIndex(work_simplex, dim-1) << "\n";
		(*output) += Term ( cellToIndex (work_simplex, dim-1 ), sign ? positive : negative );
		sign = not sign;
	}
	output -> dimension () = dim - 1;
	//std::cout << *output << "\n";
} /* SimplicialComplex::boundary */

inline void SimplicialComplex::coboundary ( Chain * output, const Index input, int dim ) const {
	//std::cout << "cbd(" << input << ", " << dim << ") = ";

	*output = coboundaries_[dim].find(input)->second;
	//std::cout << *output << "\n";

} /* SimplicialComplex::coboundary */

inline void SimplicialComplex::generateCoboundaryData ( void ) {
// should be called by finalize()
	coboundaries_ . resize ( dimension () + 1 );
	// Insert empty chains.
	for ( int dim = 0; dim <= dimension (); ++ dim ) {
		for ( Index i = 0; i < size ( dim ); ++ i ) {
			coboundaries_ [ dim ] . insert ( std::make_pair ( i, Chain () ) );
			coboundaries_ [ dim ] [ i ] . dimension () = dim + 1;
		}
	}
	// Fill in coboundary data
	for ( int dim = 0; dim <= dimension (); ++ dim ) {
		for ( Index i = 0; i < size ( dim ); ++ i ) {
			Chain bd = boundary ( i, dim );
			BOOST_FOREACH ( const Term & t, bd () ) {
				coboundaries_ [ dim - 1 ] [ t . index () ] += Term ( i,  t . coef () );
			}
		}
	}
}
  
} // namespace chomp


#endif

